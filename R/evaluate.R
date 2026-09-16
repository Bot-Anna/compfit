# ============================================================
# evaluate.R
# Solve the ODE at a given (initial_state, parms) and evaluate every data and
# dummy formula on the resulting trajectory. Self-sufficient: reads the time
# grid, solver settings, formulas and dummies from the `fit` object, so it needs
# no global variables. Used by main.R for plotting and by plot_fit().
#
#   ev <- solve_and_evaluate(fit, initial_state, parms)
#   ev$sir_out      # raw trajectory (states by time) + date column
#   ev$evaluation   # one column per data/dummy formula + date
#
# `data_dummy` is optional; pass NULL (default) to evaluate only data streams.
# ============================================================

# Recompute the model's derived quantities on an already-solved trajectory.
#
# `spec` is the `derived_spec` built by compartmentalFunction() (raw Functions,
# level membership, startpoint/cutoff). Returns a data frame with one column per
# level head count (N1..Nk, N_<level>, total_pop) and per Functions-column entry,
# each a vector over the time grid -- or NULL if there is nothing to add.
#
# Faithful to the ODE by construction: the RAW Functions (original state names,
# `time`) are replayed as plain R against the named state columns (which are
# already vectors over the grid), in declaration order, after binding the SAME
# level heads and startpoint/cutoff the generated model uses. `if_else`/`ifelse`
# and friends resolve because the eval environment's parent is the package
# namespace (where `if_else` lives).
# Parse a vector of raw `name <- rhs` Function strings into (expr, lhs) pairs,
# dropping anything that is not a valid assignment. Pure parsing -- no eval -- so
# it is safe to do ONCE and reuse across trajectories.
.derived_parse <- function(functions_raw) {
  out <- list()
  if (is.null(functions_raw)) return(out)
  for (f in functions_raw) {
    ex <- tryCatch(parse(text = f)[[1]], error = function(e) NULL)
    if (is.null(ex) || !is.call(ex) || !as.character(ex[[1]]) %in% c("<-", "=")) next
    out[[length(out) + 1L]] <- list(ex = ex, lhs = deparse(ex[[2]]))
  }
  out
}

# Attach pre-parsed Functions to a derived spec so the per-eval replay in the
# loss never re-parses. Returns the spec augmented with `$parsed_functions`;
# NULL passes through. Call once (e.g. at loss-closure build time).
.derived_prepare <- function(spec) {
  if (is.null(spec)) return(NULL)
  spec$parsed_functions <- .derived_parse(spec$functions_raw)
  spec
}

# ---- Observable-referenced Functions (Bayesian backends) -------------------
# The Julia/Stan likelihoods evaluate each data-stream Formula on the solved
# trajectory. A Formula may name a TIME-VARYING Function (e.g. observed
# diagnoses = f_R_CH(t)*(R_H+R_L)), which is neither a state nor a scalar
# parameter, so it must be recomputed as a vector over the solve grid.
#
# Returns the subset of `spec`'s Functions that (a) are state-INDEPENDENT (their
# rhs uses only parameters/time -- the realistic case for time-varying rates, and
# the only case a grid recompute can do without mapping state trajectories), and
# (b) are reachable from `formulas`, following Function-to-Function references so
# a rate defined via helper Functions still resolves. Declaration order is
# preserved, so emitting them in order is well-defined. Empty when no observable
# references a Function -- callers then emit nothing and behave exactly as before.
.derived_observable_functions <- function(spec, formulas, param_symbols = character(0)) {
  if (is.null(spec) || !length(spec$functions_raw)) return(setNames(character(0), character(0)))
  state_syms <- unique(c(spec$comp_names,
                         paste0("N", seq_along(spec$level_members)),
                         if (length(spec$level_names)) paste0("N_", spec$level_names[nzchar(spec$level_names)]),
                         "total_pop"))
  cls <- .classify_functions(spec$functions_raw, state_syms, param_symbols)
  indep <- cls$independent                      # named rhs vector, dependency order
  if (!length(indep)) return(setNames(character(0), character(0)))

  vars_of <- function(txt) {
    txt <- gsub("^\\s*(annual|cumulative)\\((.*)\\)\\s*$", "\\2", gsub("`", "", txt))
    tryCatch(all.vars(parse(text = txt)[[1]]), error = function(e) character(0))
  }
  # Seed with the names the observables mention, then close over Function refs.
  need <- unique(unlist(lapply(formulas, vars_of)))
  repeat {
    hit <- intersect(names(indep), need)
    extra <- unique(unlist(lapply(indep[hit], vars_of)))
    new <- setdiff(intersect(extra, names(indep)), need)
    if (!length(new)) break
    need <- c(need, new)
  }
  keep <- names(indep) %in% need
  indep[keep]
}

#' Replay a model's derived quantities on a solved trajectory
#'
#' Internal, but EXPORTED deliberately: [extract_code()] emits a self-contained
#' script whose loss body calls `.derived_columns()` unqualified after a plain
#' `library(compfit)`, so it must be reachable from outside the namespace. The
#' call there sits inside a `tryCatch(error = NULL)` that degrades to
#' states-and-parameters only, so un-exporting this would not raise an error in
#' the extracted script -- it would silently drop every derived column from the
#' loss. Not part of the stable user-facing API.
#'
#' @param sir_out Solved trajectory: one column per state, over the time grid.
#' @param parms Named numeric vector of parameters.
#' @param time Numeric time grid.
#' @param spec The `derived_spec` built by `compartmentalFunction()`.
#' @return A data frame of derived columns, or `NULL` if there is nothing to add.
#' @keywords internal
#' @export
.derived_columns <- function(sir_out, parms, time, spec) {
  if (is.null(spec) || !length(spec$functions_raw) && !length(spec$level_members))
    return(NULL)
  n   <- nrow(sir_out)
  env <- new.env(parent = environment(.derived_columns))

  # State trajectory columns (vectors) + every parameter as a bare name.
  for (nm in names(sir_out)) assign(nm, suppressWarnings(as.numeric(sir_out[[nm]])), envir = env)
  if (!is.null(names(parms)))
    for (nm in names(parms)) assign(nm, parms[[nm]], envir = env)
  # Grid scalars the Functions may reference.
  assign("time", as.numeric(time), envir = env)
  assign("t",    as.numeric(time), envir = env)
  assign("startpoint", spec$startpoint, envir = env)
  assign("cutoff",     spec$cutoff,     envir = env)

  # Level head counts: N1..Nk, total_pop, and the named N_<level> aliases.
  comp     <- spec$comp_names
  head_nm  <- character(0)
  total    <- numeric(n)
  for (i in seq_along(spec$level_members)) {
    idx  <- spec$level_members[[i]]
    cols <- comp[idx]
    cols <- cols[cols %in% names(sir_out)]
    Ni   <- if (length(cols)) rowSums(as.data.frame(sir_out)[, cols, drop = FALSE]) else numeric(n)
    assign(paste0("N", i), Ni, envir = env); head_nm <- c(head_nm, paste0("N", i))
    if (length(spec$level_names) >= i && nzchar(spec$level_names[i])) {
      an <- paste0("N_", spec$level_names[i])
      assign(an, Ni, envir = env); head_nm <- c(head_nm, an)
    }
    total <- total + Ni
  }
  assign("total_pop", total, envir = env)

  # Replay each Functions entry in declaration order; record its LHS name.
  # Use pre-parsed (expr, lhs) pairs when the caller supplied a spec run through
  # .derived_prepare() -- the loss does this once so the hot per-eval path never
  # re-parses; the plotting path passes a raw spec and parses inline here.
  parsed <- spec$parsed_functions
  if (is.null(parsed)) parsed <- .derived_parse(spec$functions_raw)
  fn_nm <- character(0)
  for (pf in parsed) {
    # Per-entry guard: a Function that references something the replay doesn't
    # reconstruct is skipped (as are its dependents) without losing the entries
    # that DID resolve, e.g. time-only ones like q_CH / the sigmoids.
    ok <- tryCatch({ eval(pf$ex, envir = env); TRUE }, error = function(e) FALSE)
    if (ok) fn_nm <- c(fn_nm, pf$lhs)
  }

  # Assemble output: head counts + total_pop + Function results, as grid vectors
  # (scalars, e.g. `ratio <- n_H/n_L`, are recycled to the grid length).
  out <- unique(c(head_nm, "total_pop", fn_nm))
  cols <- list()
  for (nm in out) {
    v <- tryCatch(get(nm, envir = env), error = function(e) NULL)
    if (is.null(v) || !is.numeric(v)) next
    if (length(v) == 1L) v <- rep(v, n)
    if (length(v) == n)  cols[[nm]] <- as.numeric(v)
  }
  if (!length(cols)) return(NULL)
  as.data.frame(cols, check.names = FALSE)
}

#' Solve the ODE and evaluate formulas
#'
#' Solves the model at a given `(initial_state, parms)` and evaluates every data
#' (and optional dummy) formula on the resulting trajectory. Self-sufficient:
#' reads the time grid, solver settings and formulas from the `fit` object.
#'
#' @param fit A `"compartmentalFit"` object.
#' @param initial_state Named numeric vector of initial compartment values.
#' @param parms Named numeric vector of parameters.
#' @param data_dummy Optional dummy-data data frame; `NULL` evaluates only the
#'   data streams.
#' @return A list with `sir_out` (trajectory) and `evaluation` (formula columns).
#' @examples
#' \dontrun{
#' # fit from fitCompartmentalModel(); evaluate at the fitted point
#' p  <- get_point(fit)
#' ev <- solve_and_evaluate(fit, p$initial_state, p$parms)
#' head(ev$evaluation)
#' }
#' @export
solve_and_evaluate <- function(fit, initial_state, parms, data_dummy = NULL) {
  tg        <- fit$time_grid
  time      <- tg$time
  partition <- tg$partition          # passed explicitly to evaluate_formula()
  date      <- fit$model$date
  
  X <- unlist(initial_state)
  p <- unlist(parms)
  t <- c(as.numeric(min(time)), as.numeric(max(time)))
  
  sol <- if (identical(.ode_backend(fit$solver$backend), "r")) {
    solveWithR(X, t, p, time,
               abstol = fit$solver$abstol, reltol = fit$solver$reltol,
               method = .desolve_method(fit$solver$solver))
  } else {
    solveWithJulia(X, t, p, time,
                   solver = fit$solver$solver,
                   abstol = fit$solver$abstol,
                   reltol = fit$solver$reltol)
  }
  sir_out <- as.data.frame(t(sol$matrix))
  colnames(sir_out) <- names(initial_state)
  
  # Guard: a failed/aborted ODE solve returns a trajectory whose length does not
  # match the time grid. Catch it here with a legible message rather than letting
  # the date-column assignment fail cryptically ("replacement has N rows ...").
  if (nrow(sir_out) != length(date)) {
    stop(sprintf(
      "ODE solve returned %d time points but the grid has %d. The solve likely ",
      nrow(sir_out), length(date)),
      "failed/aborted at these parameter values (often a stiff/unstable region, ",
      "e.g. plotting at a non-converged posterior mean). Check the parameter ",
      "values and the solver warnings above.")
  }
  sir_out$date <- date

  # Enrich the trajectory with the model's DERIVED quantities -- every
  # Functions-column entry (q_CH, omega_*, p_*, sigmoids, ...) and the level head
  # counts (N1..Nk, N_<level>, total_pop) -- recomputed on the solved trajectory.
  # This lets data/dummy formulas reference those names directly, exactly as the
  # ODE body does. Wrapped so an exotic model that fails the replay degrades to
  # the plain trajectory (previous behaviour) rather than breaking plotting.
  derived <- tryCatch(.derived_columns(sir_out, parms, time, fit$model$derived),
                      error = function(e) NULL)
  if (!is.null(derived) && ncol(derived)) {
    add <- setdiff(names(derived), names(sir_out))   # never clobber a state/date col
    if (length(add)) sir_out <- cbind(sir_out, derived[, add, drop = FALSE])
  }

  streams <- fit$data$names_data_points
  ev <- data.frame(date = sir_out$date)
  for (col in streams) {
    ev[[col]] <- evaluate_formula(col, sir_out, parms, time, partition)
  }
  if (!is.null(data_dummy)) {
    for (col in data_dummy$Formula) {
      ev[[col]] <- evaluate_formula(col, sir_out, parms, time, partition)
    }
  }
  
  list(sir_out = sir_out, evaluation = ev)
}