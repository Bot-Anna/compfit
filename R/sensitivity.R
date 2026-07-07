# ============================================================
# sensitivity.R
# Sensitivity analysis -- GLOBAL (variance-based Sobol over the parameter box)
# and LOCAL (derivative-based, at the fitted point); the two are complementary.
#
# ===================== GLOBAL: Sobol =========================
# Variance-based (Sobol) global sensitivity analysis.
#
# Sobol needs only the MODEL MACHINERY (a loss / observable map) and the
# parameter [lower, upper] ranges -- NOT a fit. It explores the whole parameter
# box and ignores any fitted estimates. The functions therefore accept EITHER:
#   * a `compartmentalModel` from build_compartmental_model()  (no fit needed), or
#   * a `compartmentalFit`  from fitCompartmentalModel()       (uses its machinery
#       only; the fitted point/posterior are not used).
#
# Convenience: if you pass model parameters + data directly, the functions build
# the machinery for you (a build-only step -- no optimisation or sampling).
#
#   sobol_loss(x, n = ...)      Sobol indices for the LOSS
#   sobol_streams(x, n = ...)   Sobol indices per STREAM output
#
# Both return a tidy data.frame: output, parameter, first_order, total.
#
# COST: N*(d+2) model evaluations (per stream for sobol_streams), each an ODE
# solve. Start small. Uniform sampling over the box can hit ODE-unstable regions
# (-> NA outputs); heavy NA rates make indices unreliable.
# ============================================================

# Resolve the machinery from whatever the user passed. Accepts a
# compartmentalModel, a compartmentalFit, or (modelParams, dataCombined) to
# build one. Returns an object exposing $loss, $bounds, $model, $time_grid,
# $data (the common fields both fit and model share).
.sobol_machinery <- function(x, dataCombined = NULL, solver = NULL) {
  if (inherits(x, "compartmentalModel") || inherits(x, "compartmentalFit"))
    return(x)
  # Otherwise treat x as modelParams and build (build-only; no fit).
  if (is.null(dataCombined))
    stop("Pass a compartmentalModel/compartmentalFit, or modelParams + dataCombined.")
  if (is.null(solver)) solver <- solver_control()
  build_compartmental_model(x, dataCombined, solver = solver)
}

# Tidy a sensitivity::sobol* object into first-order (S) and total (T) indices.
.sobol_tidy <- function(si, pnames, output_name) {
  S <- if (!is.null(si$S)) si$S[, 1] else rep(NA, length(pnames))
  T <- if (!is.null(si$T)) si$T[, 1] else rep(NA, length(pnames))
  data.frame(output = output_name, parameter = pnames,
             first_order = as.numeric(S), total = as.numeric(T),
             stringsAsFactors = FALSE)
}

# ---- Sobol on the LOSS -----------------------------------------------------
# The loss closure expects a NORMALISED [0,1] par vector (its bounds ARE [0,1]),
# so we sample the unit cube directly and call the loss.
#' Sobol sensitivity of the loss
#'
#' Variance-based (Sobol) global sensitivity of the loss function to the fitted
#' parameters, via [sensitivity::soboljansen()].
#'
#' @param x A `"compartmentalFit"` or `"compartmentalModel"` object.
#' @param n Base sample size for the Sobol design.
#' @param dataCombined Optional override data (defaults to the fit's data).
#' @param solver Optional solver settings (defaults to the fit's solver).
#' @param progress Show a progress bar.
#' @param ... Passed to [sensitivity::soboljansen()].
#' @return A tidy data frame of first-order and total indices per parameter.
#' @examples
#' \dontrun{
#' sobol_loss(fit, n = 500)   # fit or model from build_compartmental_model()
#' }
#' @export
sobol_loss <- function(x, n = 1000, dataCombined = NULL, solver = NULL,
                       progress = TRUE, ...) {
  if (!requireNamespace("sensitivity", quietly = TRUE))
    stop("sobol_loss() requires the 'sensitivity' package (install.packages('sensitivity')).")
  mac <- .sobol_machinery(x, dataCombined, solver)
  if (is.null(mac$loss)) stop("No loss function found in the supplied object.")

  d <- length(mac$bounds$lower)
  pnames <- names(mac$bounds$lower)
  if (is.null(pnames)) pnames <- paste0("p", seq_len(d))

  X1 <- matrix(runif(n * d), n, d, dimnames = list(NULL, pnames))
  X2 <- matrix(runif(n * d), n, d, dimnames = list(NULL, pnames))

  # Progress over the full Sobol design (soboljansen evaluates n*(d+2) points,
  # each an ODE solve via the loss -- the expensive part).
  pb <- NULL; cnt <- 0L
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)
  best_state <- new.env(); best_state$error <- Inf; best_state$par <- NULL
  model_fun <- function(X) {
    if (isTRUE(progress) && is.null(pb))
      pb <<- utils::txtProgressBar(min = 0, max = nrow(X), style = 3)
    apply(X, 1, function(par) {
      v <- tryCatch(mac$loss(as.numeric(par), best_state, FALSE),
                    error = function(e) NA_real_)
      cnt <<- cnt + 1L
      if (isTRUE(progress)) utils::setTxtProgressBar(pb, cnt)
      if (!is.finite(v)) NA_real_ else v
    })
  }

  si <- sensitivity::soboljansen(model = model_fun, X1 = X1, X2 = X2, ...)
  .sobol_tidy(si, pnames, "loss")
}

# ---- Sobol per STREAM ------------------------------------------------------
# Output = each stream's trajectory summarised to a scalar (default: mean).
# Sampling is over the fitted params/states [lower, upper] (natural scale).
#' Sobol sensitivity of model outputs
#'
#' Variance-based (Sobol) global sensitivity of a summary of each model output
#' stream to the fitted parameters.
#'
#' @param x A `"compartmentalFit"` or `"compartmentalModel"` object.
#' @param n Base sample size for the Sobol design.
#' @param summary_fun Function summarising each output trajectory to a scalar.
#' @param streams Optional subset of output streams to analyse.
#' @param dataCombined Optional override data (defaults to the fit's data).
#' @param solver Optional solver settings (defaults to the fit's solver).
#' @param progress Show a progress bar.
#' @param ... Passed to [sensitivity::soboljansen()].
#' @return A tidy data frame of first-order/total indices per output and parameter.
#' @examples
#' \dontrun{
#' sob <- sobol_streams(fit, n = 500)
#' }
#' @export
sobol_streams <- function(x, n = 500, summary_fun = mean, streams = NULL,
                          dataCombined = NULL, solver = NULL, progress = TRUE, ...) {
  if (!requireNamespace("sensitivity", quietly = TRUE))
    stop("sobol_streams() requires the 'sensitivity' package.")
  if (!exists("solve_and_evaluate"))
    stop("sobol_streams() needs solve_and_evaluate() in scope.")
  mac <- .sobol_machinery(x, dataCombined, solver)
  
  sap <- mac$model$sap
  tg  <- mac$time_grid
  lo <- c(sap$lower_states, sap$lower_params)
  hi <- c(sap$upper_states, sap$upper_params)
  pnames <- names(lo)
  d <- length(lo)
  if (is.null(streams)) streams <- mac$data$names_data_points
  
  # solve_and_evaluate() reads the time grid/solver from its first argument via
  # the same field names on both fit and model objects, so `mac` works directly.
  eval_one <- function(u) {
    par_nat <- lo + u * (hi - lo)
    names(par_nat) <- pnames
    pe <- tryCatch(.recover_solution(par_nat, sap, tg), error = function(e) NULL)
    if (is.null(pe)) return(setNames(rep(NA_real_, length(streams)), streams))
    ev <- tryCatch(solve_and_evaluate(mac, pe$initial_state, pe$parms)$evaluation,
                   error = function(e) NULL)
    if (is.null(ev)) return(setNames(rep(NA_real_, length(streams)), streams))
    vapply(streams, function(st) summary_fun(ev[[st]], na.rm = TRUE), numeric(1))
  }
  
  X1 <- matrix(runif(n * d), n, d, dimnames = list(NULL, pnames))
  X2 <- matrix(runif(n * d), n, d, dimnames = list(NULL, pnames))
  
  res <- list()
  for (i in seq_along(streams)) {
    st <- streams[i]
    if (isTRUE(progress))
      cat(sprintf("Sobol stream %d/%d: %s\n", i, length(streams), st))
    pb <- NULL; cnt <- 0L
    mf <- function(X) {
      if (isTRUE(progress) && is.null(pb))
        pb <<- utils::txtProgressBar(min = 0, max = nrow(X), style = 3)
      apply(X, 1, function(u) {
        val <- eval_one(u)[st]
        cnt <<- cnt + 1L
        if (isTRUE(progress)) utils::setTxtProgressBar(pb, cnt)
        val
      })
    }
    si <- tryCatch(sensitivity::soboljansen(model = mf, X1 = X1, X2 = X2, ...),
                   error = function(e) NULL)
    if (isTRUE(progress) && !is.null(pb)) { close(pb); cat("\n") }
    if (!is.null(si)) res[[st]] <- .sobol_tidy(si, pnames, st)
  }
  do.call(rbind, res)
}

# ---- sobol_report: inspect a sobol_streams result --------------------------
# Accepts the tidy data.frame from sobol_streams(), OR a fit/model (+ n) to
# compute it first. Prints a readable summary and returns the pieces invisibly.
#
#   sobol_report(sob)                      # from an existing result
#   sobol_report(fit, n = 500)             # compute then report
#
# Reports, per the Sobol math:
#   * VALIDITY: indices must lie in [0,1]; out-of-range or heavy NA => the
#     estimator did not converge (commonly the unstable full-box problem).
#   * Per-stream parameter ranking by total index.
#   * INTERACTION structure: gap (total - first_order) per parameter; sum(S)<=1
#     and sum(T)>=1 per stream (sum(S) near 1 => additive; well below => interactive).
#   * Global parameter ranking (mean total across streams); near-zero => candidate
#     to fix / poorly identified (cf. identifiability_report).
#   * Optional parameter x stream heatmap of the total index.
#' Summarise a Sobol sensitivity analysis
#'
#' Reports validity, per-stream top parameters, interaction structure, and a
#' global parameter ranking from a Sobol analysis. Accepts a precomputed tidy
#' table or a fit/model (computed via [sobol_streams()]).
#'
#' @param x A tidy Sobol data frame (output, parameter, first_order, total) or a
#'   `"compartmentalFit"`/`"compartmentalModel"` object.
#' @param n Base sample size when computing from a fit/model.
#' @param cor_to Optional parameter name to correlate others against.
#' @param plots If `TRUE`, build a total-index heatmap.
#' @param top_n Number of top parameters to report per stream.
#' @param ... Passed to [sobol_streams()] when computing.
#' @return Invisibly, a list with the table, stream sums, parameter ranking,
#'   not-converged flags, range/NA diagnostics, and (if `plots`) a heatmap.
#' @examples
#' \dontrun{
#' sobol_report(fit, n = 500)            # compute from a fit, or:
#' sob <- sobol_streams(fit, n = 500)
#' r   <- sobol_report(sob)              # report a precomputed table
#' }
#' @export
sobol_report <- function(x, n = 500, cor_to = NULL, plots = TRUE,
                         top_n = 5, ...) {
  sob <- if (is.data.frame(x) &&
             all(c("output","parameter","first_order","total") %in% names(x))) x
         else sobol_streams(x, n = n, ...)
  if (is.null(sob) || !nrow(sob)) { message("No Sobol results to report."); return(invisible(NULL)) }

  streams <- unique(sob$output)
  pnames  <- unique(sob$parameter)
  sob$interaction <- sob$total - sob$first_order

  oob <- function(v) sum(v < -1e-6 | v > 1 + 1e-6, na.rm = TRUE)
  n_oob   <- oob(sob$first_order) + oob(sob$total)
  na_frac <- mean(is.na(sob$total))

  cat("=== Sobol sensitivity report ===\n")
  cat(sprintf("streams: %d   parameters: %d   rows: %d\n",
              length(streams), length(pnames), nrow(sob)))

  # VALIDITY
  cat("\n[validity]\n")
  cat(sprintf("  indices outside [0,1]: %d   |   NA total indices: %.0f%%\n",
              n_oob, 100 * na_frac))
  bad_streams <- vapply(streams, function(st) {
    d <- sob[sob$output == st, ]
    oob(d$first_order) + oob(d$total) > 0 || mean(is.na(d$total)) > 0.2
  }, logical(1))
  if (any(bad_streams)) {
    cat("  NOT CONVERGED for: ", paste(streams[bad_streams], collapse = ", "), "\n")
    cat("  -> these reflect estimator breakdown (often the unstable full-box).\n",
        "     Increase n, use the posterior-range variant, or prefer local_sensitivity().\n", sep = "")
  } else cat("  all streams in range; indices look usable.\n")

  # PER-STREAM TOP PARAMETERS (by total)
  cat("\n[per-stream top parameters by total index]\n")
  for (st in streams) {
    d <- sob[sob$output == st, ]
    d <- d[order(-d$total), ]
    k <- seq_len(min(top_n, nrow(d)))
    cat(sprintf("  %s:\n", st))
    for (i in k)
      cat(sprintf("      %-14s  S=%+.3f  T=%+.3f  (interaction %+.3f)\n",
                  d$parameter[i], d$first_order[i], d$total[i], d$interaction[i]))
  }

  # INTERACTION STRUCTURE (per-stream sums)
  cat("\n[structure: per-stream index sums]\n")
  agg <- aggregate(cbind(first_order, total) ~ output, sob, sum, na.rm = TRUE)
  for (i in seq_len(nrow(agg))) {
    tag <- if (agg$first_order[i] > 0.9) "additive"
           else if (agg$first_order[i] < 0.6) "interaction-heavy" else "mixed"
    cat(sprintf("  %-16s  sum(S)=%.2f  sum(T)=%.2f   [%s]\n",
                agg$output[i], agg$first_order[i], agg$total[i], tag))
  }
  cat("  (sum(S) near 1 => additive; well below => interactions matter)\n")

  # GLOBAL PARAMETER RANKING (mean total across streams)
  cat("\n[global parameter ranking: mean total across streams]\n")
  pr <- aggregate(total ~ parameter, sob, mean, na.rm = TRUE)
  pr <- pr[order(-pr$total), ]
  for (i in seq_len(min(10, nrow(pr))))
    cat(sprintf("      %-14s  %.3f\n", pr$parameter[i], pr$total[i]))
  negligible <- pr$parameter[pr$total < 0.02]
  if (length(negligible))
    cat("  negligible everywhere (mean total < 0.02; candidates to fix): ",
        paste(negligible, collapse = ", "), "\n")

  # HEATMAP
  heat <- NULL
  if (isTRUE(plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    heat <- ggplot2::ggplot(sob, ggplot2::aes(parameter, output, fill = total)) +
      ggplot2::geom_tile(colour = "grey90") +
      ggplot2::scale_fill_viridis_c(limits = c(0, 1), oob = scales::squish) +
      ggplot2::labs(title = "Total Sobol index: parameter x stream",
                    x = NULL, y = NULL, fill = "total") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  invisible(list(table = sob,
                 stream_sums = agg,
                 parameter_ranking = pr,
                 not_converged = streams[bad_streams],
                 n_outside_range = n_oob, na_fraction = na_frac,
                 heatmap = heat))
}
# ============================================================
# ===================== LOCAL: derivative =====================
# LOCAL (derivative-based) sensitivity at the fitted point.
#
# Unlike Sobol (global, variance over the whole box), local sensitivity is the
# partial derivative of an output wrt each parameter, evaluated AT the fitted
# estimate -- "if I nudge this parameter slightly, how much does the output
# move, here?" Computed by CENTRAL finite differences on the normalised [0,1]
# scale (perturb +/- eps, re-solve, difference). Cost: ~2*d ODE solves -- cheap,
# and always feasible (you evaluate at a known-good point, no unstable-box NAs).
#
#   local_sensitivity(fit, type = "relative", output = "both",
#                     summary_fun = mean, eps = 1e-4, n_draws = NULL)
#
# Returns a list:
#   $summary    tidy df: output, parameter, sensitivity  (one number per
#               stream x parameter, from summary_fun over the trajectory)
#   $trajectory tidy df: output, parameter, date, sensitivity  (the derivative
#               as a function of time -- shows WHEN each parameter matters)
#
# type (normalisation of the derivative dg/dtheta):
#   "absolute" -> dg/dtheta            (raw units)
#   "relative" -> (dg/dtheta)*(theta/g) (elasticity; dimensionless, comparable)
#   "semi"     -> (dg/dtheta)*theta     (scaled by parameter only)
#
# Plots: plot_local_sensitivity(ls, kind = "tornado" | "trajectory").
# ============================================================

# Central-difference point set on the NORMALISED scale, evaluated via the loss
# machinery's denormalisation + solve. Returns the fitted point and the per-
# parameter perturbed evaluations.
#' Local (derivative) sensitivity at a fit
#'
#' Finite-difference sensitivity of each model output to the fitted parameters,
#' evaluated at the fitted point. `type` selects absolute derivatives, relative
#' elasticities, or semi-relative (parameter-scaled) sensitivities.
#'
#' @param fit A `"compartmentalFit"` object.
#' @param type Sensitivity type: relative (elasticity), absolute, or semi.
#' @param output What to return: both, summary, or trajectory.
#' @param summary_fun Function summarising each output trajectory to a scalar.
#' @param eps Finite-difference step on the normalised scale.
#' @param data_dummy Optional dummy-data data frame.
#' @param progress Show a progress bar.
#' @return A list with `summary` and/or `trajectory` tidy data frames.
#' @examples
#' \dontrun{
#' ls <- local_sensitivity(fit, type = "relative")   # fit from fitCompartmentalModel()
#' }
#' @export
local_sensitivity <- function(fit, type = c("relative", "absolute", "semi"),
                              output = c("both", "summary", "trajectory"),
                              summary_fun = mean, eps = 1e-4,
                              data_dummy = NULL, progress = TRUE) {
  type   <- match.arg(type); output <- match.arg(output)

  # The point to evaluate at, in NATURAL scale (named params + states), and the
  # bounds so we can perturb on the normalised [0,1] scale.
  if (inherits(fit, "compartmentalFit") && fit$method == "bayes") {
    pm    <- vapply(posterior_draws(fit), median, numeric(1), na.rm = TRUE)
    point <- .cfit_point_from_draw(fit, as.data.frame(t(pm)))
    nat   <- c(point$parms, point$initial_state)
  } else if (inherits(fit, "compartmentalFit")) {
    nat <- c(unlist(fit$point$parms), unlist(fit$point$initial_state))
    point <- list(parms = fit$point$parms, initial_state = fit$point$initial_state)
  } else {
    stop("local_sensitivity() expects a compartmentalFit.")
  }

  lo <- fit$bounds$lower; hi <- fit$bounds$upper
  fitted_names <- names(lo)                     # the FITTED quantities (others fixed)
  if (is.null(fitted_names) || !length(fitted_names))
    stop("No fitted quantities to differentiate (all parameters are fixed).")

  streams <- fit$data$names_data_points
  dummies <- if (!is.null(data_dummy)) data_dummy$Formula else character(0)
  outcols <- c(streams, dummies)

  # Evaluate the model at a natural-scale named vector -> evaluation data frame.
  eval_at <- function(natvec) {
    pe <- .ls_split_point(natvec, point)
    out <- tryCatch(solve_and_evaluate(fit, pe$initial_state, pe$parms, data_dummy),
                    error = function(e) NULL)
    if (is.null(out)) NULL else out$evaluation
  }

  base_eval <- eval_at(nat)
  if (is.null(base_eval)) stop("Baseline solve failed at the fitted point.")
  dates <- base_eval$date

  d <- length(fitted_names)
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = d, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  # For each fitted parameter: perturb +/- eps on the NORMALISED scale, map back
  # to natural scale, re-solve, central-difference each output trajectory.
  traj_rows <- list()
  for (k in seq_len(d)) {
    nm <- fitted_names[k]
    span <- hi[[nm]] - lo[[nm]]
    h_nat <- eps * span                          # natural-scale step (eps on [0,1])
    if (h_nat == 0) h_nat <- eps

    up <- nat; up[[nm]] <- nat[[nm]] + h_nat
    dn <- nat; dn[[nm]] <- nat[[nm]] - h_nat
    eu <- eval_at(up); ed <- eval_at(dn)
    if (isTRUE(progress)) utils::setTxtProgressBar(pb, k)
    if (is.null(eu) || is.null(ed)) next         # skip params whose perturbation fails

    for (st in outcols) {
      gu <- eu[[st]]; gd <- ed[[st]]; gb <- base_eval[[st]]
      if (is.null(gu) || is.null(gd) || is.null(gb)) next
      deriv <- (gu - gd) / (2 * h_nat)           # dg/dtheta (natural scale)
      sens  <- switch(type,
        absolute = deriv,
        relative = deriv * nat[[nm]] / ifelse(gb == 0, NA, gb),
        semi     = deriv * nat[[nm]])
      traj_rows[[length(traj_rows) + 1]] <-
        data.frame(output = st, parameter = nm, date = dates,
                   sensitivity = sens, stringsAsFactors = FALSE)
    }
  }
  trajectory <- if (length(traj_rows)) do.call(rbind, traj_rows) else
                  data.frame(output = character(), parameter = character(),
                             date = as.Date(character()), sensitivity = numeric())

  # Summary: collapse each (output, parameter) trajectory with summary_fun.
  summary_df <- NULL
  if (output %in% c("both", "summary") && nrow(trajectory)) {
    summary_df <- stats::aggregate(sensitivity ~ output + parameter,
                                   data = trajectory,
                                   FUN = function(v) summary_fun(v, na.rm = TRUE))
  }

  structure(list(
    summary    = if (output %in% c("both", "summary")) summary_df else NULL,
    trajectory = if (output %in% c("both", "trajectory")) trajectory else NULL,
    type = type, streams = streams, dummies = dummies,
    dataCombined = fit$data$data_combined, data_dummy = data_dummy),
    class = "localSensitivity")
}

# Replace fitted entries of a base point with values from a natural-scale named
# vector, keeping fixed params/states from the base point.
.ls_split_point <- function(natvec, base) {
  parms <- as.list(base$parms); states <- as.list(base$initial_state)
  for (nm in names(natvec)) {
    if (nm %in% names(parms))  parms[[nm]]  <- natvec[[nm]]
    if (nm %in% names(states)) states[[nm]] <- natvec[[nm]]
  }
  list(parms = unlist(parms), initial_state = unlist(states))
}

# ---- plots -----------------------------------------------------------------
# kind = "tornado": bar chart of summary sensitivity per parameter, faceted by
#                   output, sorted by magnitude.
# kind = "trajectory": sensitivity-over-time curves (one line per parameter),
#                   faceted by output.
#' Plot local sensitivity
#'
#' Visualises a [local_sensitivity()] result as a tornado plot (ranked summary
#' sensitivities) or a trajectory plot (sensitivity over time).
#'
#' @param ls A [local_sensitivity()] result.
#' @param kind Plot kind: tornado or trajectory.
#' @param top Optional number of top parameters to show.
#' @param base_size Base font size.
#' @param ncol Number of facet columns (default: automatic).
#' @return A ggplot object (or a list of them).
#' @examples
#' \dontrun{
#' plot_local_sensitivity(ls, kind = "tornado")   # ls from local_sensitivity()
#' }
#' @export
plot_local_sensitivity <- function(ls, kind = c("tornado", "trajectory"),
                                    top = NULL, base_size = 11, ncol = NULL) {
  kind <- match.arg(kind)
  if (!inherits(ls, "localSensitivity")) stop("expects a local_sensitivity() result.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("needs ggplot2.")

  if (kind == "tornado") {
    df <- ls$summary
    if (is.null(df) || !nrow(df)) stop("No summary sensitivities to plot.")
    df$absS <- abs(df$sensitivity)
    if (!is.null(top)) {
      df <- do.call(rbind, by(df, df$output, function(d)
        d[order(-d$absS), ][seq_len(min(top, nrow(d))), ]))
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(parameter, absS),
                                          y = sensitivity)) +
      ggplot2::geom_col(fill = .cfit_active_palette()$model, alpha = 0.85) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ output, scales = "free", ncol = ncol) +
      ggplot2::labs(x = NULL, y = paste0(ls$type, " local sensitivity"),
                    title = "Local sensitivity at the fitted point") +
      cfit_theme(base_size)
    return(p)
  }

  # trajectory
  df <- ls$trajectory
  if (is.null(df) || !nrow(df)) stop("No trajectory sensitivities to plot.")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = sensitivity,
                                        colour = parameter, group = parameter)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60",
                        linetype = "dotted") +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.85) +
    ggplot2::facet_wrap(~ output, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = NULL, y = paste0(ls$type, " local sensitivity"),
                  colour = NULL,
                  title = "Local sensitivity over time at the fitted point") +
    cfit_theme(base_size)
  p
}
