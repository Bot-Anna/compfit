# ============================================================
# fitCompartmentalModel.R
# Master function tying together the model-building pipeline:
#   numberOfComps -> statesAndParams -> generateExpressions
#     -> compartmentalFunction -> registerJuliaODEFunction
#     -> lossFunction -> {optim | deoptim | hypercube | bayes}
#
# Returns an object of class "compartmentalFit".
# ============================================================

# ---- Control constructors --------------------------------------------------

# Solver settings passed through to solveWithJulia(). Defaulting the solver
# here (rather than per-call) keeps fitting and plotting consistent.
#' ODE solver settings
#'
#' Bundle the Julia ODE solver and tolerances used for fitting and plotting.
#'
#' @param solver ODE solver. Its meaning is backend-specific, and the default
#'   `"AutoTsit5(Rosenbrock23())"` maps to each backend's auto-switching /
#'   non-stiff default, which suits almost every compartmental model.
#'   \itemize{
#'     \item **Julia**: the `OrdinaryDiffEq` constructor expression, passed
#'       verbatim. Common choices: `"Tsit5()"` (non-stiff), `"Vern7()"`/`"Vern9()"`
#'       (high-accuracy), `"Rosenbrock23()"`/`"Rodas5()"`/`"KenCarp4()"`/`"TRBDF2()"`
#'       (stiff), `"AutoVern7(Rodas5())"` (auto-switching). Any `OrdinaryDiffEq`
#'       solver works; see the
#'       \href{https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/}{DifferentialEquations.jl docs}.
#'     \item **R** (deSolve): a `deSolve::ode` `method` name, e.g. `"lsoda"`
#'       (default, auto-switching), `"radau"`/`"bdf"` (stiff), `"rk4"`/`"ode45"`
#'       (non-stiff). Unrecognised strings fall back to `"lsoda"`.
#'     \item **Stan**: `"rk45"` (default, non-stiff), `"bdf"` (stiff), or
#'       `"adams"` (non-stiff multistep). Set via `solver = "bdf"` etc.
#'   }
#' @param abstol Absolute tolerance (all backends).
#' @param reltol Relative tolerance (all backends).
#' @param backend Solver backend: `"julia"` (default; fast, needs a Julia
#'   session via [setup_julia()]), `"r"` (pure R via deSolve; no Julia required,
#'   but slower), or `"stan"` (Bayesian sampling via Stan/rstan NUTS -- no Julia;
#'   all R-side solves, e.g. plotting and the optional MAP init, use deSolve).
#' @return A named list of solver settings.
#' @examples
#' solver_control(abstol = 1e-6, reltol = 1e-6)
#' solver_control(backend = "r")   # Julia-free
#' @export
solver_control <- function(solver  = "AutoTsit5(Rosenbrock23())",
                           abstol  = 1e-8,
                           reltol  = 1e-8,
                           backend = c("julia", "r", "stan")) {
  backend <- match.arg(backend)
  list(solver = solver, abstol = abstol, reltol = reltol, backend = backend)
}

# Which ODE SOLVER a backend uses for R-side solving (loss, plotting, MAP init,
# predictive). "stan" has no R-side solver of its own -- its Stan program solves
# internally during sampling -- so every R-side solve falls back to deSolve, so
# a Stan user never needs Julia. "julia"/"r" map to themselves.
.ode_backend <- function(backend) if (identical(backend, "stan")) "r" else backend

# Optimiser settings. Used by method = "optim" / "deoptim".
#' Optimiser settings
#'
#' Bundle optimiser settings for the maximum-likelihood methods. Two optimisers
#' share this control: the deterministic quasi-Newton **L-BFGS-B**
#' (`method = "lbfgsb"`, gradient by finite differences) and the stochastic
#' **DEoptim** differential-evolution search (`method = "deoptim"`, a global
#' method that is more robust to multi-modal / rugged objectives but slower).
#' The `maxit`/`factr` options apply to L-BFGS-B; `itermax`/`NP_mult`/`CR`/`F`/
#' `trace`/`seed` apply to DEoptim; `progress` applies to both. Settings for the
#' other optimiser are simply ignored, so one control object works for any
#' MLE method.
#'
#' @param opt_method The `optim()` method used on the `lbfgsb` path. Normally
#'   `"L-BFGS-B"` (box-constrained); `"Brent"` is the 1-D alternative. The search
#'   runs over a normalised \[0,1\] box, so the method must support box bounds.
#' @param maxit L-BFGS-B: maximum iterations before it stops (default 1000).
#' @param factr L-BFGS-B: convergence tolerance as a multiple of machine
#'   epsilon; smaller = stricter/slower (default 1e7, i.e. ~1e-8 relative).
#' @param itermax DEoptim: maximum number of generations (default 1000).
#' @param NP_mult DEoptim: population size as a multiple of the number of fitted
#'   quantities (`NP = NP_mult * d`). The DEoptim authors recommend >= 10
#'   (the default); lower is faster but explores less.
#' @param CR DEoptim: crossover probability in \[0,1\] (default 0.9).
#' @param F DEoptim: differential weighting (mutation) factor, typically in
#'   \[0,2\] (default 0.8).
#' @param trace DEoptim: print the best value every `trace` generations
#'   (default 50). Ignored when `progress = FALSE`.
#' @param seed Optional RNG seed for the stochastic `deoptim` method, for
#'   reproducible fits. `lbfgsb`/`hypercube` are already deterministic and
#'   ignore it -- except that a `|random` parameter start (see the modelParams
#'   grammar in `vignette("compfit")`) draws its L-BFGS-B start from this seed,
#'   so setting it makes that random warm start reproducible too.
#' @param progress Show optimisation progress. `TRUE` (default) prints the
#'   running best objective each time the loss improves and (for DEoptim) its
#'   per-generation `trace`. `FALSE` silences BOTH for a quiet fit -- the best
#'   solution is still tracked and recovered. (The Bayesian sampler's progress
#'   is controlled separately, by [bayes_control()]'s `progress`.)
#' @return A named list of optimiser settings, passed as `control =` to
#'   [fitCompartmentalModel()].
#' @seealso [fitCompartmentalModel()], [bayes_control()], [hypercube_control()].
#' @examples
#' optim_control(maxit = 500)                       # L-BFGS-B, 500 iterations
#' optim_control(opt_method = "DEoptim", itermax = 200)  # DE, 200 generations
#' optim_control(progress = FALSE)                  # quiet fit (no console output)
#' optim_control(seed = 1)                          # reproducible DEoptim
#' @export
optim_control <- function(opt_method = "L-BFGS-B",
                          maxit   = 1000,
                          factr   = 1e7,
                          # DEoptim-specific:
                          itermax = 1000,
                          NP_mult = 10,
                          CR      = 0.9,
                          F       = 0.8,
                          trace   = 50,
                          seed    = NULL,
                          progress = TRUE) {
  list(opt_method = opt_method, seed = seed,
       maxit = maxit, factr = factr,
       itermax = itermax, NP_mult = NP_mult,
       CR = CR, F = F, trace = trace, progress = progress)
}

# Hypercube quasi-Monte Carlo pre-search settings.
#' Hypercube pre-search settings
#'
#' Settings for the quasi-Monte Carlo (Sobol) hypercube pre-search
#' (`method = "hypercube"`), which evaluates the loss on a low-discrepancy Sobol
#' design over the normalised \[0,1\] box and keeps the best point. Deterministic
#' and derivative-free -- useful as a stand-alone coarse fit or to seed a warm
#' start for the local optimiser.
#'
#' @param n Number of Sobol sample points to evaluate (default 10000). More
#'   points give finer coverage of the parameter box at linear cost in solves.
#'   n should approximately be of magnitude 10^m where m is the number of unknown
#'   parameters and initial states.
#' @return A named list with the sample count, passed as `hypercube =` to
#'   [fitCompartmentalModel()].
#' @seealso [fitCompartmentalModel()], [optim_control()].
#' @examples
#' hypercube_control(n = 5000)
#' @export
hypercube_control <- function(n = 10000) {
  list(n = n)
}

# Bayesian settings (Julia / Turing side).
#' Bayesian sampling settings
#'
#' Bundle the settings for `method = "bayes"`. These apply across all three
#' Bayesian backends selected by [solver_control()]: Julia/Turing NUTS
#' (`backend = "julia"`), Stan/rstan NUTS (`backend = "stan"`), and the
#' gradient-free `BayesianTools` sampler (`backend = "r"`). Priors themselves are
#' read from the model sheet's `Parameters`/`States` columns (per-quantity boxes
#' or named distributions); the settings here govern the sampler, not the priors.
#'
#' @param sampler Julia/Turing sampler expression, e.g. `"NUTS(0.65)"` (the
#'   argument is the target acceptance rate). Used by `backend = "julia"` only;
#'   Stan uses its own NUTS and the R backend uses DEzs.
#' @param chains Number of independent MCMC chains (default 4). Multiple chains
#'   enable the R-hat convergence diagnostic.
#' @param iter Total iterations per chain, **including** warmup (default 2000).
#' @param warmup Warmup / burn-in iterations per chain, discarded before
#'   inference (default 1000); the kept sample size per chain is `iter - warmup`.
#' @param seed Optional RNG seed for reproducible sampling.
#' @param init_from_optim If `TRUE` (default), run a quick MLE fit first and
#'   initialise the sampler there (a MAP-style warm start). Applies to the Julia
#'   and R backends; the Stan backend relies on Stan's own warmup adaptation and
#'   ignores it.
#' @param progress Show sampler progress. `FALSE` fully silences it for every
#'   backend: the Turing `Sampling ... ETA` bar (via `Turing.setprogress!`), and
#'   Stan's per-iteration refresh, chain messages, and auto-opened progress
#'   window. (The MLE optimiser's progress is controlled separately, by
#'   [optim_control()]'s `progress`.)
#' @param sigma_prior Prior for the Gaussian/lognormal noise SD, as a Julia
#'   expression (default a half-Normal). Used by the Julia backend.
#' @param phi_prior Prior for the negative-binomial dispersion, as a Julia
#'   expression. Used by the Julia backend.
#' @param sigma_prior_stan,phi_prior_stan The same two priors written as **Stan**
#'   expressions (defaults `"normal(0, 1)"` and `"gamma(2, 0.2)"`). Used by the
#'   Stan backend, which cannot read the Julia-syntax versions above.
#' @param adapt_delta Stan NUTS target acceptance probability (default 0.8).
#'   Raise toward 1 (e.g. 0.95) to reduce divergent transitions. Stan backend
#'   only.
#' @return A named list of Bayesian settings, passed as `bayes =` to
#'   [fitCompartmentalModel()].
#' @seealso [fitCompartmentalModel()], [solver_control()], [optim_control()].
#' @examples
#' bayes_control(chains = 4, iter = 2000, warmup = 1000)
#' bayes_control(progress = FALSE)              # silent sampling
#' bayes_control(chains = 2, iter = 1000, seed = 1)  # reproducible, lighter run
#' @export
bayes_control <- function(sampler    = "NUTS(0.65)",
                          chains     = 4,
                          iter       = 2000,
                          warmup     = 1000,
                          seed       = NULL,
                          init_from_optim = TRUE,
                          progress   = TRUE,
                          sigma_prior = "truncated(Normal(0, 1), 0, Inf)",
                          phi_prior   = "Gamma(2, 5)",
                          sigma_prior_stan = "normal(0, 1)",
                          phi_prior_stan   = "gamma(2, 0.2)",
                          adapt_delta = 0.8) {
  list(sampler = sampler,
       chains = chains, iter = iter, warmup = warmup,
       seed = seed, init_from_optim = init_from_optim,
       progress = progress,
       sigma_prior = sigma_prior, phi_prior = phi_prior,
       sigma_prior_stan = sigma_prior_stan, phi_prior_stan = phi_prior_stan,
       adapt_delta = adapt_delta)
}

# ---- Internal: build the model from a modelParams sheet --------------------

.build_model <- function(modelParams, backend = "julia") {

  # Catch entry mistakes up front with a clear message (disable with
  # options(compfit.validate = FALSE) if a valid sheet is ever wrongly rejected).
  if (isTRUE(getOption("compfit.validate", TRUE)))
    validate_modelParams(modelParams)

  compartment_structure <- numberOfComps(modelParams)

  sap <- statesAndParams(modelParams)

  expressions <- generateExpressions(
    number_of_comps  = compartment_structure$number_of_comps,
    states_fitted    = sap$states_fitted,
    states_fixed     = sap$states_fixed,
    states_functions = sap$states_functions,
    params_fitted    = sap$params_fitted,
    params_fixed     = sap$params_fixed,
    params_functions = sap$params_functions,
    conditions       = modelParams$Conditions,
    comp_names       = compartment_structure$comp_names
  )

  cf <- compartmentalFunction(
    modelParams           = modelParams,
    compartment_structure = compartment_structure,
    sir_expression        = expressions$sir,
    return_expression     = expressions$ret
  )

  # Register the ODE function for the chosen backend so the loss can solve it:
  # the Julia function (solveWithJulia) or the R closure (solveWithR). "stan"
  # solves inside its own program, so its R-side solving uses deSolve.
  if (identical(.ode_backend(backend), "r")) {
    registerRODEFunction(cf$compartmental_function)
  } else {
    registerJuliaODEFunction(cf$julia_code)
  }

  list(
    structure   = compartment_structure,
    sap         = sap,
    expressions = expressions,
    compartmental_function = cf$compartmental_function,
    julia_code  = cf$julia_code,
    stan_code   = cf$stan_code,
    date        = cf$date,
    modelParams = modelParams      # original sheet, for fill_params / write_filled_params
  )
}

# ---- Internal: time grid from a modelParams sheet --------------------------

.time_grid <- function(modelParams) {
  time_info <- modelParams$Others
  time_info <- time_info[time_info != "" & !is.na(time_info)]
  time_info <- gsub(" ", "", time_info)
  time_info <- time_info[!is.na(time_info)]

  startpoint <- extract_param(time_info, "startpoint")
  endpoint   <- extract_param(time_info, "endpoint")
  partition  <- extract_param(time_info, "partition")
  cutoff     <- extract_param(time_info, "cutoff")
  if (length(cutoff) == 0 || all(is.na(cutoff))) cutoff <- Inf   # cutoff is optional

  if (!is.numeric(partition) || length(partition) != 1 || !is.finite(partition) || partition <= 0)
    stop("'partition' must be a single positive finite number (got: ", partition, ").")
  if (!is.numeric(startpoint) || !is.finite(startpoint))
    stop("'startpoint' must be a finite number (got: ", startpoint, ").")
  if (!is.numeric(endpoint) || !is.finite(endpoint))
    stop("'endpoint' must be a finite number (got: ", endpoint, ").")
  if (startpoint > endpoint)
    stop(sprintf("'startpoint' (%g) must be <= 'endpoint' (%g).", startpoint, endpoint))

  # The time axis is ANNUAL (integer years; snapshots at year-ends). Warn if a
  # start/end entry is not a plain integer -- e.g. a date "2015-06-01" or a
  # decimal "2015.5" -- because the parser keeps only the leading number, so the
  # month/day/fraction is silently dropped.
  .warn_year <- function(key) {
    hit <- time_info[grepl(paste0("^", key, "="), time_info)]
    if (!length(hit)) return(invisible())
    raw <- sub(paste0("^", key, "="), "", hit[1])
    if (nzchar(raw) && !grepl("^[0-9]+$", raw))
      warning(sprintf(paste0(
        "'%s' should be an integer year; got '%s'. Only the leading number is ",
        "used (the time axis is annual), so any date or decimal part is dropped."),
        key, raw), call. = FALSE)
  }
  .warn_year("startpoint"); .warn_year("endpoint")

  time <- seq(0, endpoint - startpoint + 1,
              length.out = (partition * (endpoint - startpoint + 1) + 1))
  time <- as.numeric(time)

  list(time = time, startpoint = startpoint, endpoint = endpoint,
       partition = partition, cutoff = cutoff)
}

# ---- Internal: prepare the data / weight matrices --------------------------

# A well-formed ZERO-STREAM data bundle: the shape .prepare_data() returns, but
# with no observed streams. Used when a model is SIMULATED rather than fitted
# (no dataCombined). Every stream matrix/mask has n_years rows and 0 columns, so
# downstream code that iterates streams (solve_and_evaluate) simply does nothing,
# and the (unused) loss is skipped by the caller.
.empty_data <- function(tg) {
  n_years <- tg$endpoint - tg$startpoint + 1
  m0 <- matrix(numeric(0), nrow = n_years, ncol = 0)
  first_snapshot <- as.Date(sprintf("%d-12-31", tg$startpoint))
  data_points <- data.frame(
    date = seq.Date(from = first_snapshot, by = "year", length.out = n_years))
  list(
    data_combined      = data.frame(),
    data_points        = data_points,
    matrix_data_points = m0,
    obs_mask           = m0, cens_mask = m0, limit_mat = m0,
    inc_mask           = m0, lcens_mask = m0, llimit_mat = m0, linc_mask = m0,
    interval_mask      = m0, ilow_mat = m0, iupp_mat = m0,
    idev_lo_mat        = m0, idev_hi_mat = m0,
    asym_mask          = m0, asym_val_mat = m0, asym_dev_mat = m0, asym_dir_mat = m0,
    weight_matrix      = diag(numeric(0)), average_matrix = diag(numeric(0)),
    names_data_points  = character(0),
    likelihood_raw     = NULL,
    cumulative_cols    = integer(0)
  )
}

.prepare_data <- function(dataCombined, tg) {
  # No data -> simulation mode: return the zero-stream bundle. A NULL frame or a
  # frame with no rows both count as "no streams"; a non-empty frame without a
  # Formula column is still a genuine error (handled below).
  if (is.null(dataCombined) || nrow(as.data.frame(dataCombined)) == 0)
    return(.empty_data(tg))

  data_combined <- as.data.frame(dataCombined, stringsAsFactors = FALSE)
  # A missing/blank Label defaults to the stream's Formula, so plots/reports
  # always have a name to show.
  data_combined <- .default_label(data_combined)

  # Optional per-stream Likelihood column (Bayes only; MLE ignores it).
  likelihood_raw <- if ("Likelihood" %in% names(data_combined))
                      as.character(data_combined$Likelihood) else NULL

  # Weight/Average are genuine numeric inputs; fill their NAs with 0 as before.
  # (We do NOT blanket-fill the whole frame -- that would destroy the x / <L
  #  markers in the data cells, which we must parse first.)
  #
  # Type-coercion guard: if Excel stores any cell in a numeric column as text
  # (a stray text "0" is common), readxl returns the WHOLE column as character,
  # and sqrt() below would error. Coerce to numeric here. Blanks -> NA -> 0;
  # genuinely non-numeric text is a data error, so warn rather than silently 0.
  as_numeric_col <- function(x, col_name) {
    if (is.numeric(x)) return(x)
    chr  <- trimws(as.character(x))
    num  <- suppressWarnings(as.numeric(chr))
    bad  <- which(is.na(num) & nzchar(chr))   # non-blank that failed to parse
    if (length(bad))
      warning(sprintf(
        "Column '%s' has %d non-numeric text value(s) (e.g. '%s') coerced to 0.",
        col_name, length(bad), chr[bad[1]]), call. = FALSE)
    num
  }
  n_str    <- nrow(data_combined)
  meta_lab <- c("Label", "Formula", "Likelihood", "Weight", "Average")
  val_cols <- setdiff(names(data_combined), meta_lab)

  # Weight: optional per-stream importance multiplier. If the column is absent,
  # every stream gets weight 1; a blank cell within a present column also
  # defaults to 1 (with a warning, so a half-filled column is never silent).
  wt <- if ("Weight" %in% names(data_combined)) as_numeric_col(data_combined$Weight, "Weight")
        else rep(1, n_str)
  # Negative weights are invalid: the MLE loss uses diag(sqrt(wt)), so a
  # negative entry yields sqrt(<0) = NaN and silently breaks the fit.
  if (any(wt < 0, na.rm = TRUE)) {
    bad <- which(wt < 0)
    stop(sprintf("Weight must be non-negative; got %s in row%s %s of dataCombined.",
                 paste(wt[bad], collapse = ", "),
                 if (length(bad) > 1L) "s" else "",
                 paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  wt_blank <- which(is.na(wt))
  if (length(wt_blank)) {
    warning(sprintf("Weight is blank in row%s %s of dataCombined; defaulting to 1.",
                    if (length(wt_blank) > 1L) "s" else "",
                    paste(wt_blank, collapse = ", ")),
            call. = FALSE)
    wt[wt_blank] <- 1
  }

  # Average: optional per-stream residual scale for the MLE objective. Auto-
  # compute 1 / mean(observed values) per stream so streams of different
  # magnitude contribute comparably. This is used when the column is ABSENT,
  # and as the fallback for any blank cell within a present column (so a
  # half-filled column is never silently zeroed).
  auto_avg <- function(r) {
    v <- suppressWarnings(as.numeric(unlist(data_combined[r, val_cols], use.names = FALSE)))
    # Cumulative streams are differenced to annual increments before fitting,
    # so scale to the increments, not the raw cumulative totals.
    if (is_cumulative_stream(data_combined$Formula[r])) v <- diff(c(0, v))
    m <- mean(v[is.finite(v) & v != 0], na.rm = TRUE)   # plain-numeric cells only
    if (is.finite(m) && m > 0) 1 / m else 1
  }
  if ("Average" %in% names(data_combined)) {
    av <- as_numeric_col(data_combined$Average, "Average")
    av_blank <- which(is.na(av))
    if (length(av_blank)) {
      warning(sprintf("Average is blank in row%s %s of dataCombined; using the auto-computed 1 / mean scale there.",
                      if (length(av_blank) > 1L) "s" else "",
                      paste(av_blank, collapse = ", ")),
              call. = FALSE)
      av[av_blank] <- vapply(av_blank, auto_avg, numeric(1))
    }
  } else {
    av <- vapply(seq_len(n_str), auto_avg, numeric(1))
  }
  # NB: pass nrow explicitly. diag() of a length-1 vector treats the value as a
  # DIMENSION (e.g. diag(0.038) -> 0x0), which breaks single-stream fits whose
  # Weight/Average differ from 1; diag(x, nrow = length(x)) forces a diagonal.
  weight_matrix  <- diag(sqrt(wt), nrow = length(wt))
  average_matrix <- diag(sqrt(av), nrow = length(av))

  # Identify columns purely by KNOWN LABELS, independent of their order. Any
  # column whose name is one of the recognised meta labels is metadata; all
  # remaining columns are per-year data. This is robust to column reordering and
  # to extra meta columns.
  all_names   <- names(data_combined)
  meta_labels <- c("Label", "Formula", "Likelihood", "Weight", "Average")
  meta_cols   <- intersect(meta_labels, all_names)
  value_cols  <- setdiff(all_names, meta_cols)

  # The formula column drives stream identity. Require it explicitly.
  formula_col <- if ("Formula" %in% all_names) "Formula" else
    stop("dataCombined must contain a 'Formula' column.")

  formulas  <- as.character(data_combined[[formula_col]])  # one per stream (row)
  n_streams <- nrow(data_combined)
  n_years   <- length(value_cols)

  # Guard: the number of year/data columns must match the modelled horizon.
  expected_years <- tg$endpoint - tg$startpoint + 1
  if (n_years != expected_years) {
    stop(sprintf(
      "dataCombined has %d data columns but the time grid expects %d years (%d..%d).\n  Data columns detected: %s",
      n_years, expected_years, tg$startpoint, tg$endpoint,
      paste(value_cols, collapse = ", ")))
  }

  # Parse every data cell. Build, with STREAMS as rows and YEARS as cols:
  #   val_mat   : observed value (NA if not observed)
  #   obs_mat   : 1 if observed
  #   cens_mat  : 1 if LEFT-censored  (<L / <=L; truth <= L, an upper bound)
  #   lim_mat   : the upper-bound limit L (NA unless left-censored)
  #   inc_mat   : 1 if "<=" inclusive, 0 if "<" strict
  #   lcens_mat : 1 if RIGHT-censored (>L / >=L; truth >= L, a lower bound)
  #   llim_mat  : the lower-bound limit L (NA unless right-censored)
  #   linc_mat  : 1 if ">=" inclusive, 0 if ">" strict
  val_mat   <- matrix(NA_real_, n_streams, n_years)
  obs_mat   <- matrix(0L,       n_streams, n_years)
  cens_mat  <- matrix(0L,       n_streams, n_years)
  lim_mat   <- matrix(NA_real_, n_streams, n_years)
  inc_mat   <- matrix(NA,       n_streams, n_years)
  lcens_mat <- matrix(0L,       n_streams, n_years)
  llim_mat  <- matrix(NA_real_, n_streams, n_years)
  linc_mat  <- matrix(NA,       n_streams, n_years)
  # Interval [A,B]:  int_mat=1, ilow=A, iupp=B. Soft shoulders [A,B]~s carry a
  # positive scale in idlo/idhi (0 = hard edges, the CDF path).
  int_mat   <- matrix(0L,       n_streams, n_years)
  ilow_mat  <- matrix(NA_real_, n_streams, n_years)
  iupp_mat  <- matrix(NA_real_, n_streams, n_years)
  idlo_mat  <- matrix(0,        n_streams, n_years)
  idhi_mat  <- matrix(0,        n_streams, n_years)
  # Asymmetric A +/- dev:  asym_mat=1, aval=A, adev=dev, adir=+/-1.
  asym_mat  <- matrix(0L,       n_streams, n_years)
  aval_mat  <- matrix(NA_real_, n_streams, n_years)
  adev_mat  <- matrix(NA_real_, n_streams, n_years)
  adir_mat  <- matrix(NA_real_, n_streams, n_years)

  # Per-stream global asymmetric deviation (from the Likelihood column), used to
  # resolve A+/A- cells that don't carry their own deviation.
  stream_asym_dev <- if (is.null(likelihood_raw)) rep(NA_real_, n_streams) else
    vapply(likelihood_raw, function(x) parseLikelihood(x)$asym_dev, numeric(1))

  for (r in seq_len(n_streams)) {
    for (cc in seq_len(n_years)) {
      raw_cell <- data_combined[[ value_cols[cc] ]][r]
      pc <- parse_data_cell(raw_cell)
      if (pc$kind == "observed") {
        val_mat[r, cc] <- pc$value; obs_mat[r, cc] <- 1L
      } else if (pc$kind == "censored") {
        if (identical(pc$bound, "lower")) {
          lcens_mat[r, cc] <- 1L; llim_mat[r, cc] <- pc$limit
          linc_mat[r, cc]  <- if (isTRUE(pc$inclusive)) 1L else 0L
        } else {                                  # "upper" (default)
          cens_mat[r, cc] <- 1L; lim_mat[r, cc] <- pc$limit
          inc_mat[r, cc]  <- if (isTRUE(pc$inclusive)) 1L else 0L
        }
      } else if (pc$kind == "interval") {
        int_mat[r, cc] <- 1L; ilow_mat[r, cc] <- pc$limit; iupp_mat[r, cc] <- pc$upper
        idlo_mat[r, cc] <- if (is.na(pc$dev))  0 else pc$dev    # soft shoulders
        idhi_mat[r, cc] <- if (is.na(pc$dev2)) 0 else pc$dev2   # (0 = hard edges)
      } else if (pc$kind == "asym") {
        dev <- pc$dev
        if (is.na(dev)) dev <- stream_asym_dev[r]      # global deviation fallback
        if (is.na(dev))
          stop(sprintf(
            "Stream '%s' cell '%s' uses an A+/A- asymmetric value but its ",
            formulas[r], as.character(raw_cell)),
            "Likelihood column declares no 'asym=' deviation. Add e.g. ",
            "'; asym=<number>' to that stream's Likelihood, or use the ",
            "explicit 'A->B' form.", call. = FALSE)
        asym_mat[r, cc] <- 1L; aval_mat[r, cc] <- pc$value
        adev_mat[r, cc] <- dev; adir_mat[r, cc] <- pc$dir
      } # missing: leave all defaults
    }
  }

  # Transpose to YEARS x STREAMS (matching the model's observable matrices).
  matrix_data_points <- t(val_mat)         # observed values; NA for cens/missing
  obs_mask   <- t(obs_mat)
  cens_mask  <- t(cens_mat)
  limit_mat  <- t(lim_mat)
  inc_mask   <- t(inc_mat)
  lcens_mask <- t(lcens_mat)
  llimit_mat <- t(llim_mat)
  linc_mask  <- t(linc_mat)
  interval_mask <- t(int_mat);  ilow_mat_t   <- t(ilow_mat);  iupp_mat_t   <- t(iupp_mat)
  idev_lo_mat   <- t(idlo_mat); idev_hi_mat  <- t(idhi_mat)
  asym_mask     <- t(asym_mat); asym_val_mat <- t(aval_mat)
  asym_dev_mat  <- t(adev_mat); asym_dir_mat <- t(adir_mat)
  colnames(matrix_data_points) <- NULL; rownames(matrix_data_points) <- NULL

  # FITTING vs PLOTTING split for cumulative streams (unchanged logic), applied
  # to observed values only. Cumulative data differenced to annual increments.
  # cumulative_cols is stored in the returned list so callers can inspect which
  # streams were differenced and verify consistency with the loss function, which
  # independently identifies cumulative streams via is_cumulative_stream().
  cumulative_cols <- which(vapply(formulas, is_cumulative_stream, logical(1)))
  for (j in cumulative_cols) {
    matrix_data_points[, j] <- diff(c(0, matrix_data_points[, j]))
  }

  # PLOTTING frame: streams as columns, years as rows. Missing -> 0; censored ->
  # 0 (the usual-colour point); the limit L is carried separately so the plot can
  # add a second, differently-coloured marker at L. Observed -> the value.
  plot_val <- t(val_mat)
  plot_val[is.na(plot_val)] <- 0           # missing & censored both show 0
  data_points <- as.data.frame(plot_val, stringsAsFactors = FALSE)
  names(data_points) <- formulas

  first_snapshot <- as.Date(sprintf("%d-12-31", tg$startpoint))
  data_points$date <- seq.Date(from = first_snapshot, by = "year",
                               length.out = tg$endpoint - tg$startpoint + 1)

  names_data_points <- formulas

  list(
    data_combined      = data_combined,
    data_points        = data_points,        # plotting (missing/censored shown as 0)
    matrix_data_points = matrix_data_points, # fitting target (NA for cens/missing)
    obs_mask           = obs_mask,           # 1 = observed
    cens_mask          = cens_mask,          # 1 = left-censored (upper bound, <L)
    limit_mat          = limit_mat,          # upper-bound limit L (raw)
    inc_mask           = inc_mask,           # 1 = "<=" inclusive, 0 = "<" strict
    lcens_mask         = lcens_mask,         # 1 = right-censored (lower bound, >L)
    llimit_mat         = llimit_mat,         # lower-bound limit L (raw)
    linc_mask          = linc_mask,          # 1 = ">=" inclusive, 0 = ">" strict
    interval_mask      = interval_mask,      # 1 = interval-censored [A,B]
    ilow_mat           = ilow_mat_t,         # interval lower edge A
    iupp_mat           = iupp_mat_t,         # interval upper edge B
    idev_lo_mat        = idev_lo_mat,        # soft-interval lower shoulder scale (0 = hard)
    idev_hi_mat        = idev_hi_mat,        # soft-interval upper shoulder scale (0 = hard)
    asym_mask          = asym_mask,          # 1 = asymmetric A +/- dev
    asym_val_mat       = asym_val_mat,       # asymmetric anchor A
    asym_dev_mat       = asym_dev_mat,       # asymmetric deviation scale
    asym_dir_mat       = asym_dir_mat,       # +1 = soft upward, -1 = soft downward
    weight_matrix      = weight_matrix,
    average_matrix     = average_matrix,
    names_data_points  = names_data_points,
    likelihood_raw     = likelihood_raw,
    cumulative_cols    = cumulative_cols     # indices of streams that were differenced
  )
}

# ---- Internal: assemble bounds and normalised initial guesses --------------

.bounds <- function(sap) {
  upper_guesses <- c(sap$upper_states, sap$upper_params)
  lower_guesses <- c(sap$lower_states, sap$lower_params)

  # The loss looks up fitted quantities BY NAME, so non-empty bounds must be
  # named (states first, then params -- the estimation order). A FULLY-FIXED
  # sheet (every quantity *name=value, e.g. a filled counterfactual sheet) has
  # ZERO fitted quantities -> length-0 bounds, which is legitimate, not an error.
  nm <- names(lower_guesses)
  if (length(lower_guesses) > 0 && (is.null(nm) || any(!nzchar(nm))))
    stop(".bounds(): lower/upper guesses are missing names; the loss requires ",
         "named fitted quantities (states then params).")

  if (length(lower_guesses) == 0) {
    # Nothing to fit: empty named numeric bounds and empty normalised start.
    empty <- setNames(numeric(0), character(0))
    return(list(lower = empty, upper = empty, init_norm = empty,
                random_init = setNames(logical(0), character(0))))
  }

  initial_guesses_normalised <- mapply(normalise,
                                       sap$initial_guesses,
                                       lower_guesses,
                                       upper_guesses)
  names(upper_guesses)              <- nm
  names(initial_guesses_normalised) <- nm

  # '|random' mask aligned to the fitted (normalised) order. Defaults to all
  # FALSE if the sap predates this field (older serialised models).
  random_init <- c(sap$random_states, sap$random_params)
  random_init <- if (is.null(random_init)) setNames(rep(FALSE, length(nm)), nm)
                 else { ri <- as.logical(random_init[nm]); ri[is.na(ri)] <- FALSE
                        setNames(ri, nm) }

  list(lower = lower_guesses, upper = upper_guesses,
       init_norm = initial_guesses_normalised, random_init = random_init)
}

# ---- Internal: optimiser back-ends -----------------------------------------

.fit_optim <- function(loss_function, best_state, best_start,
                       lower_guesses, upper_guesses, ctrl, random_init = NULL) {
  # '|random' starts: redraw those entries uniformly from their box, which on the
  # normalised [0,1] scale the optimiser works on is simply runif(). Seeded by
  # ctrl$seed so the (otherwise deterministic) L-BFGS-B fit stays reproducible.
  rnd <- if (is.null(random_init)) character(0)
         else intersect(names(random_init)[as.logical(random_init)], names(best_start))
  if (length(rnd)) {
    if (!is.null(ctrl$seed)) set.seed(ctrl$seed)
    best_start[rnd] <- runif(length(rnd))
  }
  args <- list(
    par     = best_start,
    fn      = function(par) loss_function(par, best_state, TRUE),
    method  = ctrl$opt_method,
    control = list(maxit = ctrl$maxit, factr = ctrl$factr)
  )
  # Box constraints only apply to L-BFGS-B and Brent.
  if (ctrl$opt_method %in% c("L-BFGS-B", "Brent")) {
    args$lower <- numeric(length(lower_guesses))
    args$upper <- numeric(length(upper_guesses)) + 1
  }
  do.call(optim, args)
}

.fit_deoptim <- function(loss_function, best_state,
                         lower_guesses, upper_guesses, ctrl) {
  # DEoptim (differential evolution) is stochastic; seed it for reproducibility.
  if (!is.null(ctrl$seed)) set.seed(ctrl$seed)
  DEoptim(
    fn      = function(par) loss_function(par, best_state, TRUE),
    lower   = numeric(length(lower_guesses)),
    upper   = numeric(length(upper_guesses)) + 1,
    control = DEoptim.control(
      itermax = ctrl$itermax,
      NP      = ctrl$NP_mult * length(lower_guesses),
      CR      = ctrl$CR,
      F       = ctrl$F,
      # progress = FALSE silences BOTH the loss's running-best print (baked into
      # the loss via `verbose`) and DEoptim's own per-generation trace.
      trace   = if (isFALSE(ctrl$progress)) FALSE else ctrl$trace
    )
  )
}

.fit_hypercube <- function(loss_function, sap, bounds, hc_ctrl, best_state) {
  # Update the SHARED best_state environment (by reference) so the caller can
  # recover the natural-scale point, exactly as the optim/deoptim paths do.
  # hypercubeSampling() writes the denormalised best solution into best_state$par.
  hypercubeSampling(
    n = hc_ctrl$n,
    d = sap$hypercube_dim,
    loss_function = loss_function,
    not_sobol_states_normalised  = sap$not_sobol_states,
    not_sobol_params_normalised  = sap$not_sobol_params,
    lower_guesses_states_without = sap$lower_states_without,
    upper_guesses_states_without = sap$upper_states_without,
    lower_guesses_params_without = sap$lower_params_without,
    upper_guesses_params_without = sap$upper_params_without,
    lower_guesses = bounds$lower,
    upper_guesses = bounds$upper,
    best_state    = best_state
  )
}

# ---- Internal: Bayesian back-end (Julia / Turing) --------------------------
.fit_bayes <- function(model, data, tg, bounds, point_init, bc, modelParams) {
  # Guard reserved names: the Bayes observation model uses `sigma` (Gaussian /
  # log-normal noise scale) and `phi`/`phi<i>` (negbin dispersion) as its own
  # parameters. A model parameter or state sharing one of those names silently
  # collides -- it is sampled as the noise term and then dropped as a
  # hyperparameter, so it vanishes from the posterior. Fail early with a clear
  # message instead. (MLE is unaffected; only method = "bayes" reserves these.)
  .check_reserved_bayes_names(c(names(model$sap$params_fitted),
                                names(model$sap$params_fixed),
                                names(model$sap$states_fitted),
                                names(model$sap$states_fixed)))

  # Build the parallel prior spec from the parameter sheet (independent of the
  # MLE statesAndParams() output).
  prior_spec <- buildPriorSpec(modelParams)

  # Per-stream likelihood families from the (optional) Likelihood column.
  # Absent column -> default gaussian for every stream.
  formulas <- data$names_data_points
  if (is.null(data$likelihood_raw)) {
    like_specs <- lapply(seq_along(formulas), function(i) parseLikelihood(NA))
  } else {
    like_specs <- lapply(data$likelihood_raw, parseLikelihood)
  }

  n_years   <- tg$endpoint - tg$startpoint + 1
  partition <- tg$partition

  # --- Integer failsafe for count likelihoods (poisson/negbin) ---------------
  # Poisson/NegBin have integer support; non-integer y gives -Inf log-density
  # and the sampler can't start. Round those streams' data to nearest integer,
  # warning if rounding actually changed a value (so it's a safety net, not a
  # silent edit).
  ym <- data$matrix_data_points
  fam_of <- function(i) {
    if (is.null(data$likelihood_raw)) return("gaussian")
    v <- data$likelihood_raw[i]
    if (is.na(v) || v == "") "gaussian" else tolower(trimws(v))
  }
  for (i in seq_len(ncol(ym))) {
    if (fam_of(i) %in% c("poisson", "negbin")) {
      col <- ym[, i]
      rounded <- round(col)
      if (any(abs(col - rounded) > 1e-8, na.rm = TRUE)) {
        warning(sprintf("Stream %d uses a count likelihood (%s) but has non-integer ",
                        i, fam_of(i)),
                "data; rounding to nearest integer.", call. = FALSE)
      }
      ym[, i] <- rounded
    }
  }
  data$matrix_data_points <- ym

  # Generate + register the Turing model.
  model_code <- buildJuliaBayesModel(
    prior_spec      = prior_spec,
    like_specs      = like_specs,
    formulas        = formulas,
    n_years         = n_years,
    partition       = partition,
    number_of_comps = model$structure$number_of_comps,
    comp_names      = model$structure$comp_names,
    sigma_prior     = bc$sigma_prior,
    phi_prior       = bc$phi_prior
  )
  registerJuliaBayesModel(model_code)

  # --- Censoring: family-aware effective limits ------------------------------
  # LEFT (upper bound, <L):  contribution logcdf(dist, Lstar).
  #   "<=L" / continuous -> Lstar = L;  "<L" strict on DISCRETE -> Lstar = L-1.
  # RIGHT (lower bound, >L): contribution logccdf(dist, Lstar), where Julia's
  #   ccdf(d, x) = P(Y > x). So ">L" strict -> Lstar = L; ">=L" inclusive on
  #   DISCRETE -> Lstar = L-1 (to include L: P(Y>=L)=P(Y>L-1)). Continuous -> L.
  discrete_fam <- vapply(seq_along(formulas),
                         function(i) fam_of(i) %in% c("poisson", "negbin"),
                         logical(1))
  limit_eff <- data$limit_mat                 # upper-bound limits
  cm  <- data$cens_mask;  inc  <- data$inc_mask
  llimit_eff <- data$llimit_mat               # lower-bound limits
  lcm <- data$lcens_mask; linc <- data$linc_mask
  ilow_eff <- data$ilow_mat                   # interval lower edges
  for (j in seq_len(ncol(limit_eff))) {
    if (discrete_fam[j]) {
      # left, strict ("<L") -> L-1
      strict_j <- which(cm[, j] == 1 & (is.na(inc[, j]) | inc[, j] == 0))
      limit_eff[strict_j, j] <- limit_eff[strict_j, j] - 1
      # right, inclusive (">=L") -> L-1
      incl_j <- which(lcm[, j] == 1 & linc[, j] == 1)
      llimit_eff[incl_j, j] <- llimit_eff[incl_j, j] - 1
      # interval, lower edge on discrete counts -> A-1, so cdf(B)-cdf(A-1)=P(A<=Y<=B).
      # Only HARD intervals use the CDF; soft intervals ([A,B]~s) hinge on the
      # continuous mean and must keep the original edge A.
      int_j <- which(data$interval_mask[, j] == 1 & data$idev_lo_mat[, j] == 0)
      ilow_eff[int_j, j] <- ilow_eff[int_j, j] - 1
    }
  }

  # Sanitize matrices for Julia: y NA -> 0 (gated out by obs_mask), limit NA -> 0.
  ym_jl <- data$matrix_data_points
  ym_jl[is.na(ym_jl)] <- 0
  lim_jl <- limit_eff;  lim_jl[is.na(lim_jl)] <- 0
  llim_jl <- llimit_eff; llim_jl[is.na(llim_jl)] <- 0
  ilow_jl <- ilow_eff;          ilow_jl[is.na(ilow_jl)] <- 0
  iupp_jl <- data$iupp_mat;     iupp_jl[is.na(iupp_jl)] <- 0
  idlo_jl <- data$idev_lo_mat;  idlo_jl[is.na(idlo_jl)] <- 0
  idhi_jl <- data$idev_hi_mat;  idhi_jl[is.na(idhi_jl) | idhi_jl == 0] <- 1  # avoid /0 in soft branch
  aval_jl <- data$asym_val_mat; aval_jl[is.na(aval_jl)] <- 0
  adev_jl <- data$asym_dev_mat; adev_jl[is.na(adev_jl)] <- 1   # avoid /0 when unused
  adir_jl <- data$asym_dir_mat; adir_jl[is.na(adir_jl)] <- 0

  # Per-stream scale columns for proportional-error gaussian: each stream's mean
  # OBSERVED level (ignore censored/missing zeros) so residuals are comparable.
  obsm <- data$obs_mask
  scale_cols <- vapply(seq_len(ncol(ym_jl)), function(j) {
    vals <- ym_jl[obsm[, j] == 1, j]
    m <- mean(abs(vals))
    if (!is.finite(m) || m == 0) {
      if (sum(obsm[, j]) == 0)
        warning(sprintf(
          "Stream %d ('%s') has no observed data; its proportional-error scale defaults to 1.",
          j, formulas[j]), call. = FALSE)
      1
    } else m
  }, numeric(1))
  scale_mat <- matrix(rep(scale_cols, each = nrow(ym_jl)), nrow = nrow(ym_jl))

  # Sample.
  res <- sampleWithJulia(
    y_matrix   = ym_jl,
    tspan      = c(min(tg$time), max(tg$time)),
    n_years    = n_years,
    partition  = partition,
    scale_mat  = scale_mat,
    obs_mask   = obsm,
    cens_mask  = cm,
    limit_mat  = lim_jl,
    lcens_mask = data$lcens_mask,
    llimit_mat = llim_jl,
    interval_mask = data$interval_mask,
    ilow_mat      = ilow_jl,
    iupp_mat      = iupp_jl,
    idev_lo_mat   = idlo_jl,
    idev_hi_mat   = idhi_jl,
    asym_mask     = data$asym_mask,
    asym_val_mat  = aval_jl,
    asym_dev_mat  = adev_jl,
    asym_dir_mat  = adir_jl,
    bc         = bc
  )

  # Attach the generated code for inspection/debugging.
  res$model_code <- model_code
  res$prior_spec <- prior_spec
  res$like_specs <- like_specs
  # For posterior-predictive bands: per-stream proportional-noise scale and the
  # stream order, so the plot layer can simulate from each family correctly.
  res$scale_cols <- scale_cols
  res$stream_names <- formulas
  res
}

# ---- Internal: recover natural-scale states/params from a fit --------------

.recover_solution <- function(best_par, sap, tg) {
  actual_fit <- best_par
  # A fitted quantity is a STATE iff it is one of the fitted compartments; the
  # fitted-state names (X1..Xn in the identity case, or S/I/R/... under a named
  # registry) come straight from sap. Everything else is a fitted parameter.
  actual_fit_states <- actual_fit[names(actual_fit) %in% names(sap$states_fitted)]
  params_fitted <- actual_fit[setdiff(names(actual_fit), names(actual_fit_states))]

  # Order-contract guard: catch any future drift between the four codepaths
  # that must agree on p = [fitted_params..., fitted_states...]. A fully-fixed
  # sheet has ZERO fitted quantities, where names() may be NULL or character(0);
  # normalise both to character(0) so the legitimate no-free-parameter case
  # (e.g. a filled counterfactual sheet) passes rather than tripping the guard.
  .nm <- function(x) { n <- names(x); if (is.null(n)) character(0) else n }
  expected_params <- .nm(sap$params_fitted)
  expected_states <- .nm(sap$states_fitted)
  if (!identical(.nm(params_fitted), expected_params))
    stop(sprintf(
      ".recover_solution: fitted-parameter name/order mismatch.\n  expected: [%s]\n  got:      [%s]",
      paste(expected_params, collapse = ", "),
      paste(.nm(params_fitted), collapse = ", ")
    ))
  if (!identical(sort(.nm(actual_fit_states)), sort(expected_states)))
    stop(sprintf(
      ".recover_solution: fitted-state name mismatch.\n  expected: [%s]\n  got:      [%s]",
      paste(sort(expected_states), collapse = ", "),
      paste(sort(.nm(actual_fit_states)), collapse = ", ")
    ))

  params_fixed_final <- sapply(sap$params_fixed,
                               function(expr) eval(parse(text = expr)))
  parms <- c(params_fitted, params_fixed_final)
  names_parameters <- names(parms)
  parms <- as.numeric(parms)
  names(parms) <- names_parameters

  names_states_fixed   <- names(sap$states_fixed)
  states_fixed_final_2 <- sap$states_fixed
  parms_states <- parms
  names(parms_states) <- paste0(names(parms_states), "_0")

  states_fixed_final_2[sap$non_numeric_fixed] <- mapply(
    function(f) do.call(f, list(parms_states)),
    sap$list_states_functions
  )
  states_fixed_final_2 <- as.numeric(states_fixed_final_2)
  names(states_fixed_final_2) <- names_states_fixed

  actual_fit_states_all <- c(actual_fit_states, states_fixed_final_2)
  # Order by canonical compartment position (States-column order) so the solver
  # receives X in index order regardless of how the states were named/split. When
  # no registry is present (legacy sap), fall back to the X<i> numeric suffix.
  comp_order <- if (!is.null(sap$comp_names))
                  match(names(actual_fit_states_all), sap$comp_names)
                else as.numeric(sub("^X", "", names(actual_fit_states_all)))
  initial_state <- actual_fit_states_all[order(comp_order)]

  list(initial_state = initial_state, parms = parms)
}

# ---- Build-only entry point ------------------------------------------------
# Runs everything fitCompartmentalModel() does EXCEPT the fit dispatch: builds
# the model, time grid, prepared data, bounds and loss closure, and returns them
# as a list (the "machinery"). No optimisation or sampling happens. This lets
# tools that need only the model + bounds (e.g. Sobol sensitivity) run WITHOUT a
# prior fit. fitCompartmentalModel() keeps its own copy of these steps so its
# behaviour is unchanged.
#
# Returns a list of class "compartmentalModel":
#   model, sap, time_grid, data, bounds, loss, best_state, best_start, solver
#' Build a compartmental model (no fitting)
#'
#' Runs everything [fitCompartmentalModel()] does except the fit dispatch:
#' builds the model, time grid, prepared data, bounds and loss closure. Useful
#' for tools that need the model and bounds without a prior fit (e.g. Sobol
#' sensitivity).
#'
#' @param modelParams Model parameter sheet (data frame).
#' @param dataCombined Combined observation data (data frame). Optional: pass
#'   `NULL` (the default) to build in simulation mode with no observed streams,
#'   in which case the loss closure is `NULL` (nothing to fit against). See
#'   [simulate_model()].
#' @param solver Solver settings, see [solver_control()].
#' @param init Optional initial point (named numeric) for the loss closure.
#' @param checkpoint_file Path for the best-solution checkpoint `.rds`.
#' @return A list of class `"compartmentalModel"` (model, sap, time_grid, data,
#'   bounds, loss, best_state, best_start, solver). `loss` is `NULL` when
#'   `dataCombined` is absent.
#' @examples
#' mini <- system.file("extdata", "minimal", package = "compfit")
#' sc   <- load_scenario(mini, combined_file = "dataCombined.csv",
#'                       dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
#' # Build the machinery (no fitting). The R backend needs no Julia:
#' m <- build_compartmental_model(sc$modelParams, sc$dataCombined,
#'                                solver = solver_control(backend = "r"))
#' m$bounds
#' @export
build_compartmental_model <- function(modelParams, dataCombined = NULL,
                                      solver = solver_control(),
                                      init   = NULL,
                                      checkpoint_file = tempfile(fileext = ".rds")) {
  model <- .build_model(modelParams, backend = solver$backend)
  sap   <- model$sap
  tg    <- .time_grid(modelParams)
  dat   <- .prepare_data(dataCombined, tg)
  bounds <- .bounds(sap)

  # No observed streams (dataCombined absent) -> simulation only: the loss has
  # nothing to compare against, so skip building it. solve_and_evaluate() never
  # touches the loss, and every tool that needs it (Sobol, fits) supplies data.
  loss_function <- if (length(dat$names_data_points) == 0) NULL else lossFunction(
    names_data_points         = dat$names_data_points,
    parms_expression          = model$expressions$parms,
    initial_states_expression = model$expressions$initial_states,
    states_params_expression  = model$expressions$states_params,
    penalty_expression        = model$expressions$penalty,
    time                      = tg$time,
    startpoint                = tg$startpoint,
    endpoint                  = tg$endpoint,
    partition                 = tg$partition,
    matrix_data_points        = dat$matrix_data_points,
    weight_matrix             = dat$weight_matrix,
    average_matrix            = dat$average_matrix,
    obs_mask                  = dat$obs_mask,
    cens_mask                 = dat$cens_mask,
    limit_mat                 = dat$limit_mat,
    lcens_mask                = dat$lcens_mask,
    llimit_mat                = dat$llimit_mat,
    interval_mask             = dat$interval_mask,
    ilow_mat                  = dat$ilow_mat,
    iupp_mat                  = dat$iupp_mat,
    idev_lo_mat               = dat$idev_lo_mat,
    idev_hi_mat               = dat$idev_hi_mat,
    asym_mask                 = dat$asym_mask,
    asym_val_mat              = dat$asym_val_mat,
    asym_dev_mat              = dat$asym_dev_mat,
    asym_dir_mat              = dat$asym_dir_mat,
    lower_guesses             = bounds$lower,
    upper_guesses             = bounds$upper,
    comp_names                = sap$comp_names,
    solver                    = solver$solver,
    abstol                    = solver$abstol,
    reltol                    = solver$reltol,
    backend                   = .ode_backend(solver$backend),
    checkpoint_file           = checkpoint_file
  )

  best_state <- new.env(); best_state$error <- Inf; best_state$par <- NULL
  best_start <- if (!is.null(init)) init else bounds$init_norm

  structure(list(
    model = model, sap = sap, time_grid = tg, data = dat,
    bounds = bounds, loss = loss_function,
    best_state = best_state, best_start = best_start, solver = solver
  ), class = "compartmentalModel")
}

# ---- Master function -------------------------------------------------------

#' Fit a compartmental model
#'
#' Build a compartmental ODE model from a spreadsheet specification and fit it
#' by maximum likelihood or Bayesian sampling. By default ODE solves (and Julia
#' NUTS) run in Julia; call [setup_julia()] first or rely on the lazy
#' initialiser. A pure-R backend and a Stan backend need no Julia -- see
#' [solver_control()].
#'
#' @param modelParams Model-parameter sheet (data frame): compartments, priors,
#'   rate coefficients and the time grid. See `vignette("compfit")` for the
#'   grammar.
#' @param dataCombined Observed data, one row per stream (`Formula` + one column
#'   per year, optional `Label`/`Likelihood`/`Weight`/`Average`). Required for
#'   fitting; to run a fully-fixed model forward with no data use
#'   [simulate_model()] instead.
#' @param method Fitting method (one of):
#'   \describe{
#'     \item{`"lbfgsb"`}{deterministic local MLE (L-BFGS-B). Fast; the default.}
#'     \item{`"deoptim"`}{stochastic global MLE (differential evolution). More
#'       robust on rugged/multi-modal objectives, slower. Seed via
#'       [optim_control()].}
#'     \item{`"hypercube"`}{deterministic Sobol pre-search; a coarse global scan
#'       (see [hypercube_control()]).}
#'     \item{`"bayes"`}{full Bayesian posterior sampling (see [bayes_control()]
#'       and the `backend` in [solver_control()]).}
#'   }
#' @param solver Solver/backend settings, see [solver_control()] -- chooses
#'   `"julia"` / `"r"` / `"stan"` and the ODE tolerances.
#' @param control Optimiser settings for the MLE methods, see [optim_control()].
#' @param hypercube Pre-search settings for `method = "hypercube"`, see
#'   [hypercube_control()].
#' @param bayes Sampler settings, see [bayes_control()]; used only for
#'   `method = "bayes"` (defaults are supplied if `NULL`).
#' @param init Optional starting point as a named numeric vector on the
#'   normalised \[0,1\] scale (advanced; normally left `NULL`).
#' @param checkpoint_file Path to an `.rds` where the running best solution is
#'   written during an MLE search, so a long fit can be recovered if interrupted.
#' @param debug_env Optional environment into which intermediate build objects
#'   are stashed for debugging (no effect on the result when `NULL`).
#' @return An object of class `"compartmentalFit"`: `$point` (natural-scale
#'   estimates) for MLE, `$samples` for Bayes, plus the model, data, bounds and
#'   solver used. Inspect with [summary()], [get_point()], [posterior_report()],
#'   [plot_fit()].
#' @examples
#' mini <- system.file("extdata", "minimal", package = "compfit")
#' sc   <- load_scenario(mini, combined_file = "dataCombined.csv",
#'                       dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
#' # Julia-free maximum-likelihood fit (pure-R deSolve backend) -- runs at check:
#' fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
#'                              solver = solver_control(backend = "r"),
#'                              checkpoint_file = tempfile(fileext = ".rds"))
#' get_point(fit)
#' \dontrun{
#' # The Julia backend: faster, and required for NUTS Bayesian sampling.
#' setup_julia()
#' fit_b <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "bayes",
#'                                bayes = bayes_control(chains = 2, iter = 1000))
#' }
#' @export
fitCompartmentalModel <- function(modelParams,
                                  dataCombined,
                                  method   = c("lbfgsb", "deoptim", "hypercube", "bayes"),
                                  solver   = solver_control(),
                                  control  = optim_control(),
                                  hypercube = hypercube_control(),
                                  bayes    = NULL,
                                  init     = NULL,
                                  checkpoint_file = "checkpoint_best_solution.rds",
                                  debug_env = NULL) {

  method <- match.arg(method)

  # Stash an intermediate into debug_env (if supplied) the moment it exists,
  # so build-phase failures still leave earlier artifacts inspectable.
  # No-op when debug_env is NULL, so production runs are unaffected.
  stash <- function(name, value) {
    if (!is.null(debug_env)) assign(name, value, envir = debug_env)
    value
  }

  # --- 1. Build model (structure, expressions, ODE fns, Julia registration) ---
  model <- stash("model", .build_model(modelParams, backend = solver$backend))
  sap   <- model$sap
  stash("sap", sap)
  stash("compartmental_function", model$compartmental_function)

  # --- 2. Time grid ---
  tg <- stash("tg", .time_grid(modelParams))
  # Make `partition` etc. visible for evaluate_formula() at plot time.
  startpoint <- tg$startpoint; endpoint <- tg$endpoint
  partition  <- tg$partition;  cutoff   <- tg$cutoff
  time       <- tg$time

  # --- 3. Data / weight matrices ---
  dat <- stash("dat", .prepare_data(dataCombined, tg))

  # --- 4. Bounds (BEFORE the loss, so the closure captures the right vectors) ---
  bounds <- stash("bounds", .bounds(sap))
  lower_guesses <- bounds$lower
  upper_guesses <- bounds$upper

  # --- 5. Loss function ---
  loss_function <- stash("loss_function", lossFunction(
    names_data_points         = dat$names_data_points,
    parms_expression          = model$expressions$parms,
    initial_states_expression = model$expressions$initial_states,
    states_params_expression  = model$expressions$states_params,
    penalty_expression        = model$expressions$penalty,
    time                      = time,
    startpoint                = startpoint,
    endpoint                  = endpoint,
    partition                 = partition,
    matrix_data_points        = dat$matrix_data_points,
    weight_matrix             = dat$weight_matrix,
    average_matrix            = dat$average_matrix,
    obs_mask                  = dat$obs_mask,
    cens_mask                 = dat$cens_mask,
    limit_mat                 = dat$limit_mat,
    lcens_mask                = dat$lcens_mask,
    llimit_mat                = dat$llimit_mat,
    interval_mask             = dat$interval_mask,
    ilow_mat                  = dat$ilow_mat,
    iupp_mat                  = dat$iupp_mat,
    idev_lo_mat               = dat$idev_lo_mat,
    idev_hi_mat               = dat$idev_hi_mat,
    asym_mask                 = dat$asym_mask,
    asym_val_mat              = dat$asym_val_mat,
    asym_dev_mat              = dat$asym_dev_mat,
    asym_dir_mat              = dat$asym_dir_mat,
    lower_guesses             = lower_guesses,
    upper_guesses             = upper_guesses,
    comp_names                = sap$comp_names,
    # For a Bayesian fit the loss is only exercised by the optional MLE seeding
    # step, so its printing follows bayes_control()'s `progress` -- one knob then
    # silences the whole Bayesian fit (init + sampler). MLE fits use `control`.
    verbose                   = if (identical(method, "bayes"))
                                  isTRUE((if (is.null(bayes)) bayes_control() else bayes)$progress)
                                else isTRUE(control$progress),
    solver                    = solver$solver,
    abstol                    = solver$abstol,
    reltol                    = solver$reltol,
    backend                   = .ode_backend(solver$backend),
    checkpoint_file           = checkpoint_file
  ))

  # --- 6. Best-state tracker (environment so the loss can mutate it) ---
  best_state <- new.env()
  best_state$error <- Inf
  best_state$par   <- NULL
  stash("best_state", best_state)

  # Warm start: supplied normalised par, else normalised midpoint guesses.
  best_start <- if (!is.null(init)) init else bounds$init_norm
  stash("best_start", best_start)
  # '|random' parameters redraw their start from the box at fit time (seeded),
  # unless the caller supplied an explicit `init` (which takes precedence).
  random_start <- if (!is.null(init)) NULL else bounds$random_init

  # --- 7. Dispatch on method ---
  # Guarded so a fit-phase failure still returns a usable partial object:
  # model / loss / data / bounds are all populated and inspectable, and
  # best_state holds whatever the optimiser reached before failing.
  fit_failed  <- FALSE
  fit_message <- NULL
  fit_raw <- tryCatch(
    switch(
      method,
      "lbfgsb"    = .fit_optim(loss_function, best_state, best_start,
                               lower_guesses, upper_guesses, control, random_start),
      "deoptim"   = .fit_deoptim(loss_function, best_state,
                                 lower_guesses, upper_guesses, control),
      "hypercube" = .fit_hypercube(loss_function, sap, bounds, hypercube, best_state),
      "bayes"     = {
        if (is.null(bayes)) bayes <- bayes_control()
        if (identical(solver$backend, "stan")) {
          # Stan/rstan NUTS. Stan does its own warmup adaptation, so no MAP
          # pre-fit; the loss (deSolve) is only used later for plots/predictive.
          .fit_bayes_stan(model, dat, tg, bounds, bayes, modelParams, solver)
        } else if (identical(solver$backend, "r")) {
          # Pure-R, gradient-free path (BayesianTools DEzs); no Julia.
          # Seed the population from a quick MLE fit when requested (default).
          point_init <- NULL
          if (isTRUE(bayes$init_from_optim)) {
            .fit_optim(loss_function, best_state, best_start,
                       lower_guesses, upper_guesses, control, random_start)
            point_init <- best_state$par
          }
          .fit_bayes_r(model, dat, tg, bounds, bayes, modelParams,
                       point_init = point_init, solver = solver)
        } else {
          # Optionally MAP-initialise the sampler with a quick optim run.
          point_init <- NULL
          if (isTRUE(bayes$init_from_optim)) {
            opt <- .fit_optim(loss_function, best_state, best_start,
                              lower_guesses, upper_guesses, control, random_start)
            point_init <- best_state$par
          }
          .fit_bayes(model, dat, tg, bounds, point_init, bayes, modelParams)
        }
      }
    ),
    error = function(e) {
      fit_failed  <<- TRUE
      fit_message <<- conditionMessage(e)
      warning("Fit failed: ", fit_message,
              "\nReturning a partial compartmentalFit (model/loss/data/bounds ",
              "populated; best_state holds progress before failure).",
              call. = FALSE)
      structure(list(message = fit_message), class = "fitError")
    }
  )
  stash("fit_raw", fit_raw)

  # --- 8. Recover natural-scale solution from best_state (optim/deoptim/hypercube) ---
  point <- NULL
  if (!fit_failed &&
      method %in% c("lbfgsb", "deoptim", "hypercube") &&
      !is.null(best_state$par)) {
    point <- .recover_solution(best_state$par, sap, tg)
  }

  # --- 9. Assemble result (always returned, even on fit failure) ---
  structure(
    list(
      method      = method,
      success     = !fit_failed,    # FALSE if the fit phase errored
      error_msg   = fit_message,    # NULL unless the fit phase errored
      point       = point,          # list(initial_state, parms) on natural scale
      best_state  = best_state,     # error + best normalised par
      fit_raw     = fit_raw,        # raw optim/DEoptim/hypercube/bayes return, or fitError
      samples     = if (method == "bayes" && !fit_failed) fit_raw else NULL,
      loss        = loss_function,
      model       = model,
      data        = dat,
      time_grid   = tg,
      bounds      = bounds,
      solver      = solver,
      meta        = list(checkpoint_file = checkpoint_file)
    ),
    class = "compartmentalFit"
  )
}
