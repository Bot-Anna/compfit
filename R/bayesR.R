# ============================================================
# bayesR.R
# Gradient-free Bayesian backend (pure R, NO Julia) for solver_control(backend="r").
#
# WHY GRADIENT-FREE: the Julia/Turing path uses NUTS, which needs gradients of
# the log-posterior -- i.e. automatic differentiation THROUGH the ODE solve.
# R has no equivalent, so this backend uses a derivative-free population MCMC
# (BayesianTools' DEzs). It evaluates the same model (deSolve via solveWithR)
# and the same likelihood families/censoring as the Turing model, but mixes far
# worse and needs many more iterations.
#
#   ************************************************************************
#   * VERY INEFFICIENT. Intended only for small models / quick checks, or *
#   * where installing Julia is impossible. For real inference use         *
#   * solver_control(backend = "julia").                                   *
#   ************************************************************************
#
# It returns a `samples` object structurally identical to the Julia path
# (draws with `name_n`/`name`/`sigma`/`phi_i` columns, a summary table, and the
# prior_spec), so every posterior_*() / plot_prior_posterior() / identifiability
# helper works unchanged.
# ============================================================

# Informative-prior log-density (uniform priors are handled by the sampler's
# own uniform prior, so they are not added here).
.rprior_logd <- function(spec, x) {
  a <- suppressWarnings(as.numeric(spec$args))
  d <- switch(spec$dist,
              Normal    = stats::dnorm(x, a[1], a[2], log = TRUE),
              LogNormal = stats::dlnorm(x, a[1], a[2], log = TRUE),
              Beta      = stats::dbeta(x, a[1], a[2], log = TRUE),
              Gamma     = stats::dgamma(x, shape = a[1], scale = a[2], log = TRUE),
              # Location-scale Student-t(nu, mu, sigma): standardise then correct
              # for the scale (dt is the standard t; -log(sigma) is the Jacobian).
              StudentT  = stats::dt((x - a[2]) / a[3], df = a[1], log = TRUE) - log(a[3]),
              0)
  if (is.finite(d)) d else -Inf
}

.fit_bayes_r <- function(model, data, tg, bounds, bc, modelParams, point_init = NULL,
                         solver = solver_control(backend = "r")) {
  if (!requireNamespace("BayesianTools", quietly = TRUE))
    stop("backend = \"r\" Bayesian sampling needs the 'BayesianTools' package ",
         "(install.packages(\"BayesianTools\")).")
  message("compfit: gradient-free R Bayesian sampler (BayesianTools DEzs). This ",
          "is MUCH slower and mixes far worse than the Julia/NUTS backend -- use ",
          "it only for small models or quick checks. For production inference use ",
          "solver_control(backend = \"julia\").")

  sap        <- model$sap
  prior_spec <- buildPriorSpec(modelParams)
  formulas   <- data$names_data_points
  like_specs <- if (is.null(data$likelihood_raw))
    lapply(seq_along(formulas), function(i) parseLikelihood(NA))
  else lapply(data$likelihood_raw, parseLikelihood)

  n_years   <- tg$endpoint - tg$startpoint + 1
  partition <- tg$partition

  # Count families need integer data (same failsafe as the Julia path).
  ym   <- data$matrix_data_points
  disc <- vapply(like_specs, function(s) s$family %in% c("poisson", "negbin"), logical(1))
  for (i in seq_len(ncol(ym))) if (disc[i]) ym[, i] <- round(ym[, i])
  data$matrix_data_points <- ym

  # Effective censor limits: for COUNT families the discrete CDF needs the
  # family-aware shift so strict vs inclusive is exact (left-strict "<L" -> L-1;
  # right-inclusive ">=L" -> L-1). Continuous families use the raw limit.
  limit_eff  <- data$limit_mat
  llimit_eff <- data$llimit_mat
  cm  <- data$cens_mask;  inc  <- data$inc_mask
  lcm <- data$lcens_mask; linc <- data$linc_mask
  for (j in which(disc)) {
    strict_j <- which(cm[, j]  == 1 & (is.na(inc[, j])  | inc[, j]  == 0))   # "<L"  -> L-1
    limit_eff[strict_j, j]  <- limit_eff[strict_j, j]  - 1
    incl_j   <- which(lcm[, j] == 1 & !is.na(linc[, j]) & linc[, j] == 1)    # ">=L" -> L-1
    llimit_eff[incl_j, j]   <- llimit_eff[incl_j, j]   - 1
  }

  # Per-stream proportional-error scale = mean observed level.
  obsm <- data$obs_mask
  scale_cols <- vapply(seq_len(ncol(ym)), function(j) {
    vals <- ym[obsm[, j] == 1, j]; m <- mean(abs(vals), na.rm = TRUE)
    if (!is.finite(m) || m == 0) 1 else m
  }, numeric(1))

  # ---- Sampled-parameter layout (matches the Julia draws naming) ----------
  pf <- prior_spec$order$params_fitted
  sf <- prior_spec$order$states_fitted
  model_specs <- c(prior_spec$params[pf], prior_spec$states[sf])
  model_names <- c(pf, sf)
  is_unif     <- vapply(model_specs, function(s) identical(s$dist, "Uniform"), logical(1))
  sampled_nm  <- ifelse(is_unif, paste0(model_names, "_n"), model_names)
  n_model     <- length(model_names)

  lo <- hi <- numeric(n_model)
  for (k in seq_len(n_model)) {
    s <- model_specs[[k]]; a <- suppressWarnings(as.numeric(s$args))
    if (is_unif[k]) { lo[k] <- 0; hi[k] <- 1 }
    else {
      lo[k] <- if (is.finite(s$lower)) s$lower else a[1] - 5 * a[2]
      hi[k] <- if (is.finite(s$upper)) s$upper else a[1] + 5 * a[2]
    }
  }

  # Noise hyperparameters: sigma (gaussian/lognormal) and phi_i (dispersion).
  needs_sigma <- any(vapply(like_specs,
                            function(s) s$family %in% c("gaussian", "lognormal"), logical(1)))
  hyper_nm <- character(0); hyper_lo <- hyper_hi <- numeric(0); phi_idx <- integer(0)
  if (needs_sigma) { hyper_nm <- "sigma"; hyper_lo <- 1e-4; hyper_hi <- 5 }
  for (i in seq_along(like_specs)) if (isTRUE(like_specs[[i]]$dispersion)) {
    hyper_nm <- c(hyper_nm, sprintf("phi%d", i))
    hyper_lo <- c(hyper_lo, 1e-3); hyper_hi <- c(hyper_hi, 1e4); phi_idx <- c(phi_idx, i)
  }

  par_names <- c(sampled_nm, hyper_nm)
  low <- c(lo, hyper_lo); upp <- c(hi, hyper_hi)

  # ---- theta -> observable matrix (years x streams) -----------------------
  predict_mu <- function(theta) {
    nat <- numeric(n_model)
    for (k in seq_len(n_model)) {
      if (is_unif[k]) { a <- as.numeric(model_specs[[k]]$args); nat[k] <- a[1] + theta[k] * (a[2] - a[1]) }
      else nat[k] <- theta[k]
    }
    names(nat) <- model_names
    sol <- tryCatch(.recover_solution(nat, sap, tg), error = function(e) NULL)
    if (is.null(sol)) return(NULL)
    X <- unlist(sol$initial_state); p <- unlist(sol$parms)
    t <- c(as.numeric(min(tg$time)), as.numeric(max(tg$time)))
    r <- tryCatch(solveWithR(X, t, p, tg$time, abstol = solver$abstol,
                             reltol = solver$reltol, method = .desolve_method(solver$solver)),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    sir <- as.data.frame(t(r$matrix)); colnames(sir) <- names(sol$initial_state); sir$time <- r$t
    if (nrow(sir) != length(tg$time)) return(NULL)
    mu <- matrix(NA_real_, n_years, length(formulas))
    for (i in seq_along(formulas)) {
      v <- tryCatch(evaluate_formula(formulas[i], sir, sol$parms, tg$time, partition),
                    error = function(e) NULL)
      if (is.null(v) || length(v) != n_years) return(NULL)
      mu[, i] <- v
    }
    mu
  }

  # ---- log-posterior ------------------------------------------------------
  loglik <- function(theta) {
    mu <- predict_mu(theta)
    if (is.null(mu) || any(!is.finite(mu))) return(-Inf)
    sigma <- if (needs_sigma) theta[match("sigma", par_names)] else NA_real_
    total <- 0
    for (i in seq_along(formulas)) {
      fam <- like_specs[[i]]$family
      yi  <- data$matrix_data_points[, i]; mui <- mu[, i]; sci <- scale_cols[i]
      phi <- if (i %in% phi_idx) theta[match(sprintf("phi%d", i), par_names)] else NA_real_
      lpdf <- function(val, idx) {
        m <- pmax(mui[idx], 1e-9)
        switch(fam,
               gaussian  = stats::dnorm(val, mui[idx], sigma * sci, log = TRUE),
               lognormal = stats::dlnorm(val, log(m), sigma, log = TRUE),
               poisson   = stats::dpois(val, m, log = TRUE),
               negbin    = stats::dnbinom(val, mu = m, size = phi, log = TRUE),
               stop(sprintf("family '%s' not supported in the R backend.", fam)))
      }
      lcdf <- function(val, idx, lower) {
        m <- pmax(mui[idx], 1e-9)
        switch(fam,
               gaussian  = stats::pnorm(val, mui[idx], sigma * sci, lower.tail = lower, log.p = TRUE),
               lognormal = stats::plnorm(val, log(m), sigma, lower.tail = lower, log.p = TRUE),
               poisson   = stats::ppois(val, m, lower.tail = lower, log.p = TRUE),
               negbin    = stats::pnbinom(val, mu = m, size = phi, lower.tail = lower, log.p = TRUE),
               stop(sprintf("family '%s' not supported in the R backend.", fam)))
      }
      obs <- which(obsm[, i] == 1)
      cl  <- which(data$cens_mask[, i]  == 1)   # left-censored (<=L): P(Y <= L)
      rl  <- which(data$lcens_mask[, i] == 1)   # right-censored (>=L): P(Y >= L)
      iv  <- which(data$interval_mask[, i] == 1)
      if (length(obs)) total <- total + sum(lpdf(yi[obs], obs))
      if (length(cl))  total <- total + sum(lcdf(limit_eff[cl, i],   cl, lower = TRUE))
      if (length(rl))  total <- total + sum(lcdf(llimit_eff[rl, i],  rl, lower = FALSE))
      # Interval [A,B]: hard edges -> log P(A <= Y <= B) = log(CDF(B) - CDF(A_eff))
      # via log-diff-exp (A_eff = A-1 for a discrete family). Soft [A,B]~s -> a
      # family-general plateau penalty on the mean (no CDF), matching Julia.
      if (length(iv)) {
        disc_i <- fam %in% c("poisson", "negbin")
        for (k in iv) {
          A  <- data$ilow_mat[k, i];    B  <- data$iupp_mat[k, i]
          sl <- data$idev_lo_mat[k, i]; su <- data$idev_hi_mat[k, i]
          if (!is.na(sl) && sl > 0) {
            total <- total - (max(A - mui[k], 0) / sl + max(mui[k] - B, 0) / su)
          } else {
            lB <- lcdf(B, k, lower = TRUE)
            lA <- lcdf(if (disc_i) A - 1 else A, k, lower = TRUE)
            total <- total + if (is.finite(lB) && lB > lA) lB + log1p(-exp(lA - lB)) else -Inf
          }
        }
      }
      # Asymmetric A->B / A+/A-: a hard one-sided anchor at A (log P(Y>A) if the
      # soft direction is up, else log P(Y<=A)) plus a soft linear damping of the
      # mean's excess in that direction (scale = dev). Matches the Julia builder.
      av <- which(data$asym_mask[, i] == 1)
      if (length(av)) for (k in av) {
        A   <- data$asym_val_mat[k, i]
        dev <- data$asym_dev_mat[k, i]; if (is.na(dev) || dev <= 0) dev <- Inf
        if (data$asym_dir_mat[k, i] > 0)
          total <- total + lcdf(A, k, lower = FALSE) - max(mui[k] - A, 0) / dev
        else
          total <- total + lcdf(A, k, lower = TRUE)  - max(A - mui[k], 0) / dev
      }
      if (!is.finite(total)) return(-Inf)
    }
    # Informative priors (uniform priors handled by the sampler's prior box).
    lpri <- 0
    for (k in seq_len(n_model)) if (!is_unif[k]) lpri <- lpri + .rprior_logd(model_specs[[k]], theta[k])
    if (needs_sigma) lpri <- lpri + stats::dnorm(sigma, 0, 1, log = TRUE)              # half-Normal(0,1)
    for (i in phi_idx) lpri <- lpri + stats::dgamma(theta[match(sprintf("phi%d", i), par_names)],
                                                    shape = 2, scale = 5, log = TRUE)  # Gamma(2,5)
    val <- total + lpri
    if (is.finite(val)) val else -Inf
  }

  # ---- optional MLE seeding (parity with the Julia init_from_optim path) --
  # A good start matters a lot for a gradient-free sampler on correlated
  # posteriors. point_init is a natural-scale named vector of the fitted
  # quantities (from a quick optim run); map it into the sampled coordinates
  # and use it to initialise the DEzs population.
  start_mat <- NULL
  if (!is.null(point_init)) {
    s0 <- numeric(length(par_names)); names(s0) <- par_names
    for (k in seq_len(n_model)) {
      v <- suppressWarnings(as.numeric(point_init[[model_names[k]]]))
      if (length(v) != 1 || is.na(v)) v <- (low[k] + upp[k]) / 2
      if (is_unif[k]) {
        a <- as.numeric(model_specs[[k]]$args)
        s0[k] <- min(max((v - a[1]) / (a[2] - a[1]), 1e-3), 1 - 1e-3)
      } else s0[k] <- min(max(v, low[k]), upp[k])
    }
    if (needs_sigma) s0[match("sigma", par_names)] <- 0.5
    for (i in phi_idx) s0[match(sprintf("phi%d", i), par_names)] <- 10
    nz <- max(3L, length(par_names) + 1L)            # DEzs needs >= 3 start rows
    start_mat <- t(vapply(seq_len(nz), function(j) {
      jit <- s0 + stats::rnorm(length(s0), 0, 0.05) * (upp - low)
      pmin(pmax(jit, low + 1e-6), upp - 1e-6)
    }, numeric(length(par_names))))
  }

  # ---- sample (DEzs: differential-evolution, gradient-free) ---------------
  prior <- BayesianTools::createUniformPrior(lower = low, upper = upp)
  setup <- BayesianTools::createBayesianSetup(likelihood = loglik, prior = prior,
                                              names = par_names)
  iters <- max(3000L, as.integer(bc$iter) * max(as.integer(bc$chains), 1L))
  settings <- list(iterations = iters, burnin = floor(iters / 2),
                   message = !isFALSE(bc$progress))
  if (!is.null(start_mat)) settings$startValue <- start_mat
  out <- tryCatch(
    BayesianTools::runMCMC(setup, sampler = "DEzs", settings = settings),
    error = function(e) {                            # bad startValue -> default start
      settings$startValue <- NULL
      BayesianTools::runMCMC(setup, sampler = "DEzs", settings = settings)
    })

  draws_df <- as.data.frame(BayesianTools::getSample(out, parametersOnly = TRUE, coda = FALSE))
  names(draws_df) <- par_names

  rhat <- tryCatch(BayesianTools::gelmanDiagnostics(out)$psrf[, 1],
                   error = function(e) rep(NA_real_, length(par_names)))
  ess <- tryCatch(as.numeric(coda::effectiveSize(
    BayesianTools::getSample(out, parametersOnly = TRUE, coda = TRUE))),
    error = function(e) rep(NA_real_, length(par_names)))

  summary_df <- data.frame(
    parameter = par_names,
    mean = colMeans(draws_df),
    sd   = apply(draws_df, 2, stats::sd),
    rhat = as.numeric(rhat)[seq_along(par_names)],
    ess  = as.numeric(ess)[seq_along(par_names)],
    stringsAsFactors = FALSE)

  list(draws = draws_df, summary = summary_df,
       n_chains = 3L, iter = iters, chains = "R-DEzs (gradient-free)",
       model_code = paste0("# Gradient-free R sampler (BayesianTools DEzs); no ",
                           "Julia/Turing model is generated for backend = \"r\"."),
       prior_spec = prior_spec, like_specs = like_specs,
       scale_cols = scale_cols, stream_names = formulas)
}
