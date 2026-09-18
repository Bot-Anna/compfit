# ============================================================
# test-integration-obsfun-bayes.R   (needs Julia and/or rstan)
# A data-stream Formula that names a TIME-VARYING Function, sampled for real on
# the Julia and Stan Bayesian backends. test-unit-obsfun-backends.R checks the
# emitted strings; those checks passed while every such fit failed -- Julia with
# a MethodError (Vector*Vector), Stan with a redeclared identifier -- so this
# file samples the generated models.
#
# Data are simulated from the medium SIR fixture at beta = 0.9, gamma = 0.3 with
# 3% multiplicative noise, observing X2 and the recovery flow gamma_t*X2, where
# gamma_t <- gamma*(1+ramp*time). The gamma window is the discriminating check:
# fitting the same data with the ramp removed (what a backend that dropped
# gamma_t would effectively do) puts gamma near 0.37, well outside it.
# Each backend is its own test_that block, so one missing skips only itself.
# ============================================================

.obsfun_truth <- c(beta = 0.9, gamma = 0.3)

.obsfun_scenario <- function() {
  sc <- load_scenario(fixture_dir("medium"), combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
  mp <- sc$modelParams
  mp$Conditions[1] <- NA_character_          # Conditions are not scored on the Bayes path

  truth <- mp
  truth$Parameters[1] <- sprintf("*beta=%s",  .obsfun_truth[["beta"]])
  truth$Parameters[2] <- sprintf("*gamma=%s", .obsfun_truth[["gamma"]])
  forms <- c("X2", "gamma_t*X2")
  sim <- simulate_model(truth, data_dummy = data.frame(Label = forms, Formula = forms),
                        solver = solver_control(backend = "r"))

  # The year-end rows the loss (and both Bayesian models) read a stock at.
  tg  <- sim$time_grid
  idx <- seq(tg$partition + 1,
             tg$partition * (tg$endpoint - tg$startpoint + 1) + 1, by = tg$partition)
  yrs <- as.character(seq(tg$startpoint, tg$endpoint))

  set.seed(1)
  dc <- data.frame(Label = c("Infected", "Recovery flow"), Formula = forms,
                   Likelihood = "gaussian", stringsAsFactors = FALSE)
  for (k in seq_along(yrs))
    dc[[yrs[k]]] <- vapply(forms, function(f)
      round(sim$evaluation[[f]][idx[k]] * exp(rnorm(1, 0, 0.03)), 3), numeric(1))
  list(modelParams = mp, dataCombined = dc)
}

.obsfun_check <- function(fit, backend, code_marker) {
  chk(sprintf("%s: fit succeeded", backend), isTRUE(fit$success))
  chk(sprintf("%s: model binds the grid vector", backend),
      grepl(code_marker, fit$samples$model_code, fixed = TRUE))
  s    <- posterior_summary(fit)
  mean <- setNames(s$mean, s$parameter)
  rhat <- setNames(s$rhat, s$parameter)
  chk(sprintf("%s: beta near truth (%.4f)", backend, mean[["beta"]]),
      abs(mean[["beta"]] - .obsfun_truth[["beta"]]) < 0.05)
  chk(sprintf("%s: gamma near truth (%.4f; ~0.37 if gamma_t were ignored)", backend, mean[["gamma"]]),
      abs(mean[["gamma"]] - .obsfun_truth[["gamma"]]) < 0.03)
  chk(sprintf("%s: chains mixed (rhat beta %.3f, gamma %.3f)", backend, rhat[["beta"]], rhat[["gamma"]]),
      rhat[["beta"]] < 1.1 && rhat[["gamma"]] < 1.1)
  mean[c("beta", "gamma")]
}

.obsfun_means <- new.env()

test_that("test-integration-obsfun-bayes: Julia", {
  LBL <- "integration-obsfun-bayes (julia)"
  if (!th_have_setup()) th_skip(LBL, "Julia not available")
  sc  <- .obsfun_scenario()
  fit <- chk_ok("julia bayes fit runs", suppressWarnings(fitCompartmentalModel(
    sc$modelParams, sc$dataCombined, method = "bayes",
    solver = solver_control(backend = "julia"),
    bayes  = bayes_control(chains = 2, iter = 400, warmup = 200, seed = 7, progress = FALSE),
    checkpoint_file = tempfile(fileext = ".rds"))))
  .obsfun_means$julia <- .obsfun_check(fit, "julia", "gamma_t__grid = [x_[1]")
})

test_that("test-integration-obsfun-bayes: Stan", {
  LBL <- "integration-obsfun-bayes (stan)"
  if (!requireNamespace("rstan", quietly = TRUE)) th_skip(LBL, "rstan not installed")
  sc  <- .obsfun_scenario()
  fit <- tryCatch(suppressWarnings(fitCompartmentalModel(
    sc$modelParams, sc$dataCombined, method = "bayes",
    solver = solver_control(backend = "stan"),
    bayes  = bayes_control(chains = 2, iter = 400, warmup = 200, seed = 7, progress = FALSE),
    checkpoint_file = tempfile(fileext = ".rds"))), error = function(e) e)
  if (inherits(fit, "error"))
    th_skip(LBL, paste("Stan compile/sample unavailable:", conditionMessage(fit)))
  .obsfun_means$stan <- .obsfun_check(fit, "stan", "vector[n_grid] gamma_t__grid;")
})

test_that("test-integration-obsfun-bayes: Julia and Stan agree", {
  if (is.null(.obsfun_means$julia) || is.null(.obsfun_means$stan))
    th_skip("integration-obsfun-bayes (agreement)", "needs both the Julia and Stan fits")
  d <- abs(.obsfun_means$julia - .obsfun_means$stan)
  chk(sprintf("posterior means agree across backends (|d beta| %.4f, |d gamma| %.4f)",
              d[["beta"]], d[["gamma"]]),
      all(d < 0.02))
})
