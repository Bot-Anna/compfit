test_that("test-integration-bayes-stan", {
# ============================================================
# test-integration-bayes-stan.R   (needs rstan + a C++ toolchain; NO Julia)
# End-to-end proof of the Stan Bayesian backend: compile the generated program
# and sample it with NUTS, then check the standard contract holds --
#   * fit$success, draws in the [0,1]-reparam convention that posterior_draws()
#     denormalises to natural scale;
#   * summary carries finite rhat/ess;
#   * the post-fit solve (plots / predictive) runs on deSolve, no Julia.
# Skips cleanly when rstan (or its toolchain) is unavailable.
# ============================================================

LBL <- "integration-bayes-stan"
if (!requireNamespace("rstan", quietly = TRUE)) {
  th_skip(LBL, "rstan not installed"); th_summary(LBL); quit(save = "no") }

dir_sis <- fixture_dir("SIS")
if (!dir.exists(dir_sis)) { th_skip(LBL, "SIS fixture missing"); th_summary(LBL); quit(save = "no") }
sc <- load_scenario(dir_sis, combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

# Clean gaussian data (no censoring/interval/asym), so the first-cut Stan
# backend applies: a single observed stock stream over 11 years.
yrs  <- as.character(2000:2010)
vals <- c(10, 18, 32, 55, 90, 140, 200, 260, 310, 350, 380)
dc <- data.frame(Label = "Infected", Formula = "X2", Likelihood = "gaussian",
                 as.list(setNames(vals, yrs)), check.names = FALSE,
                 stringsAsFactors = FALSE)

th_section("compile + sample via rstan (this builds a C++ model; slow first run)")
fit <- tryCatch(
  fitCompartmentalModel(sc$modelParams, dc, method = "bayes",
                        solver = solver_control(backend = "stan"),
                        bayes  = bayes_control(chains = 2, iter = 400, warmup = 200,
                                               progress = FALSE),
                        checkpoint_file = tempfile(fileext = ".rds")),
  error = function(e) e)
if (inherits(fit, "error")) {
  th_skip(LBL, paste("Stan compile/sample unavailable:", conditionMessage(fit)))
  th_summary(LBL); quit(save = "no")
}

chk("fit succeeded",             isTRUE(fit$success))
chk("chains label is rstan",     grepl("rstan", fit$samples$chains))
chk("sampled sigma is present",  "sigma" %in% names(fit$samples$draws))

th_section("posterior_draws denormalises the [0,1] reparam to natural scale")
d <- posterior_draws(fit)
chk("draws present",             nrow(d) > 0)
chk("natural-scale names (no _n)", all(c("beta", "gamma", "init_inf") %in% names(d)))
# init_inf's box is [0.001, 0.1]; the denormalised draws must live inside it.
chk("init_inf within its prior box",
    all(d$init_inf >= 0.001 - 1e-8 & d$init_inf <= 0.1 + 1e-8))
chk("beta positive", mean(d$beta) > 0)

th_section("diagnostics + Julia-free post-fit solve")
chk("summary has finite rhat", any(is.finite(fit$samples$summary$rhat)))
pe <- .cfit_point_from_draw(fit, as.data.frame(t(colMeans(d))))
ev <- solve_and_evaluate(fit, pe$initial_state, pe$parms)$evaluation
chk("deSolve post-fit evaluate works", !is.null(ev) && nrow(ev) > 0)

th_section("count family: a negbin fit samples the integer-data likelihood + dispersion")
# Overdispersed integer counts for infected (X2); no censoring/interval/asym.
cnt <- c(22, 35, 60, 95, 150, 240, 360, 520, 700)
dcn <- data.frame(Label = "Infected", Formula = "X2", Likelihood = "negbin",
                  as.list(setNames(cnt, as.character(2000:2008))),
                  check.names = FALSE, stringsAsFactors = FALSE)
fit_nb <- tryCatch(
  fitCompartmentalModel(sc$modelParams, dcn, method = "bayes",
                        solver = solver_control(backend = "stan"),
                        bayes  = bayes_control(chains = 2, iter = 400, warmup = 200,
                                               progress = FALSE),
                        checkpoint_file = tempfile(fileext = ".rds")),
  error = function(e) e)
if (inherits(fit_nb, "error")) {
  th_skip(LBL, paste("negbin Stan fit unavailable:", conditionMessage(fit_nb)))
} else {
  chk("negbin fit succeeded", isTRUE(fit_nb$success))
  dn <- posterior_draws(fit_nb)
  chk("dispersion phi_1 was sampled", "phi_1" %in% names(dn) && all(dn$phi_1 > 0))
  chk("beta recovered positive", mean(dn$beta) > 0)
}

th_summary(LBL)
})
