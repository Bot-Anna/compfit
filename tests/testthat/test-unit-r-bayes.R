test_that("test-unit-r-bayes", {
  # Gradient-free pure-R Bayesian backend (BayesianTools DEzs): a full posterior
  # fit with NO Julia. Skipped on CRAN (it runs an MCMC + an MLE seed, so it is
  # slower than a typical example), but runs on CI and locally. Checks that the
  # samples object is structurally compatible with the Julia path so every
  # posterior_*() helper works, and that MLE-seeding lands near the MLE.
  skip_on_cran()
  skip_if_not_installed("deSolve")
  skip_if_not_installed("BayesianTools")
  skip_if_not_installed("coda")

  dir <- fixture_dir("minimal")
  sc <- load_scenario(dir, combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

  withr::local_seed(1)   # reproducible DEzs draws WITHOUT leaking the global seed
  fit <- suppressMessages(fitCompartmentalModel(
    sc$modelParams, sc$dataCombined, method = "bayes",
    solver = solver_control(backend = "r"),
    bayes  = bayes_control(iter = 500, chains = 3, progress = FALSE,
                           init_from_optim = TRUE)))

  chk("samples present (no Julia)", !is.null(fit$samples))
  chk("draws carry normalised uniform columns", "k_n" %in% names(fit$samples$draws))
  chk("draws carry the noise hyperparameter", "sigma" %in% names(fit$samples$draws))

  # posterior_*() helpers all operate on the (compatible) samples object.
  pd <- posterior_draws(fit)          # natural scale (k_n -> k, etc.)
  chk("posterior_draws maps to natural names", all(c("k","m","X1") %in% names(pd)))
  pm <- posterior_means(fit)
  chk("k within prior bounds",  pm[["k"]]  >= 0  && pm[["k"]]  <= 1)
  chk("X1 within prior bounds", pm[["X1"]] >= 50 && pm[["X1"]] <= 150)
  # MLE-seeded, the posterior mean should sit near the MLE (k ~ 0.18, X1 ~ 115).
  chk("seeded posterior near MLE for k",  abs(pm[["k"]]  - 0.178) < 0.1)
  chk("seeded posterior near MLE for X1", abs(pm[["X1"]] - 115)   < 25)

  ci <- posterior_intervals(fit, 0.9)
  chk("intervals returned", is.data.frame(ci) && nrow(ci) >= 1)

  rep <- chk_ok("posterior_report runs",
                suppressWarnings(invisible(utils::capture.output(posterior_report(fit)))))
  chk_ok("plot_prior_posterior builds", plot_prior_posterior(fit))
  chk_ok("identifiability_report runs",
         suppressWarnings(invisible(utils::capture.output(
           identifiability_report(fit, plots = FALSE)))))
})
