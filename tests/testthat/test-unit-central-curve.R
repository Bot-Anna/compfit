test_that("test-unit-central-curve", {
# ============================================================
# test-unit-central-curve.R   (pure R; no Julia)
# Which central curve plot_fit() draws for a Bayesian fit.
#
# A plug-in curve ("median"/"mean") solves the ODE once at a summary of the
# draws; a band summarises the OUTPUT of every draw at each time point. Those
# agree only if the parameter -> trajectory map is linear (means) or monotone
# (medians), and a compartmental model is neither -- so a plug-in line need not
# centre its own band. central = "pointwise" draws the ensemble's pointwise
# median, which IS that band's 50% quantile.
#
# Default is "median": it matches save_scenario(summary=) and the counterfactual
# reference, and survives a skewed/funnel posterior that drags a mean into a tail.
# ============================================================
skip_on_cran()                       # runs an MCMC
skip_if_not_installed("deSolve")
skip_if_not_installed("BayesianTools")
skip_if_not_installed("ggplot2")

sc <- load_scenario(fixture_dir("medium"), combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
withr::local_seed(1)
fit <- suppressMessages(fitCompartmentalModel(sc$modelParams, sc$dataCombined,
         method = "bayes", solver = solver_control(backend = "r"),
         bayes = bayes_control(iter = 400, chains = 2, progress = FALSE,
                               init_from_optim = TRUE)))
chk("bayes fit produced", isTRUE(fit$success))

th_section("get_central_point(summary=) picks the summary, and defaults to median")
d   <- posterior_draws(fit)
med <- get_central_point(fit, "median")$parms
mn  <- get_central_point(fit, "mean")$parms
chk_equal("median matches the draws' median", med[["beta"]],
          stats::median(d$beta), tol = 1e-8)
chk_equal("mean matches the draws' mean", mn[["beta"]], mean(d$beta), tol = 1e-8)
chk("default is the median", identical(get_central_point(fit)$parms, med))
chk_error("an unknown summary is rejected", get_central_point(fit, "mode"))

th_section("bands carry the 50% quantile of the same ensemble")
b <- suppressWarnings(compfit:::.cfit_predictive_bands(fit, n_draws = 60, band_type = "mean"))
st <- fit$data$names_data_points[1]
chk("band frame has med alongside lo/hi",
    all(c("lo", "hi", "med") %in% names(b[[st]])))
chk("med lies inside the interval", all(b[[st]]$med >= b[[st]]$lo & b[[st]]$med <= b[[st]]$hi))

th_section("central = 'pointwise' reproduces the band's median exactly")
pe <- get_central_point(fit, "median")
ev <- solve_and_evaluate(fit, pe$initial_state, pe$parms)$evaluation
withr::local_seed(2)
b2 <- suppressWarnings(compfit:::.cfit_predictive_bands(fit, n_draws = 60, band_type = "mean"))
withr::local_seed(2)
pw <- suppressWarnings(compfit:::.cfit_pointwise_central(fit, ev, n_draws = 60,
                                                         n_rep = 5, data_dummy = NULL))
hit <- match(b2[[st]]$date, pw$date)
chk_equal("pointwise curve == band med on the shared dates",
          pw[[st]][hit], b2[[st]]$med, tol = 1e-8)
chk("plug-in and pointwise curves actually differ",
    !isTRUE(all.equal(pw[[st]], ev[[st]])))

th_section("plot_fit accepts every central=, including over a spaghetti cloud")
for (ce in c("median", "mean", "pointwise"))
  chk_ok(sprintf("plot_fit(central = '%s') builds", ce),
         suppressWarnings(suppressMessages(
           plot_fit(fit, central = ce, band_type = "mean", n_draws = 40))))
chk_ok("spaghetti + pointwise builds",
       suppressWarnings(suppressMessages(
         plot_fit(fit, band_type = "spaghetti", spaghetti_draws = 20,
                  central = "pointwise", n_draws = 40))))
chk_error("an unknown central= is rejected",
          suppressMessages(plot_fit(fit, central = "typo")))

th_section("pointwise needs draws: an MLE fit warns and falls back")
mf <- suppressMessages(fitCompartmentalModel(sc$modelParams, sc$dataCombined,
        method = "lbfgsb", solver = solver_control(backend = "r")))
w <- NULL
p <- withCallingHandlers(suppressMessages(plot_fit(mf, central = "pointwise")),
       warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
chk("warns that an MLE fit has one point estimate", grepl("point estimate", w %||% ""))
chk("still returns panels", length(p$plots) > 0)
})
