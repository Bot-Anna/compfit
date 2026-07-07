test_that("test-unit-bayes-r", {
# ============================================================
# test-unit-bayes-r.R   (pure R; no Julia)
# The gradient-free R Bayesian sampler now scores every data-cell type,
# including interval and asymmetric cells. Fit the SIS scenario (which ships an
# asymmetric A->B cell) end-to-end and confirm it samples cleanly -- no
# "unsupported cell" warning, real posterior draws. Skips without BayesianTools.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R",
               "fitCompartmentalModel.R", "priorSpec.R", "bayesR.R", "evaluate.R"))

if (!requireNamespace("BayesianTools", quietly = TRUE)) {
  th_skip("bayes-r", "BayesianTools not installed"); th_summary("bayes-r"); quit(save = "no") }

sc <- load_scenario(fixture_dir("SIS"), combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

th_section("SIS (asymmetric cell) fits on the R backend, no unsupported warning")
w <- character(0)
fit <- withCallingHandlers(
  suppressMessages(fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "bayes",
    solver = solver_control(backend = "r"),
    bayes  = bayes_control(iter = 150, chains = 2, progress = FALSE),
    checkpoint_file = tempfile(fileext = ".rds"))),
  warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })

chk("fit succeeded",             isTRUE(fit$success))
chk("posterior draws produced",  !is.null(fit$samples) && nrow(fit$samples$draws) > 0)
chk("no 'unsupported/ignored cell' warning",
    !any(grepl("IGNORED|does not support|cannot use", w)))

th_summary("bayes-r")
})
