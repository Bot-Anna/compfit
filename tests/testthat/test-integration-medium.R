test_that("test-integration-medium", {
# ============================================================
# test-integration-medium.R   (needs Julia)
# Fits the MEDIUM fixture (3-compartment SIR with a quadratic infection term,
# a time-varying recovery Function, a Condition, and parameter-dependent
# function-defined initial states) end to end, confirming the ODE solves and the
# parameter-dependent initial states resolve to the expected values.
# ============================================================

LBL <- "integration-medium"
dir <- fixture_dir("medium")
if (!dir.exists(dir))  { th_skip(LBL, "medium fixture missing"); th_summary(LBL); quit(save = "no") }
if (!th_have_setup())  { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

sc <- chk_ok("load_scenario (csv fixture)",
             load_scenario(dir, combined_file = "dataCombined.csv",
                           dummy_file = "dataDummy.csv", params_file = "modelParams.csv"))

th_section("MLE fit of the SIR model")
fit <- chk_ok("fit runs", fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                                                checkpoint_file = tempfile(fileext = ".rds")))
chk("fit succeeded", isTRUE(fit$success))
chk("beta and gamma estimated", all(c("beta", "gamma") %in% names(fit$point$parms)))

th_section("parameter-dependent initial states resolve (N0=1000, init_inf=0.01)")
is0 <- fit$point$initial_state
chk_equal("X1(0) = N0*(1-init_inf) = 990", unname(is0["X1"]), 990, tol = 1e-6)
chk_equal("X2(0) = init_inf*N0 = 10",      unname(is0["X2"]), 10,  tol = 1e-6)
chk_equal("X3(0) = 0",                     unname(is0["X3"]), 0,   tol = 1e-6)

th_section("solve + plot")
ev <- chk_ok("solve_and_evaluate runs",
             solve_and_evaluate(fit, fit$point$initial_state, fit$point$parms, sc$dataDummy)$evaluation)
chk("trajectory is finite and non-negative",
    all(is.finite(ev$X2)) && all(ev$X2 >= -1e-6))
chk("conservation: S+I+R stays ~constant (closed SIR)",
    { tot <- ev$X1 + ev$X2 + ev$X3; max(abs(tot - tot[1])) < 1e-3 * tot[1] })
res <- chk_ok("plot_fit builds", plot_fit(fit, data_dummy = sc$dataDummy))
chk("per-stream plots produced", length(res$plots) > 0)

# dataCombined takes precedence: a dummy whose Formula equals a fitted stream is
# skipped (the fitted panel, with data points, is kept) and a message is emitted.
th_section("dataCombined precedence over dataDummy")
s1 <- fit$data$names_data_points[1]
dd <- data.frame(Label = paste(s1, "(dummy)"), Formula = s1, stringsAsFactors = FALSE)
saw_msg <- FALSE
res2 <- withCallingHandlers(
  plot_fit(fit, data_dummy = dd),
  message = function(m) {
    if (grepl("both dataCombined and dataDummy", conditionMessage(m))) saw_msg <<- TRUE
    invokeRestart("muffleMessage")
  })
chk("precedence message emitted", saw_msg)
chk("shared-formula panel is the fitted one (not the dummy)",
    !identical(res2$plots[[s1]]$labels$title, paste(s1, "(dummy)")))

# Every panel carries the SAME date scale. Dummy panels used to set no
# scale_x_date at all and inherited ggplot's default, so they disagreed with the
# fitted panels on both tick positions and label format.
th_section("shared x-axis scale across fitted and dummy panels")
x_breaks <- function(p) {
  b <- NULL
  for (s in p$scales$scales) if (any(c("x", "date") %in% s$aesthetics)) b <- s$breaks
  b
}
all_b <- lapply(res$plots, x_breaks)
chk("no panel falls back to the default x scale",
    !any(vapply(all_b, is.null, logical(1))))
chk("all panels share identical date breaks",
    all(vapply(all_b, function(b) identical(b, all_b[[1]]), logical(1))))

th_summary(LBL)
})
