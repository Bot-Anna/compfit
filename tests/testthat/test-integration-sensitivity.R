test_that("test-integration-sensitivity", {
# ============================================================
# test-integration-sensitivity.R   (needs Julia)
# local_sensitivity() (derivative-based, at the fit) and a small Sobol run on the
# medium SIR fixture. Confirms structure, plots, and a physical correctness
# signal: the conserved total S+I+R has ~zero sensitivity to the parameters.
# ============================================================

LBL <- "integration-sensitivity"
dir <- fixture_dir("medium")
if (!dir.exists(dir)) { th_skip(LBL, "medium fixture missing"); th_summary(LBL); quit(save = "no") }
if (!th_have_setup()) { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

sc  <- load_scenario(dir, combined_file = "dataCombined.csv",
                     dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                             checkpoint_file = tempfile(fileext = ".rds"))

th_section("local_sensitivity: structure")
ls <- chk_ok("runs", local_sensitivity(fit, type = "relative",
                                        data_dummy = sc$dataDummy, progress = FALSE))
chk("class localSensitivity", inherits(ls, "localSensitivity"))
chk("summary has output/parameter/sensitivity",
    all(c("output", "parameter", "sensitivity") %in% names(ls$summary)))
chk("one row per (output, parameter)",
    nrow(ls$summary) == length(unique(ls$summary$output)) * length(unique(ls$summary$parameter)))
chk("trajectory carries dated sensitivities",
    nrow(ls$trajectory) > 0 && "date" %in% names(ls$trajectory))

th_section("physical check: conserved total has ~zero sensitivity")
tot <- ls$summary[ls$summary$output == "X1+X2+X3", "sensitivity"]
chk("S+I+R insensitive to all parameters", length(tot) > 0 && all(abs(tot) < 1e-6))

th_section("plots build")
chk_ok("tornado plot", plot_local_sensitivity(ls, "tornado"))
chk_ok("trajectory plot", plot_local_sensitivity(ls, "trajectory"))

th_section("Sobol (small) returns indices with a progress bar")
if (requireNamespace("sensitivity", quietly = TRUE)) {
  sl <- chk_ok("sobol_loss runs", sobol_loss(fit, n = 16, progress = FALSE))
  chk("has first_order + total columns",
      all(c("first_order", "total") %in% names(sl)))
} else {
  th_skip("sobol_loss", "install the 'sensitivity' package to run")
}

th_summary(LBL)
})
