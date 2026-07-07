test_that("test-integration-fit", {
# ============================================================
# test-integration-fit.R   (needs Julia + a scenario folder)
# End-to-end: load_scenario -> MLE fit -> save_fit/load_fit -> plot_fit ->
# extract_code. Skips cleanly unless FITCM_SCENARIO_DIR is set and setup.R
# (Julia bridge) initialises. See README.md.
# ============================================================

LBL <- "integration-fit"
if (th_skip_if_no_scenario(LBL)) { th_summary(LBL); quit(save = "no") }
if (!th_have_setup())            { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

dir <- th_scenario_dir(); f <- th_scenario_files(dir)
if (!nzchar(f$dummy)) { th_skip(LBL, "no dataDummy workbook"); th_summary(LBL); quit(save = "no") }

th_section("load_scenario")
sc <- chk_ok("load_scenario reads the folder",
             load_scenario(dir, combined_file = f$combined,
                           dummy_file = f$dummy, params_file = f$params))
chk("has dataCombined", !is.null(sc$dataCombined))
chk("has modelParams",  !is.null(sc$modelParams))
chk("has plot_path",    !is.null(sc$plot_path))

th_section("fitCompartmentalModel (MLE / lbfgsb) -- may take a while")
fit <- chk_ok("MLE fit runs",
              fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                                    checkpoint_file = tempfile(fileext = ".rds")))
chk("class compartmentalFit", inherits(fit, "compartmentalFit"))
chk("fit succeeded", isTRUE(fit$success))
chk("point estimate present", !is.null(fit$point))
chk("data has censoring masks", all(c("cens_mask","lcens_mask") %in% names(fit$data)))

th_section("save_fit / load_fit round-trip")
p <- tempfile(fileext = ".rds")
chk_ok("save_fit writes", save_fit(fit, p))
fit2 <- chk_ok("load_fit restores", load_fit(p))
chk("restored is a fit", inherits(fit2, "compartmentalFit"))
chk("loss rebuilt callable", is.function(fit2$loss))
chk_equal("best error preserved", fit2$best_state$error, fit$best_state$error)

th_section("plot_fit")
res <- chk_ok("plot_fit builds", plot_fit(fit, data_dummy = sc$dataDummy))
chk("returns per-stream plots", length(res$plots) > 0)
chk("returns editable code", is.character(res$code) && nchar(res$code) > 0)
chk("returns width/height", is.numeric(res$width) && is.numeric(res$height))

th_section("extract_code")
code <- chk_ok("extract_code emits a script", extract_code(fit))
chk("non-empty script", is.character(code) && nchar(code) > 0)

th_summary(LBL)
})
