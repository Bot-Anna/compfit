test_that("test-integration-extract-code", {
# ============================================================
# test-integration-extract-code.R   (needs Julia + a scenario folder)
# Round-trips extract_code: fit a model, emit the self-contained "all" script
# (inlined data + helpers + loss + an L-BFGS-B fit), run it in a fresh
# environment, and confirm it reproduces the fitted quantities. Defaults to the
# minimal fixture; skips cleanly without a working Julia bridge.
# ============================================================

LBL <- "integration-extract-code"
if (th_skip_if_no_scenario(LBL)) { th_summary(LBL); quit(save = "no") }
if (!th_have_setup())            { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

dir <- th_scenario_dir(); f <- th_scenario_files(dir)
sc  <- load_scenario(dir, combined_file = f$combined,
                     dummy_file = if (nzchar(f$dummy)) f$dummy else f$combined,
                     params_file = f$params)
fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                             checkpoint_file = tempfile(fileext = ".rds"))

th_section("emit + execute the self-contained script")
code <- chk_ok("extract_code(what='all') emits a script", extract_code(fit, what = "all"))
tf <- tempfile(fileext = ".R"); writeLines(code, tf)
e <- new.env(parent = globalenv())
ran <- chk_ok("emitted script runs end-to-end (re-registers Julia, re-fits)",
              sys.source(tf, envir = e))
chk("script produced fitted_par", !is.null(e$fitted_par) && is.numeric(e$fitted_par))

th_section("re-fit reproduces the original fitted quantities")
orig <- c(fit$point$parms, fit$point$initial_state)         # natural scale
common <- intersect(names(e$fitted_par), names(orig))
chk("fitted quantities recovered by name", length(common) > 0)
# Simple (likely unimodal) fixture: the re-fit should land close to the original.
chk("re-fit matches original within tolerance",
    isTRUE(all.equal(unname(e$fitted_par[common]), unname(orig[common]),
                     tolerance = 1e-2)))

th_summary(LBL)
})
