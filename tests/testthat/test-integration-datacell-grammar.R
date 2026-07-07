test_that("test-integration-datacell-grammar", {
# ============================================================
# test-integration-datacell-grammar.R   (needs Julia)
# Round-trips the interval [A,B] / asymmetric A->B / A+ data-cell grammar
# through the REAL Turing/Julia Bayesian path (not just the pure-R backend
# covered by test-unit-datacell-grammar.R): a fit with mixed observed/
# interval/asym cells registers the generated @model in Julia and samples it,
# and the generated model_code carries the new branches.
# ============================================================

LBL <- "integration-datacell-grammar"
dir <- fixture_dir("minimal")
if (!dir.exists(dir))  { th_skip(LBL, "minimal fixture missing"); th_summary(LBL); quit(save = "no") }
if (!th_have_setup())  { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

sc <- chk_ok("load_scenario (csv fixture)",
             load_scenario(dir, combined_file = "dataCombined.csv",
                           dummy_file = "dataDummy.csv", params_file = "modelParams.csv"))

th_section("mixed observed / interval / asymmetric cells")
dc <- sc$dataCombined
dc$Likelihood <- c("gaussian; asym=5", "gaussian")
dc[["2016"]][1] <- "[75,90]"      # interval (hard edges, CDF)
dc[["2017"]][1] <- "67->72"      # asym, explicit target (dev 5, up)
dc[["2018"]][1] <- "55+"         # asym, global deviation (dev 5, up)
dc[["2019"]][1] <- "[70,95]~8"   # soft interval (no-CDF plateau penalty)

th_section("Bayesian fit registers + samples the Turing @model (Julia)")
fit <- chk_ok("bayes fit runs",
              suppressWarnings(fitCompartmentalModel(
                sc$modelParams, dc, method = "bayes",
                bayes = bayes_control(iter = 200, chains = 1, warmup = 100,
                                      progress = FALSE, init_from_optim = TRUE))))
chk("fit succeeded",   isTRUE(fit$success))
chk("draws produced",  !is.null(fit$samples$draws) && nrow(fit$samples$draws) > 0)
chk("model_code has the interval branch (cf_logsubexp)",
    grepl("cf_logsubexp", fit$samples$model_code, fixed = TRUE))
chk("model_code has the asymmetric branch (asym_dir_mat)",
    grepl("asym_dir_mat", fit$samples$model_code, fixed = TRUE))
chk("model_code has the soft-interval branch (idev_lo_mat, no CDF)",
    grepl("idev_lo_mat", fit$samples$model_code, fixed = TRUE))

th_summary(LBL)
})
