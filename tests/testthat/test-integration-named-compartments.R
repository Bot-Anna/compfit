test_that("test-integration-named-compartments", {
# ============================================================
# test-integration-named-compartments.R   (needs Julia)
# End-to-end proof that the naming feature works through the REAL Julia path
# (not just the pure-R backend covered by test-unit-examples.R):
#   * an MLE fit of the S/I/R-named SIR_named scenario converges to the SAME
#     natural-scale estimates as its X1..Xn twin (SIR) -- so the Julia ODE
#     codegen honours the States-column order and the name-suffixed
#     Linear<name>/Quadratic<name> columns;
#   * a Bayesian fit registers + samples the Turing @model, and the generated
#     model_code maps the named observable `I` to the correct trajectory row
#     (sol_grid[2, :]) -- proving .observable_to_julia's name->index mapping.
# Skips cleanly unless the Julia bridge initialises.
# ============================================================

LBL <- "integration-named-compartments"
dir_named <- fixture_dir("SIR_named"); dir_x <- fixture_dir("SIR")
if (!dir.exists(dir_named) || !dir.exists(dir_x)) {
  th_skip(LBL, "SIR_named/SIR fixtures missing"); th_summary(LBL); quit(save = "no") }
if (!th_have_setup()) { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

ld <- function(d) load_scenario(d, combined_file = "dataCombined.csv",
                                dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
sc_named <- chk_ok("load SIR_named", ld(dir_named))
sc_x     <- chk_ok("load SIR",       ld(dir_x))

th_section("MLE fit of named scenario matches the X1..Xn twin (Julia ODE)")
fit_named <- chk_ok("named MLE fit runs",
                    fitCompartmentalModel(sc_named$modelParams, sc_named$dataCombined,
                                          method = "lbfgsb",
                                          checkpoint_file = tempfile(fileext = ".rds")))
fit_x     <- chk_ok("X1..Xn MLE fit runs",
                    fitCompartmentalModel(sc_x$modelParams, sc_x$dataCombined,
                                          method = "lbfgsb",
                                          checkpoint_file = tempfile(fileext = ".rds")))
chk("named fit succeeded", isTRUE(fit_named$success))
chk("recovered states are named S/I/R in order",
    identical(names(fit_named$point$initial_state), c("S", "I", "R")))
# Same estimation problem, just renamed -> the fitted parameters agree.
for (p in c("beta", "gamma", "init_inf"))
  chk_equal(sprintf("%s matches the X-twin", p),
            fit_named$point$parms[[p]], fit_x$point$parms[[p]], tol = 1e-3)

th_section("Bayesian fit: Turing @model samples + name->row mapping in codegen")
fit_b <- chk_ok("named bayes fit runs",
                suppressWarnings(fitCompartmentalModel(
                  sc_named$modelParams, sc_named$dataCombined, method = "bayes",
                  bayes = bayes_control(iter = 200, chains = 1, warmup = 100,
                                        progress = FALSE, init_from_optim = TRUE))))
chk("fit succeeded",  isTRUE(fit_b$success))
chk("draws produced", !is.null(fit_b$samples$draws) && nrow(fit_b$samples$draws) > 0)
# The infected observable `I` is compartment 2, so it must map to sol_grid[2, :].
chk("named observable I maps to trajectory row 2",
    grepl("sol_grid[2, :]", fit_b$samples$model_code, fixed = TRUE))

th_summary(LBL)
})
