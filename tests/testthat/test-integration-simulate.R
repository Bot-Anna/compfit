test_that("test-integration-simulate", {
# ============================================================
# test-integration-simulate.R   (needs Julia)
# End-to-end proof that simulate_model() runs a fully-fixed sheet through the
# REAL Julia ODE path (not just the pure-R backend covered by
# test-unit-simulate.R), and that the Julia and R backends agree on the
# trajectory for the same fixed model. Skips cleanly unless Julia initialises.
# ============================================================

LBL <- "integration-simulate"
dir_sim <- fixture_dir("SEIR_sim")
if (!dir.exists(dir_sim)) {
  th_skip(LBL, "SEIR_sim fixture missing"); th_summary(LBL); quit(save = "no") }
if (!th_have_setup()) { th_skip(LBL, "setup.R/Julia not available"); th_summary(LBL); quit(save = "no") }

mp <- read_data_file(file.path(dir_sim, "modelParams.csv"))
dd <- read_data_file(file.path(dir_sim, "dataDummy.csv"))

th_section("simulate on the Julia backend, no dataCombined")
sim_jl <- chk_ok("Julia simulate runs",
                 simulate_model(mp, data_dummy = dd, solver = solver_control()))
chk("is a compartmentalSim", inherits(sim_jl, "compartmentalSim"))
chk("four compartments",     length(sim_jl$initial_state) == 4)
chk("trajectory finite",     all(is.finite(as.matrix(
    sim_jl$sir_out[names(sim_jl$initial_state)]))))
chk("attack-rate series present",
    "(X2+X3+X4)/(X1+X2+X3+X4)" %in% names(sim_jl$evaluation))

th_section("Julia and R backends agree on the fixed-model trajectory")
sim_r <- chk_ok("R simulate runs",
                simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r")))
last <- function(s) as.numeric(s$sir_out[nrow(s$sir_out), names(s$initial_state)])
for (k in seq_along(sim_jl$initial_state))
  chk_equal(sprintf("compartment %d endpoint matches R backend", k),
            last(sim_jl)[k], last(sim_r)[k], tol = 1e-3 * abs(last(sim_r)[k]) + 1e-6)

th_summary(LBL)
})
