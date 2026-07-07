test_that("test-unit-simulate", {
# ============================================================
# test-unit-simulate.R   (pure R; no Julia)
# simulate_model(): run a FULLY-FIXED sheet forward with no data and no fitting.
# Covers the zero-stream .prepare_data path, the loss-skip in
# build_compartmental_model, the fully-fixed guard, the six *_sim example
# scenarios end-to-end on the R backend, and plot_simulation().
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R",
               "fitCompartmentalModel.R", "evaluate.R", "simulate.R"))

r_backend <- solver_control(backend = "r")

load_sim <- function(name) {
  dir <- fixture_dir(name)
  chk(paste0(name, " fixture exists"), dir.exists(dir))
  list(mp = read_data_file(file.path(dir, "modelParams.csv")),
       dd = read_data_file(file.path(dir, "dataDummy.csv")))
}

th_section("zero-stream data bundle (no dataCombined)")
tg <- compfit:::.time_grid(load_sim("minimal_sim")$mp)
ed <- compfit:::.prepare_data(NULL, tg)
chk("no stream names",            length(ed$names_data_points) == 0)
chk("matrix has 0 columns",       ncol(ed$matrix_data_points) == 0)
chk("weight matrix is 0x0",       all(dim(ed$weight_matrix) == c(0, 0)))
chk("0-row frame == NULL bundle",
    identical(names(compfit:::.prepare_data(data.frame(), tg)), names(ed)))

th_section("build_compartmental_model with no data skips the loss")
m0 <- build_compartmental_model(load_sim("minimal_sim")$mp, NULL, solver = r_backend)
chk("loss is NULL when no streams", is.null(m0$loss))
chk("still a compartmentalModel",   inherits(m0, "compartmentalModel"))

th_section("simulate_model on all six *_sim scenarios (R backend)")
expected_comps <- c(minimal_sim = 2, SI_sim = 2, SIS_sim = 2,
                    SIR_sim = 3, SEIR_sim = 4, medium_sim = 3)
for (nm in names(expected_comps)) {
  ex  <- load_sim(nm)
  sim <- chk_ok(paste(nm, "simulates"),
                simulate_model(ex$mp, data_dummy = ex$dd, solver = r_backend))
  chk(paste(nm, "is a compartmentalSim"), inherits(sim, "compartmentalSim"))
  chk(paste(nm, "compartment count"),
      length(sim$initial_state) == expected_comps[[nm]])
  chk(paste(nm, "trajectory length matches grid"),
      nrow(sim$sir_out) == length(sim$time_grid$time))
  chk(paste(nm, "trajectory is finite"), all(is.finite(as.matrix(
      sim$sir_out[setdiff(names(sim$sir_out), "date")]))))
  chk(paste(nm, "compartments non-negative"),
      all(as.matrix(sim$sir_out[names(sim$initial_state)]) >= -1e-6))
  chk(paste(nm, "dummy series evaluated"),
      all(ex$dd$Formula %in% names(sim$evaluation)))
}

th_section("medium_sim conserves the closed population (SIR, no births/deaths)")
ex  <- load_sim("medium_sim")
sim <- simulate_model(ex$mp, data_dummy = ex$dd, solver = r_backend)
tot <- rowSums(sim$sir_out[c("X1", "X2", "X3")])
chk_equal("total stays at N0 = 1000", max(abs(tot - 1000)), 0, tol = 1e-3)

th_section("fully-fixed guard: a fittable sheet is rejected")
# The plain 'minimal' scenario still has [lo,hi] boxes -> not fully specified.
mp_fit <- read_data_file(file.path(fixture_dir("minimal"), "modelParams.csv"))
chk_error("fittable sheet errors in simulate_model()",
          simulate_model(mp_fit, solver = r_backend))

th_section("print + plot_simulation")
chk_ok("print.compartmentalSim runs", print(sim))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ps <- chk_ok("plot_simulation builds", plot_simulation(sim))
  # "both" = compartments + any distinctly-named series (a dummy formula that
  # coincides with a compartment name is not duplicated).
  n_expected <- length(unique(c(names(sim$initial_state),
                                setdiff(names(sim$evaluation), "date"))))
  chk("one panel per distinct compartment/series",
      length(ps$plots) == n_expected)
  chk("states-only subsets panels",
      length(plot_simulation(sim, which = "states")$plots) ==
        length(sim$initial_state))
}

th_summary("simulate")
})
