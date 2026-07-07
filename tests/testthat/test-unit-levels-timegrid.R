test_that("test-unit-levels-timegrid", {
# ============================================================
# test-unit-levels-timegrid.R   (pure R; no Julia)
# Two behaviours:
#   * levels are detected by the '_' PREFIX (not the word "Level"), and a sheet
#     with NO level column defaults to a single level holding every compartment
#     (so N1 = total_pop and quadratic terms still normalise);
#   * .time_grid() warns when startpoint/endpoint is not a plain integer year
#     (a date or decimal is silently truncated to its leading number).
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R",
               "fitCompartmentalModel.R", "evaluate.R", "simulate.R"))

sc <- load_scenario(fixture_dir("SIR"), combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

th_section("no level column -> single auto-level of all compartments")
mp0 <- sc$modelParams; mp0[["_Level1"]] <- NULL
chk("numberOfComps still counts via States", numberOfComps(mp0)$number_of_comps == 3)
chk("no '_' column detected", length(numberOfComps(mp0)$compartment_cols) == 0)
m0 <- chk_ok("quadratic model builds without a level column",
             build_compartmental_model(mp0, sc$dataCombined, solver = solver_control(backend = "r")))
chk("N1 spans every compartment",
    grepl("real N1 = X[1]+X[2]+X[3]", m0$model$stan_code, fixed = TRUE))

th_section("level column name is free (detected by '_' prefix, not \"Level\")")
mpA <- sc$modelParams; mpA[["_Age1"]] <- mpA[["_Level1"]]; mpA[["_Level1"]] <- NULL
chk("'_Age1' is recognised as a level column",
    identical(numberOfComps(mpA)$compartment_cols, "_Age1"))
chk_ok("model with _Age1 builds",
       build_compartmental_model(mpA, sc$dataCombined, solver = solver_control(backend = "r")))

th_section("no-level model simulates and conserves a closed population")
mps <- read_data_file(file.path(fixture_dir("SIR_sim"), "modelParams.csv"))
mps[["_Level1"]] <- NULL
dds <- read_data_file(file.path(fixture_dir("SIR_sim"), "dataDummy.csv"))
sim <- chk_ok("no-level SIR simulates",
              simulate_model(mps, data_dummy = dds, solver = solver_control(backend = "r")))
tot <- rowSums(sim$sir_out[c("X1", "X2", "X3")])
chk_equal("closed population conserved (N0 = 5000)", max(abs(tot - 5000)), 0, tol = 1e-2)

th_section(".time_grid warns on a non-integer year")
warns <- function(mp) {
  w <- character(0)
  withCallingHandlers(compfit:::.time_grid(mp),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  w
}
set_start <- function(v) { mp <- sc$modelParams
  mp$Others[mp$Others == "startpoint=2010"] <- paste0("startpoint=", v); mp }
chk("date startpoint warns",    any(grepl("integer year", warns(set_start("2010-06-01")))))
chk("decimal startpoint warns", any(grepl("integer year", warns(set_start("2010.5")))))
chk("clean integer startpoint is silent", length(warns(sc$modelParams)) == 0)

th_summary("levels-timegrid")
})
