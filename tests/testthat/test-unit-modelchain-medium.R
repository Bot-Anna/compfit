test_that("test-unit-modelchain-medium", {
# ============================================================
# test-unit-modelchain-medium.R   (pure R; no Julia)
# Builds the MEDIUM fixture (3-compartment SIR) through the parsing/codegen
# chain, exercising the features the minimal model lacks: a quadratic infection
# term, a time-varying Function, a Condition penalty, and parameter-dependent
# function-defined initial states.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R"))

fx <- file.path(fixture_dir("medium"), "modelParams.csv")
chk("fixture exists", file.exists(fx))
mp <- read_data_file(fx, text_cols = TRUE)

cs  <- numberOfComps(mp)
sap <- statesAndParams(mp)

th_section("structure")
chk("three compartments", cs$number_of_comps == 3)
chk("beta, gamma fitted", setequal(names(sap$params_fitted), c("beta", "gamma")))
chk("N0/init_inf/ramp fixed", all(c("N0", "init_inf", "ramp") %in% names(sap$params_fixed)))

th_section("parameter-dependent function-defined initial states")
chk("X1 and X2 are function-defined states",
    all(c("X1", "X2") %in% names(sap$list_states_functions)))
# The state functions take the params (named with a _0 suffix) and compute the
# initial condition. X1 = N0*(1-init_inf); X2 = init_inf*N0.
ps <- c(N0_0 = 1000, init_inf_0 = 0.01)
chk_equal("X1(0) = N0*(1-init_inf) = 990", sap$list_states_functions[["X1"]](ps), 990)
chk_equal("X2(0) = init_inf*N0 = 10",      sap$list_states_functions[["X2"]](ps), 10)

th_section("Condition -> penalty expression")
ex <- generateExpressions(number_of_comps = cs$number_of_comps,
  states_fitted = sap$states_fitted, states_fixed = sap$states_fixed,
  states_functions = sap$states_functions, params_fitted = sap$params_fitted,
  params_fixed = sap$params_fixed, params_functions = sap$params_functions,
  conditions = mp$Conditions)
chk("penalty enforces beta > gamma",
    any(grepl("gamma", ex$penalty)) && any(grepl("beta", ex$penalty)))

th_section("codegen: quadratic infection + time-varying recovery")
cf <- compartmentalFunction(modelParams = mp, compartment_structure = cs,
                            sir_expression = ex$sir, return_expression = ex$ret)
jc <- cf$julia_code
chk("time-varying Function spliced in", grepl("gamma_t=gamma\\*\\(1\\+ramp\\*t\\)", jc))
chk("frequency-dependent infection term", grepl("\\(-beta\\*X\\[2\\]\\*X\\[1\\]\\)/N1", jc))
chk("S loses to infection (dX[1] = secOrd)", grepl("dX\\[1\\] = secOrd_1_2", jc))
chk("I gains infection, loses to recovery", grepl("dX\\[2\\] = -gamma_t\\*X\\[2\\]-secOrd_1_2", jc))
chk("R gains recovered", grepl("dX\\[3\\] = gamma_t\\*X\\[2\\]", jc))

th_summary("modelchain-medium")
})
