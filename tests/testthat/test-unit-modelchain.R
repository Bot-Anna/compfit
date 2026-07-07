test_that("test-unit-modelchain", {
# ============================================================
# test-unit-modelchain.R   (pure R; no Julia)
# Builds the minimal fixture model through the whole parsing/codegen chain
# (numberOfComps -> statesAndParams -> generateExpressions ->
# compartmentalFunction) WITHOUT registering Julia. Validates both the fixture
# and a large slice of otherwise-untested model-building code.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R"))  # incl. buildJuliaODEFunction (no Julia call)

fx <- file.path(fixture_dir("minimal"), "modelParams.csv")
chk("fixture modelParams.csv exists", file.exists(fx))
mp <- read_data_file(fx, text_cols = TRUE)   # mirrors load_scenario()'s typing

th_section("numberOfComps")
cs <- chk_ok("runs", numberOfComps(mp))
chk("two compartments", cs$number_of_comps == 2)

th_section("statesAndParams")
sap <- chk_ok("runs", statesAndParams(mp))
chk("X1 is fitted", "X1" %in% names(sap$states_fitted))
chk("X2 is fixed (not fitted)", !("X2" %in% names(sap$states_fitted)))
chk("k, m are fitted params", all(c("k", "m") %in% names(sap$params_fitted)))

th_section("generateExpressions")
ex <- chk_ok("runs", generateExpressions(
  number_of_comps = cs$number_of_comps,
  states_fitted = sap$states_fitted, states_fixed = sap$states_fixed,
  states_functions = sap$states_functions, params_fitted = sap$params_fitted,
  params_fixed = sap$params_fixed, params_functions = sap$params_functions,
  conditions = mp$Conditions))
chk("produces a sir expression", !is.null(ex$sir))

th_section("compartmentalFunction (codegen; no Julia registration)")
cf <- chk_ok("runs", compartmentalFunction(
  modelParams = mp, compartment_structure = cs,
  sir_expression = ex$sir, return_expression = ex$ret))
chk("emits Julia source as text", is.character(cf$julia_code) && nchar(cf$julia_code) > 0)
chk("compartmental_function is an R closure", is.function(cf$compartmental_function))
chk("ODE has both state derivatives", grepl("dX\\[1\\]", cf$julia_code) && grepl("dX\\[2\\]", cf$julia_code))
chk("linear coefficients appear", grepl("k\\*X\\[1\\]", cf$julia_code) && grepl("m\\*X\\[2\\]", cf$julia_code))

th_summary("modelchain")
})
