test_that("test-unit-recover-solution", {
# ============================================================
# test-unit-recover-solution.R   (pure R; no Julia)
# .recover_solution() -- rebuilds natural-scale (initial_state, parms) from a
# best-par vector, and the order-contract guard (empty-safe, but trips on a real
# name/order mismatch).
# ============================================================
th_load_pure(c("utils.R", "fitCompartmentalModel.R"))

# Minimal sap: one fitted state (X1), one fitted param (beta), nothing fixed.
sap_ok <- list(
  params_fitted = c(beta = NA_real_),     # only names matter for the guard
  states_fitted = c(X1 = NA_real_),
  params_fixed  = setNames(list(), character(0)),
  states_fixed  = setNames(numeric(0), character(0)),
  non_numeric_fixed = integer(0),
  list_states_functions = list()
)
tg <- list(startpoint = 2013, endpoint = 2015)

th_section("recovers natural-scale solution")
rs <- chk_ok("runs on a consistent best_par",
             .recover_solution(c(X1 = 100, beta = 0.5), sap_ok, tg))
chk_equal("parms beta", rs$parms["beta"], c(beta = 0.5), check_names = TRUE)
chk_equal("initial_state X1", rs$initial_state["X1"], c(X1 = 100), check_names = TRUE)

th_section("order-contract guard")
# A param name not present in sap$params_fitted must trip the guard.
chk_error("param name mismatch errors",
          .recover_solution(c(X1 = 100, gamma = 0.5), sap_ok, tg))
# A missing fitted state must trip the guard too.
chk_error("missing fitted state errors",
          .recover_solution(c(beta = 0.5), sap_ok, tg))

th_section("fully-fixed sheet: zero FITTED, all quantities FIXED (counterfactual case)")
# A filled counterfactual sheet has nothing to fit, but every param/state is
# present as a FIXED value -- .recover_solution must rebuild them from sap.
sap_fixed <- list(
  params_fitted = setNames(numeric(0), character(0)),   # nothing fitted
  states_fitted = setNames(numeric(0), character(0)),
  params_fixed  = list(beta = "0.5"),                   # everything fixed
  states_fixed  = c(X1 = 100),
  non_numeric_fixed = integer(0),
  list_states_functions = list()
)
rf <- chk_ok("empty best_par recovers from fixed quantities",
             .recover_solution(setNames(numeric(0), character(0)), sap_fixed, tg))
chk_equal("fixed param recovered", rf$parms["beta"], c(beta = 0.5), check_names = TRUE)
chk_equal("fixed state recovered", rf$initial_state["X1"], c(X1 = 100), check_names = TRUE)

th_summary("recover-solution")
})
