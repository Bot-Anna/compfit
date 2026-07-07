test_that("test-unit-local-sensitivity", {
# ============================================================
# test-unit-local-sensitivity.R   (pure R; no Julia)
# The pure helper .ls_split_point(): overlay a perturbed natural-scale vector
# onto a base point, replacing fitted params/states and keeping the fixed ones.
# ============================================================
th_load_pure(c("sensitivity.R"))

th_section(".ls_split_point: overlay perturbed values, keep fixed ones")
base <- list(parms = c(beta = 0.5, N0 = 1000),
             initial_state = c(X1 = 990, X2 = 10))

# Perturb a param (beta) and a state (X1); leave N0 and X2 untouched.
pe <- .ls_split_point(c(beta = 0.8, X1 = 900), base)
chk_equal("param beta replaced", pe$parms["beta"], c(beta = 0.8), check_names = TRUE)
chk_equal("param N0 kept",       pe$parms["N0"],   c(N0 = 1000),  check_names = TRUE)
chk_equal("state X1 replaced",   pe$initial_state["X1"], c(X1 = 900), check_names = TRUE)
chk_equal("state X2 kept",       pe$initial_state["X2"], c(X2 = 10),  check_names = TRUE)

th_section("no-op when nothing overlaps")
pe2 <- .ls_split_point(c(unknown = 1), base)
chk_equal("params unchanged", pe2$parms, c(beta = 0.5, N0 = 1000), check_names = TRUE)
chk_equal("states unchanged", pe2$initial_state, c(X1 = 990, X2 = 10), check_names = TRUE)

th_summary("local-sensitivity")
})
