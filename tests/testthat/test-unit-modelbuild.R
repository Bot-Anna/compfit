test_that("test-unit-modelbuild", {
# ============================================================
# test-unit-modelbuild.R   (pure R; no Julia)
# numberOfComps numeric coercion of Level_ indices (the text-read fix), and
# .bounds handling of named / empty / unnamed fitted quantities.
# ============================================================
th_load_pure(c("utils.R", "numberOfComps.R", "fitCompartmentalModel.R"))

th_section("numberOfComps: Level_ indices coerced numerically")
# Simulate modelParams read all-as-text: compartment indices arrive as strings.
mp <- data.frame(`Level_1` = c("1", "2", "10", NA),
                 Other       = c("a", "b", "c", "d"),
                 check.names = FALSE, stringsAsFactors = FALSE)
nc <- numberOfComps(mp)
chk("max is 10 not lexicographic '2'", nc$number_of_comps == 10)
chk("number_of_comps is numeric", is.numeric(nc$number_of_comps))
chk("compartment_cols detected by ^Level_", identical(nc$compartment_cols, "Level_1"))

th_section(".bounds: named fitted quantities")
sap_named <- list(
  lower_states = c(X1 = 0),  upper_states = c(X1 = 10),
  lower_params = c(beta = 0), upper_params = c(beta = 1),
  initial_guesses = c(5, 0.5)            # aligned to c(states, params)
)
b <- .bounds(sap_named)
chk_equal("init_norm midpoints", b$init_norm, c(0.5, 0.5))
chk("lower names preserved", identical(names(b$lower), c("X1", "beta")))

th_section(".bounds: fully-fixed sheet (zero fitted quantities)")
sap_empty <- list(
  lower_states = numeric(0), upper_states = numeric(0),
  lower_params = numeric(0), upper_params = numeric(0),
  initial_guesses = numeric(0)
)
be <- chk_ok("empty sap does not error", .bounds(sap_empty))
chk("empty lower length 0", length(be$lower) == 0)
chk("empty init_norm length 0", length(be$init_norm) == 0)

th_section(".bounds: unnamed non-empty fitted quantities error")
sap_unnamed <- list(
  lower_states = 0, upper_states = 1,            # NOTE: no names
  lower_params = numeric(0), upper_params = numeric(0),
  initial_guesses = 0.5
)
chk_error("unnamed bounds rejected", .bounds(sap_unnamed))

th_section("statesAndParams: bare names in a State expr -> their _0 alias")
mpb <- data.frame(
  States     = c("*S=frac*Ntot", "*I=(1-frac)*Ntot", "*R=0"),
  Parameters = c("frac=[0,1]", "*Ntot=1000", "beta=[0,1]"),
  check.names = FALSE, stringsAsFactors = FALSE)
sapb <- statesAndParams(mpb)
chk("bare params rewritten to _0 in the fixed State expr",
    sapb$states_fixed[["S"]] == "frac_0*Ntot_0")
chk("second State expr rewritten too",
    sapb$states_fixed[["I"]] == "(1-frac_0)*Ntot_0")
# Already-_0 expressions are untouched (no double-suffix).
mpc <- mpb; mpc$States <- c("*S=frac_0*Ntot_0", "*I=(1-frac_0)*Ntot_0", "*R=0")
chk("existing _0 expressions are a no-op (not doubled)",
    statesAndParams(mpc)$states_fixed[["S"]] == "frac_0*Ntot_0")

th_summary("modelbuild")
})
