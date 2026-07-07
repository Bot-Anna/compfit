test_that("test-unit-modelbuild", {
# ============================================================
# test-unit-modelbuild.R   (pure R; no Julia)
# numberOfComps numeric coercion of _Level indices (the text-read fix), and
# .bounds handling of named / empty / unnamed fitted quantities.
# ============================================================
th_load_pure(c("utils.R", "numberOfComps.R", "fitCompartmentalModel.R"))

th_section("numberOfComps: _Level indices coerced numerically")
# Simulate modelParams read all-as-text: compartment indices arrive as strings.
mp <- data.frame(`_Level` = c("1", "2", "10", NA),
                 Other       = c("a", "b", "c", "d"),
                 check.names = FALSE, stringsAsFactors = FALSE)
nc <- numberOfComps(mp)
chk("max is 10 not lexicographic '2'", nc$number_of_comps == 10)
chk("number_of_comps is numeric", is.numeric(nc$number_of_comps))
chk("compartment_cols detected by ^_", identical(nc$compartment_cols, "_Level"))

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

th_summary("modelbuild")
})
