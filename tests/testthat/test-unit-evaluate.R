test_that("test-unit-evaluate", {
# ============================================================
# test-unit-evaluate.R   (pure R; no Julia)
# evaluate_formula() -- evaluates data/dummy formulas on a SOLVED trajectory.
# Used by plots, counterfactual, and stacked_plots. Tested against a synthetic
# sir_out so no ODE solve is needed.
# ============================================================
th_load_pure(c("utils.R"))

# Grid: partition = 4, two years (t in [0, 2]); 9 points, dt = 0.25.
partition <- 4
time <- seq(0, 2, length.out = partition * 2 + 1)
sir_out <- data.frame(X1 = rep(1, length(time)),     # constant flux of 1
                      X2 = as.numeric(seq_along(time)))
parms <- c(g = 3)

th_section("stock formulas (bare states + bare params)")
chk_equal("bare state X1", evaluate_formula("X1", sir_out, parms, time, partition), rep(1, length(time)))
chk_equal("arithmetic X2 + 0*X1", evaluate_formula("X2 + 0*X1", sir_out, parms, time, partition), sir_out$X2)
chk_equal("bare parameter g*X1", evaluate_formula("g*X1", sir_out, parms, time, partition), rep(3, length(time)))

th_section("annual(): rolling 1-year integral")
ann <- evaluate_formula("annual(X1)", sir_out, parms, time, partition)
chk("first `partition` points are NA", all(is.na(ann[seq_len(partition)])))
chk_equal("integral of constant 1 over 1 year = 1 (at end)", ann[length(ann)], 1)
chk_equal("integral = 1 at first full window", ann[partition + 1], 1)

th_section("cumulative(): running integral")
cum <- evaluate_formula("cumulative(X1)", sir_out, parms, time, partition)
chk_equal("starts at 0", cum[1], 0)
chk_equal("total = elapsed time (2)", cum[length(cum)], 2)
chk("monotone non-decreasing", all(diff(cum) >= -1e-9))

th_summary("evaluate")
})
