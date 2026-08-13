test_that("test-unit-derived", {
# ============================================================
# test-unit-derived.R   (pure R; no Julia)
# .derived_columns(): replaying the model's Functions + level head counts on an
# already-solved trajectory, so data/dummy formulas can reference q_CH, omega_*,
# N_<level>, etc. by name. This is the engine behind the "Functions in formula
# scope" feature used by solve_and_evaluate()/plot_fit().
# ============================================================
th_load_pure(c("utils.R", "evaluate.R"))

# A small 4-compartment, 2-level trajectory (High/Low), 3 annual grid points.
sir_out <- data.frame(
  S_H  = c(100, 90, 80),
  S_L  = c(200, 190, 180),
  I_H  = c(10, 12, 14),
  I_L  = c(5, 6, 7),
  date = as.Date(c("2014-12-31", "2015-12-31", "2016-12-31")),
  check.names = FALSE)
parms <- c(n_H = 20, n_L = 14, q0 = 0.2, q1 = 0.5, beta = 0.007, tau = 6.83)
time  <- c(0, 1, 2)                              # + startpoint -> 2015, 2016, 2017

spec <- list(
  # Declaration order = dependency order (ratio, then a time-only step, then a
  # level-head-dependent omega, then a chained foi). `bad` references an unknown
  # symbol and must be skipped WITHOUT dropping the entries that resolve.
  functions_raw = c(
    "ratio <- n_H/n_L",
    "q <- if_else(time+startpoint<2016, q0, q1)",
    "bad <- nonexistent_symbol*2",
    "omega_H <- n_H*N_H",
    "foi <- beta*tau*I_H/omega_H"),
  level_members = list(c(1L, 3L), c(2L, 4L)),    # H = {S_H,I_H}, L = {S_L,I_L}
  level_names   = c("H", "L"),
  comp_names    = c("S_H", "S_L", "I_H", "I_L"),
  startpoint    = 2015,
  cutoff        = Inf)

d <- .derived_columns(sir_out, parms, time, spec)

th_section("level head counts N_<level> / total_pop")
chk("N_H = S_H + I_H", isTRUE(all.equal(d$N_H, c(110, 102, 94))))
chk("N_L = S_L + I_L", isTRUE(all.equal(d$N_L, c(205, 196, 187))))
chk_equal("positional N1 == N_H", d$N1, d$N_H)
chk_equal("total_pop = N_H + N_L", d$total_pop, c(315, 298, 281))

th_section("time-dependent Function via startpoint (the q_CH case)")
# time+startpoint = 2015,2016,2017; `< 2016` -> only the first year takes q0.
chk_equal("q steps q0 -> q1 at the cutoff year", d$q, c(0.2, 0.5, 0.5))

th_section("scalar Function recycled to the grid")
chk_equal("ratio = n_H/n_L on every row", d$ratio, rep(20 / 14, 3))

th_section("chained Functions resolve against level heads")
chk_equal("omega_H = n_H * N_H", d$omega_H, 20 * c(110, 102, 94))
chk_equal("foi = beta*tau*I_H/omega_H",
          d$foi, 0.007 * 6.83 * c(10, 12, 14) / (20 * c(110, 102, 94)))

th_section("a Function that cannot resolve is skipped, not fatal")
chk("bad entry dropped", is.null(d[["bad"]]))
chk("resolvable entries survive the bad one", all(c("q", "omega_H", "foi") %in% names(d)))

th_summary("derived")
})
