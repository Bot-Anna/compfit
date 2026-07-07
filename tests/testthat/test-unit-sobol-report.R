test_that("test-unit-sobol-report", {
# ============================================================
# test-unit-sobol-report.R   (pure R; no Julia)
# sobol_report()'s reporting logic, driven by a synthetic tidy Sobol table
# (output, parameter, first_order, total) -- no ODE solve needed.
# ============================================================
th_load_pure(c("sensitivity.R"))

# Well-behaved indices in [0,1]; p1 dominates, p2 negligible in B.
sob <- data.frame(
  output      = c("A", "A", "B", "B"),
  parameter   = c("p1", "p2", "p1", "p2"),
  first_order = c(0.80, 0.10, 0.55, 0.00),
  total       = c(0.85, 0.15, 0.62, 0.01),
  stringsAsFactors = FALSE)

th_section("valid indices: structure + rankings")
r <- NULL; invisible(capture.output(r <- sobol_report(sob, plots = FALSE)))
chk("returns the documented pieces",
    all(c("table","stream_sums","parameter_ranking","not_converged",
          "n_outside_range","na_fraction","heatmap") %in% names(r)))
chk("nothing outside [0,1]", r$n_outside_range == 0)
chk("no NA totals", r$na_fraction == 0)
chk("no stream flagged as not-converged", length(r$not_converged) == 0)
chk("interaction column = total - first_order",
    isTRUE(all.equal(r$table$total - r$table$first_order, r$table$interaction)))
chk("global ranking orders p1 above p2", r$parameter_ranking$parameter[1] == "p1")
chk("plots = FALSE -> no heatmap", is.null(r$heatmap))

th_section("out-of-range indices flagged as not converged")
bad <- sob; bad$total[1] <- 1.6                 # impossible Sobol total
rb <- NULL; invisible(capture.output(rb <- sobol_report(bad, plots = FALSE)))
chk("counts the out-of-range index", rb$n_outside_range >= 1)
chk("flags the offending stream (A) as not converged", "A" %in% rb$not_converged)

th_section("accepts a precomputed result without recomputing")
# Passing a tidy data.frame must NOT trigger a (Julia) recompute; it just reports.
chk_ok("runs on a data.frame with no fit/model in scope",
       invisible(capture.output(sobol_report(sob, plots = FALSE))))

th_summary("sobol-report")
})
