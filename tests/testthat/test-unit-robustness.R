test_that("test-unit-robustness", {
# ============================================================
# test-unit-robustness.R   (pure R; needs deSolve)
# Optimiser robustness: the loss penalises a non-finite (blown-up) solve instead
# of scoring it as a perfect fit, and the stochastic DEoptim method is
# reproducible when seeded.
# ============================================================

th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "validate.R",
               "statesAndParams.R", "generateExpressions.R", "compartmentalFunction.R"))

pad <- function(x, n = 4) c(as.character(x), rep("", n - length(x)))
mp  <- data.frame(`_Level1` = pad(c("1", "2")),
  Others = pad(c("startpoint=2000", "endpoint=2005", "partition=4")),
  States = pad(c("*X1=10", "*X2=5")), Functions = pad(character(0)),
  Parameters = pad("m=[0.01,1]"), Conditions = pad(character(0)),
  Linear1 = pad(c("-m", "0")), Quadratic1 = pad(c("0", "0")),
  Linear2 = pad(c("0", "-m")), Quadratic2 = pad(c("0", "0")),
  check.names = FALSE, stringsAsFactors = FALSE)
mkdc <- function(formula) {
  d <- data.frame(Label = "s", Formula = formula, check.names = FALSE, stringsAsFactors = FALSE)
  for (y in 2000:2005) d[[as.character(y)]] <- "3"; d
}

th_section("a non-finite solve/observable is penalised, not rewarded")
# sqrt(X1-500) is NaN for all reachable X1 (X1 decays from 10) -> a failed
# observable; the loss must return the large penalty, not ~0.
f_nan <- suppressWarnings(fitCompartmentalModel(mp, mkdc("sqrt(X1-500)"), method = "lbfgsb",
           solver = solver_control(backend = "r"), checkpoint_file = tempfile(fileext = ".rds")))
f_ok  <- suppressWarnings(fitCompartmentalModel(mp, mkdc("X1"), method = "lbfgsb",
           solver = solver_control(backend = "r"), checkpoint_file = tempfile(fileext = ".rds")))
d <- length(f_nan$bounds$lower)
chk("non-finite observable -> large penalty",   # suppress the sqrt(-x) NaN warning
    suppressWarnings(f_nan$loss(rep(0.5, d), f_nan$best_state, FALSE)) >= 1e12)
chk("finite observable -> finite loss below the penalty",
    { v <- f_ok$loss(rep(0.5, d), f_ok$best_state, FALSE); is.finite(v) && v < 1e12 })

th_section("seeded DEoptim is reproducible")
run <- function() get_point(fitCompartmentalModel(mp, mkdc("X1"), method = "deoptim",
         control = optim_control(itermax = 30, seed = 42),
         solver = solver_control(backend = "r"), checkpoint_file = tempfile(fileext = ".rds")))$parms[["m"]]
chk("same seed -> identical estimate", isTRUE(all.equal(run(), run())))

th_summary("unit-robustness")
})
