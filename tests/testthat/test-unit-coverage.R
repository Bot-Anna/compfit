test_that("test-unit-coverage", {
# ============================================================
# test-unit-coverage.R   (pure R; no Julia)
# The coverage / scaling pattern (#6): a reported stream is a fraction of the
# true model quantity, observed = rho * true, with rho a fitted parameter (any
# prior). This needs no new cell type -- the `Formula` column already references
# parameters by name -- so here we just confirm it builds, fits, and recovers
# rho when a second stream anchors the epidemic scale.
# ============================================================

th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R"))

set.seed(7)
lin  <- matrix("0", 3, 3); lin[3, 2] <- "(gamma)"; lin[2, 2] <- "-(gamma)"
quad <- list(rep("0", 9), rep("0", 9), rep("0", 9)); quad[[1]][2] <- "*2*-beta"
others <- c("startpoint=2000", "endpoint=2005", "partition=4")
sheet <- function(params) {
  states <- c("*X1=N0_0*(1-i0_0)", "*X2=i0_0*N0_0", "*X3=0")
  nr <- 9L; pad <- function(x) c(as.character(x), rep("", nr - length(x)))
  df <- data.frame(`_Level1` = pad(1:3), Others = pad(others), States = pad(states),
                   Functions = pad(character(0)), Parameters = pad(params),
                   Conditions = pad(character(0)), check.names = FALSE, stringsAsFactors = FALSE)
  for (i in 1:3) { df[[paste0("Linear", i)]] <- pad(lin[i, ]); df[[paste0("Quadratic", i)]] <- pad(quad[[i]]) }
  df
}

th_section("a fitted rho appears in Formula and is a fitted parameter")
fit_sheet <- sheet(c("beta=[0,3]", "gamma=[0,1]", "*N0=5000", "i0=[0.0005,0.02]", "rho=Beta(2,2)"))
sap <- statesAndParams(fit_sheet)
chk("rho is a fitted parameter", "rho" %in% names(sap$params_fitted))
chk("rho is not a param-function", !("rho" %in% names(sap$params_functions)))

th_section("simulate observed = rho * true; fit recovers rho")
truth <- sheet(c("*beta=1.1", "*gamma=0.4", "*N0=5000", "*i0=0.002", "*rho=0.4"))
m <- compfit:::.build_model(truth, backend = "r"); tg <- compfit:::.time_grid(truth)
sol <- compfit:::.recover_solution(setNames(numeric(0), character(0)), m$sap, tg)
r <- solveWithR(unlist(sol$initial_state), c(min(tg$time), max(tg$time)), unlist(sol$parms), tg$time)
idx <- seq(tg$partition + 1, tg$partition * 6 + 1, by = tg$partition); I <- r$matrix[2, idx]
yrs <- 2000:2005
Itrue <- rpois(6, pmax(I, 0.1)); Irep <- rpois(6, pmax(0.4 * I, 0.1))
mkrow <- function(lab, form, vals) {
  d <- data.frame(Label = lab, Formula = form, Likelihood = "poisson", Weight = 1,
                  Average = signif(1 / mean(vals), 4), check.names = FALSE, stringsAsFactors = FALSE)
  for (k in seq_along(yrs)) d[[as.character(yrs[k])]] <- as.character(vals[k]); d
}
dc <- rbind(mkrow("True prevalence", "X2", Itrue), mkrow("Reported (rho*X2)", "rho*X2", Irep))
fit <- fitCompartmentalModel(fit_sheet, dc, method = "lbfgsb",
                             solver = solver_control(backend = "r"),
                             checkpoint_file = tempfile(fileext = ".rds"))
chk("coverage fit succeeds", isTRUE(fit$success))
chk_equal("rho recovered near truth (0.4)", get_point(fit)$parms[["rho"]], 0.4, tol = 0.15)

th_summary("unit-coverage")
})
