test_that("test-unit-examples", {
# ============================================================
# test-unit-examples.R   (pure R; no Julia)
# Guards the four shipped epidemic example scenarios (SI / SIS / SIR / SEIR)
# under inst/extdata/. For each: load_scenario -> .build_model (R backend) ->
# .prepare_data must succeed, the compartment count is as expected, and the
# special data cells each example advertises land in the right mask. Cheap
# insurance that the example sheets and their data grammar stay valid.
# ============================================================

th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R"))

load_ex <- function(name) {
  dir <- fixture_dir(name)
  chk(paste0(name, " fixture exists"), dir.exists(dir))
  load_scenario(dir, combined_file = "dataCombined.csv",
                dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
}
prep  <- function(sc) compfit:::.prepare_data(sc$dataCombined, compfit:::.time_grid(sc$modelParams))
ncomp <- function(sc) numberOfComps(sc$modelParams)$number_of_comps
nfam  <- function(sc) length(unique(tolower(trimws(sc$dataCombined$Likelihood))))
any_true <- function(x) !is.null(x) && any(x, na.rm = TRUE)

th_section("SI — 2-compartment, poisson, a missing cell")
sc <- load_ex("SI"); r <- prep(sc)
chk("SI builds + prepares",   !is.null(r$obs_mask))
chk("SI has 2 compartments",  ncomp(sc) == 2)
chk("SI has a missing cell",  any(!r$obs_mask))

th_section("SIS — time-varying Function, Condition, gaussian, <= and A->B cells")
sc <- load_ex("SIS"); r <- prep(sc)
chk("SIS has 2 compartments",            ncomp(sc) == 2)
chk("SIS has a left-censored (<=) cell", any_true(r$cens_mask))
chk("SIS has an asymmetric (A->B) cell", any_true(r$asym_mask))

th_section("SIR — two families (negbin I + gaussian R), interval + right-censored on R")
sc <- load_ex("SIR"); r <- prep(sc)
chk("SIR has 3 compartments",             ncomp(sc) == 3)
chk("SIR has an interval [A,B] cell",     any_true(r$interval_mask))
chk("SIR has a right-censored (>=) cell", any_true(r$lcens_mask))
chk("SIR declares two families",          nfam(sc) == 2)

th_section("SEIR — three families (E gaussian, I negbin, R lognormal), >= / x / [A,B]")
sc <- load_ex("SEIR"); r <- prep(sc)
chk("SEIR has 4 compartments",             ncomp(sc) == 4)
chk("SEIR has a right-censored (>=) cell", any_true(r$lcens_mask))
chk("SEIR has an interval [A,B] cell",      any_true(r$interval_mask))
chk("SEIR has a missing (x) cell",          any(!r$obs_mask))
chk("SEIR declares three families",         nfam(sc) == 3)

th_section("regression: single-stream fit with Average != 1 (diag() dimension bug)")
# diag(sqrt(av)) treats a length-1 value as a DIMENSION (diag(0.04) -> 0x0),
# which broke single-stream fits whose Average differs from 1. SI is one stream
# with Average = 1/mean, so this fit exercises the fix directly.
sc <- load_ex("SI")
chk("SI Average is 1/mean, not 1", abs(as.numeric(sc$dataCombined$Average[1]) - 1) > 1e-6)
fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                             solver = solver_control(backend = "r"),
                             checkpoint_file = tempfile(fileext = ".rds"))
chk("single-stream fit with Average != 1 succeeds", isTRUE(fit$success))

th_section("regression: sheet with n_params > n^2 builds without a recycling warning")
# SIS is a 2-compartment model (n^2 = 4) with 5 parameter rows, so the Quadratic
# column is padded to 5. The builder must slice each Quadratic column to n^2
# rather than reshape the padded column (which recycled and warned).
sc <- load_ex("SIS")
w <- character(0)
withCallingHandlers(
  compfit:::.build_model(sc$modelParams, backend = "r"),
  warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
chk("no matrix-recycling warning building SIS",
    !any(grepl("sub-multiple or multiple", w)))

th_section("missing/blank rate-coefficient columns are silently all-zero (not error)")
sc <- load_ex("minimal")
grab <- function(mp) {
  ww <- character(0)
  ok <- withCallingHandlers({ compfit:::.build_model(mp, backend = "r"); TRUE },
          warning = function(cw) { ww <<- c(ww, conditionMessage(cw)); invokeRestart("muffleWarning") })
  list(ok = isTRUE(ok), coef_warn = any(grepl("rate-coefficient", ww)))
}
chk("complete sheet: no coefficient warning", !grab(sc$modelParams)$coef_warn)
# A missing Linear<i>/Quadratic<i> column means "no such terms for this
# compartment" -- treated as all-zero, and silently (no warning, no error).
mpq <- sc$modelParams; mpq$Quadratic2 <- NULL           # omit a Quadratic column
gq <- grab(mpq)
chk("missing Quadratic column builds silently (all-zero)", gq$ok && !gq$coef_warn)
mpl <- sc$modelParams; mpl$Linear2 <- NULL              # omit a Linear column
gl <- grab(mpl)
chk("missing Linear column builds silently (all-zero)", gl$ok && !gl$coef_warn)
mpb <- sc$modelParams; mpb$Linear1[1] <- ""             # blank a single cell
gb <- grab(mpb)
chk("blank coefficient cell builds, treated as 0 (no warning)", gb$ok && !gb$coef_warn)
# A compartment left with NO terms at all (its only rate column removed) must
# still emit a well-formed zero derivative, not a malformed empty `dX[i] =`.
sc3  <- load_ex("SIR"); mp3 <- sc3$modelParams; mp3$Linear3 <- NULL   # X3's only inflow
mod3 <- suppressWarnings(compfit:::.build_model(mp3, backend = "r"))
chk("term-less compartment yields dX = 0 (well-formed)",
    any(grepl("dX[3] = 0", strsplit(mod3$stan_code, "\n")[[1]], fixed = TRUE)))

th_section("SIR_priors / SEIR_priors -- all prior distributions + families")
# The prior-showcase folders use "_Level1" (not "_Country1") and, between them,
# every prior distribution and every likelihood family.
prior_dists <- function(sc) {
  ps <- sc$modelParams$Parameters
  ps <- ps[nzchar(ps) & !grepl("^\\*", ps) & grepl("=", ps)]
  unname(vapply(ps, function(e) parsePrior(sub("^[^=]*=", "", e))$dist, character(1)))
}
sc <- load_ex("SIR_priors")
chk("SIR_priors uses _Level (not _Country)",
    any(grepl("^_Level", names(sc$modelParams))) && !any(grepl("Country", names(sc$modelParams))))
chk("SIR_priors builds", { compfit:::.build_model(sc$modelParams, backend = "r"); TRUE })
d1 <- prior_dists(sc)
chk("SIR_priors has StudentT/Normal/Beta/Uniform priors",
    all(c("StudentT", "Normal", "Beta", "Uniform") %in% d1))
chk("SIR_priors has negbin + gaussian",
    all(c("negbin", "gaussian") %in% tolower(sc$dataCombined$Likelihood)))

sc <- load_ex("SEIR_priors")
chk("SEIR_priors builds", { compfit:::.build_model(sc$modelParams, backend = "r"); TRUE })
d2 <- prior_dists(sc)
chk("SEIR_priors has Gamma/LogNormal/Beta/Uniform priors",
    all(c("Gamma", "LogNormal", "Beta", "Uniform") %in% d2))
chk("SEIR_priors has gaussian + poisson + lognormal",
    all(c("gaussian", "poisson", "lognormal") %in% tolower(sc$dataCombined$Likelihood)))
chk("together: all six prior distributions covered",
    all(c("Uniform", "Normal", "LogNormal", "Beta", "Gamma", "StudentT") %in% c(d1, d2)))

th_section("SIR_named -- named compartments (S/I/R) reproduce the X1..Xn twin")
# The naming feature: States column names the compartments AND fixes their
# order; _Level and the data formulas reference them by name; rate columns are
# Linear<name>/Quadratic<name>. Fitting must be identical to the X1..Xn SIR.
sc_named <- load_ex("SIR_named")
chk("SIR_named has 3 compartments",      ncomp(sc_named) == 3)
chk("SIR_named registry is S, I, R",
    identical(numberOfComps(sc_named$modelParams)$comp_names, c("S", "I", "R")))
chk("SIR_named validates", { validate_modelParams(sc_named$modelParams); TRUE })

sc_x <- load_ex("SIR")
m_x     <- build_compartmental_model(sc_x$modelParams,     sc_x$dataCombined,
                                     solver = solver_control(backend = "r"))
m_named <- build_compartmental_model(sc_named$modelParams, sc_named$dataCombined,
                                     solver = solver_control(backend = "r"))
chk("named model exposes comp_names S/I/R",
    identical(m_named$model$structure$comp_names, c("S", "I", "R")))
# Loss at a shared normalised point must match to numerical tolerance: same ODE,
# same data, only the compartment NAMES differ.
d  <- length(m_x$bounds$lower)
set.seed(11); pt <- runif(d)
es <- function() { e <- new.env(); e$error <- Inf; e$par <- NULL; e }
chk_equal("named-compartment loss equals X1..Xn loss",
          m_named$loss(pt, es()), m_x$loss(pt, es()))

# A full fit recovers name-ordered initial states (S, I, R in canonical order).
fit_named <- fitCompartmentalModel(sc_named$modelParams, sc_named$dataCombined,
                                   method = "lbfgsb",
                                   solver = solver_control(backend = "r"),
                                   checkpoint_file = tempfile(fileext = ".rds"))
chk("named-compartment fit succeeds", isTRUE(fit_named$success))
chk("recovered initial_state is named S/I/R in order",
    identical(names(fit_named$point$initial_state), c("S", "I", "R")))

th_summary("unit-examples")
})
