test_that("test-unit-obsfun-backends", {
# ============================================================
# test-unit-obsfun-backends.R   (pure R; no Julia / no Stan)
# Observable-referenced Functions on the BAYESIAN backends.
#
# A data-stream Formula may name a TIME-VARYING Function rather than a scalar
# parameter -- e.g. observed diagnoses = f_R(t)*(R_H+R_L) when the diagnosis rate
# scales up over time. The R/MLE loss gets these via .derived_columns(); Julia and
# Stan build their own observation blocks, so they need the same names recomputed
# as vectors over the solve grid.
#
# Guards: (1) which Functions are selected, (2) that the emitted Julia/Stan code
# defines them, and (3) that a model whose observables reference NO Function emits
# nothing (byte-identical codegen -> existing bayes fits unaffected).
#
# The vectors are named `<name>__grid` and the observables that use them are
# evaluated point by point. Both were found by running a real fit (the string
# checks alone passed): the init scope already binds the bare name as a scalar,
# which Stan rejects as a redeclaration, and `gamma_t*X2` over grid vectors is a
# MethodError in Julia and an outer product in Stan.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R",
               "fitCompartmentalModel.R", "priorSpec.R", "evaluate.R",
               "bayesJulia.R", "bayesStan.R"))

# Time-varying diagnosis rate f_R built from a helper Function (testramp), so the
# selection has to follow a Function -> Function reference, not just direct hits.
spec <- list(
  functions_raw = c(
    "testramp <- lo + (1-lo)/(1+exp(-kk*(time+startpoint-t0)))",
    "f_R <- f_R_late*testramp",
    "unused_rate <- f_R_late*0.5",
    "omega_H <- n_H*N_H",
    "foi <- beta*I_H/omega_H"),
  level_members = list(c(1L, 3L), c(2L, 4L)),
  level_names   = c("H", "L"),
  comp_names    = c("S_H", "S_L", "I_H", "I_L"),
  startpoint    = 2015,
  cutoff        = Inf)

params <- c("lo", "kk", "t0", "f_R_late", "n_H", "beta")
formulas <- c("annual(f_R*(I_H+I_L))", "S_H")

th_section("selection: observable-referenced, state-independent Functions")
of <- .derived_observable_functions(spec, formulas, param_symbols = params)
chk("f_R selected (named by the observable)", "f_R" %in% names(of))
chk("testramp selected (reached via f_R)",    "testramp" %in% names(of))
chk("unused_rate NOT selected (no observable references it)",
    !"unused_rate" %in% names(of))
chk("state-dependent omega_H NOT selected", !"omega_H" %in% names(of))
chk("state-dependent foi NOT selected",     !"foi" %in% names(of))
chk("declaration order preserved (testramp before f_R)",
    which(names(of) == "testramp") < which(names(of) == "f_R"))

th_section("no observable references a Function -> nothing selected")
of_none <- .derived_observable_functions(spec, c("S_H", "annual(beta*I_H)"),
                                         param_symbols = params)
chk("empty selection", length(of_none) == 0)

th_section("Julia: emitted block defines the names before the observables")
jl <- .derived_functions_to_julia(of, startpoint = spec$startpoint, cutoff = spec$cutoff)
chk("emits a non-empty block", length(jl) > 0)
jl_txt <- paste(jl, collapse = "\n")
chk("defines f_R__grid as a grid vector",      grepl("f_R__grid = \\[x_\\[", jl_txt))
chk("defines testramp__grid as a grid vector", grepl("testramp__grid = \\[x_\\[", jl_txt))
chk("never rebinds a bare Function name outside the comprehension",
    !any(grepl("^    (f_R|testramp) = ", jl)))
chk("comprehension assignments are local (no capture of the outer scalar)",
    grepl("local f_R = ", jl_txt) && grepl("local testramp = ", jl_txt))
chk("iterates the solve grid",           grepl("for t_grid_ in grid_t", jl_txt))
chk("binds time to the grid point",      grepl("local time = t_grid_", jl_txt))
chk("startpoint injected as a literal (no bare symbol left)",
    grepl("2015", jl_txt) && !grepl("\\bstartpoint\\b", jl_txt))
chk("R's if_else/ifelse translated away", !grepl("if_else\\(|\\bifelse\\(", jl_txt))
chk("no R assignment arrows survive",     !grepl("<-", jl_txt, fixed = TRUE))
chk("empty selection -> no Julia lines",
    length(.derived_functions_to_julia(of_none, 2015, Inf)) == 0)

th_section("Stan: emitted block defines the names before the observables")
st <- .derived_functions_to_stan(of, startpoint = spec$startpoint, cutoff = spec$cutoff)
chk("emits a non-empty block", length(st) > 0)
st_txt <- paste(st, collapse = "\n")
chk("declares f_R__grid over the grid",      grepl("vector[n_grid] f_R__grid;", st_txt, fixed = TRUE))
chk("declares testramp__grid over the grid", grepl("vector[n_grid] testramp__grid;", st_txt, fixed = TRUE))
chk("never declares a bare Function name (init scope owns it)",
    !grepl("vector[n_grid] f_R;", st_txt, fixed = TRUE))
chk("helper chain reads the earlier grid vector at g",
    grepl("f_R__grid[g] = f_R_late*testramp__grid[g];", st_txt, fixed = TRUE))
chk("startpoint injected as a literal",
    grepl("2015", st_txt) && !grepl("\\bstartpoint\\b", st_txt))
chk("uses a Stan for-loop over the grid", grepl("for \\(", st_txt))
chk("no R assignment arrows survive", !grepl("<-", st_txt, fixed = TRUE))
chk("every statement is ;-terminated",
    all(grepl(";\\s*$", grep("=", st, value = TRUE))))
chk("empty selection -> no Stan lines",
    length(.derived_functions_to_stan(of_none, 2015, Inf)) == 0)

th_section("observables: per-point over the grid only when a Function is used")
cn <- c("X1", "X2", "X3"); gf <- "gamma_t"
jl_obs <- function(f) compfit:::.observable_to_julia(f, 6, 4, cn, grid_funs = gf)
chk("Julia stock reads each operand at the grid point",
    identical(jl_obs("gamma_t*X2"),
              "([ (gamma_t__grid[g_]*sol_grid[2, g_]) for g_ in eachindex(grid_t) ])[annual_idx]"))
chk("Julia annual() integrates the per-point flux",
    identical(jl_obs("annual(gamma_t*X2)"),
              "annual_integral([ (gamma_t__grid[g_]*sol_grid[2, g_]) for g_ in eachindex(grid_t) ], grid_t, 4, 6)"))
chk("Julia scalar-minus-Function stays scalar per point (1-v is a MethodError on vectors)",
    grepl("(1-gamma_t__grid[g_])*sol_grid[2, g_]", jl_obs("(1-gamma_t)*X2"), fixed = TRUE))
st_obs <- function(f) compfit:::.observable_lines_stan(f, 3, cn, grid_funs = gf)
chk("Stan fills obsflux point by point",
    identical(st_obs("gamma_t*X2"),
              c("  vector[n_grid] obsflux3;",
                "  for (g in 1:n_grid) obsflux3[g] = gamma_t__grid[g]*traj[2, g];",
                "  vector[n_years] mu3 = obsflux3[annual_idx];")))
chk("Stan annual() integrates obsflux",
    identical(st_obs("annual(gamma_t*X2)")[3],
              "  vector[n_years] mu3 = annual_integral(obsflux3, tgrid, partition, n_years);"))
for (f in c("X2", "rho*X2", "annual(X2)", "X2/(X1+X2+X3)")) {
  chk(sprintf("Julia '%s' (no Function) byte-identical to the legacy form", f),
      identical(jl_obs(f), compfit:::.observable_to_julia(f, 6, 4, cn)))
  chk(sprintf("Stan '%s' (no Function) byte-identical to the legacy form", f),
      identical(st_obs(f),
                sprintf("  vector[n_years] mu3 = %s;", compfit:::.observable_to_stan(f, cn))))
}

th_section("whole program: medium fixture, observable uses its time-varying gamma_t")
# gamma_t is ALSO an init-scope Function here, so its bare name is already bound
# in the model body -- the collision the per-name `__grid` suffix exists for.
sc  <- load_scenario(fixture_dir("medium"), combined_file = "dataCombined.csv",
                     dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
mod <- compfit:::.build_model(sc$modelParams, backend = "r")
ps  <- buildPriorSpec(sc$modelParams)
tg  <- compfit:::.time_grid(sc$modelParams); ny <- tg$endpoint - tg$startpoint + 1
lk  <- rep(list(parseLikelihood("gaussian")), 2)
fm  <- c("X2", "gamma_t*X2")
chk("fixture binds gamma_t in init scope (precondition)",
    "gamma_t" %in% mod$init_fun_defs$name)

jl_prog <- compfit:::buildJuliaBayesModel(ps, lk, fm, ny, tg$partition,
             mod$structure$number_of_comps, mod$structure$comp_names,
             init_fun_defs = mod$init_fun_defs, derived_spec = mod$derived)
chk("Julia program binds gamma_t__grid",   grepl("gamma_t__grid = [x_[1]", jl_prog, fixed = TRUE))
chk("Julia mu2 uses the per-point form",   grepl("mu2 = ([ (gamma_t__grid[g_]*sol_grid[2, g_])", jl_prog, fixed = TRUE))
chk("Julia mu1 (no Function) unchanged",   grepl("mu1 = (sol_grid[2, :])[annual_idx]", jl_prog, fixed = TRUE))

st_prog <- compfit:::buildStanModel(ps, lk, fm, ny, tg$partition,
             mod$structure$number_of_comps, mod$structure$comp_names, mod$stan_code,
             init_fun_defs = mod$init_fun_defs, derived_spec = mod$derived)
chk("Stan program declares gamma_t__grid", grepl("vector[n_grid] gamma_t__grid;", st_prog, fixed = TRUE))
chk("Stan program declares the bare gamma_t exactly once (init scope)",
    length(gregexpr("real gamma_t = ", st_prog, fixed = TRUE)[[1]]) == 1 &&
      !grepl("vector[n_grid] gamma_t;", st_prog, fixed = TRUE))
if (requireNamespace("rstan", quietly = TRUE)) {
  st_ok <- tryCatch({ rstan::stanc(model_code = st_prog, model_name = "obsfun"); TRUE },
                    error = function(e) conditionMessage(e))
  chk("Stan program passes stanc (no redeclaration, no outer product)", isTRUE(st_ok))
}
})
