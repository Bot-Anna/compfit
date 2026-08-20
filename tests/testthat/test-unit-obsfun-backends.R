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
# ============================================================
th_load_pure(c("utils.R", "evaluate.R", "bayesJulia.R", "bayesStan.R"))

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
chk("defines f_R as a grid vector",     grepl("f_R = \\[x_\\[", jl_txt))
chk("defines testramp as a grid vector", grepl("testramp = \\[x_\\[", jl_txt))
chk("iterates the solve grid",           grepl("for t_grid_ in grid_t", jl_txt))
chk("binds time to the grid point",      grepl("time = t_grid_", jl_txt))
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
chk("declares f_R over the grid",      grepl("f_R\\[", st_txt))
chk("declares testramp over the grid", grepl("testramp\\[", st_txt))
chk("startpoint injected as a literal",
    grepl("2015", st_txt) && !grepl("\\bstartpoint\\b", st_txt))
chk("uses a Stan for-loop over the grid", grepl("for \\(", st_txt))
chk("no R assignment arrows survive", !grepl("<-", st_txt, fixed = TRUE))
chk("every statement is ;-terminated",
    all(grepl(";\\s*$", grep("=", st, value = TRUE))))
chk("empty selection -> no Stan lines",
    length(.derived_functions_to_stan(of_none, 2015, Inf)) == 0)
})
