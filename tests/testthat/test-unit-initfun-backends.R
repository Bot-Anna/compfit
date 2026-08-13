test_that("test-unit-initfun-backends", {
# ============================================================
# test-unit-initfun-backends.R   (pure R; no Julia/Stan runtime)
# The Julia + Stan emitters put the state-independent Functions into the init
# (`_0`) scope as `bare` + `_0` (matching how parameters are emitted), so the
# derived-state expressions that reference the bare name resolve. And validation
# rejects a State cell that references a STATE-DEPENDENT Function.
# ============================================================
th_load_pure(c("utils.R", "statesAndParams.R", "generateExpressions.R", "numberOfComps.R",
               "compartmentalFunction.R", "priorSpec.R", "bayesJulia.R", "bayesStan.R",
               "validate.R", "fitCompartmentalModel.R"))

mk <- function(states, funs) data.frame(
  Others     = c("startpoint=2015", "endpoint=2018", "partition=4", rep("", max(0, length(states) - 3))),
  States     = states,
  Parameters = c("*Ntot=1000", "*a=2", "*b=3", "beta=[0,1]",
                 rep("", max(0, length(states) - 4))),
  Functions  = c(funs, rep("", max(0, length(states) - length(funs)))),
  check.names = FALSE, stringsAsFactors = FALSE)

th_section("Julia + Stan emit bare + _0 Function defs; derived states resolve")
old <- options(compfit.validate = FALSE)
mp <- mk(c("*S=frac*Ntot", "*I=(1-frac)*Ntot", "*Rr=g*Ntot", ""),
         c("frac <- a/(a+b)", "g <- if_else(time+startpoint<2016, 0.1, 0.9)"))
m  <- .build_model(mp, backend = "r")
ps <- buildPriorSpec(mp)
jl <- buildJuliaBayesModel(ps, list(), character(0), 4, 4, 3, c("S", "I", "Rr"),
                           init_fun_defs = m$init_fun_defs)
st <- buildStanModel(ps, list(), character(0), 4, 4, 3, c("S", "I", "Rr"),
                     ode_function = m$stan_code, init_fun_defs = m$init_fun_defs)
options(old)

chk("Julia defines bare frac", any(grepl("    frac = a_0/(a_0+b_0)", jl, fixed = TRUE)))
chk("Julia aliases frac_0 = frac", any(grepl("    frac_0 = frac", jl, fixed = TRUE)))
chk("Julia g via native ifelse at t=-1",
    any(grepl("    g = ifelse((-1)+2015<2016", jl, fixed = TRUE)))
chk("Stan defines bare frac", any(grepl("  real frac = a_0/(a_0+b_0);", st, fixed = TRUE)))
chk("Stan aliases frac_0 = frac", any(grepl("  real frac_0 = frac;", st, fixed = TRUE)))
chk("Stan g via ternary at t=-1",
    any(grepl("  real g = ((-1)+2015<2016 ? 0.1 : 0.9);", st, fixed = TRUE)))

th_section("validation rejects a State that references a state-dependent Function")
mp_bad <- mk(c("*S=frac*Ntot", "*I=pool*0.1", "*Rr=0", ""),
             c("frac <- a/(a+b)", "pool <- S+I"))   # pool uses states -> state-dependent
chk_error("State referencing a state-dependent Function is rejected",
          validate_modelParams(mp_bad))
chk_ok("State referencing only a state-independent Function is accepted",
       validate_modelParams(mk(c("*S=frac*Ntot", "*I=0", "*Rr=0", ""),
                               c("frac <- a/(a+b)"))))

th_summary("initfun-backends")
})
