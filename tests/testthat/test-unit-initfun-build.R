test_that("test-unit-initfun-build", {
# ============================================================
# test-unit-initfun-build.R   (pure R; no Julia)
# End-to-end (R backend, no solve): an initial-state cell may reference a
# state-independent Function -- an algebraic one (frac) and a time-dependent one
# evaluated at t=-1 (g). Exercises BOTH R init paths: the loss's states_params
# block (checked as generated code) and the .recover_solution() closures (checked
# as the recovered initial state).
# ============================================================
th_load_pure(c("utils.R", "statesAndParams.R", "generateExpressions.R",
               "numberOfComps.R", "compartmentalFunction.R", "fitCompartmentalModel.R"))
old <- options(compfit.validate = FALSE); on.exit(options(old))

mp <- data.frame(
  Others     = c("startpoint=2015", "endpoint=2018", "partition=4", "", ""),
  States     = c("*S=frac*Ntot", "*I=(1-frac)*Ntot", "*Rr=g*Ntot", "", ""),
  Parameters = c("*Ntot=1000", "*a=2", "*b=3", "beta=[0,1]", ""),
  Functions  = c("frac <- a/(a+b)",                          # indep (params)
                 "g <- if_else(time+startpoint<2016, 0.1, 0.9)",  # indep (time): 0.1 at t=-1
                 "tot <- S+I+Rr",                            # DEP (states) -- unused at init
                 "", ""),
  check.names = FALSE, stringsAsFactors = FALSE)

m <- .build_model(mp, backend = "r")

th_section("classification picked the right init Functions")
chk("frac + g are init-usable", setequal(names(m$init_funs$independent), c("frac", "g")))
chk("tot is state-dependent", "tot" %in% m$init_funs$dependent)

th_section("loss states_params block defines the Function _0 aliases (t=-1)")
sp <- m$expressions$states_params
chk("frac_0 defined from params", grepl("frac_0 <- a_0/\\(a_0\\+b_0\\)", sp))
chk("g_0 defined with time->-1", grepl("g_0 <- if_else\\(\\(-1\\)\\+2015<2016", sp))
chk("state cell rewritten to _0 alias", grepl("Rr_0<-g_0\\*Ntot_0", gsub(" ", "", sp)))

th_section(".recover_solution resolves the initial state via the closures")
sol <- .recover_solution(c(beta = 0.5), m$sap, m$time_grid)
X <- sol$initial_state
chk_equal("S = frac*Ntot = (2/5)*1000", X[["S"]], 400)
chk_equal("I = (1-frac)*Ntot", X[["I"]], 600)
chk_equal("Rr = g(t=-1)*Ntot = 0.1*1000", X[["Rr"]], 100)

th_summary("initfun-build")
})
