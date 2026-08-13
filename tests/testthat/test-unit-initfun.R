test_that("test-unit-initfun", {
# ============================================================
# test-unit-initfun.R   (pure R; no Julia)
# .classify_functions() (state-independent vs state-dependent Functions) and
# .init_function_lines() (their `<name>_0 <- ...` definitions evaluated at the
# init time t=-1). This is the engine behind "initial states may reference
# state-independent Functions".
# ============================================================
th_load_pure(c("utils.R"))

state_symbols <- c("S_H", "I_H", "N1", "N_High", "total_pop")   # compartments + aliases
param_symbols <- c("Pi_CH", "total_switch", "n_L", "q14", "q15", "beta")

fns <- c(
  "ntilde <- n_L",                                  # indep: parameter only
  "sigma_HL <- (1-Pi_CH)*total_switch",             # indep: parameters
  "q <- if_else(time+startpoint<2015, q14, q15)",   # indep: params + time/startpoint
  "omega <- ntilde*N1",                             # DEP:   N1 is a state alias
  "foi <- beta*I_H/omega",                          # DEP:   I_H (state) + omega (dep)
  "bad <- nonexistent*2",                           # DEP:   unknown symbol (fail-safe)
  "chain <- sigma_HL + ntilde")                     # indep: refs independent Functions

cl <- .classify_functions(fns, state_symbols, param_symbols)

th_section("classification: independent set (transitive, order-preserving)")
chk("independent names", identical(names(cl$independent),
                                   c("ntilde", "sigma_HL", "q", "chain")))
chk("dependent names", setequal(cl$dependent, c("omega", "foi", "bad")))
chk("a state alias (N1) makes it dependent", "omega" %in% cl$dependent)
chk("depending on a dependent Function propagates", "foi" %in% cl$dependent)
chk("unknown symbol is fail-safe dependent", "bad" %in% cl$dependent)
chk("referencing only independent Functions stays independent",
    "chain" %in% names(cl$independent))

th_section("init lines: `_0` rewrite, time -> -1, startpoint literalised")
lines <- .init_function_lines(cl$independent, param_symbols, startpoint = 2015, cutoff = Inf)
chk("sigma_HL uses _0 params",
    any(grepl("sigma_HL_0 <- \\(1-Pi_CH_0\\)\\*total_switch_0", lines)))
chk("q substitutes time and startpoint literals",
    any(grepl("q_0 <- if_else\\(\\(-1\\)\\+2015<2015, q14_0, q15_0\\)", lines)))
chk("chain references the _0 aliases of other Functions",
    any(grepl("chain_0 <- sigma_HL_0 \\+ ntilde_0", lines)))

th_section("evaluating the lines at t=-1 gives the pre-window values")
e <- list2env(list(n_L_0 = 14, Pi_CH_0 = 0.525, total_switch_0 = 2,
                   q14_0 = 0.24, q15_0 = 0.28, beta_0 = 0.007))
for (ln in lines) eval(parse(text = ln), envir = e)
chk_equal("ntilde_0 = n_L", get("ntilde_0", e), 14)
chk_equal("sigma_HL_0 = (1-Pi)*switch", get("sigma_HL_0", e), 0.95)
chk_equal("q_0 picks q14 at t=-1 (year 2014)", get("q_0", e), 0.24)
chk_equal("chain_0 = sigma_HL + ntilde", get("chain_0", e), 14.95)

th_summary("initfun")
})
