test_that("test-unit-validate", {
# ============================================================
# test-unit-validate.R   (pure R; no Julia)
# validate_modelParams(): every shipped fixture passes (no false positives), and
# each class of bad entry raises a clear error naming the column/cell and the
# problem -- the checks that turn a cryptic builder failure into a useful message.
# ============================================================

th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "validate.R"))

mk <- function(params, states = c("*X1=990", "*X2=10"),
               q1 = c("0", "*2*-beta"), linear1 = c("0", "0"),
               others = c("startpoint=2000", "endpoint=2005", "partition=4"),
               functions = c("", "")) {
  cols <- list(`_Level1` = c("1", "2"), Others = others, States = states,
               Functions = functions, Parameters = params, Conditions = "",
               Linear1 = linear1, Quadratic1 = q1,
               Linear2 = c("0", "-gamma"), Quadratic2 = c("0", "0"))
  nr <- max(lengths(cols))
  cols <- lapply(cols, function(x) c(as.character(x), rep("", nr - length(x))))
  do.call(data.frame, c(cols, check.names = FALSE, stringsAsFactors = FALSE))
}
errmsg <- function(mp) tryCatch({ validate_modelParams(mp); NA_character_ },
                                error = function(e) conditionMessage(e))

th_section("valid sheets pass (no false positives)")
chk_ok("a clean sheet validates", validate_modelParams(mk(c("beta=[0,2]", "gamma=[0,1]"))))
# param-dependent init + a time-varying function must NOT be flagged
ok2 <- mk(c("beta=[0,2]", "gamma=[0,1]", "*N0=1000", "*init_inf=0.01", "*ramp=0.05"),
          states = c("*X1=N0_0*(1-init_inf_0)", "*X2=init_inf_0*N0_0"),
          functions = c("gamma_t<-gamma*(1+ramp*time)", ""))
ok2$Linear2 <- c("0", "-(gamma_t)", rep("", nrow(ok2) - 2))
chk_ok("param-dependent init + time function validate", validate_modelParams(ok2))
chk_ok("every distribution prior validates",
       validate_modelParams(mk(c("beta=StudentT(4,1,0.4)[0,3]", "gamma=Normal(0.4,0.1)[0,1]",
                                 "a=LogNormal(-1,0.5)", "b=Beta(2,2)", "c=Gamma(2,3)"))))

th_section("non-numeric values are caught")
chk("fixed param non-numeric",  grepl("abc", errmsg(mk(c("*beta=abc", "gamma=[0,1]")))))
chk("fixed state non-numeric",  grepl("States cell", errmsg(mk(c("beta=[0,1]", "gamma=[0,1]"),
                                                               states = c("*X1=foo", "*X2=10")))))
chk("text inside a box",        grepl("two numbers", errmsg(mk(c("beta=[0,foo]", "gamma=[0,1]")))))
chk("non-numeric prior arg",    grepl("non-numeric argument", errmsg(mk(c("beta=Normal(1,foo)", "gamma=[0,1]")))))
chk("non-numeric partition",    grepl("partition.*numeric", errmsg(mk(c("beta=[0,1]", "gamma=[0,1]"),
                                                                       others = c("startpoint=2000", "endpoint=2005", "partition=abc")))))

th_section("structural / symbol mistakes are caught")
chk("reversed box lo>=hi",      grepl("lower < upper", errmsg(mk(c("beta=[5,3]", "gamma=[0,1]")))))
chk("undefined symbol in coeff", grepl("bta", errmsg(mk(c("beta=[0,1]", "gamma=[0,1]"), q1 = c("0", "*2*-bta")))))
chk("Quadratic target out of range",
    grepl("out of range", errmsg(mk(c("beta=[0,1]", "gamma=[0,1]"), q1 = c("0", "*9*-beta")))))
chk("missing endpoint",         grepl("missing required .endpoint", errmsg(mk(c("beta=[0,1]", "gamma=[0,1]"),
                                                                               others = c("startpoint=2000", "partition=4")))))
# a coefficient references 'beta' but it is not declared -- and 'beta' is a base
# R function, so exists('beta') is TRUE; it must STILL be flagged (not masked).
chk("undeclared base-fn-named symbol ('beta') is caught",
    grepl("'beta'.*not a declared", errmsg(mk(c("gamma=[0,1]", "*m=0.1"),
                                              q1 = c("0", "0"),
                                              linear1 = c("-beta", "0")))))

th_section("reserved names + column groups")
chk("Julia keyword 'end' is a reserved name",
    grepl("reserved", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "end=[0,1]")))))
chk("codegen variable 't' is a reserved name",
    grepl("reserved", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "t=[0,1]")))))
# Codegen-generated variables/functions: level pops, sums, per-term helpers.
for (nm in c("N1=[0,1]", "N2=[0,1]", "total_pop=[0,1]", "f12=[0,1]", "g23=[0,1]",
             "cst1=[0,1]", "secOrd_1_2=[0,1]", "time=[0,1]"))
  chk(paste("collision name", sub("=.*", "", nm), "rejected"),
      grepl("reserved", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", nm)))))
# But close look-alikes that are NOT generated must still pass.
chk("N0 (fixed initial pop) is allowed",
    is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "*N0=1000")))))
chk("single-digit f1 / g2 are allowed",
    is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "f1=[0,1]", "g2=[0,1]")))))
# A STATE named after a reserved word is now caught too (was a gap).
chk("reserved word as a STATE name is caught",
    grepl("reserved", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]"),
                                states = c("*end=990", "*X2=10")))))
# Declaring both `beta` and its `_0` alias clashes.
chk("param clashing with a declared quantity's _0 alias is caught",
    grepl("_0", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "beta_0=[0,1]")))))
chk("a lone _0-suffixed name (no matching base) is fine",
    is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "kappa_0=[0,1]")))))
chk("ordinary names are fine", is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]")))))
chk("missing 'States' column is an error", {
  m <- mk(c("beta=[0,2]", "gamma=[0,1]")); m$States <- NULL
  grepl("missing 'States'", errmsg(m)) })
chk("unrecognised column warns", {
  w <- character(0); m <- mk(c("beta=[0,2]", "gamma=[0,1]")); m$Prameters <- ""
  withCallingHandlers(validate_modelParams(m),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  any(grepl("unrecognised", w)) })

th_section("all shipped fixtures pass")
for (nm in c("minimal", "medium", "SI", "SIS", "SIR", "SEIR", "SIR_priors", "SEIR_priors")) {
  dir <- fixture_dir(nm)
  if (!dir.exists(dir)) next
  sc <- load_scenario(dir, combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
  chk(paste(nm, "fixture validates"), isTRUE(validate_modelParams(sc$modelParams)))
}

th_summary("unit-validate")
})
