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
  cols <- list(`Level_1` = c("1", "2"), Others = others, States = states,
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
    grepl("not a valid compartment", errmsg(mk(c("beta=[0,1]", "gamma=[0,1]"), q1 = c("0", "*9*-beta")))))
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
# Codegen-generated variables: level pops, weighted pools, sums, term temporaries.
for (nm in c("N1=[0,1]", "N2=[0,1]", "Nw1=[0,1]", "Nw2=[0,1]", "total_pop=[0,1]",
             "secOrd_1_2=[0,1]", "time=[0,1]"))
  chk(paste("collision name", sub("=.*", "", nm), "rejected"),
      grepl("reserved", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", nm)))))
# But close look-alikes that are NOT generated must still pass.
chk("N0 (fixed initial pop) is allowed",
    is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "*N0=1000")))))
# The old inline-$ helper functions f<ij>/g<ij>/cst<i> are no longer generated
# (the $-function feature was removed), so those names are free again.
chk("former $-helper names f12 / g23 / cst1 are allowed",
    is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]", "f12=[0,1]", "g23=[0,1]", "cst1=[0,1]")))))
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

th_section("removed inline $-functions are rejected with a migration hint")
chk("$ in a Quadratic coefficient is rejected",
    grepl("\\$.*removed|removed.*Functions", errmsg(mk(c("beta=[0,2]", "gamma=[0,1]"),
                                                        q1 = c("0", "*2*$beta*exp(-0.1*time)")))))
chk("$ in a Constant cell is rejected", {
  m <- mk(c("beta=[0,2]", "gamma=[0,1]")); m$Constant <- c("$0.3*time", "0", rep("", nrow(m) - 2))
  grepl("removed", errmsg(m)) })

th_section("Mixing_<level> columns (full-pool denominator weights)")
add_mix <- function(mp, col, vals) { mp[[col]] <- c(as.character(vals), rep("", nrow(mp) - length(vals))); mp }
chk("valid numeric Mixing column passes",
    is.na(errmsg(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Mixing_1", c("1", "0.3")))))
chk("Mixing weight referencing a declared parameter passes",
    is.na(errmsg(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]", "w=[0,1]")), "Mixing_1", c("1", "w")))))
chk("negative Mixing weight is rejected",
    grepl("must be >= 0", errmsg(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Mixing_1", c("1", "-0.5")))))
chk("undeclared symbol in a Mixing weight is caught",
    grepl("not a declared", errmsg(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Mixing_1", c("1", "wat")))))

# A Mixing_<name> column REQUIRES a Level_<name> column with the same <name>
# (the name may be a number or a string). The fixture's level is Level_1.
chk("numeric-id Mixing_1 matches Level_1 -> passes",
    is.na(errmsg(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Mixing_1", c("1", "1")))))
chk("numeric-id Mixing_2 without Level_2 is rejected",
    grepl("no matching level column", errmsg(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Mixing_2", c("1", "1")))))
# String-named level: rename the fixture level to Level_grpA.
mk_named <- function(params, mixcol, vals) {
  m <- mk(params); m[["Level_grpA"]] <- m[["Level_1"]]; m[["Level_1"]] <- NULL
  add_mix(m, mixcol, vals)
}
chk("string-id Mixing_grpA matches Level_grpA -> passes",
    is.na(errmsg(mk_named(c("beta=[0,2]", "gamma=[0,1]"), "Mixing_grpA", c("1", "0.4")))))
chk("string-id Mixing_grpB without Level_grpB is rejected",
    grepl("no matching level column", errmsg(mk_named(c("beta=[0,2]", "gamma=[0,1]"), "Mixing_grpB", c("1", "1")))))

th_section("Pool_<level> (denominator as a function of the level populations)")
perr <- function(mp) suppressMessages(errmsg(mp))   # Pool columns emit an info note
chk("Pool expression referencing the level head counts passes",
    is.na(perr(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Pool_1", c("N1 + 0.5*total_pop")))))
chk("undeclared symbol in a Pool expression is caught",
    grepl("not a declared", perr(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Pool_1", c("N1 + wat")))))
chk("Pool without a matching level column is rejected",
    grepl("no matching level column", perr(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Pool_Ghost", c("N1")))))
chk("Pool and Mixing on the same level is rejected", {
  m <- add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Pool_1", c("N1"))
  m <- add_mix(m, "Mixing_1", c("1", "1"))
  grepl("both a", perr(m)) })
chk("a Pool column emits the divide-by-zero note", {
  msgs <- character(0)
  withCallingHandlers(
    validate_modelParams(add_mix(mk(c("beta=[0,2]", "gamma=[0,1]")), "Pool_1", c("N1"))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  any(grepl("floored at a small positive value", msgs)) })

th_section("a parameter/function may not share a State name (state-slot rewrite)")
chk("parameter named after a State is rejected",
    grepl("both as a State", errmsg(mk(c("X1=[0,1]", "gamma=[0,1]")))))   # X1 is a compartment name
# The level is `Level_1`, so its head-count alias is N_1; a param may not take it.
chk("parameter named after a level alias N_<level> is rejected",
    grepl("head-count alias", errmsg(mk(c("N_1=[0,1]", "gamma=[0,1]")))))
chk("a Function may reference the level alias N_<level>",
    is.na(errmsg({ m <- mk(c("beta=[0,2]", "gamma=[0,1]"))
                   m$Functions[1] <- "frac<-N_1/total_pop"; m })))

th_section("operator characters in a name are rejected (legal: letters/digits/_)")
# Names must be legal identifiers on EVERY backend; operator characters are parsed
# as operators, never as part of a name, so `f^R` reads as `f ^ R`. (Learned the
# hard way aligning the two-country HIV model: `f^R_CH`, `Phiw_(H,CH)`, etc.)
for (bad in c("f^R", "Phiw_(H,CH)", "p_H,H", "beta-x", "rate/yr", "a b")) {
  msg <- errmsg(mk(c(paste0(bad, "=[0,1]"), "beta=[0,1]", "gamma=[0,1]")))
  chk(paste("illegal name", sQuote(bad), "rejected"),
      grepl("not valid names|parsed as operators", msg))
}
chk("a legal underscore name (f_R_CH, like the fixed HIV model) is accepted",
    is.na(errmsg(mk(c("f_R_CH=[0,1]", "beta=[0,1]", "gamma=[0,1]")))))
# The superscript-style names the HIV sheet uses are fine once written with '_'.
chk("underscore superscript names (S_HighCH_nP style) are accepted",
    is.na(errmsg(mk(c("beta=[0,2]", "gamma=[0,1]"),
                    states = c("*S_HighCH_nP=990", "*S_HighCH_P=10")))))
# An operator character used inside an EXPRESSION (not a declaration) is caught too,
# via the existing unknown-symbol check (`beta^R_CH` -> free var `R_CH`).
chk("a caret inside a coefficient expression is caught", {
  m <- mk(c("beta=[0,2]", "gamma=[0,1]"))
  m$Linear1 <- c("beta^R_CH", "0", rep("", nrow(m) - 2))
  !is.na(errmsg(m)) })

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
