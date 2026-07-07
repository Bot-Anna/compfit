test_that("test-unit-utils", {
# ============================================================
# test-unit-utils.R   (pure R; no Julia, no data)
# Grammar parsers and numeric helpers in R/utils.R.
# ============================================================
th_load_pure(c("utils.R"))

th_section("parse_data_cell: cell grammar")
o <- parse_data_cell("142");   chk("observed value",      o$kind == "observed" && o$value == 142)
o <- parse_data_cell("0");     chk("observed zero",        o$kind == "observed" && o$value == 0)
o <- parse_data_cell("x");     chk("'x' -> missing",       o$kind == "missing")
o <- parse_data_cell("");      chk("blank -> missing",     o$kind == "missing")
o <- parse_data_cell(NA);      chk("NA -> missing",        o$kind == "missing")
o <- parse_data_cell("<5");    chk("'<5' left strict",     o$kind=="censored" && o$bound=="upper" && o$inclusive==FALSE && o$limit==5)
o <- parse_data_cell("<=5");   chk("'<=5' left inclusive", o$kind=="censored" && o$bound=="upper" && o$inclusive==TRUE)
o <- parse_data_cell(">3");    chk("'>3' right strict",    o$kind=="censored" && o$bound=="lower" && o$inclusive==FALSE && o$limit==3)
o <- parse_data_cell(">=3");   chk("'>=3' right inclusive",o$kind=="censored" && o$bound=="lower" && o$inclusive==TRUE)
o <- parse_data_cell(" >= 3 ");chk("whitespace tolerated", o$kind=="censored" && o$bound=="lower" && o$limit==3)
chk_error("unparseable cell errors", parse_data_cell("abc"))

th_section("parse_data_cell: interval + asymmetric grammar")
o <- parse_data_cell("[3,7]");    chk("interval [3,7]",      o$kind=="interval" && o$limit==3 && o$upper==7)
o <- parse_data_cell("120->124"); chk("asym A->B upward",    o$kind=="asym" && o$value==120 && o$dev==4  && o$dir==1)
o <- parse_data_cell("120->110"); chk("asym A->B downward",  o$kind=="asym" && o$value==120 && o$dev==10 && o$dir==-1)
o <- parse_data_cell("120+");     chk("asym A+ (global up)",  o$kind=="asym" && o$value==120 && is.na(o$dev) && o$dir==1)
o <- parse_data_cell("120-");     chk("asym A- (global down)",o$kind=="asym" && o$value==120 && is.na(o$dev) && o$dir==-1)
chk_error("interval A>=B errors", parse_data_cell("[7,3]"))
o <- parse_data_cell("[110,130]~10");   chk("soft interval ~s",     o$kind=="interval" && o$limit==110 && o$upper==130 && o$dev==10 && o$dev2==10)
o <- parse_data_cell("[110,130]~10,25"); chk("soft interval ~sl,su", o$kind=="interval" && o$dev==10 && o$dev2==25)
o <- parse_data_cell("[110,130]");      chk("hard interval has no shoulders", o$kind=="interval" && is.na(o$dev) && is.na(o$dev2))
chk_error("soft shoulder must be positive", parse_data_cell("[1,2]~-3"))
chk_error("soft shoulder must be numeric",  parse_data_cell("[1,2]~x"))

th_section("parseLikelihood: families + dispersion flag")
chk("blank -> gaussian default", parseLikelihood(NA)$family == "gaussian" && parseLikelihood(NA)$dispersion == FALSE)
chk("'pois' alias",     parseLikelihood("pois")$family == "poisson")
chk("'nb' -> negbin + dispersion", { s <- parseLikelihood("nb"); s$family=="negbin" && s$dispersion })
chk("'bb' -> betabinom + dispersion", { s <- parseLikelihood("bb"); s$family=="betabinom" && s$dispersion })
chk("binomial no dispersion", parseLikelihood("binomial")$dispersion == FALSE)
chk_error("unknown family errors", parseLikelihood("weibull"))
chk("'asym=' deviation parsed", { s <- parseLikelihood("gaussian; asym=4"); s$family=="gaussian" && s$asym_dev==4 })
chk("no asym -> NA",            is.na(parseLikelihood("poisson")$asym_dev))
chk_error("negative asym errors", parseLikelihood("gaussian; asym=-1"))

th_section("parsePrior: prior grammar")
p <- parsePrior("[0,1]");          chk("box -> Uniform", p$kind=="estimated" && p$dist=="Uniform" && p$lower==0 && p$upper==1)
p <- parsePrior("Normal(2,0.5)");  chk("Normal unbounded", p$dist=="Normal" && p$args[1]==2 && is.infinite(p$lower))
p <- parsePrior("Normal(2,0.5)[0,5]"); chk("Normal truncated", p$dist=="Normal" && p$lower==0 && p$upper==5)
p <- parsePrior("3.5");            chk("plain number -> fixed", p$kind=="fixed" && p$value==3.5)
p <- parsePrior("StudentT(3,0.5,0.2)"); chk("StudentT parsed", p$dist=="StudentT" && all(p$args==c(3,0.5,0.2)) && p$start==0.5)
p <- parsePrior("StudentT(4,1,0.3)[0,2]"); chk("StudentT truncated", p$dist=="StudentT" && p$lower==0 && p$upper==2)
chk_error("StudentT needs 3 args",  parsePrior("StudentT(1,2)"))
chk_error("StudentT needs nu>0",    parsePrior("StudentT(-1,0,1)"))
chk_error("StudentT needs sigma>0", parsePrior("StudentT(3,0,-1)"))
chk_equal("StudentT MLE box = mu +/- 6 sigma",
          compfit:::.mle_box(parsePrior("StudentT(3,0.5,0.2)")), c(-0.7, 1.7), tol = 1e-8)
chk("StudentT R-backend log-density finite",
    is.finite(compfit:::.rprior_logd(parsePrior("StudentT(3,0.5,0.2)"), 0.6)))

th_section(".is_prior_entry / .mle_box: prior routing + finite L-BFGS-B box")
chk("box is a prior",           compfit:::.is_prior_entry("beta=[0,2]"))
chk("Normal is a prior",        compfit:::.is_prior_entry("beta=Normal(0.5,0.2)"))
chk("truncated Normal is a prior", compfit:::.is_prior_entry("beta=Normal(0.5,0.2)[0,2]"))
chk("expression is NOT a prior", !compfit:::.is_prior_entry("beta=alpha*gamma"))
chk("bare constant is NOT a prior (fixed)", !compfit:::.is_prior_entry("beta=0.5"))
chk_equal("box box = [lo,hi]",         compfit:::.mle_box(parsePrior("[0,2]")),        c(0, 2))
chk_equal("Normal box = mean +/- 4sd", compfit:::.mle_box(parsePrior("Normal(0.5,0.2)")), c(-0.3, 1.3), tol = 1e-8)
chk_equal("truncation honoured",       compfit:::.mle_box(parsePrior("Normal(0.5,0.2)[0,2]")), c(0, 2))
chk("Beta box = [0,1]",  isTRUE(all.equal(compfit:::.mle_box(parsePrior("Beta(2,2)")), c(0, 1))))
chk("Gamma box floored at 0", compfit:::.mle_box(parsePrior("Gamma(2,3)"))[1] == 0)

th_section(".check_reserved_bayes_names: sigma / phi are reserved for Bayes")
chk_error("'sigma' is reserved",  compfit:::.check_reserved_bayes_names(c("beta", "sigma")))
chk_error("'phi' is reserved",    compfit:::.check_reserved_bayes_names("phi"))
chk_error("'phi2' is reserved",   compfit:::.check_reserved_bayes_names(c("gamma", "phi2")))
chk_ok("'sigma_ei' is allowed",   compfit:::.check_reserved_bayes_names(c("beta", "sigma_ei", "gamma")))
chk_ok("'sigmoid'/'alpha' allowed", compfit:::.check_reserved_bayes_names(c("sigmoid", "alpha", "phix")))

th_section("statesAndParams: distributional prior becomes a fitted param (not a function)")
mp_np <- data.frame(
  `_Level1` = c("1","2"), Others = c("startpoint=2000","endpoint=2005"),
  States = c("*X1=990","*X2=10"), Functions = c("",""),
  Parameters = c("beta=Normal(0.5,0.2)[0,2]","gamma=[0,1]"), Conditions = c("",""),
  Linear1=c("0","0"), Quadratic1=c("0","*2*-beta"),
  Linear2=c("0","-gamma"), Quadratic2=c("0","0"),
  check.names = FALSE, stringsAsFactors = FALSE)
sap_np <- statesAndParams(mp_np)
chk("beta is a fitted parameter",   "beta" %in% names(sap_np$params_fitted))
chk("beta box honours truncation",  sap_np$lower_params[["beta"]] == 0 && sap_np$upper_params[["beta"]] == 2)
chk("beta is NOT a param-function", !("beta" %in% names(sap_np$params_functions)))
# States share the same parser: a distributional STATE prior is a fitted state,
# not a mis-parsed state-function (the latent bug the shared .parse_entry_group fixes).
mp_st <- mp_np; mp_st$States <- c("X1=Beta(2,2)", "*X2=10")
sap_st <- statesAndParams(mp_st)
chk("X1 fitted state from a distribution prior", "X1" %in% names(sap_st$states_fitted))
chk("X1 box from Beta(2,2) = [0,1]", sap_st$lower_states[["X1"]] == 0 && sap_st$upper_states[["X1"]] == 1)
chk("X1 NOT a state-function", !("X1" %in% names(sap_st$states_functions)))
# the shared parser classifies entries the same for states and params
grp <- compfit:::.parse_entry_group(c("*a=1", "b=[0,2]", "c=Normal(0,1)", "d=e*2"))
chk(".parse_entry_group: fixed/prior/function split",
    identical(names(grp$fixed), "a") && all(c("b", "c") %in% grp$without_names) &&
      "d" %in% names(grp$functions))

th_section("family registry constants")
chk(".DISCRETE_FAMILIES", identical(.DISCRETE_FAMILIES, c("poisson","negbin")))
chk(".NEEDS_PHI_FAMILIES has betabinom", "betabinom" %in% .NEEDS_PHI_FAMILIES)
chk(".ALL_FAMILIES has 6", length(.ALL_FAMILIES) == 6)

th_section("numeric helpers")
chk_equal("normalise midpoint", normalise(5, 0, 10), 0.5)
chk_equal("denormalise inverse", denormalise(normalise(7.3, 2, 9), 2, 9), 7.3)
chk("is_cumulative_stream TRUE",  is_cumulative_stream("cumulative(g*X1)"))
chk("is_cumulative_stream FALSE", !is_cumulative_stream("annual(g*X1)"))
chk_equal("extract_param startpoint", extract_param(c("startpoint=2013","endpoint=2020"), "startpoint"), 2013)
chk_equal("extract_numbers [0,1]", extract_numbers("[0,1]"), c(0,1))
chk_equal("snap_to_step", snap_to_step(0,10,2.3,1), 2)
chk("remove_trailing_plus", remove_trailing_plus("a+b+ ") == "a+b")
chk("reduce_expression collapses signs", reduce_expression("a--b") == "a+b")

th_summary("utils")
})
