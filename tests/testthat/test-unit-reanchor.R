test_that("test-unit-reanchor", {
# ============================================================
# test-unit-reanchor.R   (pure R; no Julia)
# reanchor_scenario(): move a fit's startpoint EARLIER, freeze its parameter
# estimates, respecify the initial state at the new start year.
#
# Re-anchoring rather than backward integration: the model is always solved
# FORWARD, so the earlier initial state is a supplied (or estimated) assumption.
# The new years are padded with BLANK data columns, which the data-cell grammar
# reads as missing -- no data is invented and nothing enters the loss.
# ============================================================
skip_if_not_installed("deSolve")

sc <- load_scenario(fixture_dir("medium"), combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
fit <- suppressMessages(fitCompartmentalModel(sc$modelParams, sc$dataCombined,
         method = "lbfgsb", solver = solver_control(backend = "r")))
chk("MLE fit produced", isTRUE(fit$success))
est <- get_central_point(fit)$parms

th_section("argument checks")
chk_error("a later startpoint is refused",  reanchor_scenario(fit, 2018))
chk_error("the same startpoint is refused", reanchor_scenario(fit, 2015))
chk_error("a fractional year is refused",   reanchor_scenario(fit, 2010.5))
chk_error("a vector is refused",            reanchor_scenario(fit, c(2010, 2011)))
chk("the refusal points at endpoint= for the other direction",
    grepl("endpoint", tryCatch(reanchor_scenario(fit, 2018),
                               error = function(e) conditionMessage(e)), fixed = TRUE))
chk_error("an unknown compartment in states= is refused",
          reanchor_scenario(fit, 2010, states = c(NOPE = "[0,1]")))
chk_error("an unnamed states= is refused",
          reanchor_scenario(fit, 2010, states = "[0,1]"))

th_section("inheriting the initial state warns (it describes the OLD start year)")
w <- NULL
re0 <- withCallingHandlers(reanchor_scenario(fit, 2010),
         warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
chk("warns, naming both years", grepl("2015", w %||% "") && grepl("2010", w %||% ""))

th_section("the rewritten sheet: start year, frozen parameters, new states")
re <- reanchor_scenario(fit, 2010,
        states = c(X1 = "[900,1000]", X2 = "[0,50]", X3 = "0"))
tg <- compfit:::.time_grid(re$modelParams)
chk_equal("startpoint moved", tg$startpoint, 2010)
chk_equal("endpoint untouched", tg$endpoint,
          compfit:::.time_grid(sc$modelParams)$endpoint)
chk_equal("partition untouched", tg$partition,
          compfit:::.time_grid(sc$modelParams)$partition)
chk("other Others entries survive (cutoff)",
    any(grepl("cutoff=", re$modelParams$Others)))
pcells <- gsub(" ", "", stats::na.omit(re$modelParams$Parameters))
chk("every parameter is frozen with the '*' fixed marker",
    all(grepl("^\\*", pcells)))
chk("beta frozen at its estimate",
    any(grepl(sprintf("^\\*beta=%s$", format(est[["beta"]])), pcells)))
scells <- gsub(" ", "", stats::na.omit(re$modelParams$States))
chk("a box state is fitted (no '*')",  "X1=[900,1000]" %in% scells)
chk("a bare number state is fixed",    "X3=0" %in% scells)
chk("validate_modelParams accepts the result",
    is.null(tryCatch({ validate_modelParams(re$modelParams); NULL },
                     error = function(e) e)))

th_section("data padded with BLANK columns for the new years")
meta <- intersect(c("Label","Formula","Likelihood","Weight","Average","Asym"),
                  names(re$dataCombined))
yrs <- setdiff(names(re$dataCombined), meta)
chk_equal("one column per year of the new grid", yrs, as.character(2010:2020))
chk("the new years are blank",
    all(vapply(as.character(2010:2014),
               function(y) all(is.na(re$dataCombined[[y]])), logical(1))))
chk("the original years keep their values",
    identical(as.character(re$dataCombined[["2015"]]),
              as.character(sc$dataCombined[["2015"]])))
d <- compfit:::.prepare_data(re$dataCombined, tg)
chk_equal("only the original cells are observed", sum(d$obs_mask == 1),
          sum(compfit:::.prepare_data(sc$dataCombined,
                compfit:::.time_grid(sc$modelParams))$obs_mask == 1))

th_section("the re-anchored sheet refits, forward, from the earlier year")
f2 <- suppressMessages(fitCompartmentalModel(re$modelParams, re$dataCombined,
        method = "lbfgsb", solver = solver_control(backend = "r")))
chk("refit succeeded", isTRUE(f2$success))
p2 <- get_central_point(f2)
chk_equal("frozen beta is unchanged by the refit", p2$parms[["beta"]], est[["beta"]],
          tol = 1e-8)
chk("the earlier initial state was estimated", "X1" %in% names(p2$initial_state))
ev <- solve_and_evaluate(f2, p2$initial_state, p2$parms)$evaluation
chk_equal("the solve spans the longer grid", nrow(ev), tg$partition *
          (tg$endpoint - tg$startpoint + 1) + 1)
chk("trajectory stays non-negative (forward integration, unlike a backward solve)",
    all(vapply(setdiff(names(ev), "date"),
               function(c0) all(ev[[c0]] >= -1e-8, na.rm = TRUE), logical(1))))
})
