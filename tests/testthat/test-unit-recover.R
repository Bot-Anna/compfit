test_that("test-unit-recover", {
# ============================================================
# test-unit-recover.R   (pure R; no Julia)
# .recover_solution() must accept fitted parameters in whatever order the fit
# returns them (the Julia backend uses sheet order; sap groups with-init params
# first), reordering by NAME to sap's canonical order -- and still reject a
# genuine name-set mismatch.
# ============================================================
th_load_pure(c("utils.R", "statesAndParams.R", "generateExpressions.R", "numberOfComps.R",
               "compartmentalFunction.R", "fitCompartmentalModel.R"))

# `b` (plain box) precedes `a` (`|init`) in the sheet, so sheet order is [b, a]
# while sap groups the with-init `a` first -> [a, b]: the interleave that used to
# trip the order guard.
mp <- data.frame(
  Others     = c("startpoint=2015", "endpoint=2018", "partition=4"),
  States     = c("*S=Ntot", "*I=0", ""),
  Parameters = c("*Ntot=1000", "b=[0,10]", "a=[0,10]|3"),
  check.names = FALSE, stringsAsFactors = FALSE)
m <- .build_model(mp, backend = "r")

th_section("sap groups with-init first; the fit returns sheet order")
chk("sap params_fitted = [a, b] (with-init first)",
    identical(names(m$sap$params_fitted), c("a", "b")))

th_section("recover reorders the fit's [b, a] to sap's [a, b] by name")
sol <- .recover_solution(c(b = 5, a = 3), m$sap, m$time_grid)   # fit/sheet order
chk_equal("a recovered by name (not position)", sol$parms[["a"]], 3)
chk_equal("b recovered by name", sol$parms[["b"]], 5)
chk_equal("fixed param carried through", sol$parms[["Ntot"]], 1000)
chk("fitted block is in sap's canonical order",
    identical(names(sol$parms)[1:2], c("a", "b")))

th_section("a genuine name-set mismatch is still rejected")
chk_error("unknown fitted-parameter name is rejected",
          .recover_solution(c(b = 5, z = 3), m$sap, m$time_grid))

th_summary("recover")
})
