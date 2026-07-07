test_that("test-unit-fillparams", {
# ============================================================
# test-unit-fillparams.R   (pure R; no Julia)
# Writing fitted estimates back into a parameter sheet: .fp_rewrite_cell and
# fill_params on a constructed MLE fit. xlsx round-trip (write + verify) is
# tested only when writexl + readxl are installed.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "fillParams.R"))

th_section(".fp_rewrite_cell: per-cell rewriting")
est <- list(beta = 0.5, X1 = 10)
chk("fitted cell -> *name=value", .fp_rewrite_cell("beta=[0,1]", est, NULL) == "*beta=0.5")
chk("fixed cell left unchanged",  .fp_rewrite_cell("*p=0.6", est, NULL) == "*p=0.6")
chk("unknown name left unchanged",.fp_rewrite_cell("q=[0,1]", est, NULL) == "q=[0,1]")
chk("NA left unchanged",          is.na(.fp_rewrite_cell(NA, est, NULL)))
chk("digits rounding applied",    .fp_rewrite_cell("beta=[0,1]", list(beta = 0.123456), 3) == "*beta=0.123")

th_section("fill_params: rewrite a whole sheet from an MLE fit")
fit <- structure(list(method = "lbfgsb",
                      point = list(parms = c(beta = 0.5),
                                   initial_state = c(X1 = 10))),
                 class = "compartmentalFit")
mp <- data.frame(
  States     = c("X1=[0,100]"),
  Parameters = c("beta=[0,1]"),
  `_Level1`  = c("1"),
  check.names = FALSE, stringsAsFactors = FALSE)
filled <- fill_params(fit, mp)
chk("state cell filled",  filled$States[1] == "*X1=10")
chk("param cell filled",  filled$Parameters[1] == "*beta=0.5")
chk("structural column untouched", filled$`_Level1`[1] == "1")
chk("column count preserved", ncol(filled) == ncol(mp))

th_section("write_filled_params + verify_filled round-trip (needs writexl+readxl)")
if (requireNamespace("writexl", quietly = TRUE) && requireNamespace("readxl", quietly = TRUE)) {
  d <- file.path(tempdir(), paste0("fp_", as.integer(runif(1, 1, 1e6)))); dir.create(d)
  p <- chk_ok("write_filled_params writes a file",
              write_filled_params(fit, d, modelParams = mp, out_file = "filled.xlsx"))
  chk("file exists", file.exists(file.path(d, "filled.xlsx")))
  chk_ok("verify_filled passes round-trip",
         verify_filled(file.path(d, "filled.xlsx"), mp))
} else {
  th_skip("xlsx round-trip", "install writexl + readxl to run")
}

th_summary("fillparams")
})
