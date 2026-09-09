test_that("test-unit-asym-column", {
# ============================================================
# test-unit-asym-column.R   (pure R; no Julia)
# The dedicated `Asym` column of dataCombined: it supplies the stream-global
# deviation for A+/A- cells, takes precedence over the deprecated
# `Likelihood: ...; asym=` clause, is validated, and -- critically -- is treated
# as a META column rather than as a year of observations.
# ============================================================
th_load_pure(c("utils.R", "fitCompartmentalModel.R"))

mk <- function(..., formula = "X3") {
  extra <- list(...)
  dc <- data.frame(Label = "A", Formula = formula,
                   check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(extra)) dc[[nm]] <- extra[[nm]]
  dc[["2013"]] <- "10"
  dc[["2014"]] <- "20-"      # A- : soft downward, needs a deviation
  dc[["2015"]] <- "30"
  dc
}
tg <- list(startpoint = 2013, endpoint = 2015)

th_section("Asym column supplies the deviation")
dat <- .prepare_data(mk(Asym = 3), tg)
chk("the A- cell is flagged asymmetric", dat$asym_mask[2, 1] == 1)
chk_equal("anchor recorded", dat$asym_val_mat[2, 1], 20)
chk_equal("direction is soft-downward", dat$asym_dir_mat[2, 1], -1)
chk_equal("deviation taken from Asym", dat$asym_dev_mat[2, 1], 3)

th_section("Asym is metadata, not a year of data")
# Three year columns are declared and the horizon is three years; if `Asym` were
# mistaken for a data column the year-count guard would fire.
chk("three observation years survive", nrow(dat$obs_mask) == 3)
chk("one stream", ncol(dat$obs_mask) == 1)
chk("plain cells still observed", dat$obs_mask[1, 1] == 1 && dat$obs_mask[3, 1] == 1)

th_section("precedence: Asym beats the deprecated Likelihood clause")
dat2 <- .prepare_data(mk(Likelihood = "negbin; asym=25", Asym = 4), tg)
chk_equal("Asym wins over asym=", dat2$asym_dev_mat[2, 1], 4)

th_section("backward compatibility: the Likelihood clause still works")
dat3 <- .prepare_data(mk(Likelihood = "negbin; asym=25"), tg)
chk_equal("asym= still honoured when Asym absent", dat3$asym_dev_mat[2, 1], 25)

th_section("a blank Asym cell falls back to the clause")
dat4 <- .prepare_data(mk(Likelihood = "negbin; asym=7", Asym = NA), tg)
chk_equal("blank Asym -> clause used", dat4$asym_dev_mat[2, 1], 7)

th_section("per-cell A->B still overrides both")
dc5 <- mk(Likelihood = "negbin; asym=25", Asym = 4)
dc5[["2014"]] <- "20->17"                    # dev = |17-20| = 3, dir = -1
dat5 <- .prepare_data(dc5, tg)
chk_equal("per-cell deviation wins", dat5$asym_dev_mat[2, 1], 3)

th_section("validation")
chk("non-positive Asym is rejected",
    inherits(try(.prepare_data(mk(Asym = 0), tg), silent = TRUE), "try-error"))
chk("negative Asym is rejected",
    inherits(try(.prepare_data(mk(Asym = -2), tg), silent = TRUE), "try-error"))
chk("A- with no deviation anywhere still errors",
    inherits(try(.prepare_data(mk(), tg), silent = TRUE), "try-error"))

th_section("symmetric case: Asym = 1 gives equal slopes")
dat6 <- .prepare_data(mk(Asym = 1), tg)
chk_equal("dev 1 recorded", dat6$asym_dev_mat[2, 1], 1)

th_summary("test-unit-asym-column")
})
