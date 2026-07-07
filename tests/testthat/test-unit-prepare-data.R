test_that("test-unit-prepare-data", {
# ============================================================
# test-unit-prepare-data.R   (pure R; no Julia)
# .prepare_data() -- the data contract feeding the loss and the Bayes path:
# censoring matrices, cumulative differencing, Weight/Average coercion, and the
# column-count-vs-horizon guard.
# ============================================================
th_load_pure(c("utils.R", "fitCompartmentalModel.R"))

# Three streams x three years (2013..2015):
#   X3              : observed / <5 (left strict) / observed
#   X4              : observed / >=10 (right inclusive) / observed
#   cumulative(X1)  : fully observed -> differenced to annual increments
dc <- data.frame(Label = c("A", "B", "C"),
                 Formula = c("X3", "X4", "cumulative(X1)"),
                 Weight = c("1", "0", "2"),     # "0" stored as text -> coercion
                 Average = c(1, 1, 1),
                 check.names = FALSE, stringsAsFactors = FALSE)
dc[["2013"]] <- c("10", "8",  "5")
dc[["2014"]] <- c("<5", ">=10", "10")
dc[["2015"]] <- c("20", "12", "18")
tg <- list(startpoint = 2013, endpoint = 2015)

dat <- .prepare_data(dc, tg)

th_section("orientation + identity")
chk("matrices are years x streams", nrow(dat$obs_mask) == 3 && ncol(dat$obs_mask) == 3)
chk("names_data_points == formulas", identical(dat$names_data_points, dc$Formula))
chk("data_points has 3 dated rows", nrow(dat$data_points) == 3 && !is.null(dat$data_points$date))

th_section("left-censoring (<5 on stream 1, year 2)")
chk("exactly one left-censored cell", sum(dat$cens_mask) == 1)
chk("cens at [year2, stream1]", dat$cens_mask[2, 1] == 1)
chk_equal("limit recorded", dat$limit_mat[2, 1], 5)
chk("strict -> inc_mask 0", dat$inc_mask[2, 1] == 0)

th_section("right-censoring (>=10 on stream 2, year 2)")
chk("exactly one right-censored cell", sum(dat$lcens_mask) == 1)
chk("rcens at [year2, stream2]", dat$lcens_mask[2, 2] == 1)
chk_equal("lower limit recorded", dat$llimit_mat[2, 2], 10)
chk("inclusive -> linc_mask 1", dat$linc_mask[2, 2] == 1)

th_section("observation mask")
chk("censored cells are not 'observed'", dat$obs_mask[2, 1] == 0 && dat$obs_mask[2, 2] == 0)
chk("observed cells flagged", dat$obs_mask[1, 1] == 1 && dat$obs_mask[3, 3] == 1)

th_section("cumulative differencing")
chk("cumulative_cols identifies stream 3", identical(as.integer(dat$cumulative_cols), 3L))
chk_equal("annual increments diff(c(0,5,10,18))", dat$matrix_data_points[, 3], c(5, 5, 8))
chk("non-cumulative stream keeps observed value", dat$matrix_data_points[1, 1] == 10)

th_section("Weight/Average numeric coercion + matrices")
chk("weight_matrix is 3x3", all(dim(dat$weight_matrix) == c(3, 3)))
chk("text '0' weight coerced (no NA)", !any(is.na(dat$weight_matrix)))
chk_equal("weight diag = sqrt(c(1,0,2))", diag(dat$weight_matrix), sqrt(c(1, 0, 2)))

th_section("blank/non-numeric Weight warns and defaults to 1 (not dropped)")
dc2 <- dc; dc2$Weight <- c("abc", "1", "1")     # "abc" -> NA -> default 1 (warned)
warned <- FALSE
d2 <- withCallingHandlers(.prepare_data(dc2, tg),
                          warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
chk("warning emitted for non-numeric Weight", warned)
chk("blank/unparseable weight defaults to 1", d2$weight_matrix[1, 1] == 1)

th_section("negative Weight is rejected")
dc4 <- dc; dc4$Weight <- c("1", "-2", "1")
chk("negative weight errors",
    isTRUE(tryCatch({ .prepare_data(dc4, tg); FALSE },
                    error = function(e) grepl("non-negative", conditionMessage(e)))))

th_section("blank Average cell falls back to the auto-computed scale (warned)")
dc5 <- dc; if (!"Average" %in% names(dc5)) dc5$Average <- rep(1, nrow(dc5))
dc5$Average[1] <- NA                              # blank one cell of a present column
awarned <- FALSE
d5 <- withCallingHandlers(.prepare_data(dc5, tg),
                          warning = function(w) { awarned <<- TRUE; invokeRestart("muffleWarning") })
chk("warning emitted for blank Average cell", awarned)
chk("blank Average not zeroed (auto-computed instead)", d5$average_matrix[1, 1] != 0)

th_section("Weight/Average are optional (auto-computed when absent)")
dc3 <- dc; dc3$Weight <- NULL; dc3$Average <- NULL
d3  <- .prepare_data(dc3, tg)
chk("builds without Weight/Average columns", all(dim(d3$average_matrix) == c(3, 3)))
chk("absent Weight defaults to 1", all(diag(d3$weight_matrix) == 1))
# Average auto = 1/mean(observed) per stream; stream 1 observed {10,.,20} -> mean 15
# (year-2 is '<5' censored -> NA -> ignored), so average diag = sqrt(1/mean).
chk("absent Average auto-computed (not 1)", all(diag(d3$average_matrix) != 1))
chk_equal("stream-1 Average = 1/mean(10,20)=1/15", diag(d3$average_matrix)[1], sqrt(1/15), tol = 1e-8)
# cumulative streams are fit on their annual increments, so auto-Average scales
# to the increments (diff(c(0, 5,10,18)) = 5,5,8 -> mean 6), not the raw totals.
dcC <- data.frame(Label = "C", Formula = "cumulative(X1)", check.names = FALSE, stringsAsFactors = FALSE)
dcC[["2013"]] <- "5"; dcC[["2014"]] <- "10"; dcC[["2015"]] <- "18"
dC <- .prepare_data(dcC, tg)
chk_equal("cumulative Average = 1/mean(increments)=1/6", diag(dC$average_matrix)[1], sqrt(1/6), tol = 1e-8)

th_section("column-count vs horizon guard")
tg_bad <- list(startpoint = 2013, endpoint = 2016)   # expects 4 cols, sheet has 3
chk_error("mismatched horizon errors", .prepare_data(dc, tg_bad))

th_summary("prepare-data")
})
