# ============================================================
# helper-compfit.R
# Bridges the original zero-dependency harness API (chk*, th_*) onto testthat,
# so the migrated test bodies run unchanged. Loaded automatically by testthat.
#
#   chk / chk_equal / chk_error / chk_ok  -> testthat expectations
#   th_skip                               -> testthat::skip() (throws; any
#                                            `quit()`/return after it is dead code)
#   th_load_pure                          -> no-op (the package is already loaded)
#   th_have_setup                         -> initialises Julia (skips on failure)
# ============================================================

# ---- Project root (testthat runs with wd = tests/testthat) ------------------
TH_root <- local({
  if (nzchar(Sys.getenv("FITCM_ROOT"))) return(Sys.getenv("FITCM_ROOT"))
  cand <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  if (dir.exists(file.path(cand, "R"))) return(cand)
  if (dir.exists("R")) return(normalizePath("."))
  cand
})

# ---- Assertions (testthat-backed) ------------------------------------------
chk <- function(label, cond) testthat::expect_true(isTRUE(cond), info = label)

chk_equal <- function(label, actual, expected, tol = 1e-8, check_names = FALSE) {
  a <- actual; e <- expected
  if (!check_names) { a <- unname(a); e <- unname(e) }
  testthat::expect_equal(a, e, tolerance = tol, info = label)
}

chk_error <- function(label, expr) testthat::expect_error(expr, info = label)

chk_ok <- function(label, expr) {
  res <- NULL
  err <- tryCatch({ res <- expr; NULL }, error = function(e) e)
  if (is.null(err)) testthat::succeed(label)
  else testthat::fail(sprintf("%s: %s", label, conditionMessage(err)))
  invisible(res)
}

# ---- Flow / reporting (no-ops under testthat) ------------------------------
th_skip    <- function(label, reason) testthat::skip(sprintf("%s -- %s", label, reason))
th_section <- function(name) invisible(NULL)
th_summary <- function(title = "") invisible(TRUE)

# ---- Loading project code --------------------------------------------------
# The package is already loaded (devtools::test / R CMD check), so sourcing
# individual R/ files is unnecessary. Kept as a no-op for body compatibility.
th_load_pure <- function(files) invisible(NULL)

# Initialise Julia for integration tests; skip the test if it is unavailable.
th_have_setup <- function() {
  if (isTRUE(getOption("th_setup_done"))) return(TRUE)
  ok <- tryCatch({ setup_julia(); TRUE },
                 error = function(e) { message("setup_julia() failed: ",
                                               conditionMessage(e)); FALSE })
  options(th_setup_done = ok)
  ok
}

# ---- Integration data config (committed fixtures; env override) ------------
# Fixtures ship as package data under inst/extdata; resolve via system.file()
# (works once installed / under load_all), with a source-tree fallback.
fixture_dir <- function(name) {
  d <- system.file("extdata", name, package = "compfit")
  if (nzchar(d) && dir.exists(d)) return(d)
  file.path(TH_root, "inst", "extdata", name)
}

th_scenario_dir <- function() {
  d <- Sys.getenv("FITCM_SCENARIO_DIR", unset = "")
  if (nzchar(d)) return(d)
  fx <- fixture_dir("minimal")
  if (dir.exists(fx)) fx else ""
}

th_pick <- function(dir, candidates) {
  hit <- candidates[file.exists(file.path(dir, candidates))]
  if (length(hit)) hit[1] else ""
}

th_scenario_files <- function(dir) {
  list(
    combined = th_pick(dir, c(Sys.getenv("FITCM_COMBINED", "dataCombined.xlsx"),
                              "dataCombined.xlsx", "dataCombined.csv")),
    dummy    = th_pick(dir, c(Sys.getenv("FITCM_DUMMY", "dataDummy.xlsx"),
                              "dataDummy.xlsx", "dataDummy.csv")),
    params   = th_pick(dir, c(Sys.getenv("FITCM_PARAMS", "model.xlsx"),
                              "model.xlsx", "modelParams.xlsx",
                              "modelParams.csv", "model.csv"))
  )
}

th_skip_if_no_scenario <- function(label) {
  d <- th_scenario_dir()
  if (!nzchar(d) || !dir.exists(d)) {
    th_skip(label, "set FITCM_SCENARIO_DIR to a scenario folder"); return(TRUE)
  }
  f <- th_scenario_files(d)
  if (!nzchar(f$combined) || !nzchar(f$params)) {
    th_skip(label, "scenario folder missing dataCombined/params workbook"); return(TRUE)
  }
  FALSE
}
