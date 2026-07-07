test_that("test-unit-extract-code", {
# ============================================================
# test-unit-extract-code.R   (pure R; no Julia)
# Guards extract_code's inlined-objects list against lossFunction's formals.
# This is the test that would have caught the lcens_mask/llimit_mat omission:
# whenever lossFunction gains a captured (data) argument, .xc_loss_captured must
# list it, or the emitted self-contained script references an undefined object.
# ============================================================
th_load_pure(c("lossFunction.R", "extract_code.R"))

th_section(".xc_loss_captured stays in sync with lossFunction() formals")
# Every lossFunction formal is either:
#   - an EXPRESSION arg spliced into the body (`*_expression`, names_data_points),
#   - a solver setting baked into the solve call (solver/abstol/reltol), or
#   - a DATA/build object the closure CAPTURES -> must be inlined by extract_code.
fmls <- names(formals(lossFunction))
# comp_names is resolved into the body at build time (it rewrites state
# references in the data formulas, like the *_expression args) -- not captured
# by the returned closure, so it is build-time, not an inlined data object.
# (verbose IS captured by the closure -> inlined, so it stays in .xc_loss_captured.)
spliced  <- c(fmls[grepl("_expression$", fmls)], "names_data_points", "comp_names")
baked    <- c("solver", "abstol", "reltol", "backend")
expected <- setdiff(fmls, c(spliced, baked))

chk("captured set equals lossFunction's data formals",
    setequal(expected, .xc_loss_captured))

missing_from_capture <- setdiff(expected, .xc_loss_captured)
chk("no data formal is missing from .xc_loss_captured (would break the script)",
    length(missing_from_capture) == 0)
if (length(missing_from_capture))
  cat("    MISSING:", paste(missing_from_capture, collapse = ", "), "\n")

stale_in_capture <- setdiff(.xc_loss_captured, expected)
chk("no stale entry in .xc_loss_captured",
    length(stale_in_capture) == 0)
if (length(stale_in_capture))
  cat("    STALE:", paste(stale_in_capture, collapse = ", "), "\n")

# Sanity: the right-censoring args (the ones that were once missing) are present.
chk("lcens_mask is captured", "lcens_mask" %in% .xc_loss_captured)
chk("llimit_mat is captured", "llimit_mat" %in% .xc_loss_captured)

th_summary("extract-code")
})
