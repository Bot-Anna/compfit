# ============================================================
# fillParams.R
# Write a fitted model back into a modelParams sheet: each FITTED entry
# (name=[lo,hi] or name=[lo,hi]|(init)) is replaced by its estimate in FIXED
# form (*name=value), so the resulting sheet encodes a concrete, fully-specified
# model with no free parameters. Fixed entries (*...) and blanks are untouched.
#
# This is the basis of a counterfactual workflow: fill in the fit, then edit a
# single cell in the new sheet (e.g. change *p=0.6 to *p=0.9) and re-run.
#
#   filled <- fill_params(fit, modelParams)         # returns the filled data frame
#   write_filled_params(fit, scenario_dir, ...)     # writes it to the scenario folder
#
# Estimates: posterior means for a bayes fit, the point estimate for MLE.
# Matching is BY NAME against the fit's parameters and initial states.
# ============================================================

# Extract the natural-scale named estimates (params + states) from a fit.
# summary = "median" (default; robust to skew, transformation-invariant) or
# "mean" (posterior mean). MLE fits ignore this (they have a single point).
.fp_estimates <- function(fit, summary = c("median", "mean")) {
  summary <- match.arg(summary)
  if (!inherits(fit, "compartmentalFit"))
    stop("fill_params() expects a compartmentalFit.")
  if (fit$method == "bayes") {
    if (is.null(fit$samples)) stop("Bayes fit has no samples.")
    # Summarise the posterior draws per parameter. Drop likelihood
    # hyperparameters (sigma, phi*) -- they are not model quantities.
    draws <- posterior_draws(fit)
    keep  <- !grepl("^sigma$|^phi[0-9]*$", names(draws))
    draws <- draws[, keep, drop = FALSE]
    est <- if (summary == "median")
             vapply(draws, median, numeric(1), na.rm = TRUE)
           else
             vapply(draws, mean, numeric(1), na.rm = TRUE)
    as.list(est)
  } else {
    if (is.null(fit$point)) stop("MLE-type fit has no point estimate.")
    c(as.list(unlist(fit$point$parms)),
      as.list(unlist(fit$point$initial_state)))
  }
}

# Rewrite a single cell string. If it is a FITTED entry (no leading '*', of the
# form name=...), and `name` has an estimate, return "*name=<value>". Otherwise
# return the cell unchanged (fixed entries, blanks, functions, unmatched names).
.fp_rewrite_cell <- function(cell, estimates, digits = NULL) {
  if (is.na(cell)) return(cell)
  s <- gsub(" ", "", as.character(cell))
  if (s == "") return(cell)
  if (startsWith(s, "*")) return(cell)              # already fixed -> leave

  nm <- sub("=.*", "", s)                           # name before '='
  if (!nzchar(nm) || is.null(estimates[[nm]])) return(cell)  # not fitted/known

  val <- estimates[[nm]]
  val_str <- if (is.null(digits)) format(val, scientific = FALSE, trim = TRUE)
             else formatC(val, format = "g", digits = digits)
  paste0("*", nm, "=", val_str)
}

# Fill a modelParams data frame: rewrite the States and Parameters columns,
# leaving ALL other columns (the _Level* compartment columns, Others, Functions,
# Conditions, Linear1..n, Quadratic1..n, ...) exactly as-is so the filled sheet still
# rebuilds the identical compartmental model. check.names = FALSE is essential:
# the compartment columns are detected by a leading underscore (^_), which R's
# default name-mangling would destroy.
#' Write fitted estimates back into a model sheet
#'
#' Returns a copy of the model parameter sheet with fitted States/Parameters
#' replaced by their point estimates (fixing them as `*name=value`); structural
#' columns are left intact.
#'
#' @param fit A `"compartmentalFit"` object.
#' @param modelParams The model parameter sheet (data frame) to fill.
#' @param digits Optional rounding for written values.
#' @param summary Which point estimate to write (`"median"` or `"mean"`).
#' @return The filled model parameter data frame.
#' @examples
#' \dontrun{
#' # fit + modelParams from a fitted scenario
#' filled <- fill_params(fit, sc$modelParams, digits = 4)
#' }
#' @export
fill_params <- function(fit, modelParams, digits = NULL,
                        summary = c("median", "mean")) {
  estimates <- .fp_estimates(fit, summary = match.arg(summary))
  mp <- as.data.frame(modelParams, stringsAsFactors = FALSE,
                      check.names = FALSE)

  for (col in intersect(c("States", "Parameters"), names(mp))) {
    mp[[col]] <- vapply(mp[[col]],
                        function(cell) .fp_rewrite_cell(cell, estimates, digits),
                        character(1))
  }
  mp
}

# Fill and WRITE the sheet into the scenario folder. Returns the path.
#   fit          a compartmentalFit
#   scenario_dir folder to write into (e.g. sc$dir)
#   out_file     output filename (default derived with a _filled suffix)
#   digits       optional significant digits for the written values
#   summary      "median" (default) or "mean" for the bayes point estimate
#' Write a filled model sheet to a workbook
#'
#' Fills the sheet via [fill_params()] and writes it to an `.xlsx` in the
#' scenario directory.
#'
#' @param fit A `"compartmentalFit"` object.
#' @param scenario_dir Directory to write the workbook into.
#' @param modelParams The model parameter sheet (defaults to the fit's stored
#'   sheet).
#' @param out_file Output workbook name.
#' @param digits Optional rounding for written values.
#' @param summary Which point estimate to write (`"median"` or `"mean"`).
#' @return The output path, invisibly.
#' @examples
#' \dontrun{
#' write_filled_params(fit, tempdir(), modelParams = sc$modelParams)
#' }
#' @export
write_filled_params <- function(fit, scenario_dir,
                                modelParams = fit$model$modelParams,
                                out_file = "modelParams_filled.xlsx",
                                digits = NULL,
                                summary = c("median", "mean")) {
  if (is.null(modelParams))
    stop("No modelParams available. Pass modelParams = sc$modelParams explicitly ",
         "(older fits may not store the sheet).")

  filled <- fill_params(fit, modelParams, digits = digits,
                        summary = match.arg(summary))
  out_path <- file.path(scenario_dir, out_file)

  ext <- tolower(tools::file_ext(out_file))
  if (ext %in% c("xlsx", "xls")) {
    # writexl preserves column order and names exactly (no name mangling) and
    # writes NA as an empty cell -- matching the blank-cell semantics that the
    # Linearj / Quadraticj / _Level* columns rely on. Preferred over openxlsx for
    # this fidelity.
    if (requireNamespace("writexl", quietly = TRUE)) {
      writexl::write_xlsx(filled, out_path)
    } else if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(filled, out_path, keepNA = FALSE)
    } else {
      stop("Writing xlsx needs the 'writexl' (preferred) or 'openxlsx' package.")
    }
  } else if (ext == "csv") {
    utils::write.csv(filled, out_path, row.names = FALSE, na = "")
  } else if (ext %in% c("tsv", "txt")) {
    utils::write.table(filled, out_path, sep = "\t", row.names = FALSE, na = "")
  } else if (ext == "rds") {
    saveRDS(filled, out_path)        # rds round-trips perfectly (no Excel quirks)
  } else {
    stop("Unsupported output extension: ", ext)
  }

  message("Wrote filled parameter sheet to ", out_path,
          "  (", ncol(filled), " columns preserved)")
  invisible(out_path)
}

# Verify a written filled sheet rebuilds the SAME model structure as the original
# (column names, compartment structure, first/second-order coefficient matrices,
# time grid). Reads the file back and compares against the original modelParams.
# Returns TRUE invisibly if all checks pass; otherwise stops with the mismatch.
#' Verify a filled model sheet
#'
#' Re-reads a filled workbook and checks it rebuilds the same model structure as
#' the original sheet (catches accidental structural edits).
#'
#' @param filled_path Path to the filled `.xlsx`.
#' @param modelParams The original model parameter sheet (data frame).
#' @return `TRUE`, invisibly, if all checks pass; otherwise an error.
#' @examples
#' \dontrun{
#' verify_filled(file.path(tempdir(), "modelParams_filled.xlsx"), sc$modelParams)
#' }
#' @export
verify_filled <- function(filled_path, modelParams) {
  # Read as text to match how load_scenario() reads modelParams, so the
  # structural comparison below isn't confused by a column the reader happens to
  # type numerically on one side and as text on the other.
  reread <- read_data_file(filled_path, text_cols = TRUE)

  # 1. Same columns, same order (critical: ^_ Level cols, Quadratic numbering).
  if (!identical(names(reread), names(modelParams)))
    stop("Column names/order differ after round-trip:\n  original: ",
         paste(names(modelParams), collapse = ", "), "\n  reread:   ",
         paste(names(reread), collapse = ", "))

  # 2. Same number of rows.
  if (nrow(reread) != nrow(modelParams))
    stop(sprintf("Row count differs: original %d, reread %d.",
                 nrow(modelParams), nrow(reread)))

  # 3. Structural columns identical (everything EXCEPT States/Parameters, which
  #    we intentionally changed). Compare as trimmed character, NA-aware.
  norm <- function(x) { v <- trimws(as.character(x)); v[is.na(v)] <- ""; v }
  struct_cols <- setdiff(names(modelParams), c("States", "Parameters"))
  for (col in struct_cols) {
    if (!identical(norm(reread[[col]]), norm(modelParams[[col]])))
      stop("Structural column '", col, "' changed after round-trip; the rebuilt ",
           "model would differ. Check the xlsx writer's blank-cell handling.")
  }

  message("verify_filled: OK -- ", length(struct_cols),
          " structural columns identical; the filled sheet rebuilds the same model.")
  invisible(TRUE)
}
