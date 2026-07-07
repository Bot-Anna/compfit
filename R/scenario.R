# ============================================================
# scenario.R
# A "scenario" is a single folder containing the three main input
# files (dataCombined, dataDummy, modelParams). Plots and other
# outputs are written to subfolders created relative to that folder.
#
#   load_scenario("/path/to/Data Switzerland and Austria")
#     -> reads the three input files (any supported format)
#     -> ensures an output subfolder exists (created if missing)
#     -> returns everything fitCompartmentalModel() needs + paths
#
# Supported input formats (auto-detected by extension):
#   .xlsx/.xls  Excel            (readxl)
#   .csv        comma-separated  (readr if available, else base)
#   .tsv/.txt   tab-separated    (readr if available, else base)
#   .rds        native R object  (readRDS)
#   .parquet    Parquet          (arrow, optional)
#   .feather    Feather          (arrow, optional)
#   .json       JSON tabular     (jsonlite, optional)
# ============================================================

# Read a single tabular data file, dispatching on its extension. Returns a
# data.frame / tibble. Gives a clear message if an optional package is needed.
#
# text_cols = TRUE forces EVERY column to be read as character, bypassing the
# reader's per-column type guessing. This is essential for the data sheets
# (dataCombined / dataDummy): their cells mix numbers with markers like `x`,
# `<L`, `>=L`, and -- critically -- numbers that Excel may have stored as TEXT
# (a stray text "0" is common). With type guessing, readxl turns such a text
# cell in an otherwise-numeric column into NA, silently dropping the value
# BEFORE any downstream parser sees it (the "works for some, not others"
# symptom). Reading everything as text is lossless: parse_data_cell() and the
# Weight/Average coercion convert each cell to a number afterwards.
#' Read a data file (csv or xlsx)
#'
#' Reads a `.csv` or `.xlsx` data file, optionally forcing every column to be
#' read as text so numbers stored as text round-trip unchanged.
#'
#' @param path Path to a `.csv` or `.xlsx` file.
#' @param text_cols If `TRUE`, read all columns as character.
#' @param ... Passed to the underlying reader.
#' @return A data frame.
#' @examples
#' f <- system.file("extdata", "minimal", "modelParams.csv", package = "compfit")
#' read_data_file(f, text_cols = TRUE)
#' @export
read_data_file <- function(path, text_cols = FALSE, ...) {
  if (!file.exists(path))
    stop("File does not exist:\n  ", path)

  ext <- tolower(tools::file_ext(path))

  need <- function(pkg, fmt)
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Reading %s files needs the '%s' package. Install it with install.packages('%s').",
                   fmt, pkg, pkg))

  # Per-reader argument that forces all-character reads, when text_cols = TRUE.
  has_readr <- requireNamespace("readr", quietly = TRUE)
  xl_types  <- if (text_cols) "text" else NULL
  readr_types <- if (text_cols && has_readr) readr::cols(.default = readr::col_character()) else NULL
  base_classes <- if (text_cols) "character" else NA

  do_read <- function() switch(ext,
         xlsx = { need("readxl", "Excel"); readxl::read_excel(path, col_types = xl_types, ...) },
         xls  = { need("readxl", "Excel"); readxl::read_excel(path, col_types = xl_types, ...) },
         csv  = if (has_readr)
           readr::read_csv(path, show_col_types = FALSE, col_types = readr_types, ...)
         else read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = base_classes, ...),
         tsv  = if (has_readr)
           readr::read_tsv(path, show_col_types = FALSE, col_types = readr_types, ...)
         else read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                         colClasses = base_classes, ...),
         txt  = if (has_readr)
           readr::read_tsv(path, show_col_types = FALSE, col_types = readr_types, ...)
         else read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                         colClasses = base_classes, ...),
         rds  = readRDS(path),
         parquet = { need("arrow", "Parquet"); as.data.frame(arrow::read_parquet(path, ...)) },
         feather = { need("arrow", "Feather"); as.data.frame(arrow::read_feather(path, ...)) },
         json = { need("jsonlite", "JSON");
           as.data.frame(jsonlite::fromJSON(path, ...), stringsAsFactors = FALSE) },
         stop(sprintf("Unsupported file extension '.%s' for:\n  %s\n",
                      ext, path),
              "Supported: xlsx, xls, csv, tsv, txt, rds, parquet, feather, json.")
  )

  # Run the reader; on failure, name the FILE and the reason so the user knows
  # exactly which input failed and why (rather than a bare, path-less reader error).
  out <- tryCatch(
    withCallingHandlers(
      do_read(),
      warning = function(w)                         # readr's generic parse warning
        if (grepl("parsing (issue|problem)", conditionMessage(w), ignore.case = TRUE))
          invokeRestart("muffleWarning")),          # re-surfaced below WITH the location
    error = function(e)
      stop(sprintf("Could not read the data file\n  %s\n  (parsed as a '.%s' file) -- %s",
                   path, ext, conditionMessage(e)), call. = FALSE))

  # readr keeps recoverable parse problems on the returned object; surface the
  # first one WITH its row/column so the user knows where in the file to look.
  if (has_readr && ext %in% c("csv", "tsv", "txt")) {
    pr <- tryCatch(readr::problems(out), error = function(e) NULL)
    if (!is.null(pr) && nrow(pr) > 0) {
      p1 <- pr[1, ]
      warning(sprintf(
        "%s: %d parsing problem(s); first at row %s, column %s (expected %s, got '%s'). Call readr::problems() on the loaded object for the full list.",
        basename(path), nrow(pr), p1[["row"]], p1[["col"]], p1[["expected"]], p1[["actual"]]),
        call. = FALSE)
    }
  }
  out
}

# Resolve a scenario file: if `file` has an extension, use it as given; if it is
# a bare basename (no extension), search the folder for any supported extension.
.resolve_scenario_file <- function(scenario_dir, file, role) {
  # If an extension is present, take the path literally.
  if (nzchar(tools::file_ext(file))) {
    p <- file.path(scenario_dir, file)
    if (!file.exists(p))
      stop(sprintf("Missing %s file:\n  %s\n(pass %s_file= if the name differs.)",
                   role, p, role))
    return(p)
  }
  # No extension: search for basename.<supported-ext>.
  exts <- c("xlsx", "xls", "csv", "tsv", "txt", "rds", "parquet", "feather", "json")
  cand <- file.path(scenario_dir, paste0(file, ".", exts))
  hit  <- cand[file.exists(cand)]
  if (length(hit) == 0)
    stop(sprintf("No %s file found for basename '%s' in:\n  %s\n(tried: %s)",
                 role, file, scenario_dir, paste(exts, collapse = ", ")))
  if (length(hit) > 1)
    warning(sprintf("Multiple %s files match '%s'; using:\n  %s",
                    role, file, hit[1]))
  hit[1]
}

#' Load a scenario folder
#'
#' Reads the input workbooks (combined data, optional dummy data, model
#' parameter sheet) from a scenario directory and assembles the inputs needed by
#' [fitCompartmentalModel()] / [build_compartmental_model()], plus an output
#' plot path.
#'
#' @param scenario_dir Folder containing the scenario inputs.
#' @param combined_file Combined-data workbook name.
#' @param dummy_file Optional dummy-data workbook name.
#' @param params_file Model parameter sheet name.
#' @param plots_subdir Sub-folder for output plots.
#' @param plot_file Default output plot file name.
#' @return A list with the loaded `dataCombined`, `dataDummy`, `modelParams`,
#'   and the resolved plot path.
#' @examples
#' # The package ships a minimal example scenario as csv:
#' mini <- system.file("extdata", "minimal", package = "compfit")
#' sc <- load_scenario(mini,
#'                     combined_file = "dataCombined.csv",
#'                     dummy_file    = "dataDummy.csv",
#'                     params_file   = "modelParams.csv")
#' str(sc, max.level = 1)
#' @export
load_scenario <- function(scenario_dir,
                          combined_file = "dataCombined.xlsx",
                          dummy_file    = "dataDummy.xlsx",
                          params_file   = "modelParams.xlsx",
                          plots_subdir  = "Plots",
                          plot_file     = "Plot.pdf") {
  
  if (!dir.exists(scenario_dir)) {
    stop(sprintf("Scenario folder does not exist:\n  %s", scenario_dir))
  }
  
  # --- Resolve the three main input files (extension-aware) ---
  combined_path <- .resolve_scenario_file(scenario_dir, combined_file, "combined")
  dummy_path    <- .resolve_scenario_file(scenario_dir, dummy_file,    "dummy")
  params_path   <- .resolve_scenario_file(scenario_dir, params_file,   "params")
  
  # --- Read them (format auto-detected per file; they may even differ) ---
  # dataCombined AND modelParams are read all-as-text (text_cols = TRUE) so that
  # numbers Excel may have stored as text survive the read instead of being
  # silently dropped to NA by the reader's per-column type guessing.
  #
  #  * dataCombined: a stray text "0" in a data cell would otherwise vanish.
  #    Cells are converted to numbers downstream by parse_data_cell() and the
  #    Weight/Average coercion in .prepare_data().
  #  * modelParams: the Linearj / Quadraticj coefficient columns are STRING
  #    expressions (e.g. `g_HA`, `*5*$beta`, or `0`). When such a column mixes
  #    text coefficients with 0s, the guesser may type it numeric and turn the
  #    text coefficients into NA -- which compartmentalFunction() then treats as
  #    0, silently DROPPING model terms (the intermittent "works for some,
  #    not others" bug). Reading as text preserves them; the only columns that
  #    must be numeric are the `_Level*` compartment INDICES, which
  #    numberOfComps() and compartmentalFunction() coerce explicitly.
  dataCombined <- read_data_file(combined_path, text_cols = TRUE)
  dataDummy    <- read_data_file(dummy_path)
  modelParams  <- read_data_file(params_path, text_cols = TRUE)
  
  # --- Resolve the output plot path (do NOT create anything here) ---
  # The Plots folder is created lazily by save_grid()/next_plot_path() only when
  # a plot is actually written, so merely loading a scenario -- including a
  # read-only one such as the bundled fixture -- has no filesystem side effects.
  plot_dir  <- file.path(scenario_dir, plots_subdir)
  plot_path <- file.path(plot_dir, plot_file)
  
  list(
    dir          = scenario_dir,
    dataCombined = dataCombined,
    dataDummy    = dataDummy,
    modelParams  = modelParams,
    plot_dir     = plot_dir,
    plot_path    = plot_path
  )
}

# Auto-incrementing numbered plot path. Given an intended path like
# ".../Plots/Plot.pdf", returns ".../Plots/Plot_(N+1).pdf" where N is the highest
# existing Plot_<number> in that folder (so successive saves never overwrite).
# The first save becomes Plot_1.pdf. A legacy un-numbered "Plot.pdf" is ignored
# for the numbering (only Plot_<digits> files are scanned).
#' Next numbered plot path
#'
#' Given a plot path, returns the next available numbered variant (e.g.
#' `Plot_1.pdf`, `Plot_2.pdf`, ...) so existing plots are never overwritten.
#'
#' @param path A plot file path.
#' @return A character path with the next free numeric suffix.
#' @examples
#' next_plot_path(file.path(tempdir(), "Plot.pdf"))
#' @export
next_plot_path <- function(path) {
  dir  <- dirname(path)
  base <- basename(path)
  ext  <- tools::file_ext(base)                       # e.g. "pdf"
  stem <- tools::file_path_sans_ext(base)             # e.g. "Plot"
  ext_dot <- if (nzchar(ext)) paste0(".", ext) else ""
  
  # Escape any regex-special characters in the stem and extension robustly.
  esc <- function(s) gsub("([][{}().^$*+?\\\\|])", "\\\\\\1", s)
  pat <- sprintf("^%s_(\\d+)%s$",
                 esc(stem),
                 if (nzchar(ext)) paste0("\\.", esc(ext)) else "")
  
  existing <- list.files(dir, pattern = pat)
  next_n <- if (length(existing)) {
    nums <- as.integer(sub(pat, "\\1", existing))
    max(nums, na.rm = TRUE) + 1L
  } else 1L
  
  file.path(dir, sprintf("%s_%d%s", stem, next_n, ext_dot))
}

# Save a plot result (anything with $grid/$width/$height, e.g. from plot_fit,
# plot_counterfactual, plot_counterfactual_effect, stacked_plots) to a NEW
# numbered file under `path` (Plot_1, Plot_2, ...), so previous plots are never
# overwritten. Closes any open graphics devices first. No-op with a message if
# patchwork produced no $grid (the individual panels are in res$plots).
#' Save a plot grid
#'
#' Writes the assembled `$grid` of a plot result to disk. No-op (with a message)
#' when `$grid` is `NULL` (e.g. patchwork unavailable).
#'
#' @param res A plot result list containing a `$grid` element.
#' @param path Output file path.
#' @param limitsize Passed to [ggplot2::ggsave()].
#' @param ... Passed to [ggplot2::ggsave()].
#' @return The output path, invisibly (or `NULL` if there is no grid).
#' @examples
#' \dontrun{
#' res <- plot_fit(fit)            # fit from fitCompartmentalModel()
#' save_grid(res, file.path(tempdir(), "Plot.pdf"))
#' }
#' @export
save_grid <- function(res, path, limitsize = FALSE, ...) {
  if (is.null(res$grid)) {
    message("No $grid (patchwork missing?); individual panels are in res$plots.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("save_grid() needs ggplot2.")
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  out_path <- next_plot_path(path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)  # ensure folder exists
  ggplot2::ggsave(out_path, res$grid, width = res$width, height = res$height,
                  limitsize = limitsize, ...)
  message("Saved plot to ", out_path)
  invisible(out_path)
}