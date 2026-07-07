# ============================================================
# tools/migrate_linear_columns.R
# One-off migration: convert a modelParams sheet from the OLD single column-major
# `Linear` column to the NEW one-column-per-compartment layout Linear1..Linear<n>
# (mirroring Quadratic<j>). Output columns are interleaved per compartment:
#   <meta columns...>, Linear1, Quadratic1, Linear2, Quadratic2, ...
# The unused `Description` column is dropped.
#
# The old `Linear` column listed the n*n first-order coefficients column-major;
# Linear<j> is the j-th block of n consecutive entries (= column j of the
# coefficient matrix). Behaviour-preserving: the rebuilt ODE is identical (the
# model reads columns by NAME, so the interleaving is purely cosmetic).
#
# CELL FILL COLOURS are carried over for .xlsx inputs (needs 'openxlsx'): every
# cell keeps its background fill -- including THEME colours with tint, resolved
# against the workbook's theme palette -- and each new Linear<j> cell inherits
# the colour of the old Linear cell whose value it now holds. Other formatting
# (fonts, borders, widths, number formats) is NOT preserved. Disable with
# preserve_colours = FALSE.
#
# USAGE
#   source("tools/migrate_linear_columns.R")
#   migrate_linear_columns("model.xlsx")                 # -> model_LinearJ.xlsx
#   migrate_linear_columns("model.xlsx", "model_new.xlsx")
#   Rscript tools/migrate_linear_columns.R model.xlsx [output.(xlsx|csv)]
#
# Handles .xlsx/.xls and .csv. Reads values as text. Migrates the FIRST sheet.
# ============================================================

# Read a sheet as an all-character data frame, empties as NA, names verbatim.
.mig_read <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading Excel needs the 'readxl' package.")
    as.data.frame(readxl::read_excel(path, col_types = "text"), check.names = FALSE)
  } else if (ext == "csv") {
    read.csv(path, colClasses = "character", check.names = FALSE, na.strings = "")
  } else {
    stop("Unsupported input extension '.", ext, "' (use .xlsx, .xls, or .csv).")
  }
}

# Plain (no-colour) writer.
.mig_write <- function(df, path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (requireNamespace("writexl", quietly = TRUE)) {
      writexl::write_xlsx(df, path)
    } else if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(df, path, keepNA = FALSE)
    } else stop("Writing Excel needs 'writexl' (preferred) or 'openxlsx'.")
  } else if (ext == "csv") {
    utils::write.csv(df, path, row.names = FALSE, na = "")
  } else {
    stop("Unsupported output extension '.", ext, "' (use .xlsx or .csv).")
  }
}

# ---- colour helpers --------------------------------------------------------
# ARGB ("FFFF0000") or RGB ("FF0000") -> "#RRGGBB"; NA/invalid -> NA.
.mig_argb_to_hex <- function(argb) {
  if (is.null(argb) || is.na(argb) || !nzchar(argb)) return(NA_character_)
  s <- toupper(gsub("[^0-9A-F]", "", toupper(argb)))
  if (nchar(s) == 8) s <- substr(s, 3, 8)      # drop alpha
  if (nchar(s) != 6) return(NA_character_)
  paste0("#", s)
}

# Theme colour palette indexed by Excel's theme attribute (0..11), resolved from
# the workbook's xl/theme/theme*.xml. Returns 6-hex strings, or NULL if absent.
.mig_theme_palette <- function(path) {
  tryCatch({
    z  <- utils::unzip(path, list = TRUE)
    tf <- grep("xl/theme/theme.*\\.xml$", z$Name, value = TRUE)
    if (!length(tf)) return(NULL)
    ex <- file.path(tempdir(), paste0("themx_", as.integer(runif(1, 1, 1e8))))
    dir.create(ex); on.exit(unlink(ex, recursive = TRUE), add = TRUE)
    utils::unzip(path, files = tf[1], exdir = ex)
    xml   <- paste(readLines(file.path(ex, tf[1]), warn = FALSE), collapse = "")
    block <- sub('.*<a:clrScheme[^>]*>(.*?)</a:clrScheme>.*', '\\1', xml, perl = TRUE)
    if (identical(block, xml)) return(NULL)
    m   <- regmatches(block, gregexpr('(?:val|lastClr)="[0-9A-Fa-f]{6}"', block, perl = TRUE))[[1]]
    clr <- toupper(gsub('.*"([0-9A-Fa-f]{6})".*', '\\1', m))   # clrScheme order
    if (length(clr) < 12) return(NULL)
    # clrScheme order: dk1,lt1,dk2,lt2,accent1..6,hlink,folHlink.
    # Theme index maps with the first two pairs swapped (0<->lt1, 1<->dk1, ...).
    c(clr[2], clr[1], clr[4], clr[3], clr[5:12])
  }, error = function(e) NULL)
}

.mig_rgb2hsl <- function(r, g, b) {
  mx <- max(r, g, b); mn <- min(r, g, b); l <- (mx + mn) / 2
  if (mx == mn) return(c(0, 0, l))
  d <- mx - mn
  s <- if (l > 0.5) d / (2 - mx - mn) else d / (mx + mn)
  h <- if (mx == r) (g - b) / d + (if (g < b) 6 else 0)
       else if (mx == g) (b - r) / d + 2 else (r - g) / d + 4
  c(h / 6, s, l)
}
.mig_hsl2rgb <- function(h, s, l) {
  if (s == 0) return(c(l, l, l))
  hue <- function(p, q, t) {
    if (t < 0) t <- t + 1; if (t > 1) t <- t - 1
    if (t < 1/6) p + (q - p) * 6 * t
    else if (t < 1/2) q
    else if (t < 2/3) p + (q - p) * (2/3 - t) * 6
    else p
  }
  q <- if (l < 0.5) l * (1 + s) else l + s - l * s
  p <- 2 * l - q
  c(hue(p, q, h + 1/3), hue(p, q, h), hue(p, q, h - 1/3))
}
# Apply an Excel tint to a 6-hex colour (lighten if tint>0, darken if <0).
.mig_apply_tint <- function(hex6, tint) {
  if (is.na(tint) || tint == 0) return(hex6)
  rgb <- strtoi(substring(hex6, c(1, 3, 5), c(2, 4, 6)), 16L) / 255
  hsl <- .mig_rgb2hsl(rgb[1], rgb[2], rgb[3]); l <- hsl[3]
  l   <- if (tint < 0) l * (1 + tint) else l * (1 - tint) + tint
  out <- .mig_hsl2rgb(hsl[1], hsl[2], max(0, min(1, l)))
  sprintf("%02X%02X%02X", round(out[1] * 255), round(out[2] * 255), round(out[3] * 255))
}
# Resolve an openxlsx fillFg (named: rgb=, or theme=+tint=) to "#RRGGBB" or NA.
# NOTE: when openxlsx loads an Excel-authored file the attribute names carry a
# leading space (" rgb", " theme", " tint"), so match on TRIMMED names.
.mig_resolve_fill <- function(fg, palette) {
  if (is.null(fg) || !length(fg)) return(NA_character_)
  nm  <- trimws(names(fg))
  val <- function(key) { i <- which(nm == key); if (length(i)) unname(fg[[i[1]]]) else NA_character_ }
  rgb <- val("rgb")
  if (!is.na(rgb)) return(.mig_argb_to_hex(rgb))
  thm <- val("theme")
  if (!is.na(thm)) {
    idx <- suppressWarnings(as.integer(thm))
    if (is.na(idx) || idx < 0 || idx > 11 || is.null(palette)) return(NA_character_)
    base <- palette[idx + 1L]; if (is.na(base)) return(NA_character_)
    tnt  <- val("tint"); tint <- if (!is.na(tnt)) suppressWarnings(as.numeric(tnt)) else 0
    return(paste0("#", .mig_apply_tint(base, tint)))
  }
  NA_character_
}

# "A12" -> column number (12 -> 12). Vectorised over addresses.
.mig_addr_rc <- function(addr) {
  letters_part <- toupper(gsub("[0-9]+$", "", addr))
  rows <- as.integer(gsub("^[A-Z]+", "", addr))
  cols <- vapply(letters_part, function(s) {
    chars <- utf8ToInt(s) - utf8ToInt("A") + 1L
    Reduce(function(a, b) a * 26L + b, chars, accumulate = FALSE)
  }, integer(1))
  list(row = rows, col = unname(cols))
}

# Legacy 56-colour indexed palette (for the rare <fgColor indexed="N"/>).
.mig_indexed_palette <- c(
  "000000","FFFFFF","FF0000","00FF00","0000FF","FFFF00","FF00FF","00FFFF",
  "000000","FFFFFF","FF0000","00FF00","0000FF","FFFF00","FF00FF","00FFFF",
  "800000","008000","000080","808000","800080","008080","C0C0C0","808080",
  "9999FF","993366","FFFFCC","CCFFFF","660066","FF8080","0066CC","CCCCFF",
  "000080","FF00FF","FFFF00","00FFFF","800080","800000","008080","0000FF",
  "00CCFF","CCFFFF","CCFFCC","FFFF99","99CCFF","FF99CC","CC99FF","FFCC99",
  "3366FF","33CCCC","99CC00","FFCC00","FF9900","FF6600","666699","969696",
  "003366","339966","003300","333300","993300","993366","333399","333333")

# Resolve one <fgColor .../> (attributes rgb / theme+tint / indexed) to "#RRGGBB".
.mig_fgcolor_hex <- function(rgb, theme, tint, indexed, palette) {
  if (!is.na(rgb)) return(.mig_argb_to_hex(rgb))
  if (!is.na(theme)) {
    idx <- suppressWarnings(as.integer(theme))
    if (is.na(idx) || idx < 0 || idx > 11 || is.null(palette)) return(NA_character_)
    base <- palette[idx + 1L]; if (is.na(base)) return(NA_character_)
    tn <- if (!is.na(tint)) suppressWarnings(as.numeric(tint)) else 0
    return(paste0("#", .mig_apply_tint(base, tn)))
  }
  if (!is.na(indexed)) {
    i <- suppressWarnings(as.integer(indexed))
    if (!is.na(i) && i >= 0 && i < length(.mig_indexed_palette))
      return(paste0("#", .mig_indexed_palette[i + 1L]))
  }
  NA_character_
}

# Per-cell fill colours of the first worksheet as a (sheet_row x sheet_col)
# matrix of "#RRGGBB" (NA where unfilled). Row 1 is the header. Reads the raw
# xlsx XML directly (authoritative for EVERY cell) when xml2 is available --
# handling RGB, theme+tint, and indexed fills -- and falls back to openxlsx.
.mig_fill_matrix <- function(path, n_rows, n_cols) {
  sheet1 <- tryCatch(openxlsx::sheets(openxlsx::loadWorkbook(path))[1],
                     error = function(e) "Sheet1")
  M <- matrix(NA_character_, n_rows, n_cols)

  ok <- requireNamespace("xml2", quietly = TRUE)
  if (ok) ok <- tryCatch({
    ns  <- c(d = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
    ex  <- file.path(tempdir(), paste0("xlsxx_", as.integer(runif(1, 1, 1e8))))
    dir.create(ex); on.exit(unlink(ex, recursive = TRUE), add = TRUE)
    utils::unzip(path, exdir = ex)
    palette <- .mig_theme_palette(path)

    sty   <- xml2::read_xml(file.path(ex, "xl/styles.xml"))
    # fillId -> resolved hex (solid fills only)
    fills <- xml2::xml_find_all(sty, ".//d:fills/d:fill", ns)
    fill_hex <- vapply(fills, function(f) {
      pf <- xml2::xml_find_first(f, "./d:patternFill", ns)
      if (is.na(pf) || !identical(xml2::xml_attr(pf, "patternType"), "solid"))
        return(NA_character_)
      fg <- xml2::xml_find_first(pf, "./d:fgColor", ns)
      if (is.na(fg)) return(NA_character_)
      .mig_fgcolor_hex(xml2::xml_attr(fg, "rgb"), xml2::xml_attr(fg, "theme"),
                       xml2::xml_attr(fg, "tint"), xml2::xml_attr(fg, "indexed"),
                       palette)
    }, character(1))
    # xf index -> fillId
    xfs    <- xml2::xml_find_all(sty, ".//d:cellXfs/d:xf", ns)
    xf_fid <- as.integer(xml2::xml_attr(xfs, "fillId"))
    xf_hex <- ifelse(is.na(xf_fid), NA_character_, fill_hex[xf_fid + 1L])

    sh    <- xml2::read_xml(file.path(ex, "xl/worksheets/sheet1.xml"))
    cells <- xml2::xml_find_all(sh, ".//d:c", ns)
    s     <- as.integer(xml2::xml_attr(cells, "s")); s[is.na(s)] <- 0L
    hex   <- xf_hex[s + 1L]
    keep  <- which(!is.na(hex))
    if (length(keep)) {
      rc <- .mig_addr_rc(xml2::xml_attr(cells[keep], "r"))
      for (k in seq_along(keep)) {
        r <- rc$row[k]; cc <- rc$col[k]
        if (r >= 1 && r <= n_rows && cc >= 1 && cc <= n_cols) M[r, cc] <- hex[keep[k]]
      }
    }
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {                                   # fallback: openxlsx styleObjects
    wb      <- openxlsx::loadWorkbook(path)
    palette <- .mig_theme_palette(path)
    for (so in wb$styleObjects) {
      if (!identical(so$sheet, sheet1)) next
      hx <- .mig_resolve_fill(so$style$fill$fillFg, palette)
      if (is.na(hx)) next
      for (k in seq_along(so$rows))
        if (so$rows[k] <= n_rows && so$cols[k] <= n_cols) M[so$rows[k], so$cols[k]] <- hx
    }
  }
  list(M = M, sheet = sheet1)
}

# Colour-preserving writer: out_df + column provenance + the original fill matrix.
.mig_write_coloured <- function(out_df, prov, M, sheet_name, n, n_data, path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, 1, out_df, keepNA = FALSE)

  by_colour <- list()
  add_cell <- function(hex, row, col) {
    if (is.na(hex)) return(invisible())
    by_colour[[hex]] <<- rbind(by_colour[[hex]], c(row, col))
  }
  for (oc in seq_along(prov)) {
    p <- prov[[oc]]; src <- p$src_col
    add_cell(M[1L, src], 1L, oc)                          # header colour
    for (r in seq_len(n_data)) {
      hex <- if (identical(p$kind, "linear")) {
               if (r <= n) M[(p$j - 1L) * n + r + 1L, src] else NA_character_
             } else M[r + 1L, src]                        # copy: same data row
      add_cell(hex, r + 1L, oc)
    }
  }
  for (hex in names(by_colour)) {
    rc <- by_colour[[hex]]
    openxlsx::addStyle(wb, 1, openxlsx::createStyle(fgFill = hex),
                       rows = rc[, 1], cols = rc[, 2], gridExpand = FALSE, stack = TRUE)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

migrate_linear_columns <- function(in_path, out_path = NULL,
                                   drop_description = TRUE,
                                   preserve_colours = TRUE) {
  mp <- .mig_read(in_path)

  if (!"Linear" %in% names(mp)) {
    if (any(grepl("^Linear[0-9]+$", names(mp))))
      stop("This sheet already uses Linear<j> columns; nothing to migrate.")
    stop("No single 'Linear' column found in: ", in_path)
  }

  lin      <- as.character(mp$Linear)
  lin_vals <- lin[!is.na(lin) & trimws(lin) != ""]
  cnt      <- length(lin_vals)
  if (cnt == 0L) stop("The 'Linear' column is empty.")
  n <- sqrt(cnt)
  if (n != floor(n))
    stop(sprintf(paste0("The 'Linear' column has %d non-empty entries, not a ",
                        "perfect square; cannot infer the number of compartments."),
                 cnt))
  n <- as.integer(n)

  lev <- grep("^_", names(mp), value = TRUE, ignore.case = TRUE)
  if (length(lev)) {
    idx  <- suppressWarnings(as.numeric(unlist(lapply(mp[lev], as.character))))
    nmax <- suppressWarnings(max(idx, na.rm = TRUE))
    if (is.finite(nmax) && nmax != n)
      warning(sprintf(paste0("Inferred n = %d from the Linear column, but the ",
                             "_Level columns imply %d compartments. Proceeding ",
                             "with n = %d -- check the sheet."), n, nmax, n))
  }

  nrow_sheet <- nrow(mp)
  if (nrow_sheet < n)
    stop(sprintf("Sheet has only %d rows but the model has %d compartments.",
                 nrow_sheet, n))
  new_lin <- lapply(seq_len(n), function(j) {
    chunk <- lin_vals[((j - 1L) * n + 1L):(j * n)]
    v <- rep(NA_character_, nrow_sheet); v[seq_len(n)] <- chunk; v
  })

  linear_idx   <- match("Linear", names(mp))
  quad_names   <- grep("^Quadratic[0-9]+$", names(mp), value = TRUE)
  quad_names   <- quad_names[order(as.integer(sub("^Quadratic", "", quad_names)))]
  dropped_desc <- isTRUE(drop_description) && "Description" %in% names(mp)

  # Meta columns = everything except the old Linear, the Quadratic<j> columns,
  # and (if dropping) Description -- in their original order.
  exclude    <- c("Linear", quad_names, if (dropped_desc) "Description")
  meta_names <- setdiff(names(mp), exclude)

  out <- list(); prov <- list()
  push <- function(name, value, p) { out[[name]] <<- value; prov[[length(prov) + 1L]] <<- p }
  for (cn in meta_names) push(cn, mp[[cn]], list(kind = "copy", src_col = match(cn, names(mp))))
  # Interleaved per compartment: Linear1, Quadratic1, Linear2, Quadratic2, ...
  for (i in seq_len(n)) {
    push(paste0("Linear", i), new_lin[[i]], list(kind = "linear", j = i, src_col = linear_idx))
    qn <- paste0("Quadratic", i)
    if (qn %in% names(mp)) push(qn, mp[[qn]], list(kind = "copy", src_col = match(qn, names(mp))))
  }
  # Any extra Quadratic columns beyond n (unusual) keep their place at the end.
  for (qn in quad_names) {
    qi <- as.integer(sub("^Quadratic", "", qn))
    if (qi > n) push(qn, mp[[qn]], list(kind = "copy", src_col = match(qn, names(mp))))
  }
  out_df <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)

  if (is.null(out_path)) {
    ext <- tools::file_ext(in_path)
    out_path <- file.path(dirname(in_path),
                          paste0(tools::file_path_sans_ext(basename(in_path)),
                                 "_LinearJ.", ext))
  }

  in_xlsx  <- tolower(tools::file_ext(in_path))  %in% c("xlsx", "xls")
  out_xlsx <- tolower(tools::file_ext(out_path)) %in% c("xlsx", "xls")
  have_ox  <- requireNamespace("openxlsx", quietly = TRUE)

  coloured <- FALSE
  if (isTRUE(preserve_colours) && in_xlsx && out_xlsx && have_ox) {
    fm <- .mig_fill_matrix(in_path, n_rows = nrow_sheet + 1L, n_cols = ncol(mp))
    .mig_write_coloured(out_df, prov, fm$M, fm$sheet, n, nrow_sheet, out_path)
    coloured <- TRUE
  } else {
    .mig_write(out_df, out_path)
    if (isTRUE(preserve_colours) && in_xlsx && out_xlsx && !have_ox)
      message("(Install 'openxlsx' to carry over cell fill colours.)")
    else if (isTRUE(preserve_colours) && !(in_xlsx && out_xlsx))
      message("(Colour preservation applies only to .xlsx -> .xlsx migrations.)")
  }

  message(sprintf("Migrated a %d-compartment model: split 'Linear' into %s, interleaved as Linear<i>/Quadratic<i>%s%s.",
                  n, paste0("Linear1..Linear", n),
                  if (dropped_desc) "; dropped 'Description'" else "",
                  if (coloured) "; cell fill colours carried over" else ""))
  message("Wrote: ", out_path)
  invisible(out_path)
}

# ---- CLI entry point (only when run directly via Rscript) ------------------
if (identical(environment(), globalenv()) && sys.nframe() == 0L) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (!length(.args)) {
    cat("Usage: Rscript tools/migrate_linear_columns.R <input.(xlsx|csv)> [output.(xlsx|csv)]\n")
    quit(status = 1L)
  }
  migrate_linear_columns(.args[1], if (length(.args) >= 2L) .args[2] else NULL)
}
