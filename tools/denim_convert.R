# ============================================================
# denimConvert.R
# Convert between denim's DSL transition syntax and a compfit modelParams sheet.
#
#   denim_to_modelParams(transitions, initialValues, parameters, ...)  # DSL  -> sheet
#   modelParams_to_denim(modelParams)                                  # sheet -> DSL
#
# This is a FORMAT converter for the ODE-representable subset of denim, i.e. the
# transitions that map onto compfit's first-order (Linear) / mass-action
# (Quadratic) rate matrices:
#
#   from -> to = <rate> * from            (linear)      <-> Linear columns
#   from -> to = d_exponential(rate)      (linear)      <-> Linear columns
#   from -> to = <coeff> * from*other/N   (mass action) <-> Quadratic columns
#
# The non-exponential dwell-time distributions (d_gamma, d_weibull, d_lognormal,
# nonparametric, ...) have no plain rate-matrix equivalent, so they are reported
# as unsupported rather than silently mis-converted.
#
# Encoding (verified against inst/extdata/medium):
#   Linear<i>[j]     = coefficient of X_j in dX_i
#   Quadratic<i>[j]  = "*goes_to*coeff" = coeff * X_i * X_j / N added to dX_i and
#                      subtracted from dX_goes_to  (compfit divides by the total
#                      population N automatically).
# ============================================================

# ---- DSL parsing helpers ---------------------------------------------------

# Normalise the accepted DSL inputs to a list of list(from, to, rhs) triples.
# Accepts: a character vector of "from -> to = rhs" lines, OR a named list/vector
# whose names are "from -> to" and values are the rhs strings.
.dsl_transitions <- function(transitions) {
  out <- list()
  add <- function(lhs, rhs) {
    fromto <- strsplit(lhs, "->", fixed = TRUE)[[1]]
    if (length(fromto) != 2)
      stop("Transition '", lhs, "' is not of the form 'from -> to'.")
    out[[length(out) + 1L]] <<- list(from = trimws(fromto[1]),
                                     to   = trimws(fromto[2]),
                                     rhs  = trimws(rhs))
  }
  if (!is.null(names(transitions)) && all(nzchar(names(transitions)))) {
    for (nm in names(transitions)) add(nm, as.character(transitions[[nm]]))
  } else {
    for (line in transitions) {
      line <- trimws(sub("#.*$", "", line))          # strip comments
      if (!nzchar(line)) next
      parts <- strsplit(line, "=", fixed = TRUE)[[1]] # split on first '='
      if (length(parts) < 2) stop("Transition '", line, "' has no '=' RHS.")
      add(parts[1], paste(parts[-1], collapse = "="))
    }
  }
  out
}

# Flatten a multiplicative/division expression into numerator + denominator
# factors (each a deparsed string).
.dsl_flatten <- function(e) {
  if (is.call(e)) {
    op <- as.character(e[[1]])
    if (op == "*") {
      a <- .dsl_flatten(e[[2]]); b <- .dsl_flatten(e[[3]])
      return(list(num = c(a$num, b$num), den = c(a$den, b$den)))
    }
    if (op == "/") {
      a <- .dsl_flatten(e[[2]]); b <- .dsl_flatten(e[[3]])
      return(list(num = c(a$num, b$den), den = c(a$den, b$num)))
    }
    if (op == "(") return(.dsl_flatten(e[[2]]))
  }
  list(num = deparse(e), den = character(0))
}

# Multiply a set of factor strings into one expression string ("1" if empty).
.dsl_join <- function(factors) if (length(factors)) paste(factors, collapse = " * ") else "1"

# ---- DSL -> modelParams ----------------------------------------------------

#' Convert denim DSL transitions to a compfit model sheet
#'
#' Translates a model written in denim's DSL transition syntax into a compfit
#' `modelParams` data frame (optionally written to CSV). Handles the
#' ODE-representable subset: linear rates (`rate * from` or
#' `d_exponential(rate)`) and normalised mass action (`coeff * from*other/N`).
#' Non-exponential dwell-time distributions are reported as unsupported.
#'
#' @param transitions denim transitions: a character vector of
#'   `"from -> to = rhs"` lines, or a named list/vector with names `"from -> to"`
#'   and rhs string values.
#' @param initialValues Named numeric vector of initial compartment values (the
#'   names define the compartments and their order).
#' @param parameters Optional named numeric vector of parameter values; written
#'   as fixed (`*name=value`) — edit to `name=[lo,hi]` to fit.
#' @param startpoint,endpoint,partition,cutoff `Others` settings for the sheet.
#' @param file Optional path; if given, the sheet is written there as CSV.
#' @return A `modelParams` data frame (invisibly if `file` is written), with a
#'   `"compartments"` attribute mapping denim names to `X1, X2, ...`.
#' @export
denim_to_modelParams <- function(transitions, initialValues, parameters = NULL,
                                 startpoint = 0, endpoint = 10, partition = 1,
                                 cutoff = NULL, file = NULL) {
  comps <- names(initialValues)
  if (is.null(comps) || !all(nzchar(comps)))
    stop("initialValues must be a named numeric vector (names = compartments).")
  n   <- length(comps)
  idx <- setNames(seq_len(n), comps)                 # name -> X index
  Xof <- function(nm) sprintf("X%d", idx[[nm]])

  tr <- .dsl_transitions(transitions)

  # Accumulators for the rate matrices (as expression strings).
  lin  <- matrix("0", n, n)   # lin[i, j] = coeff of X_j in dX_i (Linear<i> = row i)
  # Quadratic<i> is read by compfit as matrix(column, n, n): an n x n matrix,
  # so the column MUST hold n^2 entries (column-major). The pair (eq i, other j)
  # lives at matrix[j, 1] -> column-major position j. Anything shorter would be
  # recycled by matrix() and fabricate spurious interaction terms.
  quad <- vector("list", n)
  for (i in seq_len(n)) quad[[i]] <- rep("0", n * n)

  add_lin <- function(i, j, coeff) {
    lin[i, j] <<- if (lin[i, j] == "0") coeff else paste(lin[i, j], "+", coeff)
  }

  for (t in tr) {
    if (!t$from %in% comps) stop("Unknown source compartment '", t$from, "'.")
    if (!t$to   %in% comps) stop("Unknown target compartment '", t$to, "'.")
    a <- idx[[t$from]]; b <- idx[[t$to]]
    e <- parse(text = t$rhs)[[1]]

    # d_exponential(rate) -> linear rate; other d_*/nonparametric -> unsupported.
    if (is.call(e) && as.character(e[[1]]) == "d_exponential") {
      args <- as.list(e)[-1]
      rate <- deparse(if (!is.null(args$rate)) args$rate else args[[1]])
      add_lin(a, a, paste0("-(", rate, ")"))
      add_lin(b, a, rate)
      next
    }
    if (is.call(e) && grepl("^(d_gamma|d_weibull|d_lognormal|nonparametric|constant|multinomial|transprob)$",
                            as.character(e[[1]]))) {
      stop("Transition '", t$from, " -> ", t$to, " = ", t$rhs,
           "': non-exponential dwell distributions have no plain rate-matrix ",
           "equivalent and are not handled by the converter (they need ",
           "stage-expansion / a non-Markovian engine).")
    }

    # Math expression: classify as linear (in `from`) or mass action.
    fl     <- .dsl_flatten(e)
    states_num <- intersect(fl$num, comps)
    if (!(t$from %in% states_num))
      stop("Transition '", t$from, " -> ", t$to, " = ", t$rhs,
           "': the flow must be proportional to the source compartment '",
           t$from, "'.")

    if (length(states_num) == 1L) {
      # linear: rate = (num without `from`) / den
      rate_num <- fl$num[fl$num != t$from]
      rate <- .dsl_join(rate_num)
      if (length(fl$den)) rate <- paste0("(", rate, ") / (", .dsl_join(fl$den), ")")
      add_lin(a, a, paste0("-(", rate, ")"))
      add_lin(b, a, rate)
    } else if (length(states_num) == 2L) {
      # mass action: coeff * from * other / N  (N = a symbol in the denominator)
      other <- setdiff(states_num, t$from)
      c_idx <- idx[[other]]
      num_rest <- fl$num[!(fl$num %in% c(t$from, other))]   # drop the two states
      # Drop ONE population normaliser from the denominator (compfit divides by N).
      den_rest <- fl$den
      if (length(den_rest) >= 1L) den_rest <- den_rest[-1L] else
        stop("Mass-action transition '", t$rhs, "' must be normalised by the ",
             "population, e.g. beta * from * other / N.")
      coeff <- .dsl_join(num_rest)
      if (length(den_rest)) coeff <- paste0("(", coeff, ") / (", .dsl_join(den_rest), ")")
      if (b <= a)
        stop("Mass-action transition '", t$from, " -> ", t$to,
             "': compfit's quadratic encoding needs the target index (", b,
             ") greater than the source index (", a, "). Reorder the compartments.")
      quad[[a]][c_idx] <- sprintf("*%d*-(%s)", b, coeff)   # matrix[c_idx, 1] (col-major pos c_idx)
    } else {
      stop("Transition '", t$rhs, "' involves more than two compartments; only ",
           "linear and pairwise mass-action flows are convertible.")
    }
  }

  # ---- assemble the sheet ----
  states_col <- vapply(seq_len(n),
                       function(i) sprintf("*X%d=%s", i, format(initialValues[[i]])),
                       character(1))
  params_col <- if (length(parameters))
    vapply(names(parameters), function(p) sprintf("*%s=%s", p, format(parameters[[p]])),
           character(1)) else character(0)
  others <- c(sprintf("startpoint=%s", startpoint),
              sprintf("endpoint=%s", endpoint),
              sprintf("partition=%s", partition))
  if (!is.null(cutoff)) others <- c(others, sprintf("cutoff=%s", cutoff))  # optional

  # The Quadratic<i> columns need n^2 rows (compfit reshapes them to n x n), so
  # that sets the sheet's row count.
  nrows <- max(n * n, length(others), length(params_col))
  pad <- function(x, fill = "") c(x, rep(fill, nrows - length(x)))

  df <- data.frame(
    `_Level1` = pad(as.character(seq_len(n))),
    Others      = pad(others),
    States      = pad(states_col),
    Functions   = pad(character(0)),
    Parameters  = pad(params_col),
    Conditions  = pad(character(0)),
    check.names = FALSE, stringsAsFactors = FALSE)

  for (i in seq_len(n)) {
    df[[sprintf("Linear%d", i)]]    <- pad(lin[i, ], fill = "")
    df[[sprintf("Quadratic%d", i)]] <- pad(quad[[i]], fill = "")
  }
  attr(df, "compartments") <- idx

  if (!is.null(file)) {
    utils::write.csv(df, file, row.names = FALSE)
    message("Wrote model sheet to ", file,
            "  (compartments: ", paste(sprintf("%s=X%d", comps, idx), collapse = ", "), ")")
    return(invisible(df))
  }
  df
}

# ---- modelParams -> DSL ----------------------------------------------------

#' Convert a compfit model sheet to denim DSL transitions
#'
#' The inverse of [denim_to_modelParams()]: reads a compfit `modelParams` sheet
#' and emits the equivalent denim DSL transitions, plus the `initialValues` and
#' `parameters`. Time-varying `Functions` are inlined into the expressions.
#'
#' @param modelParams A compfit `modelParams` data frame (or a path to its CSV).
#' @param file Optional path to write the DSL block (as text).
#' @return A list with `dsl` (character transitions), `initialValues`, and
#'   `parameters` (invisibly if `file` is written).
#' @export
modelParams_to_denim <- function(modelParams, file = NULL) {
  if (is.character(modelParams) && length(modelParams) == 1L && file.exists(modelParams))
    modelParams <- utils::read.csv(modelParams, colClasses = "character",
                                   check.names = FALSE)
  mp <- as.data.frame(lapply(modelParams, as.character), stringsAsFactors = FALSE,
                      check.names = FALSE)
  nz <- function(col) { v <- col[!is.na(col) & nzchar(trimws(col))]; trimws(v) }

  n <- sum(grepl("^Linear[0-9]+$", names(mp)))
  get_col <- function(nm) if (nm %in% names(mp)) mp[[nm]][seq_len(n)] else rep("0", n)
  cell <- function(x) { x <- trimws(x); if (is.na(x) || !nzchar(x)) "0" else x }

  # Inline Function definitions (name<-expr) into coefficients.
  fns <- nz(mp$Functions)
  subst <- function(expr) {
    for (f in fns) {
      kv <- strsplit(f, "<-", fixed = TRUE)[[1]]
      if (length(kv) == 2)
        expr <- gsub(paste0("\\b", trimws(kv[1]), "\\b"),
                     paste0("(", trimws(kv[2]), ")"), expr)
    }
    expr
  }

  dsl <- character(0)
  # Linear transitions: off-diagonal positive entry Linear<i>[j] (i != j) = inflow
  # to i from j at that rate -> "Xj -> Xi = rate * Xj".
  for (i in seq_len(n)) {
    col <- get_col(sprintf("Linear%d", i))
    for (j in seq_len(n)) {
      if (i == j) next
      c_ij <- cell(col[j])
      if (c_ij == "0" || startsWith(c_ij, "-")) next    # skip zero / outflow
      dsl <- c(dsl, sprintf("X%d -> X%d = %s * X%d", j, i, subst(c_ij), j))
    }
  }
  # Mass-action: Quadratic<i> is a column reshaped (column-major) into an n x n
  # matrix; a "*goes_to*coeff" at matrix[j, k] means coeff*X[i]*X[j]/N (source i,
  # other j) flowing to goes_to  ->  "Xi -> X<goes_to> = coeff * Xi*Xj/N".
  get_full <- function(nm) {
    v <- if (nm %in% names(mp)) as.character(mp[[nm]]) else character(0)
    length(v) <- n * n          # take/pad to exactly n^2 (matches matrix(.,n,n))
    v
  }
  for (i in seq_len(n)) {
    M <- matrix(get_full(sprintf("Quadratic%d", i)), n, n)
    for (j in seq_len(n)) for (k in seq_len(n)) {
      c_jk <- cell(M[j, k])
      if (!grepl("^\\*[0-9]+\\*", c_jk)) next
      m <- regmatches(c_jk, regexec("^\\*([0-9]+)\\*(.*)$", c_jk))[[1]]
      goes_to <- as.integer(m[2]); coeff <- m[3]
      coeff <- sub("^-\\((.*)\\)$", "\\1", coeff)       # -(beta) -> beta (flow magnitude)
      coeff <- sub("^-", "", coeff)
      dsl <- c(dsl, sprintf("X%d -> X%d = %s * X%d*X%d/N",
                            i, goes_to, subst(coeff), i, j))
    }
  }

  # initialValues + parameters from States / Parameters columns.
  parse_kv <- function(entries) {
    entries <- sub("^\\*", "", entries)                # drop fixed marker
    nm  <- sub("=.*$", "", entries)
    rhs <- sub("^[^=]*=", "", entries)
    setNames(rhs, nm)
  }
  initialValues <- parse_kv(nz(mp$States))
  parameters    <- parse_kv(nz(mp$Parameters))

  res <- list(dsl = dsl, initialValues = initialValues, parameters = parameters)
  if (!is.null(file)) {
    writeLines(c("# denim transitions (converted from a compfit modelParams sheet)",
                 dsl, "",
                 "# initialValues: " , paste(names(initialValues), initialValues, sep = " = "),
                 "# parameters:    ",  paste(names(parameters), parameters, sep = " = ")),
               file)
    message("Wrote denim DSL to ", file)
    return(invisible(res))
  }
  res
}
