# ============================================================
# validate.R -- upfront validity checks for a modelParams sheet.
# Turns entry mistakes (a stray letter in a numeric value, a reversed box, a
# typo'd symbol in a rate coefficient) into a single clear error naming the
# offending column/cell and the problem, instead of a cryptic failure deep in
# the builder. Called at the top of .build_model() (gate: options(compfit.validate=)).
# ============================================================

# Names that would break the generated model: Julia keywords, or the ODE/loss
# codegen's own variables (a parameter named `end` -> `end = p[1]`, or `t`/`X`/
# `du` shadowing an internal). Checked against parameter/function names.
.JULIA_RESERVED <- c(
  "baremodule", "begin", "break", "catch", "const", "continue", "do", "else",
  "elseif", "end", "export", "false", "finally", "for", "function", "global",
  "if", "import", "in", "isa", "let", "local", "macro", "module", "quote",
  "return", "struct", "true", "try", "using", "where", "while", "abstract",
  "mutable", "primitive", "type",
  "t", "p", "X", "du", "dX", "parms", "initial_states")

#' Validate a modelParams sheet
#'
#' Check a model-parameter sheet for common entry errors before it is built.
#' On any problem it raises a single error listing every issue found, each naming
#' the column/cell and what is wrong. Checks:
#' \itemize{
#'   \item fixed \code{States}/\code{Parameters} values and \code{Linear}/
#'     \code{Quadratic}/\code{Constant}/\code{Functions}/\code{Conditions}
#'     expressions reference only declared symbols (parameters, \code{_0}
#'     aliases, the declared states -- \code{X1..Xn} or named, e.g.
#'     \code{S}/\code{I}/\code{R} --, functions, \code{time}), so a typo like
#'     \code{*2*-bta} is caught;
#'   \item box priors \code{[lo,hi]} are two numbers with \code{lo < hi};
#'   \item distribution priors have numeric arguments;
#'   \item no parameter/state/function name collides with a Julia keyword or a
#'     codegen variable/function (\code{t}/\code{p}/\code{X}/\code{du}/\code{dX}/
#'     \code{parms}, \code{N1..}, \code{total_pop}, \code{f<ij>}, \code{g<ij>},
#'     \code{cst<i>}, \code{secOrd_i_j}, or another quantity's \code{_0} alias);
#'   \item \code{Quadratic} cells are \code{*goto*coeff} with a target in
#'     \code{1..n};
#'   \item \code{Others} has numeric \code{startpoint}/\code{endpoint}/\code{partition}.
#' }
#'
#' @param modelParams The model-parameter data frame (as read by
#'   \code{load_scenario()}).
#' @return Invisibly \code{TRUE} when valid; otherwise stops with the collected
#'   problems.
#' @examples
#' mini <- system.file("extdata", "minimal", package = "compfit")
#' sc <- load_scenario(mini, combined_file = "dataCombined.csv",
#'                     dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
#' validate_modelParams(sc$modelParams)   # TRUE (the fixture is valid)
#' @export
validate_modelParams <- function(modelParams) {
  errs <- character(0)
  add  <- function(fmt, ...) errs[[length(errs) + 1L]] <<- sprintf(fmt, ...)

  cs <- tryCatch(numberOfComps(modelParams), error = function(e) NULL)
  if (is.null(cs) || !is.finite(cs$number_of_comps) || cs$number_of_comps < 1)
    stop("validate_modelParams: could not determine the number of compartments ",
         "from the '_'-prefixed index column(s).", call. = FALSE)
  n <- cs$number_of_comps
  state_names <- .compartments(modelParams)          # canonical names (X1..Xn or S,I,R,...)
  if (length(state_names) < n) state_names <- union(state_names, paste0("X", seq_len(n)))

  nz <- function(x) { x <- as.character(x); trimws(x[!is.na(x) & nzchar(trimws(x))]) }
  P  <- nz(modelParams$Parameters)
  param_names <- unique(sub("^\\*?\\s*([A-Za-z.][A-Za-z0-9_.]*)\\s*=.*", "\\1", P))
  Fn <- nz(modelParams$Functions)
  func_names <- unique(sub("^\\s*([A-Za-z.][A-Za-z0-9_.]*)\\s*<-.*", "\\1", Fn))
  allowed <- unique(c(param_names, paste0(param_names, "_0"), state_names,
                      func_names, "time", "t", "N", "pi"))

  # Column groups: a missing/misspelled group (e.g. 'Prameters') otherwise reads
  # as empty and fails cryptically later. Warn on any unrecognised column; a
  # missing States column is fatal.
  # Rate columns may be suffixed by index (Linear1..) or compartment name (LinearS..).
  lin_quad_ok <- grepl("^(Linear|Quadratic)[0-9]+$", names(modelParams)) |
    names(modelParams) %in% c(paste0("Linear", state_names),
                              paste0("Quadratic", state_names))
  recognised <- grepl("^_", names(modelParams)) | lin_quad_ok |
    names(modelParams) %in% c("Others", "States", "Functions", "Parameters", "Conditions", "Constant")
  if (any(!recognised))
    warning(sprintf(
      paste0("modelParams has unrecognised column(s): %s -- a typo? Expected ",
             "States / Parameters / Others / Functions / Conditions / Constant / ",
             "Linear<j> / Quadratic<j> and a '_'-prefixed index column."),
      paste(names(modelParams)[!recognised], collapse = ", ")), call. = FALSE)
  if (!"States" %in% names(modelParams))     add("missing 'States' column.")
  if (!"Parameters" %in% names(modelParams))
    warning("modelParams has no 'Parameters' column.", call. = FALSE)

  # Reserved names: Julia keywords, codegen scalar/temporary variables, and the
  # per-term helper functions the generator emits. Any of these used as a
  # parameter / STATE / function name would shadow an internal and silently break
  # the generated model, so they are rejected up front. States are checked too --
  # they also become codegen variables -- closing a gap where a state named `end`
  # or `t` slipped past. The generated names guarded against:
  #   N1, N2, ...            level populations (N0 stays free -- it is a common
  #                          fixed initial-population parameter)
  #   total_pop              sum of level populations
  #   f<ij> / g<ij>          first/second-order time-varying coefficient functions
  #   cst<i>                 Constant-column time-functions
  #   secOrd_<i>_<j>         second-order term temporaries
  #   time                   rewritten to the codegen time variable `t`
  .codegen_reserved <- function(nm) {
    nm %in% c(.JULIA_RESERVED, "time", "total_pop") |
      grepl("^N[1-9][0-9]*$", nm) |
      grepl("^[fg][0-9]{2,}$", nm) |
      grepl("^cst[0-9]+$", nm) |
      grepl("^secOrd_[0-9]+_[0-9]+$", nm)
  }
  nm_all <- unique(c(param_names, state_names, func_names))
  reserved_hit <- unique(nm_all[.codegen_reserved(nm_all)])
  if (length(reserved_hit))
    add(paste0("reserved name(s) %s -- these collide with a Julia keyword or an ",
               "internal codegen variable/function (t, p, X, du, dX, parms, ",
               "N<level>, total_pop, f<ij>, g<ij>, cst<i>, secOrd_i_j) and would ",
               "break the generated model; rename them."),
        paste(sQuote(reserved_hit), collapse = ", "))

  # A parameter/state whose name is ALSO another declared quantity's `_0`
  # initial-value alias (declaring both `beta` and `beta_0`, say) shadows that
  # alias. Only a genuine pair -- `<name>` and `<name>_0` both declared -- clashes;
  # a lone `_0` reference inside an expression is fine.
  has0        <- grep("_0$", nm_all, value = TRUE)
  alias_clash <- unique(has0[sub("_0$", "", has0) %in% nm_all])
  if (length(alias_clash))
    add(paste0("name(s) %s clash with the auto-generated `_0` initial-value alias ",
               "of another declared quantity; rename them."),
        paste(sQuote(alias_clash), collapse = ", "))

  is_num <- function(x) !is.na(suppressWarnings(as.numeric(x)))

  # A coefficient / value expression must be a number or a parseable expression
  # whose free variables are all declared (or resolvable base-R objects).
  check_expr <- function(expr, where) {
    expr <- trimws(expr)
    if (!nzchar(expr) || expr == "0" || is_num(expr)) return(invisible())
    vars <- tryCatch(all.vars(parse(text = expr)), error = function(e) NULL)
    if (is.null(vars)) {
      add("%s: '%s' is not a number or a valid expression.", where, expr)
      return(invisible())
    }
    bad <- setdiff(vars, allowed)
    # Allow a resolvable base *value* (e.g. `pi`), but NOT a base *function* used
    # as a variable: beta(), gamma(), c(), t() exist as functions and would
    # otherwise mask an undeclared parameter named beta/gamma/c/t.
    bad <- bad[!vapply(bad, function(v) {
      o <- get0(v, inherits = TRUE); !is.null(o) && !is.function(o)
    }, logical(1))]
    if (length(bad))
      add("%s: unknown symbol(s) %s in '%s' -- not a declared parameter/state/function.",
          where, paste(sprintf("'%s'", bad), collapse = ", "), expr)
    invisible()
  }
  check_box <- function(rhs, where) {
    inner <- sub("^\\[\\s*(.*?)\\s*\\]$", "\\1", sub("\\|.*$", "", trimws(rhs)))
    parts <- strsplit(inner, ",", fixed = TRUE)[[1]]
    nums  <- suppressWarnings(as.numeric(trimws(parts)))
    if (length(nums) != 2 || any(is.na(nums)))
      add("%s: box '%s' must be two numbers '[lo,hi]'.", where, rhs)
    else if (!(nums[1] < nums[2]))
      add("%s: box '%s' needs lower < upper.", where, rhs)
    invisible()
  }

  ## ---- States ----
  for (s in nz(modelParams$States)) {
    rhs <- sub("^[^=]*=", "", s)
    if (grepl("^\\*", s))                         check_expr(rhs, sprintf("States cell '%s'", s))
    else if (grepl("^\\[", trimws(rhs)))          check_box(rhs, sprintf("States cell '%s'", s))
  }
  ## ---- Parameters ----
  for (p in P) {
    rhs <- sub("^[^=]*=", "", p)
    if (grepl("^\\*", p)) { check_expr(rhs, sprintf("Parameters cell '%s'", p)); next }
    spec <- tryCatch(suppressWarnings(parsePrior(rhs)), error = function(e) e)
    if (inherits(spec, "error"))                  check_expr(rhs, sprintf("Parameters cell '%s'", p))
    else if (identical(spec$dist, "Uniform"))     check_box(rhs, sprintf("Parameters cell '%s'", p))
    else if (any(is.na(spec$args)))               add("Parameters cell '%s': prior has non-numeric argument(s).", p)
  }
  ## ---- Functions ----
  for (f in Fn) check_expr(sub("^[^<]*<-", "", f), sprintf("Functions cell '%s'", f))
  ## ---- Conditions ---- (comparators replaced so the expression parses)
  for (cnd in nz(modelParams$Conditions))
    check_expr(gsub("[<>=]+", "-", cnd), sprintf("Conditions cell '%s'", cnd))
  ## ---- Constant ---- (one value per compartment; $-prefix = time-function)
  if ("Constant" %in% names(modelParams))
    for (v in as.character(modelParams$Constant)[seq_len(n)]) {
      v <- trimws(v)
      if (is.na(v) || !nzchar(v) || v == "0") next
      check_expr(sub("^\\$", "", v), sprintf("Constant cell '%s'", v))
    }

  ## ---- Linear<j> / Quadratic<j> ---- (j = index, or LinearS.. by comp name)
  resolve_col <- function(prefix, i) {
    cand <- paste0(prefix, c(i, state_names[i]))
    hit  <- cand[cand %in% names(modelParams)]
    if (length(hit)) hit[1] else cand[1]
  }
  for (j in seq_len(n)) {
    lc <- resolve_col("Linear", j); qc <- resolve_col("Quadratic", j)
    if (lc %in% names(modelParams))
      for (v in as.character(modelParams[[lc]])[seq_len(n)])
        if (!is.na(v)) check_expr(v, sprintf("%s cell '%s'", lc, v))
    if (qc %in% names(modelParams))
      for (v in as.character(modelParams[[qc]])[seq_len(n * n)]) {
        v <- trimws(v)
        if (is.na(v) || !nzchar(v) || v == "0") next
        mm <- regmatches(v, regexec("^\\*(\\d+)\\*(.+)$", v))[[1]]
        if (length(mm) != 3L) { add("%s cell '%s' must be '*goto*coeff' or 0.", qc, v); next }
        goto <- suppressWarnings(as.integer(mm[2]))
        if (is.na(goto) || goto < 1L || goto > n)
          add("%s cell '%s': target compartment %s is out of range 1..%d.", qc, v, mm[2], n)
        check_expr(mm[3], sprintf("%s cell '%s' coefficient", qc, v))
      }
  }
  ## ---- Others (time grid) ----
  o <- nz(modelParams$Others)
  for (key in c("startpoint", "endpoint", "partition")) {
    hit <- grep(sprintf("^%s=", key), o, value = TRUE)
    if (length(hit) == 0L) { add("Others: missing required '%s='.", key); next }
    if (!is.finite(suppressWarnings(as.numeric(sub("^[^=]*=", "", hit[1])))))
      add("Others: '%s' must be numeric (got '%s').", key, hit[1])
  }

  if (length(errs))
    stop("Invalid modelParams:\n  - ", paste(errs, collapse = "\n  - "), call. = FALSE)
  invisible(TRUE)
}
