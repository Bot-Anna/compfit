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
#'     \code{Quadratic}/\code{Constant}/\code{Functions}/\code{Conditions}/
#'     \code{Mixing_<level>}/\code{Pool_<level>} expressions reference only
#'     declared symbols (parameters, \code{_0} aliases, the declared states --
#'     \code{X1..Xn} or named, e.g. \code{S}/\code{I}/\code{R}, which are rewritten
#'     to their state slot --, functions, \code{time}, \code{N<level>}), so a typo
#'     like \code{*2*-bta} is caught; a parameter/function may not share a State
#'     name;
#'   \item box priors \code{[lo,hi]} are two numbers with \code{lo < hi};
#'   \item distribution priors have numeric arguments;
#'   \item no parameter/state/function name collides with a Julia keyword or a
#'     codegen variable (\code{t}/\code{p}/\code{X}/\code{du}/\code{dX}/
#'     \code{parms}, \code{N1..}, \code{Nw1..}, \code{total_pop},
#'     \code{secOrd_i_j}, or another quantity's \code{_0} alias);
#'   \item every declared name (Parameters/Functions/States) is a legal identifier
#'     -- letters, digits, \code{_}, \code{.} -- so operator characters (\code{^},
#'     \code{()}, \code{,}, \code{-}, \code{/}, ...) that would be parsed as
#'     operators rather than a name are rejected;
#'   \item \code{Quadratic} cells are \code{*goto*coeff} with a target given as a
#'     compartment index (\code{1..n}) or a State name (e.g. \code{*I*});
#'   \item each \code{Mixing_<level>} / \code{Pool_<level>} column has a matching
#'     \code{_<level>} index column (\code{Mixing} holds non-negative per-compartment
#'     weights; \code{Pool} holds a single pool expression, and the two are mutually
#'     exclusive for a level);
#'   \item \code{Others} has numeric \code{startpoint}/\code{endpoint}/\code{partition};
#'   \item the \code{Functions} column can be put in a working evaluation order.
#'     Row order itself does not matter -- a Function may be written above the
#'     entries it references, and the builder sorts the column into dependency
#'     order -- but a name \strong{defined more than once}, a
#'     \strong{self-reference} (\code{x <- x + 1}; a Function is a definition, not
#'     an update), and a \strong{reference cycle} have no valid order and are
#'     rejected, the cycle reported as the loop it forms.
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
         "from the 'Level_<id>' index column(s).", call. = FALSE)
  n <- cs$number_of_comps
  state_names <- .compartments(modelParams)          # canonical names (X1..Xn or S,I,R,...)
  if (length(state_names) < n) state_names <- union(state_names, paste0("X", seq_len(n)))

  nz <- function(x) { x <- as.character(x); trimws(x[!is.na(x) & nzchar(trimws(x))]) }
  P  <- nz(modelParams$Parameters)
  param_names <- unique(sub("^\\*?\\s*([A-Za-z.][A-Za-z0-9_.]*)\\s*=.*", "\\1", P))
  Fn <- nz(modelParams$Functions)
  func_names <- unique(sub("^\\s*([A-Za-z.][A-Za-z0-9_.]*)\\s*<-.*", "\\1", Fn))
  # Level head counts N1..Nk, their named aliases N_<level>, and total_pop are
  # codegen quantities in scope for any coefficient / Pool_<level> expression, so
  # allow them as referenceable symbols. `startpoint`/`cutoff` are the time-grid
  # constants the generator injects as literals (Julia/Stan) / closes over (R), so
  # an expression may use them too -- e.g. a calendar-year step `time+startpoint`.
  n_levels     <- max(1L, length(cs$compartment_cols))
  n_alias_names <- if (length(cs$compartment_cols))
                     paste0("N_", sub("^Level_", "", cs$compartment_cols, ignore.case = TRUE)) else character(0)
  allowed <- unique(c(param_names, paste0(param_names, "_0"), state_names,
                      func_names, "time", "t", "N", "pi", "startpoint", "cutoff",
                      paste0("N", seq_len(n_levels)), n_alias_names, "total_pop"))

  # A Function may be used in an INITIAL-STATE expression only if it is state-
  # INDEPENDENT (its value at t=-1 is fixed by parameters + time). A state-
  # DEPENDENT Function (omega/p_*/anything transitively using a compartment or a
  # level head N_<level>) is circular at init. Classify so a State cell that
  # references a dependent Function is rejected with a pointed message.
  state_syms_init <- unique(c(state_names, "total_pop",
                              paste0("N", seq_len(n_levels)), n_alias_names,
                              if (length(cs$compartment_cols))
                                paste0("Nw", seq_along(cs$compartment_cols))))
  fn_dependent <- .classify_functions(Fn, state_syms_init, param_names)$dependent

  # Column groups: a missing/misspelled group (e.g. 'Prameters') otherwise reads
  # as empty and fails cryptically later. Warn on any unrecognised column; a
  # missing States column is fatal.
  # Rate columns may be suffixed by index (Linear1..) or compartment name (LinearS..).
  lin_quad_ok <- grepl("^(Linear|Quadratic)[0-9]+$", names(modelParams)) |
    names(modelParams) %in% c(paste0("Linear", state_names),
                              paste0("Quadratic", state_names))
  recognised <- grepl("^Level_", names(modelParams), ignore.case = TRUE) | lin_quad_ok |
    grepl("^(Mixing|Pool)_", names(modelParams)) |
    names(modelParams) %in% c("Others", "States", "Functions", "Parameters", "Conditions", "Constant")
  if (any(!recognised))
    warning(sprintf(
      paste0("modelParams has unrecognised column(s): %s -- a typo? Expected ",
             "States / Parameters / Others / Functions / Conditions / Constant / ",
             "Linear<j> / Quadratic<j> / Mixing_<level> / Pool_<level> and a ",
             "'Level_<id>' index column."),
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
  #   Nw1, Nw2, ...          weighted mixing pools (Mixing_<level> denominators)
  #   total_pop              sum of level populations
  #   secOrd_<i>_<j>         second-order term temporaries
  #   time                   rewritten to the codegen time variable `t`
  .codegen_reserved <- function(nm) {
    nm %in% c(.JULIA_RESERVED, "time", "total_pop") |
      grepl("^N[1-9][0-9]*$", nm) |
      grepl("^Nw[1-9][0-9]*$", nm) |
      grepl("^secOrd_[0-9]+_[0-9]+$", nm)
  }
  nm_all <- unique(c(param_names, state_names, func_names))
  reserved_hit <- unique(nm_all[.codegen_reserved(nm_all)])
  if (length(reserved_hit))
    add(paste0("reserved name(s) %s -- these collide with a Julia keyword or an ",
               "internal codegen variable (t, p, X, du, dX, parms, N<level>, ",
               "Nw<level>, total_pop, secOrd_i_j) and would break the generated ",
               "model; rename them."),
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

  # A parameter/function may not share a State (compartment) name: a compartment
  # name is rewritten to its state slot X[k] in every coefficient/Functions/Pool
  # expression, so such a token always resolves to the state -- a like-named
  # parameter/function would be silently shadowed.
  name_state_clash <- unique(intersect(c(param_names, func_names), state_names))
  if (length(name_state_clash))
    add(paste0("name(s) %s are used both as a State (compartment) and as a ",
               "parameter/function; a compartment name resolves to its state in ",
               "expressions, so rename the parameter/function."),
        paste(sQuote(name_state_clash), collapse = ", "))

  # A parameter/state/function may not take a generated level head-count alias
  # name N_<level> -- it would be shadowed by (or shadow) that codegen variable.
  alias_clash2 <- unique(intersect(nm_all, n_alias_names))
  if (length(alias_clash2))
    add(paste0("name(s) %s collide with the level head-count alias N_<level>; ",
               "rename them."),
        paste(sQuote(alias_clash2), collapse = ", "))

  # A declared name (Parameters / Functions / States) must be a legal identifier:
  # letters, digits, `_`, `.`, starting with a letter or dot. Operator characters
  # -- `^` `(` `)` `,` `-` `/` `*` `+` spaces, ... -- do NOT work as names, on any
  # backend: an expression parses them as operators, so e.g. `f^R_CH` reads as
  # `f ^ R_CH` (two names and a power), never one name. Caught here so such a sheet
  # fails with a clear message instead of silently mis-parsing deep in codegen.
  raw_decl_names <- trimws(c(
    sub("=.*$",  "", sub("^\\*", "", P)),                      # Parameters: name before '='
    sub("(<-|=).*$", "", Fn),                                  # Functions:  name before '<-'/'='
    sub("=.*$",  "", sub("^\\*", "", nz(modelParams$States)))  # States:     name before '='
  ))
  raw_decl_names <- unique(raw_decl_names[nzchar(raw_decl_names)])
  illegal_names  <- raw_decl_names[!grepl("^[A-Za-z.][A-Za-z0-9_.]*$", raw_decl_names)]
  if (length(illegal_names))
    add(paste0("name(s) %s are not valid names -- use only letters, digits and '_' ",
               "(starting with a letter). Operator characters such as ^ ( ) , - / * ",
               "are parsed as operators, not part of a name, so they cannot be used."),
        paste(sQuote(illegal_names), collapse = ", "))

  is_num <- function(x) !is.na(suppressWarnings(as.numeric(x)))

  # A coefficient / value expression must be a number or a parseable expression
  # whose free variables are all declared (or resolvable base-R objects).
  check_expr <- function(expr, where) {
    expr <- trimws(expr)
    if (grepl("^\\$", expr)) {
      add(paste0("%s: inline '$' time-functions have been removed -- define the ",
                 "expression in the Functions column and reference it by name."),
          where)
      return(invisible())
    }
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

  check_state_fun_dep <- function(rhs, where) {
    if (!length(fn_dependent)) return(invisible())
    vars <- tryCatch(all.vars(parse(text = rhs)), error = function(e) character(0))
    hit  <- intersect(vars, fn_dependent)
    if (length(hit))
      add(paste0("%s references state-dependent Function(s) %s, which cannot be ",
                 "evaluated in an initial-state expression -- they depend on the ",
                 "compartments being initialised. Reference only parameters or ",
                 "state-independent Functions (sigma/q/sigmoid/...), or inline the value."),
          where, paste(sprintf("'%s'", hit), collapse = ", "))
    invisible()
  }

  ## ---- States ----
  for (s in nz(modelParams$States)) {
    rhs <- sub("^[^=]*=", "", s)
    if (grepl("^\\*", s)) {
      check_expr(rhs, sprintf("States cell '%s'", s))
      check_state_fun_dep(rhs, sprintf("States cell '%s'", s))
    } else if (grepl("^\\[", trimws(rhs)))        check_box(rhs, sprintf("States cell '%s'", s))
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
  # Functions are reordered into dependency order before codegen, so writing one
  # above the entries it uses is fine. What CANNOT be resolved by reordering is a
  # duplicate definition, a self-reference, or a reference cycle -- reported here
  # so they surface with the rest of the sheet errors, before any backend runs.
  for (e in .function_order_errors(.function_order(Fn))) add("%s", e)
  ## ---- Conditions ---- (comparators replaced so the expression parses)
  if ("Conditions" %in% names(modelParams))
    for (cnd in nz(modelParams$Conditions))
      check_expr(gsub("[<>=]+", "-", cnd), sprintf("Conditions cell '%s'", cnd))
  ## ---- Constant ---- (one value per compartment: number/parameter/expression)
  if ("Constant" %in% names(modelParams))
    for (v in as.character(modelParams$Constant)[seq_len(n)]) {
      v <- trimws(v)
      if (is.na(v) || !nzchar(v) || v == "0") next
      check_expr(v, sprintf("Constant cell '%s'", v))
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
        mm <- regmatches(v, regexec("^\\*([^*]+)\\*(.+)$", v))[[1]]
        if (length(mm) != 3L) { add("%s cell '%s' must be '*goto*coeff' or 0.", qc, v); next }
        goto <- .comp_index(mm[2], state_names)   # target by index (*5*) or name (*I*)
        if (is.na(goto) || goto < 1L || goto > n)
          add("%s cell '%s': target compartment '%s' is not a valid compartment (1..%d or a State name).",
              qc, v, mm[2], n)
        check_expr(mm[3], sprintf("%s cell '%s' coefficient", qc, v))
      }
  }
  ## ---- Mixing_<level> ---- (full-pool denominator weights; one per compartment)
  # Each needs a matching '_<level>' index column; cells are non-negative
  # numbers or coefficient expressions (blank -> membership default).
  mix_cols  <- grep("^Mixing_", names(modelParams), value = TRUE)
  lvl_names <- sub("^Level_", "", cs$compartment_cols, ignore.case = TRUE)
  for (mc in mix_cols) {
    lname <- sub("^Mixing_", "", mc)
    if (!(lname %in% lvl_names))
      add("column '%s' has no matching level column '_%s'.", mc, lname)
    for (v in as.character(modelParams[[mc]])[seq_len(n)]) {
      v <- trimws(v)
      if (is.na(v) || !nzchar(v)) next
      if (is_num(v)) {
        if (as.numeric(v) < 0)
          add("%s cell '%s': mixing weight must be >= 0.", mc, v)
      } else {
        check_expr(v, sprintf("%s cell '%s'", mc, v))
      }
    }
  }

  ## ---- Pool_<level> ---- (a whole mixing pool as one expression of the levels)
  # Function of the level head counts N1..Nk / total_pop (+ params / functions /
  # time); mutually exclusive with Mixing_<level>.
  pool_cols <- grep("^Pool_", names(modelParams), value = TRUE)
  for (pc in pool_cols) {
    lname <- sub("^Pool_", "", pc)
    if (!(lname %in% lvl_names))
      add("column '%s' has no matching level column '_%s'.", pc, lname)
    if (paste0("Mixing_", lname) %in% names(modelParams))
      add(paste0("level '%s' has both a 'Mixing_%s' and a 'Pool_%s' column -- use ",
                 "one (per-compartment weights OR a single pool expression)."),
          lname, lname, lname)
    pv <- trimws(as.character(modelParams[[pc]]))
    pv <- pv[!is.na(pv) & nzchar(pv)]
    if (length(pv)) check_expr(pv[1], sprintf("%s cell '%s'", pc, pv[1]))
  }
  # Note: a Pool expression is the transmission denominator, so it must stay
  # positive over the whole solve; it is floored at a small positive value to
  # avoid divide-by-zero, but a pool that reaches 0 or goes negative is a
  # modelling error (spurious transmission spikes). Informational, not an error.
  if (length(pool_cols))
    message("validate_modelParams: Pool_<level> denominators are floored at a ",
            "small positive value to guard against divide-by-zero; keep each pool ",
            "expression strictly positive across the time span to avoid spurious ",
            "transmission when it approaches the floor.")

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
