# ============================================================
# utils.R
# Pure, reusable helper functions with no side effects.
# Sourced once; used by the model-building subscripts and by
# fitCompartmentalModel().
# ============================================================

# Null-coalescing operator: `a %||% b` is `a` unless it is NULL, then `b`.
# Package-internal; handy for defaulting optional control-list fields.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Parse the Conditions column into "violation" expression strings: each is >= 0
# exactly when the constraint is VIOLATED, so a soft penalty on max(0, viol)
# enforces it. `beta>gamma` -> "(gamma)-(beta)". Used by both the MLE loss
# (already) and the Bayesian samplers, which add -.COND_PENALTY * max(0,viol)^2
# to the log-posterior. Comparators: >, >=, <, <=.
.COND_PENALTY <- 1000
.parse_conditions <- function(modelParams) {
  cv <- as.character(modelParams$Conditions)
  cv <- gsub(" ", "", cv[!is.na(cv) & trimws(cv) != ""])
  out <- character(0)
  for (cnd in cv) {
    op <- if (grepl(">=", cnd)) ">=" else if (grepl("<=", cnd)) "<=" else
          if (grepl(">",  cnd)) ">"  else if (grepl("<",  cnd)) "<"  else NA_character_
    if (is.na(op)) next
    ab <- strsplit(cnd, op, fixed = TRUE)[[1]]
    if (length(ab) != 2L || !all(nzchar(ab))) next
    out <- c(out, if (op %in% c(">", ">="))         # want lhs > rhs
                    sprintf("(%s)-(%s)", ab[2], ab[1])   # violated if rhs-lhs > 0
                  else sprintf("(%s)-(%s)", ab[1], ab[2]))  # want lhs < rhs
  }
  out
}

# Functions the Julia/Stan observable translators can carry over verbatim:
# arithmetic operators and elementwise math common to R, Stan and Julia. NOT
# max/min/pmax/pmin/sum/ifelse (they are scalar-reducing or non-existent in Stan,
# so they would silently mistranslate a vectorised observable).
.TRANSLATABLE_FUNS <- c("+", "-", "*", "/", "^", "(",
  "exp", "log", "log2", "log10", "sqrt", "abs", "expm1", "log1p",
  "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh")

# Error clearly when a data-stream Formula uses a function the codegen cannot
# translate, instead of emitting subtly wrong Stan/Julia. Only the elementwise
# set above is safe on a trajectory vector.
.assert_translatable <- function(formula_str, where = "Formula") {
  f <- gsub("`", "", as.character(formula_str))
  f <- sub("^(annual|cumulative)\\((.*)\\)$", "\\2", f)   # drop the reducer wrapper
  e <- tryCatch(parse(text = f)[[1]], error = function(err) NULL)
  if (is.null(e)) return(invisible())                    # unparseable -> let it fail downstream
  funs <- character(0)
  walk <- function(x) if (is.call(x)) {
    if (is.symbol(x[[1]])) funs <<- c(funs, as.character(x[[1]]))
    for (a in as.list(x)[-1]) walk(a)
  }
  walk(e)
  bad <- setdiff(unique(funs), .TRANSLATABLE_FUNS)
  if (length(bad))
    stop(sprintf(paste0("%s '%s' uses %s, which the Julia/Stan code generator ",
      "cannot translate. Use only arithmetic and elementwise math (exp, log, ",
      "sqrt, ...); for anything else, precompute it as a model State or Function."),
      where, formula_str, paste(sQuote(bad), collapse = ", ")), call. = FALSE)
  invisible()
}

# ---- Timing ----
# Custom toc() function. Reads .start_time from the calling/global scope.
toc_fmt <- function() {
  secs <- as.numeric((proc.time() - .start_time)[["elapsed"]])
  if (secs >= 3600) {
    cat(round(secs / 3600, 2), "hours elapsed\n")
  } else if (secs >= 60) {
    cat(round(secs / 60, 2), "min elapsed\n")
  } else {
    cat(round(secs, 3), "sec elapsed\n")
  }
  beepr::beep(1)
}

# ---- Parameter / string parsing ----

# Get the value of a named parameter (e.g. "startpoint=2013") from a vector of
# "name=value" strings.
extract_param <- function(x, name) {
  as.numeric(sub(paste0(".*", name, "=([0-9.]+).*"), "\\1", x[grep(name, x)]))
}

# Determine when the cutoff year takes place -- time is 0, 1, 2, ... not years.
transform_time <- function(x, startpoint) {
  return(x + 1 - startpoint)
}

# Nearest step value, vectorised with recycling.
snap_to_step <- function(min, max, value, step) {
  n <- max(length(min), length(max), length(value), length(step))
  min   <- rep(min,   length.out = n)
  max   <- rep(max,   length.out = n)
  value <- rep(value, length.out = n)
  step  <- rep(step,  length.out = n)
  snapped <- round((value - min) / step) * step + min
  snapped <- pmin(pmax(snapped, min), max)
  return(snapped)
}

# Insert an expression into an existing function body.
# (Convert a string into an expression first via parse(text=...)!)
funins <- function(f, expr, after) {
  body(f) <- as.call(append(as.list(body(f)), expr, after = after))
  f
}

# Remove a trailing "+" (and any trailing spaces before it).
remove_trailing_plus <- function(input_string) {
  while (grepl("\\ $", input_string)) {
    input_string <- sub("\\ $", "", input_string)
  }
  if (grepl("\\+$", input_string)) {
    input_string <- sub("\\+$", "", input_string)
  }
  return(input_string)
}

# Leading asterisk: indicates the term also goes to another equation.
has_leading_asterisk <- function(text) {
  substr(text, 1, 1) == "*"
}

# Leading dollar sign: indicates the value is a parameter to be fitted.
has_leading_dollar_sign <- function(text) {
  substr(text, 1, 1) == "$"
}

# A data stream is "cumulative" iff its formula is wrapped in cumulative(...).
# Used by BOTH the loss (to difference the model output to annual increments
# for FITTING) and .prepare_data (to difference the data the same way), so they
# must share one definition. Plotting is unaffected: evaluate_formula() still
# renders cumulative() as a running total.
is_cumulative_stream <- function(formula_str) {
  grepl("^cumulative\\(", gsub("`", "", formula_str))
}

# Classify one data-value cell into observed / censored / missing.
# Grammar (case-insensitive, whitespace-tolerant):
#   "x", "", NA            -> missing
#   "<=5", "<= 5"          -> left-censored  (upper bound), limit 5, inclusive
#   "<5",  "< 5"           -> left-censored  (upper bound), limit 5, strict
#   ">=5", ">= 5"          -> right-censored (lower bound), limit 5, inclusive
#   ">5",  "> 5"           -> right-censored (lower bound), limit 5, strict
#   "0", "142", "-3.1"     -> observed value
# Returns list(kind, value, limit, inclusive, bound).
#   kind     : "observed" | "censored" | "missing"
#   value    : numeric observed value (NA unless observed)
#   limit    : numeric detection / bound limit (NA unless censored)
#   inclusive: TRUE for "<=" / ">=", FALSE for "<" / ">" (NA unless censored)
#   bound    : "upper" (left-censored, truth <= L) | "lower" (right-censored,
#              truth >= L) | NA
# NOTE: two-character operators (<=, >=) MUST be tested before one-character.
parse_data_cell <- function(cell) {
  # Uniform return shape across every kind. New fields upper/dev/dir default to
  # NA so existing consumers of kind/value/limit/inclusive/bound are unaffected.
  mk <- function(kind, value = NA_real_, limit = NA_real_, inclusive = NA,
                 bound = NA_character_, upper = NA_real_, dev = NA_real_,
                 dir = NA_real_, dev2 = NA_real_)
    list(kind = kind, value = value, limit = limit, inclusive = inclusive,
         bound = bound, upper = upper, dev = dev, dir = dir, dev2 = dev2)

  s <- trimws(as.character(cell))
  if (is.na(cell) || s == "" || tolower(s) == "x") return(mk("missing"))

  num1 <- function(x, what) {
    v <- suppressWarnings(as.numeric(trimws(x)))
    if (is.na(v)) stop(sprintf("Unparseable %s in data cell: '%s'", what, s))
    v
  }

  # Precedence: most specific / longest first.
  # 1. INTERVAL  [A,B]        -- truth in [A,B], flat within (hard edges, CDF).
  #    SOFT INTERVAL [A,B]~s  -- plateau [A,B] with exponential soft shoulders of
  #    scale s (or ~sl,su for asymmetric shoulders): a family-general soft prior
  #    on the model mean, log g(mu), evaluated WITHOUT the CDF -- so it also
  #    works for the discrete families under NUTS. Requires A < B.
  if (grepl("^\\[", s)) {
    body <- s; soft <- NA_character_
    if (grepl("\\]~", s)) {                      # split off a ~shoulder suffix
      sp <- strsplit(s, "~", fixed = TRUE)[[1]]
      if (length(sp) != 2L) stop(sprintf("Malformed soft-interval cell: '%s'", s))
      body <- sp[1]; soft <- sp[2]
    }
    inner <- sub("^\\[\\s*(.*?)\\s*\\]$", "\\1", body)
    if (identical(inner, body)) stop(sprintf("Malformed interval cell: '%s'", s))
    ab <- strsplit(inner, ",", fixed = TRUE)[[1]]
    if (length(ab) != 2L) stop(sprintf("Interval cell must be '[A,B]': '%s'", s))
    A <- num1(ab[1], "interval lower"); B <- num1(ab[2], "interval upper")
    if (!(A < B)) stop(sprintf("Interval cell '[A,B]' requires A < B: '%s'", s))
    if (is.na(soft)) return(mk("interval", limit = A, upper = B))
    sd <- suppressWarnings(as.numeric(strsplit(soft, ",", fixed = TRUE)[[1]]))
    if (length(sd) < 1L || length(sd) > 2L || any(is.na(sd)) || any(sd <= 0))
      stop(sprintf("Soft-interval shoulder(s) '~s' (or '~sl,su') must be positive: '%s'", s))
    dlo <- sd[1]; dhi <- if (length(sd) >= 2L) sd[2] else sd[1]
    return(mk("interval", limit = A, upper = B, dev = dlo, dev2 = dhi))
  }
  # 2. ASYMMETRIC with explicit target  A->B  -- centred A, soft toward B.
  if (grepl("->", s, fixed = TRUE)) {
    ab <- strsplit(s, "->", fixed = TRUE)[[1]]
    if (length(ab) != 2L) stop(sprintf("Asymmetric cell must be 'A->B': '%s'", s))
    A <- num1(ab[1], "asym value"); B <- num1(ab[2], "asym target")
    return(mk("asym", value = A, dev = abs(B - A), dir = sign(B - A)))
  }
  # 3. Existing censoring (<= before <, >= before >).
  if (grepl("^<=", s)) return(mk("censored", limit = num1(sub("^<=", "", s), "censored"), inclusive = TRUE,  bound = "upper"))
  if (grepl("^>=", s)) return(mk("censored", limit = num1(sub("^>=", "", s), "censored"), inclusive = TRUE,  bound = "lower"))
  if (grepl("^<",  s)) return(mk("censored", limit = num1(sub("^<",  "", s), "censored"), inclusive = FALSE, bound = "upper"))
  if (grepl("^>",  s)) return(mk("censored", limit = num1(sub("^>",  "", s), "censored"), inclusive = FALSE, bound = "lower"))
  # 4. ASYMMETRIC with GLOBAL deviation  A+ / A-  (dev resolved later from the stream).
  if (grepl("[+-]$", s)) {
    dir <- if (endsWith(s, "+")) 1 else -1
    A   <- num1(substr(s, 1L, nchar(s) - 1L), "asym value")
    return(mk("asym", value = A, dev = NA_real_, dir = dir))
  }
  # 5. Plain observed number.
  mk("observed", value = num1(s, "data"))
}

# Repeatedly collapse sign/operator patterns until stable.
reduce_expression <- function(vec) {
  while (TRUE) {
    new_vec <- gsub("\\+\\+", "+", vec)
    new_vec <- gsub("\\+-", "-", new_vec)
    new_vec <- gsub("-\\+", "-", new_vec)
    new_vec <- gsub("--", "+", new_vec)
    new_vec <- gsub(" ", "", new_vec)
    if (identical(new_vec, vec)) {
      break
    }
    vec <- new_vec
  }
  return(vec)
}

# Extract the numbers from a "[x,y]" limit specification.
extract_numbers <- function(text) {
  bracket_content <- regmatches(text, regexpr("\\[.*?\\]", text))[[1]]
  matches <- gregexpr("-?\\d*\\.?\\d+", bracket_content)
  numbers <- as.numeric(unlist(regmatches(bracket_content, matches)))
  return(numbers)
}

# ---- Likelihood family registry -------------------------------------------
# Single source of truth used by bayesJulia.R, plots.R, and
# fitCompartmentalModel.R. Adding a new family means editing only this block
# (plus the actual distribution logic in each file).
.DISCRETE_FAMILIES    <- c("poisson", "negbin")
.NEEDS_SIGMA_FAMILIES <- c("gaussian", "lognormal")
.NEEDS_PHI_FAMILIES   <- c("negbin", "betabinom")
.ALL_FAMILIES         <- c("gaussian", "lognormal", "poisson", "negbin", "binomial", "betabinom")

# ---- Bayesian grammar parsers ----------------------------------------------
# These are PURE and used only on the Bayesian path. They never run for MLE,
# so the MLE behaviour is completely unaffected by their presence.

# parsePrior: interpret the right-hand side of a fitted entry (the part after
# "name=") as a PRIOR specification. The MLE path continues to read bounds via
# extract_numbers as before; this is a parallel reader used only for Bayes.
#
# Recognised forms (rhs, after stripping spaces and any "|initial"):
#   [a,b]                 -> Uniform(a, b)                      (the existing box)
#   Normal(m,s)           -> Normal(m, s)                       (unbounded)
#   Normal(m,s)[a,b]      -> truncated(Normal(m,s), a, b)       (prior + box)
#   LogNormal(m,s)        -> LogNormal(m, s)
#   Beta(a,b)             -> Beta(a, b)
#   Gamma(a,b)            -> Gamma(a, b)
#   StudentT(nu,m,s)      -> location-scale Student-t: m + s * TDist(nu)
#                            (heavy-tailed; nu = df, m = location, s = scale)
#   <plain number>        -> fixed (not estimated; handled upstream, but we
#                            return kind="fixed" defensively)
#
# Returns a list with at least: kind ("estimated"/"fixed"), dist (Julia
# distribution name or "Uniform"), args (numeric vector), lower, upper
# (truncation bounds; -Inf/Inf if none), and start (a sensible init value).
parsePrior <- function(rhs) {
  rhs <- gsub("\\s", "", rhs)
  rhs <- sub("\\|.*$", "", rhs)   # drop any "|initial" suffix (Bayes ignores it)
  
  # Pure box form: [a,b] -> Uniform(a,b)
  if (grepl("^\\[.*\\]$", rhs)) {
    nums <- extract_numbers(rhs)
    return(list(kind = "estimated", dist = "Uniform",
                args = nums, lower = nums[1], upper = nums[2],
                start = midpoint(nums[1], nums[2])))
  }
  
  # Distribution form, optionally followed by a [a,b] truncation
  m <- regexec("^([A-Za-z]+)\\(([^)]*)\\)(\\[.*\\])?$", rhs)
  g <- regmatches(rhs, m)[[1]]
  if (length(g) >= 3 && nzchar(g[2])) {
    dist <- g[2]
    args <- as.numeric(strsplit(g[3], ",")[[1]])
    trunc <- if (length(g) >= 4 && nzchar(g[4])) extract_numbers(g[4]) else c(-Inf, Inf)
    # StudentT(nu, mu, sigma): location-scale Student-t (df nu, location mu,
    # scale sigma) -- takes THREE positive-df/scale arguments.
    if (dist == "StudentT") {
      if (length(args) != 3L || any(is.na(args)) || args[1] <= 0 || args[3] <= 0)
        stop(sprintf(
          "StudentT prior needs 'StudentT(nu, mu, sigma)' with nu>0, sigma>0: '%s'", rhs))
    }
    # A sensible start: the distribution's "centre" where obvious.
    start <- switch(dist,
                    Normal    = args[1],
                    LogNormal = exp(args[1]),
                    Beta      = args[1] / (args[1] + args[2]),
                    Gamma     = args[1] * args[2],
                    StudentT  = args[2],   # location
                    args[1])
    # Respect truncation for the start value.
    if (is.finite(trunc[1])) start <- max(start, trunc[1])
    if (is.finite(trunc[2])) start <- min(start, trunc[2])
    return(list(kind = "estimated", dist = dist, args = args,
                lower = trunc[1], upper = trunc[2], start = start))
  }
  
  # Plain number => fixed (not estimated).
  num <- suppressWarnings(as.numeric(rhs))
  if (!is.na(num)) {
    return(list(kind = "fixed", value = num))
  }
  
  stop(sprintf("parsePrior: could not interpret prior specification '%s'", rhs))
}

# .compartments: the ordered compartment names. The canonical order and index
# (compartment i = the i-th entry) come from the States column, so compartments
# can be named freely (S, I, R) instead of X1..Xn. The name before '=' is used
# (leading '*' stripped). A sheet written with X1..Xn yields c("X1", ...) -- the
# identity case, so existing sheets behave exactly as before.
.compartments <- function(modelParams) {
  s <- as.character(modelParams$States)
  s <- s[!is.na(s) & trimws(s) != ""]
  s <- gsub(" ", "", s)
  trimws(sub("^\\*?([^=]+)=.*", "\\1", s))
}

# .comp_index: map a compartment reference -- a name from .compartments() OR a
# 1-based integer index -- to its integer index. Vectorised; NA for anything
# unrecognised. Lets _Level and Linear<>/Quadratic<> accept names or numbers.
.comp_index <- function(x, comp_names) {
  x   <- trimws(as.character(x))
  n   <- length(comp_names)
  by_name <- match(x, comp_names)
  by_num  <- suppressWarnings(as.integer(x))
  # Accept a positive integer within range; when no registry is available (legacy
  # sheet with no States column, n == 0) accept any positive integer -- the
  # index IS the compartment there.
  ok_num  <- !is.na(by_num) & by_num >= 1L & (n == 0L | by_num <= n)
  ifelse(is.na(by_name) & ok_num, by_num, by_name)
}

# .default_label: ensure a data frame (dataCombined / dataDummy) has a Label
# column, filling any missing/blank entry from the row's Formula so labels are
# always available for plots/reports. A frame with no Formula column is returned
# unchanged (nothing to default from).
.default_label <- function(df) {
  if (is.null(df) || !("Formula" %in% names(df))) return(df)
  form <- as.character(df$Formula)
  lab  <- if ("Label" %in% names(df)) as.character(df$Label) else rep(NA_character_, length(form))
  miss <- is.na(lab) | !nzchar(trimws(lab))
  lab[miss] <- form[miss]
  df$Label <- lab
  df
}

# .is_prior_entry: TRUE when a fitted Parameters/States entry's right-hand side
# is a PRIOR that parsePrior can interpret as an estimated quantity -- a box
# [lo,hi] or a named distribution (Normal/LogNormal/Beta/Gamma, optionally
# truncated). Anything parsePrior can't read as estimated (an expression of
# other parameters, a bare constant) is NOT a prior and is left to the
# parameter-function path. Used to route entries in statesAndParams().
.is_prior_entry <- function(entry) {
  s <- tryCatch(parsePrior(sub("^[^=]*=", "", entry)), error = function(e) NULL)
  !is.null(s) && identical(s$kind, "estimated")
}

# Names the Bayes observation model reserves for its own parameters: `sigma`
# (Gaussian/log-normal noise scale) and `phi`/`phi<i>` (negbin dispersion).
.cf_reserved_bayes_pattern <- "^sigma$|^phi[0-9]*$"

# .check_reserved_bayes_names: error if any model parameter/state name collides
# with a reserved Bayes hyperparameter name. Called on the method="bayes" path.
.check_reserved_bayes_names <- function(names) {
  names <- names[!is.na(names) & nzchar(names)]
  bad   <- unique(names[grepl(.cf_reserved_bayes_pattern, names)])
  if (length(bad))
    stop(sprintf(
      paste0("Reserved parameter name(s) for a Bayesian fit: %s.\n",
             "The observation model uses 'sigma' (Gaussian/log-normal noise scale) ",
             "and 'phi'/'phi<i>' (negbin dispersion) internally, so a model ",
             "parameter or state cannot share those names. Rename it (e.g. ",
             "'sigma' -> 'sigma_ei', 'rate_ei'); MLE fits are unaffected."),
      paste(sQuote(bad), collapse = ", ")), call. = FALSE)
}

# .mle_box: finite (lower, upper) search box for L-BFGS-B from a parsePrior spec.
# Box priors pass straight through; an unbounded named distribution gets a wide
# but finite range derived from the family (mean +/- 4 sd, or the natural
# support), and any explicit truncation is honoured. This only sets the MLE
# search region -- the Bayesian path applies the actual distribution as the
# prior via buildPriorSpec()/dist_to_julia().
.mle_box <- function(spec) {
  lo <- spec$lower; hi <- spec$upper; a <- spec$args
  rng <- switch(spec$dist,
                Uniform   = c(lo, hi),
                Normal    = c(a[1] - 4 * a[2], a[1] + 4 * a[2]),
                LogNormal = c(exp(a[1] - 4 * a[2]), exp(a[1] + 4 * a[2])),
                Beta      = c(0, 1),
                Gamma     = { m <- a[1] * a[2]; s <- sqrt(a[1]) * a[2]; c(0, m + 4 * s) },
                # StudentT(nu, mu, sigma): mu +/- 6 sigma (wider than Normal for
                # the heavy tails) as a finite MLE search box.
                StudentT  = c(a[2] - 6 * a[3], a[2] + 6 * a[3]),
                c(lo, hi))
  if (is.finite(lo)) rng[1] <- lo    # honour an explicit truncation [lo,hi]
  if (is.finite(hi)) rng[2] <- hi
  if (spec$dist %in% c("LogNormal", "Gamma", "Beta")) rng[1] <- max(rng[1], 0)
  rng
}

# parseLikelihood: interpret one cell of the dataCombined "Likelihood" column
# into a canonical family name + whether it carries a dispersion parameter.
# A blank/NA cell means the DEFAULT: proportional-error Gaussian (i.e. exactly
# the current weighted-least-squares behaviour). Case-insensitive; aliases
# accepted. Only ever called on the Bayesian path.
parseLikelihood <- function(cell) {
  if (is.null(cell) || is.na(cell) || !nzchar(trimws(as.character(cell)))) {
    return(list(family = "gaussian", dispersion = FALSE, asym_dev = NA_real_))  # default
  }
  raw <- trimws(as.character(cell))

  # Optional global asymmetric deviation: "...; asym=<number>" (absolute units).
  # Used by A+/A- data cells that don't carry their own deviation.
  asym_dev <- NA_real_
  m <- regmatches(raw, regexec("asym\\s*=\\s*([0-9.eE+-]+)", raw))[[1]]
  if (length(m) == 2L) {
    asym_dev <- suppressWarnings(as.numeric(m[2]))
    if (is.na(asym_dev) || asym_dev <= 0)
      stop(sprintf("parseLikelihood: 'asym=' must be a positive number in '%s'", cell))
  }

  # Family = the token before any ';' or 'asym=' clause (empty -> gaussian).
  fam_part <- trimws(sub("asym\\s*=.*$", "", sub(";.*$", "", raw)))
  key <- tolower(fam_part)
  if (!nzchar(key)) key <- "gaussian"
  fam <- switch(key,
                gaussian = , normal = , prop = "gaussian",
                lognormal = , lnorm = "lognormal",
                poisson = , pois = "poisson",
                negbin = , nbinom = , nb = "negbin",
                binomial = , binom = "binomial",
                betabinom = , bb = "betabinom",
                stop(sprintf("parseLikelihood: unknown family '%s'", fam_part)))
  list(family = fam,
       dispersion = fam %in% .NEEDS_PHI_FAMILIES,
       asym_dev = asym_dev)
}

# ---- Numeric helpers ----

midpoint <- function(a, b) {
  return((a + b) / 2)
}

normalise <- function(params, min_vals, max_vals) {
  (params - min_vals) / (max_vals - min_vals)
}

denormalise <- function(params_norm, min_vals, max_vals) {
  params_norm * (max_vals - min_vals) + min_vals
}

# ---- Trajectory evaluation for plotting / assessment ----
# Evaluates a data/dummy formula on a solved trajectory. `partition` is an
# explicit argument (used by the annual rolling integral); pass tg$partition.
evaluate_formula <- function(col, sir_out, parms, time, partition) {
  
  # One evaluation environment shared by all formula types: state columns
  # (whatever they are named -- X1, X2, ... or S, I, R, ...), the `parms` vector
  # (so parms["g_HA"] still works), AND each
  # parameter as a BARE variable (so formulas can be written `g_HA*...`). This
  # mirrors how the loss resolves bare names, keeping plotting and fitting
  # consistent.
  build_env <- function() {
    sir_list <- lapply(as.list(sir_out), function(x) unname(as.numeric(x)))
    e <- list2env(c(sir_list, list(parms = parms)))
    if (!is.null(names(parms)))
      for (nm in names(parms)) assign(nm, parms[[nm]], envir = e)
    e
  }
  eval_env <- build_env()
  
  get_flux <- function(expr_str) {
    inner_expr <- parse(text = expr_str)[[1]]
    as.numeric(eval(inner_expr, envir = eval_env))
  }
  
  integrate_rolling <- function(flux) {
    dt <- 1 / partition
    result <- numeric(length(time))
    for (i in seq_along(time)) {
      if (i <= partition) {
        result[i] <- NA
      } else {
        idx <- (i - partition):i
        w <- rep(dt, length(idx))
        w[1] <- dt / 2
        w[length(w)] <- dt / 2
        result[i] <- sum(w * flux[idx])
      }
    }
    result
  }
  
  integrate_cumulative <- function(flux) {
    dt <- diff(time)
    c(0, cumsum(0.5 * (flux[-length(flux)] + flux[-1]) * dt))
  }
  
  if (grepl("^annual\\(", col)) {
    inner_str <- sub("^annual\\((.*)\\)$", "\\1", col)
    integrate_rolling(get_flux(inner_str))
  } else if (grepl("^cumulative\\(", col)) {
    inner_str <- sub("^cumulative\\((.*)\\)$", "\\1", col)
    integrate_cumulative(get_flux(inner_str))
  } else {
    # Stock formula: evaluate in the SAME enriched environment (states + parms +
    # bare param names), not just sir_out -- so stock formulas may also use
    # bare parameter names.
    expr <- parse(text = col)[[1]]
    as.numeric(eval(expr, envir = eval_env))
  }
}