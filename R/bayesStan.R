# ============================================================
# bayesStan.R
# Stan Bayesian backend (solver_control(backend = "stan")). Parallels the Julia
# backend in bayesJulia.R, but emits a self-contained Stan program and samples it
# with rstan. Everything upstream -- the sheet grammar, buildPriorSpec(),
# parseLikelihood(), the data/censoring prep -- is shared; only the model
# codegen + sampler differ.
#
#   buildStanODEFunction(...)  -> the ODE as a Stan `functions{}` definition
#                                 (a third emitter alongside the Julia one and
#                                 the R closure), from the SAME ODE IR.
#   buildStanModel(...)        -> the full .stan program string.
#   .fit_bayes_stan(...)       -> compile + sample via rstan, return the standard
#                                 `samples` contract (draws/summary/...).
#
# This file is sourced always, but its functions run ONLY on
# solver_control(backend = "stan"). Nothing here touches the MLE or Julia paths.
# ============================================================


# ---- R-expression -> Stan-expression helpers -------------------------------
# The ODE IR (vec_main, sir_expression, functions_expression) is in R syntax.
# Stan arithmetic is close to R's: + - * / ^ and exp/log/sqrt/max/min all carry
# over. The differences we translate: `<-` assignment, the if_else/ifelse ternary
# (Stan uses the same `cond ? a : b`), and R's TRUE/FALSE (Stan uses 1/0).

# Extract the three comma-separated args of a parenthesised call whose opening
# "(" is at position `pos`; returns list(arg1, arg2, arg3, pos_after_close) or
# NULL if the parse fails. Depth-aware, so nested parens/commas are respected.
.stan_three_args <- function(s, pos) {
  chars <- strsplit(s, "")[[1]]
  if (chars[pos] != "(") return(NULL)
  depth <- 1L; i <- pos + 1L; arg_i <- 1L; buf <- ""; args <- character(3)
  while (i <= length(chars) && arg_i <= 3L) {
    ch <- chars[i]
    if (ch == "(") { depth <- depth + 1L; buf <- paste0(buf, ch) }
    else if (ch == ")") {
      depth <- depth - 1L
      if (depth == 0L) { args[arg_i] <- trimws(buf); return(list(args[1], args[2], args[3], i + 1L)) }
      buf <- paste0(buf, ch)
    } else if (ch == "," && depth == 1L) { args[arg_i] <- trimws(buf); arg_i <- arg_i + 1L; buf <- "" }
    else buf <- paste0(buf, ch)
    i <- i + 1L
  }
  NULL
}

# Replace the outermost `fname(cond, a, b)` with Stan `(cond ? a : b)`.
.stan_ternary_once <- function(s, fname) {
  m <- regexpr(paste0("\\b", fname, "\\("), s)
  if (m == -1) return(s)
  open <- m[1] + attr(m, "match.length") - 1L         # position of "("
  a <- .stan_three_args(s, open)
  if (is.null(a)) return(s)
  paste0(substr(s, 1, m[1] - 1L),
         "(", a[[1]], " ? ", a[[2]], " : ", a[[3]], ")",
         substr(s, a[[4]], nchar(s)))
}

.r_to_stan <- function(s) {
  it <- 0L
  while (grepl("if_else\\(", s) && it < 200L) { s <- .stan_ternary_once(s, "if_else"); it <- it + 1L }
  it <- 0L
  while (grepl("\\bifelse\\(", s) && it < 200L) { s <- .stan_ternary_once(s, "ifelse"); it <- it + 1L }
  s <- gsub("\\bTRUE\\b", "1", s)
  s <- gsub("\\bFALSE\\b", "0", s)
  s <- gsub("<-", "=", s, fixed = TRUE)
  trimws(s)
}

# Turn a scalar R assignment "lhs = rhs" into a typed Stan declaration
# "real lhs = <rhs>;". Assignments into the dX vector keep their subscript.
.stan_decl_line <- function(line) {
  line <- .r_to_stan(line)
  if (!nzchar(line)) return(NULL)
  if (grepl("^dX\\[", line)) return(paste0("    ", line, ";"))
  lhs <- sub("^([A-Za-z._][A-Za-z0-9._]*)\\s*=.*", "\\1", line)
  if (identical(lhs, line)) return(NULL)   # no '=' -> not an assignment
  paste0("    real ", line, ";")
}

# Split a "\n"-joined block of R assignments into individual typed Stan lines.
.stan_decl_block <- function(block) {
  lines <- strsplit(block, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  out <- vapply(lines, function(l) {
    d <- .stan_decl_line(l)
    if (is.null(d)) NA_character_ else d
  }, character(1))
  out[!is.na(out)]
}


# ------------------------------------------------------------
# buildStanODEFunction: the ODE as a Stan functions{}-block definition.
# Signature mirrors buildJuliaODEFunction() and consumes the SAME IR, so the
# States-column order and the Linear/Quadratic decoding are honoured identically.
# Returns the `vector cf_ode(...) { ... }` definition (buildStanModel wraps it in
# the functions{} block together with annual_integral()).
# ------------------------------------------------------------
buildStanODEFunction <- function(sir_expression,
                                 functions_expression,
                                 vec_help_expressions_second_order,
                                 vec_main,
                                 number_of_comps,
                                 level_compartments,
                                 cutoff,
                                 startpoint) {
  stmts <- character(0)

  # 1. Population sums per level: N1 = X[1]+X[2]+...; total_pop = N1+N2+...
  for (i in seq_along(level_compartments)) {
    idx <- level_compartments[[i]]
    terms <- paste0("X[", idx, "]", collapse = "+")
    stmts <- c(stmts, sprintf("    real N%d = %s;", i, terms))
  }
  all_N <- paste0("N", seq_along(level_compartments), collapse = "+")
  stmts <- c(stmts, sprintf("    real total_pop = %s;", all_N))

  # 2. Parameter unpacking (sir_expression: `beta = p[1]` lines).
  stmts <- c(stmts, .stan_decl_block(sir_expression))

  # 3. Time-varying functions. cutoff/startpoint are R-side constants the ODE
  #    references, injected as literals exactly as the Julia emitter does.
  fe <- functions_expression
  fe <- gsub("\\bcutoff\\b",     as.character(cutoff),     fe)
  fe <- gsub("\\bstartpoint\\b", as.character(startpoint), fe)
  stmts <- c(stmts, .stan_decl_block(fe))

  # 4. Second-order initialisation temporaries.
  for (expr in vec_help_expressions_second_order) {
    if (!is.na(expr) && nzchar(trimws(expr))) {
      d <- .stan_decl_line(expr)
      if (!is.null(d)) stmts <- c(stmts, d)
    }
  }

  # 5. The derivatives themselves, into the declared dX vector.
  dstmts <- sprintf("    dX[%d] = %s;", seq_len(number_of_comps),
                    vapply(vec_main[seq_len(number_of_comps)], .r_to_stan, character(1)))

  paste(c(
    "  vector cf_ode(real t, vector X, array[] real p) {",
    stmts,
    sprintf("    vector[%d] dX;", number_of_comps),
    dstmts,
    "    return dX;",
    "  }"
  ), collapse = "\n")
}


# ---- Compiled-model cache --------------------------------------------------
# rstan::stan_model() runs a ~1-2 min C++ build. Cache compiled models within
# the session, keyed by an md5 of the program text, so refitting the same model
# (new data, more iterations, a tweaked prior that leaves the code unchanged)
# reuses the binary. Keyed on the code, so any change to the program recompiles.
.stan_model_cache <- new.env(parent = emptyenv())

.stan_compile <- function(code, quiet) {
  f <- tempfile(fileext = ".stan"); on.exit(unlink(f))
  writeLines(code, f)
  key <- unname(tools::md5sum(f))
  hit <- get0(key, envir = .stan_model_cache, inherits = FALSE)
  if (!is.null(hit)) return(hit)
  sm <- if (quiet) suppressMessages(rstan::stan_model(model_code = code))
        else rstan::stan_model(model_code = code)
  assign(key, sm, envir = .stan_model_cache)
  sm
}

# ---- Families the Stan backend supports ------------------------------------
# Continuous families use the real data matrix; count families use the integer
# matrices (Stan's discrete *_lpmf need integer support). .fit_bayes_stan()
# rejects anything outside this set with a clear message.
.STAN_FAMILIES <- c("gaussian", "lognormal", "poisson", "negbin")


# ------------------------------------------------------------
# .observable_to_stan: one stream's R observable formula -> a Stan expression
# giving vector[n_years] mu. Mirrors .observable_to_julia: compartments map to
# rows of `traj` (K x n_grid), stock formulas are read at annual snapshot
# indices, annual()/cumulative() fluxes go through annual_integral().
# ------------------------------------------------------------
.observable_to_stan <- function(formula_str, comp_names) {
  .assert_translatable(formula_str, "data-stream Formula")
  f <- gsub("`", "", formula_str)
  to_state <- function(s) {
    for (k in order(-nchar(comp_names)))
      s <- gsub(paste0("\\b\\Q", comp_names[k], "\\E\\b"),
                sprintf("traj[%d, ]", k), s, perl = TRUE)
    s
  }
  if (grepl("^annual\\(", f) || is_cumulative_stream(f)) {
    inner <- sub("^(annual|cumulative)\\((.*)\\)$", "\\2", f)
    sprintf("annual_integral(to_vector(%s), tgrid, partition, n_years)", to_state(inner))
  } else {
    sprintf("to_vector((%s)[annual_idx])", to_state(f))
  }
}


# ------------------------------------------------------------
# .stan_dist: an informative parsePrior spec -> Stan constructor + the parameter
# support bounds + the model `~` line (with T[lo,hi] when truncated).
# ------------------------------------------------------------
.stan_dist <- function(nm, spec) {
  d <- spec$dist; a <- spec$args
  core <- switch(d,
                 Normal    = sprintf("normal(%s, %s)", a[1], a[2]),
                 LogNormal = sprintf("lognormal(%s, %s)", a[1], a[2]),
                 Beta      = sprintf("beta(%s, %s)", a[1], a[2]),
                 Gamma     = sprintf("gamma(%s, %s)", a[1], a[2]),
                 StudentT  = sprintf("student_t(%s, %s, %s)", a[1], a[2], a[3]),
                 stop(sprintf("Prior distribution '%s' not supported in the Stan backend.", d)))
  nat_lo <- if (d %in% c("LogNormal", "Gamma", "Beta")) 0 else -Inf
  nat_hi <- if (d == "Beta") 1 else Inf
  tl <- if (is.finite(spec$lower)) spec$lower else NA_real_
  th <- if (is.finite(spec$upper)) spec$upper else NA_real_
  decl_lo <- max(nat_lo, if (!is.na(tl)) tl else -Inf)
  decl_hi <- min(nat_hi, if (!is.na(th)) th else Inf)
  bnd <- if (is.finite(decl_lo) && is.finite(decl_hi))
           sprintf("<lower=%s, upper=%s>", decl_lo, decl_hi)
         else if (is.finite(decl_lo)) sprintf("<lower=%s>", decl_lo)
         else if (is.finite(decl_hi)) sprintf("<upper=%s>", decl_hi) else ""
  trunc <- if (!is.na(tl) || !is.na(th))
             sprintf(" T[%s, %s]", if (!is.na(tl)) tl else "", if (!is.na(th)) th else "") else ""
  list(bound = bnd, model_line = sprintf("  %s ~ %s%s;", nm, core, trunc))
}


# ------------------------------------------------------------
# .likelihood_to_stan: per-stream masked likelihood block.
#   obs   -> lpdf/lpmf(y | ...);   left-censored -> lcdf(limit | ...);
#   right-censored -> lccdf(llimit | ...);   missing -> skip.
# CONTINUOUS families (gaussian/lognormal) read the real data matrix `y` and the
# real censor limits, scaled by `sigma`. COUNT families (poisson/negbin) read the
# INTEGER matrices `y_int` / `lim_int` / `llim_int` (Stan's discrete *_lpmf /
# *_lcdf need integer support), and negbin adds its per-stream dispersion phi_<i>.
# ------------------------------------------------------------
.likelihood_to_stan <- function(family, i) {
  mu  <- sprintf("mu%d", i)
  mup <- sprintf("fmax(%s[k], 1e-9)", mu)          # nonneg mean

  # Soft interval [A,B]~s: a family-general plateau penalty on the mean (no CDF).
  # Shared by every family.
  soft <- sprintf(paste0("      target += -(fmax(ilow_mat[k, %d] - %s[k], 0) / idev_lo_mat[k, %d]",
                         " + fmax(%s[k] - iupp_mat[k, %d], 0) / idev_hi_mat[k, %d]);"),
                  i, mu, i, mu, i, i)

  # Asymmetric anchor at A: hard one-sided CDF term + soft damping of the mean's
  # excess in the soft direction. `edge` is the CDF argument (int for counts).
  asym_block <- function(fam, argfmt, edge) c(
    sprintf("    else if (asym_mask[k, %d] == 1) {", i),
    sprintf("      if (asym_dir_mat[k, %d] > 0)", i),
    sprintf("        target += %s_lccdf(%s | %s) - fmax(%s[k] - asym_val_mat[k, %d], 0) / asym_dev_mat[k, %d];",
            fam, edge, argfmt, mu, i, i),
    "      else",
    sprintf("        target += %s_lcdf(%s | %s) - fmax(asym_val_mat[k, %d] - %s[k], 0) / asym_dev_mat[k, %d];",
            fam, edge, argfmt, i, mu, i),
    "    }")

  if (family %in% c("poisson", "negbin")) {
    fam  <- if (family == "poisson") "poisson" else "neg_binomial_2"
    args <- if (family == "poisson") mup else sprintf("%s, phi_%d", mup, i)
    # Hard interval on counts: log(P(A<=Y<=B)) on integer edges (ilow_int = A-1).
    hard <- sprintf("      target += log_diff_exp(%s_lcdf(iupp_int[k, %d] | %s), %s_lcdf(ilow_int[k, %d] | %s));",
                    fam, i, args, fam, i, args)
    return(paste(c(
      "  for (k in 1:n_years) {",
      sprintf("    if (obs_mask[k, %d] == 1)", i),
      sprintf("      target += %s_lpmf(y_int[k, %d] | %s);", fam, i, args),
      sprintf("    else if (cens_mask[k, %d] == 1)", i),
      sprintf("      target += %s_lcdf(lim_int[k, %d] | %s);", fam, i, args),
      sprintf("    else if (lcens_mask[k, %d] == 1)", i),
      sprintf("      target += %s_lccdf(llim_int[k, %d] | %s);", fam, i, args),
      sprintf("    else if (interval_mask[k, %d] == 1) {", i),
      sprintf("      if (idev_lo_mat[k, %d] > 0)", i), soft, "      else", hard, "    }",
      asym_block(fam, args, sprintf("asym_val_int[k, %d]", i)),
      "  }"), collapse = "\n"))
  }

  ls <- switch(family,
    gaussian  = list(loc = sprintf("%s[k]", mu),
                     scl = sprintf("sigma * scale_mat[k, %d]", i), fam = "normal"),
    lognormal = list(loc = sprintf("log(%s)", mup), scl = "sigma", fam = "lognormal"),
    stop(sprintf("Likelihood family '%s' not supported in the Stan backend.", family)))
  # Hard interval on a continuous stream: log(CDF(B) - CDF(A)) via log_diff_exp.
  hard <- sprintf("      target += log_diff_exp(%s_lcdf(iupp_mat[k, %d] | %s, %s), %s_lcdf(ilow_mat[k, %d] | %s, %s));",
                  ls$fam, i, ls$loc, ls$scl, ls$fam, i, ls$loc, ls$scl)
  scl2 <- sprintf("%s, %s", ls$loc, ls$scl)
  paste(c(
    "  for (k in 1:n_years) {",
    sprintf("    if (obs_mask[k, %d] == 1)", i),
    sprintf("      target += %s_lpdf(y[k, %d] | %s, %s);", ls$fam, i, ls$loc, ls$scl),
    sprintf("    else if (cens_mask[k, %d] == 1)", i),
    sprintf("      target += %s_lcdf(limit_mat[k, %d] | %s, %s);", ls$fam, i, ls$loc, ls$scl),
    sprintf("    else if (lcens_mask[k, %d] == 1)", i),
    sprintf("      target += %s_lccdf(llimit_mat[k, %d] | %s, %s);", ls$fam, i, ls$loc, ls$scl),
    sprintf("    else if (interval_mask[k, %d] == 1) {", i),
    sprintf("      if (idev_lo_mat[k, %d] > 0)", i), soft, "      else", hard, "    }",
    asym_block(ls$fam, scl2, sprintf("asym_val_mat[k, %d]", i)),
    "  }"), collapse = "\n")
}


# ------------------------------------------------------------
# buildStanModel: assemble the full .stan program.
#   ode_function : the `vector cf_ode(...)` definition (model$stan_code).
# Mirrors buildJuliaBayesModel's structure and ordering contract:
#   p = [fitted params..., fixed params...]; X0 in compartment order.
# Uniform priors are sampled on [0,1] and mapped back (matches posterior_draws);
# informative priors are sampled on the natural scale.
# ------------------------------------------------------------
buildStanModel <- function(prior_spec, like_specs, formulas, n_years, partition,
                           number_of_comps, comp_names, ode_function,
                           sigma_prior_stan = "normal(0, 1)",
                           phi_prior_stan   = "gamma(2, 0.2)",
                           ode_solver       = "rk45") {
  if (is.null(ode_function))
    stop("buildStanModel(): no Stan ODE was generated for this model.")

  ## --- parameters + their priors ---
  param_decls <- character(0)
  tparam_denorm <- character(0)
  model_priors  <- character(0)
  emit_prior <- function(nm, spec) {
    if (spec$dist == "Uniform") {
      lo <- spec$args[1]; hi <- spec$args[2]
      param_decls   <<- c(param_decls, sprintf("  real<lower=0, upper=1> %s_n;", nm))
      tparam_denorm <<- c(tparam_denorm,
                          sprintf("  real %s = %s + %s_n * (%s - %s);", nm, lo, nm, hi, lo))
    } else {
      dd <- .stan_dist(nm, spec)
      param_decls  <<- c(param_decls, sprintf("  real%s %s;", dd$bound, nm))
      model_priors <<- c(model_priors, dd$model_line)
    }
  }
  for (nm in prior_spec$order$params_fitted) emit_prior(nm, prior_spec$params[[nm]])
  for (nm in prior_spec$order$states_fitted) emit_prior(nm, prior_spec$states[[nm]])

  ## --- noise / dispersion hyperparameters ---
  fam_of      <- vapply(like_specs, function(s) s$family, character(1))
  needs_sigma <- any(fam_of %in% c("gaussian", "lognormal"))  # continuous noise scale
  negbin_i    <- which(fam_of == "negbin")                    # per-stream dispersion
  has_counts  <- any(fam_of %in% c("poisson", "negbin"))
  if (needs_sigma) {
    param_decls  <- c(param_decls, "  real<lower=0> sigma;")
    model_priors <- c(model_priors, sprintf("  sigma ~ %s;", sigma_prior_stan))
  }
  for (i in negbin_i) {
    param_decls  <- c(param_decls, sprintf("  real<lower=0> phi_%d;", i))
    model_priors <- c(model_priors, sprintf("  phi_%d ~ %s;", i, phi_prior_stan))
  }

  ## --- definitions: fixed params, _0 aliases, derived states ---
  def_lines <- tparam_denorm
  is_num <- function(s) !is.na(suppressWarnings(as.numeric(s)))
  for (nm in prior_spec$order$params_fixed) {
    val <- prior_spec$fixed_params[[nm]]
    rhs <- if (is_num(val)) val else val   # literal or expression of other names
    def_lines <- c(def_lines, sprintf("  real %s = %s;", nm, rhs),
                              sprintf("  real %s_0 = %s;", nm, nm))
  }
  for (nm in prior_spec$order$params_fitted)
    def_lines <- c(def_lines, sprintf("  real %s_0 = %s;", nm, nm))
  derived_state_names <- character(0)
  for (nm in prior_spec$order$states_fixed) {
    rhs <- prior_spec$fixed_states[[nm]]
    if (!is_num(rhs)) {
      def_lines <- c(def_lines, sprintf("  real %s = %s;", nm, rhs))
      derived_state_names <- c(derived_state_names, nm)
    }
  }

  ## --- assemble p (ODE order) and X0 (compartment order) ---
  p_terms <- c(prior_spec$order$params_fitted, prior_spec$order$params_fixed)
  p_line  <- sprintf("  array[%d] real p = {%s};", length(p_terms),
                     paste(p_terms, collapse = ", "))
  x_terms <- vapply(seq_len(number_of_comps), function(k) {
    nm <- comp_names[k]
    if (nm %in% prior_spec$order$states_fitted) nm
    else if (nm %in% derived_state_names)       nm
    else if (nm %in% prior_spec$order$states_fixed) prior_spec$fixed_states[[nm]]
    else "0.0"
  }, character(1))
  x_line <- sprintf("  vector[%d] X0 = [%s]';", number_of_comps,
                    paste(x_terms, collapse = ", "))

  ## --- solve + trajectory matrix (col 1 = X0, cols 2.. = ODE solver output) ---
  # ode_<solver>_tol takes the tolerances + step cap as data, so solver_control's
  # abstol/reltol reach Stan; ode_solver is "rk45" (non-stiff) or "bdf" (stiff).
  solve_lines <- c(
    sprintf(paste0("  array[n_grid - 1] vector[%d] sol_raw = ",
                   "ode_%s_tol(cf_ode, X0, t0, ts, rel_tol, abs_tol, max_num_steps, p);"),
            number_of_comps, ode_solver),
    sprintf("  matrix[%d, n_grid] traj;", number_of_comps),
    "  traj[, 1] = X0;",
    "  for (g in 2:n_grid) traj[, g] = sol_raw[g - 1];")

  ## --- observables ---
  obs_lines <- vapply(seq_along(formulas), function(i)
    sprintf("  vector[n_years] mu%d = %s;", i, .observable_to_stan(formulas[i], comp_names)),
    character(1))

  ## --- likelihood ---
  like_lines <- vapply(seq_along(formulas), function(i)
    .likelihood_to_stan(like_specs[[i]]$family, i), character(1))

  ## --- data block (integer matrices only when a count family is present) ---
  data_lines <- c(
    "  int<lower=1> n_years;",
    "  int<lower=1> partition;",
    "  int<lower=1> n_grid;",
    "  int<lower=1> n_streams;",
    "  real<lower=0> rel_tol;",
    "  real<lower=0> abs_tol;",
    "  int<lower=1> max_num_steps;",
    "  vector[n_grid] tgrid;",
    "  matrix[n_years, n_streams] y;",
    "  matrix[n_years, n_streams] obs_mask;",
    "  matrix[n_years, n_streams] cens_mask;",
    "  matrix[n_years, n_streams] limit_mat;",
    "  matrix[n_years, n_streams] lcens_mask;",
    "  matrix[n_years, n_streams] llimit_mat;",
    "  matrix[n_years, n_streams] scale_mat;",
    "  matrix[n_years, n_streams] interval_mask;",  # 1 = interval [A,B] cell
    "  matrix[n_years, n_streams] ilow_mat;",       # A (real edge; soft + continuous hard)
    "  matrix[n_years, n_streams] iupp_mat;",       # B
    "  matrix[n_years, n_streams] idev_lo_mat;",    # soft lower shoulder (0 = hard)
    "  matrix[n_years, n_streams] idev_hi_mat;",    # soft upper shoulder
    "  matrix[n_years, n_streams] asym_mask;",      # 1 = asymmetric cell
    "  matrix[n_years, n_streams] asym_val_mat;",   # anchor A (real)
    "  matrix[n_years, n_streams] asym_dev_mat;",   # soft deviation scale
    "  matrix[n_years, n_streams] asym_dir_mat;")   # +1 soft-up, -1 soft-down
  if (has_counts)
    data_lines <- c(data_lines,
      "  array[n_years, n_streams] int y_int;",     # count observations
      "  array[n_years, n_streams] int lim_int;",   # left-censor effective limit
      "  array[n_years, n_streams] int llim_int;",  # right-censor effective limit
      "  array[n_years, n_streams] int ilow_int;",  # discrete hard-interval A-1
      "  array[n_years, n_streams] int iupp_int;",  # discrete hard-interval B
      "  array[n_years, n_streams] int asym_val_int;")  # discrete asymmetric anchor A

  ## --- assemble the program ---
  join <- function(x) if (length(x)) paste(x, collapse = "\n") else ""
  paste(
    "functions {",
    ode_function,
    "  vector annual_integral(vector flux, vector tgrid, int partition, int n_years) {",
    "    vector[n_years] out;",
    "    for (yr in 1:n_years) {",
    "      int lo = (yr - 1) * partition + 1;",
    "      int hi = yr * partition + 1;",
    "      real s = 0;",
    "      for (k in lo:(hi - 1))",
    "        s += 0.5 * (flux[k] + flux[k + 1]) * (tgrid[k + 1] - tgrid[k]);",
    "      out[yr] = s;",
    "    }",
    "    return out;",
    "  }",
    "}",
    "data {",
    join(data_lines),
    "}",
    "transformed data {",
    "  real t0 = tgrid[1];",
    "  array[n_grid - 1] real ts;",
    "  for (g in 2:n_grid) ts[g - 1] = tgrid[g];",
    "  array[n_years] int annual_idx;",
    "  for (yr in 1:n_years) annual_idx[yr] = yr * partition + 1;",
    "}",
    "parameters {",
    join(param_decls),
    "}",
    "transformed parameters {",
    join(def_lines),
    p_line,
    x_line,
    join(solve_lines),
    join(obs_lines),
    "}",
    "model {",
    join(model_priors),
    join(like_lines),
    "}",
    sep = "\n")
}


# ------------------------------------------------------------
# .fit_bayes_stan: build the Stan program, compile + sample with rstan, and
# return the standard `samples` contract (same shape as .fit_bayes / .fit_bayes_r
# so posterior_draws / posterior_report / plotting all work unchanged).
# ------------------------------------------------------------
.fit_bayes_stan <- function(model, dat, tg, bounds, bc, modelParams,
                            solver = solver_control(backend = "stan")) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("backend = \"stan\" needs the 'rstan' package (install.packages(\"rstan\")).",
         call. = FALSE)
  if (is.null(model$stan_code))
    stop("No Stan ODE was generated for this model (stan_code is NULL); ",
         "the Stan backend cannot build the program.", call. = FALSE)

  prior_spec <- buildPriorSpec(modelParams)
  formulas   <- dat$names_data_points
  like_specs <- if (is.null(dat$likelihood_raw))
    lapply(seq_along(formulas), function(i) parseLikelihood(NA))
  else lapply(dat$likelihood_raw, parseLikelihood)

  # --- capability guards: report what the Stan backend can't do yet ---
  fams <- vapply(like_specs, function(s) s$family, character(1))
  bad_fam <- setdiff(unique(fams), .STAN_FAMILIES)
  if (length(bad_fam))
    stop(sprintf(paste0("The Stan backend supports the {%s} families; got '%s'. ",
                        "Use solver_control(backend = \"julia\") for that one."),
                 paste(.STAN_FAMILIES, collapse = ", "),
                 paste(bad_fam, collapse = ", ")), call. = FALSE)
  discrete_fam <- fams %in% c("poisson", "negbin")
  has_counts   <- any(discrete_fam)

  n_years   <- tg$endpoint - tg$startpoint + 1
  partition <- tg$partition
  n_comps   <- model$structure$number_of_comps
  comp_names <- model$structure$comp_names

  # Per-stream proportional-error scale = mean observed level (matches .fit_bayes).
  ym   <- dat$matrix_data_points
  obsm <- dat$obs_mask
  scale_cols <- vapply(seq_len(ncol(ym)), function(j) {
    vals <- ym[obsm[, j] == 1, j]
    m <- mean(abs(vals), na.rm = TRUE)
    if (!is.finite(m) || m == 0) 1 else m
  }, numeric(1))
  scale_mat <- matrix(rep(scale_cols, each = n_years), nrow = n_years)

  z  <- function(M) { M <- as.matrix(M); M[is.na(M)] <- 0; M }            # real, NA -> 0
  zi <- function(M) { M <- round(as.matrix(M)); M[is.na(M)] <- 0L         # integer, NA -> 0
                      storage.mode(M) <- "integer"; M }
  zdev <- function(M) { M <- as.matrix(M); M[is.na(M) | M <= 0] <- 1; M } # deviation, no /0
  tspan  <- c(min(tg$time), max(tg$time))
  n_grid <- partition * n_years + 1
  tgrid  <- seq(tspan[1], tspan[2], length.out = n_grid)
  stan_data <- list(
    n_years = n_years, partition = partition, n_grid = n_grid,
    n_streams = ncol(ym),
    rel_tol = solver$reltol, abs_tol = solver$abstol, max_num_steps = 100000L,
    tgrid = tgrid,
    y = z(ym), obs_mask = z(obsm), cens_mask = z(dat$cens_mask),
    limit_mat = z(dat$limit_mat), lcens_mask = z(dat$lcens_mask),
    llimit_mat = z(dat$llimit_mat), scale_mat = scale_mat,
    interval_mask = z(dat$interval_mask),
    ilow_mat = z(dat$ilow_mat), iupp_mat = z(dat$iupp_mat),
    idev_lo_mat = z(dat$idev_lo_mat), idev_hi_mat = z(dat$idev_hi_mat),
    asym_mask = z(dat$asym_mask), asym_val_mat = z(dat$asym_val_mat),
    asym_dev_mat = zdev(dat$asym_dev_mat), asym_dir_mat = z(dat$asym_dir_mat))

  # Count families: integer observations + integer censor limits. The effective
  # limits carry the family-aware L-1 shift (left-strict "<L" -> L-1; right-
  # inclusive ">=L" -> L-1), matching the discrete Bayesian CDF in the Julia
  # backend. Non-integer count data is rounded, with a warning.
  if (has_counts) {
    ym_int <- ym
    for (j in which(discrete_fam)) {
      col <- ym[, j]; r <- round(col)
      if (any(abs(col - r) > 1e-8, na.rm = TRUE))
        warning(sprintf(paste0("Stream %d uses a count family but has non-integer ",
                               "data; rounding to nearest integer."), j), call. = FALSE)
      ym_int[, j] <- r
    }
    lim_eff  <- dat$limit_mat
    llim_eff <- dat$llimit_mat
    cm  <- dat$cens_mask;  inc  <- dat$inc_mask
    lcm <- dat$lcens_mask; linc <- dat$linc_mask
    for (j in which(discrete_fam)) {
      strict_j <- which(cm[, j]  == 1 & (is.na(inc[, j])  | inc[, j]  == 0))   # "<L"  -> L-1
      lim_eff[strict_j, j]  <- lim_eff[strict_j, j]  - 1
      incl_j   <- which(lcm[, j] == 1 & !is.na(linc[, j]) & linc[, j] == 1)    # ">=L" -> L-1
      llim_eff[incl_j, j]   <- llim_eff[incl_j, j]   - 1
    }
    stan_data$y_int    <- zi(ym_int)
    stan_data$lim_int  <- zi(lim_eff)
    stan_data$llim_int <- zi(llim_eff)
    # Discrete hard-interval edges: lower shifted A-1 so lcdf(B)-lcdf(A-1) = P(A<=Y<=B).
    stan_data$ilow_int <- zi(dat$ilow_mat - 1)
    stan_data$iupp_int <- zi(dat$iupp_mat)
    stan_data$asym_val_int <- zi(dat$asym_val_mat)   # integer asymmetric anchor
  }

  # Stan's three ODE solvers: rk45 (non-stiff, default), bdf (stiff), adams
  # (non-stiff multistep). Selected by solver_control(solver = "bdf"/"adams").
  ss <- tolower(solver$solver %||% "")
  ode_solver <- if (grepl("bdf", ss)) "bdf" else if (grepl("adams", ss)) "adams" else "rk45"
  code <- buildStanModel(prior_spec, like_specs, formulas, n_years, partition,
                         n_comps, comp_names, model$stan_code,
                         sigma_prior_stan = bc$sigma_prior_stan %||% "normal(0, 1)",
                         phi_prior_stan   = bc$phi_prior_stan   %||% "gamma(2, 0.2)",
                         ode_solver       = ode_solver)

  # Compile (cached by program text; a ~1-2 min C++ build the first time) + sample.
  # progress = FALSE silences the sampler's per-iteration output, the chain
  # start/finish messages, and the auto-opened progress window (the C++
  # compiler's own stderr on first build is outside R's control).
  quiet <- isFALSE(bc$progress)
  sm   <- .stan_compile(code, quiet)
  seed <- if (is.null(bc$seed)) sample.int(.Machine$integer.max, 1L) else bc$seed
  refresh <- if (quiet) 0L else max(1L, floor(bc$iter / 10))
  adapt_delta <- bc$adapt_delta %||% 0.8
  do_sample <- function()
    rstan::sampling(sm, data = stan_data, chains = bc$chains,
                    iter = bc$iter, warmup = bc$warmup, seed = seed,
                    refresh = refresh, open_progress = !quiet,
                    show_messages = !quiet, verbose = FALSE,
                    control = list(adapt_delta = adapt_delta))
  # progress = FALSE also swallows rstan's post-sampling console WARNINGS (R-hat
  # NA, low ESS): those are warning()s, not messages, so suppressMessages alone
  # misses them. The diagnostics are still available via posterior_report() and
  # fit$samples$n_divergent, so nothing is lost -- the fit is just silent.
  sfit <- if (quiet) suppressWarnings(suppressMessages(do_sample())) else do_sample()

  # Record sampler pathologies. Warn about divergences only when not quiet (they
  # are always kept on the fit and shown by posterior_report()).
  n_divergent <- tryCatch(rstan::get_num_divergent(sfit), error = function(e) NA_integer_)
  if (!quiet && !is.na(n_divergent) && n_divergent > 0)
    warning(sprintf(paste0("Stan sampling had %d divergent transition(s): treat the ",
                           "posterior with caution. Try bayes_control(adapt_delta = 0.95) ",
                           "or reparameterise."), n_divergent), call. = FALSE)

  # Sampled quantities, in the estimation order, using the SAME reparam
  # convention as the Julia backend (uniform priors sampled as <name>_n on
  # [0,1], denormalised by posterior_draws()); plus the noise/dispersion
  # hyperparameters actually present (sigma for continuous families, phi_<i>
  # per negbin stream).
  samp_name <- function(nm, spec) if (spec$dist == "Uniform") paste0(nm, "_n") else nm
  sampled <- c(
    vapply(prior_spec$order$params_fitted,
           function(nm) samp_name(nm, prior_spec$params[[nm]]), character(1)),
    vapply(prior_spec$order$states_fitted,
           function(nm) samp_name(nm, prior_spec$states[[nm]]), character(1)))
  if (any(fams %in% c("gaussian", "lognormal"))) sampled <- c(sampled, "sigma")
  for (i in which(fams == "negbin")) sampled <- c(sampled, sprintf("phi_%d", i))

  draws_df <- as.data.frame(sfit, pars = sampled)
  smry <- rstan::summary(sfit, pars = sampled)$summary
  summary_df <- data.frame(
    parameter = rownames(smry), mean = smry[, "mean"], sd = smry[, "sd"],
    rhat = smry[, "Rhat"], ess = smry[, "n_eff"],
    row.names = NULL, stringsAsFactors = FALSE)

  list(draws = draws_df, summary = summary_df,
       n_chains = bc$chains, iter = bc$iter, chains = "rstan (NUTS)",
       n_divergent = n_divergent,
       model_code = code, prior_spec = prior_spec, like_specs = like_specs,
       scale_cols = scale_cols, stream_names = formulas)
}
