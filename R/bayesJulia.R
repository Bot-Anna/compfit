# ============================================================
# bayesJulia.R
# The Bayesian inference layer. Mirrors the Julia ODE bridge in compartmentalFunction.R:
#   buildJuliaBayesModel(...)    -> emits a Turing @model as a Julia string
#   registerJuliaBayesModel(...) -> defines that model in the Julia session
#   sampleWithJulia(...)         -> pushes data in once, runs NUTS, pulls draws
#
# DESIGN
# ------
# The @model reuses the already-registered pure-Julia ODE
# `compartmental_function_jl` (from registerJuliaODEFunction). Per iteration it:
#   1. draws each estimated quantity from its prior,
#   2. reassembles p = [fitted params..., fixed params...] and X = [X1_0, ...]
#      in the EXACT order the ODE expects (matching generateExpressions.R),
#   3. solves the ODE in-Julia (gradients flow through via AD),
#   4. forms each stream's observable (stock / annual / cumulative->annual),
#   5. adds one likelihood term per stream from its declared family.
#
# Everything per-iteration stays in Julia: data + save times cross the R<->Julia
# boundary ONCE (in sampleWithJulia), and only the draws come back.
#
# This file is sourced always, but its functions run ONLY on method="bayes".
# Nothing here is touched by the MLE path.
# ============================================================


# ------------------------------------------------------------
# Translate one stream's R observable formula into a Julia expression over the
# solved trajectory. `sol` in Julia is the ODE solution; we evaluate observables
# on a dense grid (sol(tt)) and reduce to annual values.
#
# Mirrors lossFunction.R's three cases, with cumulative folded to annual:
#   stock  "X3"               -> value at annual snapshot indices
#   annual "annual(flux)"     -> per-year integral of flux
#   cumul. "cumulative(flux)" -> per-year integral (DIFFERENCED), i.e. annual
#
# Returns a Julia expression (string) producing a length-(endpoint-startpoint+1)
# vector mu for that stream.
# ------------------------------------------------------------
.observable_to_julia <- function(formula_str, n_years, partition, comp_names = NULL) {
  .assert_translatable(formula_str, "data-stream Formula")
  f <- gsub("`", "", formula_str)

  # Convert a compartment reference -> sol_grid[k, :]  (state k across the dense
  # grid). A compartment may be referenced by X<k> position (identity case) or by
  # its declared name (S, I, R, ...); either resolves to its canonical row k.
  # NOTE: parameter names need NO translation -- every parameter (fitted,
  # fixed, derived) is bound as a Julia variable in the definitions block, so a
  # bare `g_HA` or `epsilonA` in the formula resolves directly. We only map the
  # state symbols to the solved trajectory rows.
  to_state <- function(s) {
    if (is.null(comp_names) || !length(comp_names))
      return(gsub("\\bX(\\d+)\\b", "sol_grid[\\1, :]", s))
    for (k in order(-nchar(comp_names)))   # longest names first
      s <- gsub(paste0("\\b\\Q", comp_names[k], "\\E\\b"),
                sprintf("sol_grid[%d, :]", k), s, perl = TRUE)
    s
  }
  
  if (grepl("^annual\\(", f) || is_cumulative_stream(f)) {
    inner <- sub("^(annual|cumulative)\\((.*)\\)$", "\\2", f)
    inner_jl <- to_state(inner)
    # Per-year trapezoidal integral of the flux over each one-year window.
    sprintf("annual_integral(%s, grid_t, %d, %d)", inner_jl, partition, n_years)
  } else {
    # Stock: value at the annual snapshot indices.
    inner_jl <- to_state(f)
    sprintf("(%s)[annual_idx]", inner_jl)
  }
}


# ------------------------------------------------------------
# Emit one Julia likelihood statement for a stream, given its family and the
# Julia variable `mu_i` holding its predicted observable, `y_i` the data column,
# and (optionally) a dispersion variable name.
#
# Families:
#   gaussian   : proportional-error Normal. sigma_i = sigma * scale_i (scale =
#                stream level), matching the average-matrix logic. Uses a shared
#                global noise scale `sigma`.
#   lognormal  : log y ~ Normal(log mu, sigma)
#   poisson    : y ~ Poisson(mu)              (mu must be a nonneg count)
#   negbin     : y ~ NegativeBinomial2(mu, phi)  (mean-dispersion param.)
#   (binomial/betabinom deferred until a Denominator column exists)
# ------------------------------------------------------------
.likelihood_to_julia <- function(family, i, disp_var = NULL) {
  # Three-way masked loop per observation:
  #   obs_mask[k,i]==1  -> observed:  y[k,i] ~ dist        (adds logpdf)
  #   cens_mask[k,i]==1 -> censored:  @addlogprob! logcdf(dist, limit_mat[k,i])
  #   else              -> missing:   skip entirely
  # The distribution `dist_i` is built once per k so logpdf and logcdf share it.
  # limit_mat already holds the family-aware effective limit (L or L-1), computed
  # in .fit_bayes, so logcdf(dist_i, limit_mat[k,i]) is correct as-is.
  yk  <- sprintf("y[k, %d]", i)
  ni  <- "size(y, 1)"
  mui <- sprintf("mu%d", i)
  sci <- sprintf("scale%d", i)
  
  dist_expr <- switch(family,
                      gaussian  = sprintf("Normal(%s[k], sigma * %s[k])", mui, sci),
                      lognormal = sprintf("LogNormal(log(max(%s[k], 1e-9)), sigma)", mui),
                      poisson   = sprintf("Poisson(max(%s[k], 1e-9))", mui),
                      negbin    = sprintf("NegativeBinomial2(max(%s[k], 1e-9), %s)", mui, disp_var),
                      stop(sprintf("Likelihood family '%s' not yet supported in Julia builder.", family))
  )
  
  paste0(
    sprintf("for k in 1:%s\n", ni),
    sprintf("        dist_i = %s\n", dist_expr),
    sprintf("        if obs_mask[k, %d] == 1\n", i),
    sprintf("            %s ~ dist_i\n", yk),
    sprintf("        elseif cens_mask[k, %d] == 1\n", i),
    sprintf("            Turing.@addlogprob! logcdf(dist_i, limit_mat[k, %d])\n", i),
    sprintf("        elseif lcens_mask[k, %d] == 1\n", i),
    sprintf("            Turing.@addlogprob! logccdf(dist_i, llimit_mat[k, %d])\n", i),
    # INTERVAL [A,B]: hard edges give P(A <= Y <= B) via a stable log-difference
    # of the CDF (ilow_mat carries the family-aware A-1 shift for discrete
    # families). SOFT shoulders [A,B]~s (idev_lo_mat[k,i] > 0) instead add a
    # family-general plateau penalty log g(mu) = -(under/sl + over/su) with NO
    # CDF, so a soft interval also works for the discrete families under NUTS.
    sprintf("        elseif interval_mask[k, %d] == 1\n", i),
    sprintf("            if idev_lo_mat[k, %d] > 0\n", i),
    sprintf("                Turing.@addlogprob! -(max(ilow_mat[k, %d] - %s[k], 0.0) / idev_lo_mat[k, %d] + max(%s[k] - iupp_mat[k, %d], 0.0) / idev_hi_mat[k, %d])\n", i, mui, i, mui, i, i),
    "            else\n",
    sprintf("                Turing.@addlogprob! cf_logsubexp(logcdf(dist_i, iupp_mat[k, %d]), logcdf(dist_i, ilow_mat[k, %d]))\n", i, i),
    "            end\n",
    # ASYMMETRIC: hard one-sided anchor at A (logccdf/logcdf) + soft linear damping
    # of the model mean's excess in the soft direction (scale = dev).
    sprintf("        elseif asym_mask[k, %d] == 1\n", i),
    sprintf("            if asym_dir_mat[k, %d] > 0\n", i),
    sprintf("                Turing.@addlogprob! logccdf(dist_i, asym_val_mat[k, %d]) - (max(%s[k] - asym_val_mat[k, %d], 0.0) / asym_dev_mat[k, %d])\n", i, mui, i, i),
    "            else\n",
    sprintf("                Turing.@addlogprob! logcdf(dist_i, asym_val_mat[k, %d]) - (max(asym_val_mat[k, %d] - %s[k], 0.0) / asym_dev_mat[k, %d])\n", i, i, mui, i),
    "            end\n",
    "        end\n",
    "    end"
  )
}


# ------------------------------------------------------------
# buildJuliaBayesModel
# Assemble the full Turing @model string.
#
# Args:
#   prior_spec     : buildPriorSpec(modelParams) output
#   like_specs     : list of parseLikelihood() results, one per stream, in the
#                    same column order as the data matrix y
#   formulas       : character vector of stream formulas (same order as y cols)
#   n_years        : endpoint - startpoint + 1
#   partition      : time partition (steps per year)
#   sir_expression : the parameter-unpacking expression (from generateExpressions),
#                    so the ODE sees p in the right order  [reserved for future use]
#   number_of_comps: total compartments (length of X)
# ------------------------------------------------------------
buildJuliaBayesModel <- function(prior_spec,
                                 like_specs,
                                 formulas,
                                 n_years,
                                 partition,
                                 number_of_comps,
                                 comp_names  = paste0("X", seq_len(number_of_comps)),
                                 sigma_prior = "truncated(Normal(0, 1), 0, Inf)",
                                 phi_prior   = "Gamma(2, 5)") {
  
  ## --- Prior block: one ~ per estimated quantity, in ODE order ---
  # Map an R/parsePrior dist spec to a Julia Distributions.jl constructor.
  dist_to_julia <- function(spec) {
    d <- spec$dist
    a <- spec$args
    core <- switch(d,
                   Uniform   = sprintf("Uniform(%s, %s)", a[1], a[2]),
                   Normal    = sprintf("Normal(%s, %s)", a[1], a[2]),
                   LogNormal = sprintf("LogNormal(%s, %s)", a[1], a[2]),
                   Beta      = sprintf("Beta(%s, %s)", a[1], a[2]),
                   Gamma     = sprintf("Gamma(%s, %s)", a[1], a[2]),
                   # Location-scale Student-t: mu + sigma * TDist(nu). Distributions
                   # has no 3-arg TDist, so build it by the affine transform.
                   StudentT  = sprintf("(%s + %s * TDist(%s))", a[2], a[3], a[1]),
                   stop(sprintf("Prior distribution '%s' not supported.", d)))
    # Apply truncation if finite bounds were given (and not already Uniform).
    if (d != "Uniform" && (is.finite(spec$lower) || is.finite(spec$upper))) {
      lo <- if (is.finite(spec$lower)) spec$lower else "-Inf"
      hi <- if (is.finite(spec$upper)) spec$upper else "Inf"
      core <- sprintf("truncated(%s, %s, %s)", core, lo, hi)
    }
    core
  }
  
  prior_lines  <- character(0)
  denorm_lines <- character(0)   # affine maps back to natural scale (definitions block)
  
  # For UNIFORM priors we sample on [0,1] (`name_n ~ Uniform(0,1)`) and map back
  # to natural scale (`name = lo + name_n*(hi-lo)`). This conditions the NUTS
  # geometry (all sampled dims on a common scale) and is STATISTICALLY EXACT for
  # uniform priors -- a pure reparameterization, the implied prior on `name` is
  # still Uniform(lo,hi). Mirrors the MLE normalise/denormalise.
  # For NON-uniform (informative) priors we sample on the natural scale directly
  # (a [0,1] stretch would distort the intended distribution), so `name` is the
  # sampled variable itself and no denorm line is needed.
  emit_prior <- function(nm, spec) {
    if (spec$dist == "Uniform") {
      lo <- spec$args[1]; hi <- spec$args[2]
      prior_lines  <<- c(prior_lines,  sprintf("    %s_n ~ Uniform(0, 1)", nm))
      denorm_lines <<- c(denorm_lines,
                         sprintf("    %s = %s + %s_n * (%s - %s)", nm, lo, nm, hi, lo))
    } else {
      prior_lines <<- c(prior_lines, sprintf("    %s ~ %s", nm, dist_to_julia(spec)))
    }
  }
  
  # Fitted params first, then fitted states -- the estimation order.
  for (nm in prior_spec$order$params_fitted) emit_prior(nm, prior_spec$params[[nm]])
  for (nm in prior_spec$order$states_fitted) emit_prior(nm, prior_spec$states[[nm]])
  
  ## --- Dispersion / noise hyperparameters with default priors ---
  hyper_lines <- character(0)
  needs_sigma <- any(vapply(like_specs, function(s) s$family %in% c("gaussian", "lognormal"),
                            logical(1)))
  if (needs_sigma) {
    hyper_lines <- c(hyper_lines, sprintf("    sigma ~ %s", sigma_prior))
  }
  disp_vars <- rep(list(NULL), length(like_specs))
  for (i in seq_along(like_specs)) {
    if (isTRUE(like_specs[[i]]$dispersion)) {
      dv <- sprintf("phi%d", i)
      disp_vars[[i]] <- dv
      hyper_lines <- c(hyper_lines, sprintf("    %s ~ %s", dv, phi_prior))
    }
  }
  
  ## --- Named-definitions block --------------------------------------------
  # Every parameter is bound as a Julia variable so formulas (observables AND
  # derived states) can reference bare names with no bucket logic.
  #   - fitted params: already bound via the prior `~` lines above.
  #   - fixed params : bound here as `name = value`.
  # The MLE loss uses an `name_0` form for parameters inside state expressions
  # (e.g. prop_lowriskCH_0). We mirror that: bind BOTH `name` and `name_0` for
  # every parameter, so derived-state expressions written with `_0` resolve.
  def_lines <- character(0)
  
  # Denormalize sampled-on-[0,1] quantities back to natural scale FIRST, so
  # every later reference (fixed-param _0 aliases, derived states, p, X0) sees
  # the natural-scale `name`. Empty when there are no uniform priors.
  def_lines <- c(def_lines, denorm_lines)
  
  # Fixed params -> `name = value`  and  `name_0 = value`.
  for (nm in prior_spec$order$params_fixed) {
    val <- prior_spec$fixed_params[[nm]]
    num <- suppressWarnings(as.numeric(val))
    rhs <- if (!is.na(num)) as.character(num) else val   # literal or expression
    def_lines <- c(def_lines,
                   sprintf("    %s = %s", nm, rhs),
                   sprintf("    %s_0 = %s", nm, nm))
  }
  # Fitted params also get a `_0` alias (the bare name is the sampled variable).
  for (nm in prior_spec$order$params_fitted) {
    def_lines <- c(def_lines, sprintf("    %s_0 = %s", nm, nm))
  }
  
  # Derived states: a fixed state whose RHS is an expression (not a plain
  # number) becomes a definition `Xk = <expr>`. Plain-number fixed states and
  # sampled states are NOT defined here (they go straight into X0).
  # Translate Xk references inside the expression to the just-defined state vars
  # (they are scalars here, not trajectories, so no sol_grid mapping).
  is_numeric_str <- function(s) !is.na(suppressWarnings(as.numeric(s)))
  derived_state_names <- character(0)
  for (nm in prior_spec$order$states_fixed) {
    rhs <- prior_spec$fixed_states[[nm]]
    if (!is_numeric_str(rhs)) {
      # expression in params (using `_0` names) -> emit as scalar definition
      def_lines <- c(def_lines, sprintf("    %s = %s", nm, rhs))
      derived_state_names <- c(derived_state_names, nm)
    }
  }
  
  ## --- Assemble p (ODE order) ---------------------------------------------
  # p = [fitted params..., fixed params...]; every term is a bound variable now
  # (fitted = sampled var, fixed = defined var), so we emit names, not literals.
  p_terms <- c(prior_spec$order$params_fitted, prior_spec$order$params_fixed)
  p_line  <- sprintf("    p = [%s]", paste(p_terms, collapse = ", "))
  
  ## --- Assemble X0 (compartment order) ------------------------------------
  # For each compartment k: sampled state -> its variable; derived state -> its
  # defined variable; plain-number fixed state -> the literal.
  x_terms <- character(number_of_comps)
  for (k in seq_len(number_of_comps)) {
    nm <- comp_names[k]          # canonical name of compartment k (X<k> or S/I/R/...)
    if (nm %in% prior_spec$order$states_fitted) {
      x_terms[k] <- nm                                   # sampled variable
    } else if (nm %in% derived_state_names) {
      x_terms[k] <- nm                                   # defined variable
    } else if (nm %in% prior_spec$order$states_fixed) {
      x_terms[k] <- prior_spec$fixed_states[[nm]]        # plain-number literal
    } else {
      x_terms[k] <- "0.0"
    }
  }
  x_line <- sprintf("    X0 = [%s]", paste(x_terms, collapse = ", "))
  
  ## --- Observable + likelihood blocks ---
  obs_lines  <- character(0)
  like_lines <- character(0)
  for (i in seq_along(formulas)) {
    mu_expr <- .observable_to_julia(formulas[i], n_years, partition, comp_names)
    obs_lines  <- c(obs_lines, sprintf("    mu%d = %s", i, mu_expr))
    like_lines <- c(like_lines, paste0("    ",
                                       .likelihood_to_julia(like_specs[[i]]$family, i, disp_vars[[i]])))
  }
  
  ## --- Model preamble: solve + grid ---
  # annual_idx mirrors the R loss EXACTLY: seq(partition+1, partition*n_years+1,
  # by=partition) in R is 1-based; Julia is also 1-based, so the same formula.
  # This gives n_years snapshot indices starting one full year in (t=0 initial
  # condition is skipped), matching the differenced/annual data matrix.
  preamble <- paste(
    "    # dense evaluation grid over the integration span",
    "    grid_t = range(tspan[1], tspan[2]; length = n_grid)",
    "    prob = ODEProblem(compartmental_function_jl, X0, (tspan[1], tspan[2]), p)",
    "    sol  = solve(prob, AutoTsit5(Rosenbrock23()); saveat = grid_t,",
    "                 abstol = 1e-6, reltol = 1e-6)",
    "    # If the solve fails (a region NUTS wandered into where the ODE is",
    "    # pathological), reject this proposal -- contribute -Inf to the",
    "    # log-density and return -- rather than letting the observables error.",
    "    if !SciMLBase.successful_retcode(sol)",
    "        Turing.@addlogprob! -Inf",
    "        return nothing",
    "    end",
    "    sol_grid = Array(sol)            # (n_comps x n_grid)",
    "    annual_idx = collect((partition + 1):partition:(partition * n_years + 1))",
    sep = "\n")
  
  scale_block <- paste(
    sprintf("    scale%d = scale_mat[:, %d]", seq_along(formulas), seq_along(formulas)),
    collapse = "\n")
  
  # Join the optional blocks with explicit "\n", never paste(x, '\n') which
  # would inject a stray space and break Julia's parser.
  join_block <- function(lines) if (length(lines)) paste(lines, collapse = "\n") else ""
  
  ## --- Full @model ---
  parts <- c(
    "@model function bayes_fit_model(y, tspan, n_years, partition, n_grid, scale_mat, obs_mask, cens_mask, limit_mat, lcens_mask, llimit_mat, interval_mask, ilow_mat, iupp_mat, idev_lo_mat, idev_hi_mat, asym_mask, asym_val_mat, asym_dev_mat, asym_dir_mat)",
    join_block(prior_lines),
    join_block(hyper_lines),
    "    # ---- named definitions: fixed params, _0 aliases, derived states ----",
    join_block(def_lines),
    p_line,
    x_line,
    preamble,
    "    # ---- per-stream scale columns ----",
    scale_block,
    "    # ---- observables ----",
    join_block(obs_lines),
    "    # ---- likelihood ----",
    join_block(like_lines),
    "    return nothing",
    "end"
  )
  parts <- parts[parts != ""]          # drop any empty optional blocks
  model_code <- paste(parts, collapse = "\n")
  
  model_code
}


# ------------------------------------------------------------
# Helper Julia code defining annual_integral and NegativeBinomial2, registered
# once alongside the model. annual_integral mirrors integrate_for_loss in R.
# ------------------------------------------------------------
.bayes_julia_helpers <- function() {
  paste(
    "function annual_integral(flux, tgrid, partition, n_years)",
    "    out = Vector{eltype(flux)}(undef, n_years)",
    "    for yr in 1:n_years",
    "        lo = (yr - 1) * partition + 1",
    "        hi = yr * partition + 1",
    "        s = 0.0",
    "        for k in lo:(hi - 1)",
    "            s += 0.5 * (flux[k] + flux[k + 1]) * (tgrid[k + 1] - tgrid[k])",
    "        end",
    "        out[yr] = s",
    "    end",
    "    return out",
    "end",
    "",
    "# Mean-dispersion Negative Binomial: mean mu, dispersion phi (Var = mu + mu^2/phi)",
    "function NegativeBinomial2(mu, phi)",
    "    p = phi / (phi + mu)",
    "    r = phi",
    "    return NegativeBinomial(r, p)",
    "end",
    "",
    "# Numerically stable log(1 - exp(x)) for x <= 0, and log(exp(a) - exp(b)) for",
    "# a >= b, used by interval-censored [A,B] likelihood contributions.",
    "cf_log1mexp(x) = x < -log(2.0) ? log1p(-exp(x)) : log(-expm1(x))",
    "cf_logsubexp(a, b) = a + cf_log1mexp(b - a)",
    sep = "\n")
}


# ------------------------------------------------------------
# registerJuliaBayesModel: ensure Turing + helpers + the model are defined.
# ------------------------------------------------------------
registerJuliaBayesModel <- function(model_code) {
  .compfit_ensure_julia()
  JuliaCall::julia_command("using Turing")
  JuliaCall::julia_command("using Distributions")
  
  # IMPORTANT: JuliaCall's string-eval path (julia_eval / julia_command) parses
  # only a SINGLE top-level expression and throws "extra token after end of
  # expression" on multi-definition code (e.g. two `function ... end` blocks).
  # include() parses a file statement-by-statement, so we write to temp files
  # and include them -- this mirrors the only method that reliably works.
  helpers_file <- tempfile(fileext = ".jl")
  model_file   <- tempfile(fileext = ".jl")
  writeLines(.bayes_julia_helpers(), helpers_file)
  writeLines(model_code, model_file)
  
  # Forward-slash paths for Julia even on Windows.
  hf <- gsub("\\\\", "/", helpers_file)
  mf <- gsub("\\\\", "/", model_file)
  JuliaCall::julia_command(sprintf('include("%s")', hf))
  JuliaCall::julia_command(sprintf('include("%s")', mf))
  
  message("Julia Bayesian model 'bayes_fit_model' registered.")
}


# ------------------------------------------------------------
# sampleWithJulia: push data once, run NUTS, return draws + diagnostics.
#
# Args:
#   y_matrix   : data matrix (rows = years, cols = streams), already differenced
#                for cumulative streams (matches the loss / .prepare_data).
#   tspan      : c(t0, t1) integration span.
#   n_years, partition : as in the model.
#   scale_mat  : per-stream scale columns (e.g. stream means) for proportional
#                gaussian; pass a matrix of ones if unused.
#   bc         : bayes_control() list (chains, iter, warmup, sampler, seed).
#   n_grid     : dense grid length for the solve (default partition*n_years+1).
#
# Returns: list(draws = data.frame, summary = data.frame, chains = <julia obj name>)
# ------------------------------------------------------------
sampleWithJulia <- function(y_matrix, tspan, n_years, partition, scale_mat, bc,
                            obs_mask = NULL, cens_mask = NULL, limit_mat = NULL,
                            lcens_mask = NULL, llimit_mat = NULL,
                            interval_mask = NULL, ilow_mat = NULL, iupp_mat = NULL,
                            idev_lo_mat = NULL, idev_hi_mat = NULL,
                            asym_mask = NULL, asym_val_mat = NULL,
                            asym_dev_mat = NULL, asym_dir_mat = NULL,
                            n_grid = NULL) {
  .compfit_ensure_julia()
  if (is.null(n_grid)) n_grid <- partition * n_years + 1

  # Default masks: everything observed, nothing censored/interval/asym (back-compat).
  z0  <- function() matrix(0, nrow(y_matrix), ncol(y_matrix))
  if (is.null(obs_mask))   obs_mask   <- matrix(1L, nrow(y_matrix), ncol(y_matrix))
  if (is.null(cens_mask))  cens_mask  <- matrix(0L, nrow(y_matrix), ncol(y_matrix))
  if (is.null(limit_mat))  limit_mat  <- z0()
  if (is.null(lcens_mask)) lcens_mask <- matrix(0L, nrow(y_matrix), ncol(y_matrix))
  if (is.null(llimit_mat)) llimit_mat <- z0()
  if (is.null(interval_mask)) interval_mask <- matrix(0L, nrow(y_matrix), ncol(y_matrix))
  if (is.null(ilow_mat))      ilow_mat      <- z0()
  if (is.null(iupp_mat))      iupp_mat      <- z0()
  if (is.null(idev_lo_mat))   idev_lo_mat   <- z0()      # 0 = hard interval edges
  if (is.null(idev_hi_mat))   idev_hi_mat   <- z0() + 1  # avoid /0 in the soft branch
  if (is.null(asym_mask))     asym_mask     <- matrix(0L, nrow(y_matrix), ncol(y_matrix))
  if (is.null(asym_val_mat))  asym_val_mat  <- z0()
  if (is.null(asym_dev_mat))  asym_dev_mat  <- z0() + 1   # avoid /0 when unused
  if (is.null(asym_dir_mat))  asym_dir_mat  <- z0()

  JuliaCall::julia_assign("_y_jl",       as.matrix(y_matrix))
  JuliaCall::julia_assign("_tspan_jl",   as.numeric(tspan))
  JuliaCall::julia_assign("_scale_jl",   as.matrix(scale_mat))
  JuliaCall::julia_assign("_n_years_jl", as.integer(n_years))
  JuliaCall::julia_assign("_part_jl",    as.integer(partition))
  JuliaCall::julia_assign("_ngrid_jl",   as.integer(n_grid))
  JuliaCall::julia_assign("_obs_jl",     matrix(as.integer(obs_mask),  nrow(obs_mask)))
  JuliaCall::julia_assign("_cens_jl",    matrix(as.integer(cens_mask), nrow(cens_mask)))
  JuliaCall::julia_assign("_limit_jl",   as.matrix(limit_mat))
  JuliaCall::julia_assign("_lcens_jl",   matrix(as.integer(lcens_mask), nrow(lcens_mask)))
  JuliaCall::julia_assign("_llimit_jl",  as.matrix(llimit_mat))
  JuliaCall::julia_assign("_int_jl",     matrix(as.integer(interval_mask), nrow(interval_mask)))
  JuliaCall::julia_assign("_ilow_jl",    as.matrix(ilow_mat))
  JuliaCall::julia_assign("_iupp_jl",    as.matrix(iupp_mat))
  JuliaCall::julia_assign("_idlo_jl",    as.matrix(idev_lo_mat))
  JuliaCall::julia_assign("_idhi_jl",    as.matrix(idev_hi_mat))
  JuliaCall::julia_assign("_asym_jl",    matrix(as.integer(asym_mask), nrow(asym_mask)))
  JuliaCall::julia_assign("_aval_jl",    as.matrix(asym_val_mat))
  JuliaCall::julia_assign("_adev_jl",    as.matrix(asym_dev_mat))
  JuliaCall::julia_assign("_adir_jl",    as.matrix(asym_dir_mat))

  seed_cmd <- if (!is.null(bc$seed)) sprintf("using Random; Random.seed!(%d); ", bc$seed) else ""

  build_cmd <- paste0(
    seed_cmd,
    "_model_jl = bayes_fit_model(_y_jl, _tspan_jl, _n_years_jl, _part_jl, _ngrid_jl, _scale_jl, _obs_jl, _cens_jl, _limit_jl, _lcens_jl, _llimit_jl, _int_jl, _ilow_jl, _iupp_jl, _idlo_jl, _idhi_jl, _asym_jl, _aval_jl, _adev_jl, _adir_jl); ",
    "nothing"
  )
  JuliaCall::julia_command(build_cmd)

  prog <- if (isFALSE(bc$progress)) "false" else "true"
  # progress = false on the sample() call can still leak the per-iteration ETA
  # bar under MCMCThreads(); Turing.setprogress!(false) is the global off-switch.
  # Wrap it so its own one-line info notice is silenced too.
  if (isFALSE(bc$progress))
    JuliaCall::julia_command(
      "Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do; Turing.setprogress!(false); end;")
  sample_cmd <- sprintf(
    "_chain_jl = sample(_model_jl, %s, MCMCThreads(), %d, %d; progress = %s); nothing",
    bc$sampler, bc$iter, bc$chains, prog
  )
  JuliaCall::julia_command(sample_cmd)
  
  # Pull draws back to R (iter*chains x params, chains flattened).
  draws   <- JuliaCall::julia_eval("Array(_chain_jl)")
  pnames  <- JuliaCall::julia_eval("string.(names(_chain_jl, :parameters))")
  
  # Summary (mean/sd/rhat/ess). MCMCChains column names vary across versions
  # (:ess vs :ess_bulk), so build the table defensively and surface any error
  # instead of silently returning NULL.
  summ <- tryCatch(
    JuliaCall::julia_eval('
      let
        s  = summarystats(_chain_jl)
        nt = s.nt
        getcol(n) = haskey(nt, n) ? collect(getproperty(nt, n)) : fill(NaN, length(nt.parameters))
        ess = haskey(nt, :ess_bulk) ? collect(nt.ess_bulk) :
              (haskey(nt, :ess) ? collect(nt.ess) : fill(NaN, length(nt.parameters)))
        ( string.(nt.parameters), getcol(:mean), getcol(:std), getcol(:rhat), ess )
      end'),
    error = function(e) {
      message("Summary extraction failed: ", conditionMessage(e),
              "\n(draws are still available; diagnostics can be computed in R).")
      NULL
    }
  )
  
  draws_df <- as.data.frame(draws)
  if (length(pnames) == ncol(draws_df)) names(draws_df) <- pnames
  
  summary_df <- NULL
  if (!is.null(summ)) {
    summary_df <- data.frame(parameter = summ[[1]], mean = summ[[2]],
                             sd = summ[[3]], rhat = summ[[4]], ess = summ[[5]],
                             stringsAsFactors = FALSE)
  }
  
  list(draws = draws_df, summary = summary_df,
       n_chains = bc$chains, iter = bc$iter, chains = "_chain_jl")
}