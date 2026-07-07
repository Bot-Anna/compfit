# ============================================================
# bayesPost.R
# Post-processing helpers for a Bayesian compartmentalFit.
#
# The sampler works on NORMALIZED parameters (uniform priors sampled on [0,1]
# as `name_n`), so raw draws/summary are on [0,1]. These helpers map back to
# natural scale using the bounds in prior_spec, and provide convenience
# extractors so everything is reachable from `fit` directly:
#
#   posterior_summary(fit)          -> natural-scale mean/sd/rhat/ess table
#   posterior_draws(fit)            -> natural-scale draws data frame
#   posterior_means(fit)            -> named natural-scale posterior means
#   posterior_intervals(fit, p)     -> natural-scale credible intervals
#   compare_to_mle(fit, fit_mle)    -> side-by-side Bayes-mean vs MLE (correctness check)
#
# All are no-ops / informative errors for non-Bayesian fits.
# ============================================================

# ---- internal: bounds lookup for uniform-prior (normalized) quantities ----
.bayes_bounds <- function(prior_spec) {
  bounds <- list()
  add <- function(lst) for (nm in names(lst)) {
    sp <- lst[[nm]]
    if (!is.null(sp$dist) && sp$dist == "Uniform")
      bounds[[nm]] <<- as.numeric(sp$args[1:2])
  }
  add(prior_spec$params)
  add(prior_spec$states)
  bounds
}

.check_bayes <- function(fit) {
  if (is.null(fit$samples))
    stop("This fit has no posterior samples (not a method='bayes' fit, or it failed).")
  invisible(TRUE)
}

# ---- Natural-scale draws -------------------------------------------------
# Maps each `name_n` column on [0,1] back to `name` on [lo,hi]. Non-normalized
# columns (e.g. sigma, phi) pass through unchanged.
#' Posterior draws (natural scale)
#'
#' Returns the posterior draws with each normalised `name_n` column mapped back
#' to its natural scale `name` on `[lo, hi]`; non-normalised columns (e.g.
#' `sigma`, `phi`) pass through unchanged.
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @return A data frame of natural-scale draws.
#' @examples
#' \dontrun{
#' d <- posterior_draws(fit_b)   # fit_b a method = "bayes" fit
#' head(d)
#' }
#' @export
posterior_draws <- function(fit) {
  .check_bayes(fit)
  draws  <- fit$samples$draws
  bounds <- .bayes_bounds(fit$samples$prior_spec)
  
  out <- draws
  for (cn in names(out)) {
    if (grepl("_n$", cn)) {
      base <- sub("_n$", "", cn)
      b <- bounds[[base]]
      if (!is.null(b)) {
        out[[cn]] <- b[1] + out[[cn]] * (b[2] - b[1])
        names(out)[names(out) == cn] <- base
      }
    }
  }
  out
}

# ---- Natural-scale summary ----------------------------------------------
# mean maps affinely; sd scales by width; rhat/ess are scale-invariant.
#' Posterior summary (natural scale)
#'
#' Per-parameter posterior summary (mean, sd, rhat, ess) with means/sds mapped
#' to the natural scale.
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @return A data frame, one row per parameter.
#' @examples
#' \dontrun{
#' posterior_summary(fit_b)   # fit_b a method = "bayes" fit
#' }
#' @export
posterior_summary <- function(fit) {
  .check_bayes(fit)
  s <- fit$samples$summary
  if (is.null(s))
    stop("fit$samples$summary is NULL (summary extraction failed during sampling). ",
         "Re-run with the updated sampleWithJulia, or compute from posterior_draws(fit).")
  bounds <- .bayes_bounds(fit$samples$prior_spec)
  
  for (i in seq_len(nrow(s))) {
    pn <- s$parameter[i]
    if (grepl("_n$", pn)) {
      base <- sub("_n$", "", pn)
      b <- bounds[[base]]
      if (!is.null(b)) {
        lo <- b[1]; w <- b[2] - b[1]
        s$parameter[i] <- base
        s$mean[i] <- lo + s$mean[i] * w
        s$sd[i]   <- s$sd[i] * w
        # rhat, ess unchanged (scale-invariant)
      }
    }
  }
  s
}

# ---- Posterior means (named, natural scale) ------------------------------
#' Posterior means (natural scale)
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @return A named numeric vector of posterior means.
#' @examples
#' \dontrun{
#' posterior_means(fit_b)   # fit_b a method = "bayes" fit
#' }
#' @export
posterior_means <- function(fit) {
  d <- posterior_draws(fit)
  colMeans(d, na.rm = TRUE)
}

# ---- Credible intervals (natural scale) ----------------------------------
# p = central mass (default 0.95 -> 2.5% / 97.5% quantiles).
#' Posterior credible intervals (natural scale)
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @param p Central credible mass (default `0.95`).
#' @return A data frame with lower/upper interval bounds per parameter.
#' @examples
#' \dontrun{
#' posterior_intervals(fit_b, p = 0.9)   # fit_b a method = "bayes" fit
#' }
#' @export
posterior_intervals <- function(fit, p = 0.95) {
  d <- posterior_draws(fit)
  a <- (1 - p) / 2
  qs <- t(apply(d, 2, quantile, probs = c(a, 0.5, 1 - a), na.rm = TRUE))
  colnames(qs) <- c("lower", "median", "upper")
  as.data.frame(qs)
}

# ---- Correctness check: Bayes posterior summaries vs MLE point -----------
# If the p-ordering and wiring are correct, the natural-scale posterior centre
# (mean AND median) of the fitted parameters should land near the MLE
# estimates. Large discrepancies flag a wiring/ordering bug rather than a
# modelling difference. Both summaries are shown: the mean is sensitive to a
# skewed posterior, the median more robust -- agreement of both with the MLE is
# the strongest signal.
#' Compare a Bayesian fit to an MLE fit
#'
#' Tabulates posterior mean and median against the MLE point estimate, with
#' relative differences for each.
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @param fit_mle A maximum-likelihood `"compartmentalFit"` object.
#' @return A data frame with `parameter`, `bayes_mean`, `bayes_median`, `mle`,
#'   and relative-difference columns.
#' @examples
#' \dontrun{
#' compare_to_mle(fit_b, fit_mle)   # Bayesian vs maximum-likelihood fit
#' }
#' @export
compare_to_mle <- function(fit, fit_mle) {
  .check_bayes(fit)
  if (is.null(fit_mle$point))
    stop("fit_mle has no point estimate (did the MLE fit succeed?).")

  bmeans <- posterior_means(fit)
  draws  <- posterior_draws(fit)
  bmed   <- apply(draws, 2, median, na.rm = TRUE)[names(bmeans)]   # align to means
  mle    <- c(fit_mle$point$parms,
              setNames(as.numeric(fit_mle$point$initial_state),
                       names(fit_mle$point$initial_state)))

  common <- intersect(names(bmeans), names(mle))
  if (length(common) == 0)
    warning("No overlapping parameter names between Bayes summaries and MLE point; ",
            "check naming (e.g. states X1.. vs params).")

  rel <- function(b) as.numeric((b - mle[common]) / pmax(abs(mle[common]), 1e-8))
  data.frame(
    parameter       = common,
    bayes_mean      = as.numeric(bmeans[common]),
    bayes_median    = as.numeric(bmed[common]),
    mle             = as.numeric(mle[common]),
    rel_diff_mean   = rel(bmeans[common]),
    rel_diff_median = rel(bmed[common]),
    stringsAsFactors = FALSE
  )
}

# ---- One-call report -----------------------------------------------------
# Bundles everything you typically want from a Bayesian fit into one list, and
# prints a short convergence verdict. Everything is on the NATURAL scale.
#
#   rep <- posterior_report(fit)
#   rep$summary     # mean/sd/rhat/ess (natural scale)
#   rep$means       # named posterior means
#   rep$intervals   # 95% credible intervals
#   rep$draws       # full natural-scale draws
#   rep$diagnostics # max_rhat, min_ess, n flagged, etc.
#' Report on a Bayesian posterior
#'
#' Prints a posterior summary and convergence diagnostics (R-hat / ESS) with
#' threshold flags, and returns the summary invisibly.
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @param p Central credible mass (default `0.95`).
#' @param rhat_thresh R-hat threshold above which a parameter is flagged.
#' @param ess_thresh Effective-sample-size threshold below which a parameter is
#'   flagged.
#' @return The posterior summary data frame, invisibly.
#' @examples
#' \dontrun{
#' posterior_report(fit_b)   # fit_b a method = "bayes" fit
#' }
#' @export
posterior_report <- function(fit, p = 0.95, rhat_thresh = 1.01, ess_thresh = 400) {
  .check_bayes(fit)
  
  draws     <- posterior_draws(fit)
  means     <- colMeans(draws, na.rm = TRUE)
  intervals <- posterior_intervals(fit, p = p)
  
  summary <- tryCatch(posterior_summary(fit), error = function(e) NULL)
  
  diagnostics <- NULL
  if (!is.null(summary)) {
    # R-backend (or any) fits may have un-computable diagnostics -> all NA.
    has_rhat <- !all(is.na(summary$rhat))
    has_ess  <- !all(is.na(summary$ess))
    max_rhat <- if (has_rhat) max(summary$rhat, na.rm = TRUE) else NA_real_
    min_ess  <- if (has_ess)  min(summary$ess,  na.rm = TRUE) else NA_real_
    bad_rhat <- if (has_rhat) summary$parameter[which(summary$rhat > rhat_thresh)] else character(0)
    low_ess  <- if (has_ess)  summary$parameter[which(summary$ess  < ess_thresh)]  else character(0)
    diagnostics <- list(max_rhat = max_rhat, min_ess = min_ess,
                        bad_rhat = bad_rhat, low_ess = low_ess,
                        rhat_thresh = rhat_thresh, ess_thresh = ess_thresh)

    cat("Posterior convergence:\n")
    if (has_rhat)
      cat(sprintf("  max R-hat: %.3f  (target <= %.2f)  %s\n", max_rhat, rhat_thresh,
                  if (max_rhat <= rhat_thresh) "OK" else "** above threshold **"))
    else
      cat("  max R-hat: not available (could not be computed for this fit)\n")
    if (has_ess)
      cat(sprintf("  min ESS:   %.0f   (target >= %d)   %s\n", min_ess, ess_thresh,
                  if (min_ess >= ess_thresh) "OK" else "** below threshold **"))
    else
      cat("  min ESS:   not available\n")
    if (length(bad_rhat)) cat("  high R-hat:", paste(bad_rhat, collapse = ", "), "\n")
    if (length(low_ess))  cat("  low ESS:   ", paste(low_ess,  collapse = ", "), "\n")
  } else {
    cat("Posterior summary unavailable (summary extraction returned NULL).\n",
        "Means and intervals are still computed from the draws.\n")
  }
  
  invisible(list(summary = summary, means = means, intervals = intervals,
                 draws = draws, diagnostics = diagnostics))
}

# ---- print method --------------------------------------------------------
#' Print a compartmental fit
#'
#' @param x A `"compartmentalFit"` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' \dontrun{
#' fit   # fit from fitCompartmentalModel(); auto-prints via this method
#' }
#' @exportS3Method print compartmentalFit
print.compartmentalFit <- function(x, ...) {
  cat(sprintf("<compartmentalFit>  method = %s   success = %s\n",
              x$method, x$success))
  if (!isTRUE(x$success) && !is.null(x$error_msg))
    cat("  error:", x$error_msg, "\n")
  
  if (x$method == "bayes" && !is.null(x$samples)) {
    d <- x$samples$draws
    cat(sprintf("  posterior: %d draws x %d sampled quantities (%d chains, %d iter)\n",
                nrow(d), ncol(d),
                x$samples$n_chains %||% NA, x$samples$iter %||% NA))
    s <- tryCatch(posterior_summary(x), error = function(e) NULL)
    if (!is.null(s)) {
      # Guard the all-NA case (a degenerate fit) so max()/min() do not warn.
      rhat <- if (any(is.finite(s$rhat))) sprintf("%.3f", max(s$rhat, na.rm = TRUE)) else "NA"
      ess  <- if (any(is.finite(s$ess)))  sprintf("%.0f", min(s$ess,  na.rm = TRUE)) else "NA"
      cat(sprintf("  max R-hat = %s   min ESS = %s\n", rhat, ess))
    }
    nd <- x$samples$n_divergent
    if (!is.null(nd) && !is.na(nd) && nd > 0)
      cat(sprintf("  WARNING: %d divergent transition(s) -- raise adapt_delta or reparameterise\n", nd))
    cat("  -> posterior_report(fit), posterior_summary(fit), posterior_intervals(fit)\n")
  } else if (!is.null(x$point)) {
    cat("  point estimate available: fit$point$parms, fit$point$initial_state\n")
  }
  invisible(x)
}