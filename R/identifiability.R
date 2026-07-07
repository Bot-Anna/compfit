# ============================================================
# identifiability.R
# Practical identifiability from the POSTERIOR (Bayes only).
#
#   identifiability_report(fit, components = ...)   # BAYES ONLY
#     "narrowing"   (a) prior -> posterior width shrinkage per parameter
#     "coverage"    (b) does the 95% CrI still span most of the prior?
#     "correlation" (c) posterior correlations + flagged trade-offs (+ plots)
#
# (a)/(b)/(c) all read the posterior, so this errors on an MLE fit.
#
# For variance-based (Sobol) sensitivity -- which is a property of the MODEL,
# not the posterior, and works for any fit -- see sensitivity.R
# (sobol_loss(), sobol_streams()). The two are deliberately separate: this file
# answers "given the data, what did we learn?"; Sobol answers "across the
# parameter space, what drives the output?".
# ============================================================

# ---- helpers ---------------------------------------------------------------

.ident_bounds <- function(prior_spec) {
  b <- list()
  add <- function(lst) for (nm in names(lst)) {
    sp <- lst[[nm]]
    if (!is.null(sp$dist) && sp$dist == "Uniform") b[[nm]] <<- as.numeric(sp$args[1:2])
  }
  add(prior_spec$params); add(prior_spec$states)
  b
}

# =====================================================================
#  identifiability_report  (BAYES ONLY)
# =====================================================================
#' Practical identifiability report
#'
#' Diagnoses practical identifiability of a Bayesian fit: prior-to-posterior
#' narrowing, posterior coverage of the prior range, and posterior parameter
#' correlations, flagging weakly-identified or correlated parameters.
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @param components Which diagnostics to run: narrowing, coverage, correlation.
#' @param cor_threshold Absolute correlation above which a pair is flagged.
#' @param coverage_threshold Posterior/prior coverage above which a parameter is
#'   flagged as poorly narrowed.
#' @param cri Central credible mass used for the diagnostics (default `0.95`).
#' @param plots If `TRUE`, include diagnostic plots.
#' @return Invisibly, a list of the diagnostic tables (and plots).
#' @examples
#' \dontrun{
#' identifiability_report(fit_b)   # fit_b a method = "bayes" fit
#' }
#' @export
identifiability_report <- function(fit,
                                   components = c("narrowing", "coverage", "correlation"),
                                   cor_threshold = 0.8,
                                   coverage_threshold = 0.7,
                                   cri = 0.95,
                                   plots = TRUE) {
  if (!inherits(fit, "compartmentalFit"))
    stop("identifiability_report() expects a 'compartmentalFit'.")
  if (is.null(fit$samples))
    stop("identifiability_report() works only for a method='bayes' fit ",
         "(it reads the posterior). For an MLE fit, there is no posterior to ",
         "assess; use sobol_loss()/sobol_streams() for sensitivity instead.")
  
  draws  <- posterior_draws(fit)
  ps     <- fit$samples$prior_spec
  bounds <- .ident_bounds(ps)
  out    <- list()
  alpha  <- (1 - cri) / 2
  
  # Parameters that have a uniform prior (so narrowing/coverage are defined).
  pars <- intersect(names(draws), names(bounds))
  
  # ---- (a) narrowing: posterior width / prior width -----------------------
  if ("narrowing" %in% components) {
    nar <- lapply(pars, function(nm) {
      lo <- bounds[[nm]][1]; hi <- bounds[[nm]][2]
      prior_w <- hi - lo
      q <- quantile(draws[[nm]], c(alpha, 1 - alpha), na.rm = TRUE)
      post_w <- as.numeric(q[2] - q[1])
      data.frame(parameter = nm,
                 prior_width = prior_w,
                 post_width  = post_w,
                 width_ratio = post_w / prior_w,         # ~1 = uninformative
                 informative = (post_w / prior_w) < 0.5) # rough flag
    })
    nar <- do.call(rbind, nar)
    nar <- nar[order(-nar$width_ratio), ]                # least-informed first
    out$narrowing <- nar
  }
  
  # ---- (b) coverage: CrI span as a fraction of the prior range ------------
  if ("coverage" %in% components) {
    cov <- lapply(pars, function(nm) {
      lo <- bounds[[nm]][1]; hi <- bounds[[nm]][2]
      q <- quantile(draws[[nm]], c(alpha, 1 - alpha), na.rm = TRUE)
      frac <- as.numeric((q[2] - q[1]) / (hi - lo))
      data.frame(parameter = nm,
                 cri_lower = as.numeric(q[1]), cri_upper = as.numeric(q[2]),
                 prior_lower = lo, prior_upper = hi,
                 cri_frac_of_prior = frac,
                 poorly_constrained = frac > coverage_threshold)
    })
    cov <- do.call(rbind, cov)
    cov <- cov[order(-cov$cri_frac_of_prior), ]
    out$coverage <- cov
  }
  
  # ---- (c) posterior correlations + trade-off flags -----------------------
  if ("correlation" %in% components) {
    # Use the sampled quantities (drop hyperparameters for clarity of the
    # parameter trade-off picture, but keep them available via $cor_full).
    M <- as.matrix(draws)
    keep <- !grepl("^sigma$|^phi[0-9]*$", colnames(M))
    Cp <- suppressWarnings(cor(M[, keep, drop = FALSE]))
    out$cor_matrix <- Cp
    
    # Strongest off-diagonal pairs above threshold.
    cn <- colnames(Cp)
    pairs <- which(upper.tri(Cp), arr.ind = TRUE)
    pr <- data.frame(
      p1 = cn[pairs[, 1]], p2 = cn[pairs[, 2]],
      cor = Cp[upper.tri(Cp)], stringsAsFactors = FALSE)
    pr <- pr[order(-abs(pr$cor)), ]
    out$top_correlations <- pr[abs(pr$cor) >= cor_threshold, , drop = FALSE]
    
    if (isTRUE(plots) && requireNamespace("ggplot2", quietly = TRUE)) {
      out$cor_plot <- .ident_cor_heatmap(Cp)
    }
  }
  
  # ---- console summary -----------------------------------------------------
  cat("Identifiability report (practical; from the posterior)\n")
  if (!is.null(out$narrowing)) {
    flagged <- out$narrowing$parameter[out$narrowing$width_ratio > 0.8]
    cat(sprintf("  (a) narrowing: %d/%d parameters barely informed (width_ratio > 0.8)\n",
                length(flagged), nrow(out$narrowing)))
    if (length(flagged)) cat("      ", paste(flagged, collapse = ", "), "\n")
  }
  if (!is.null(out$coverage)) {
    pc <- out$coverage$parameter[out$coverage$poorly_constrained]
    cat(sprintf("  (b) coverage: %d parameters with CrI spanning > %.0f%% of prior\n",
                length(pc), 100 * coverage_threshold))
    if (length(pc)) cat("      ", paste(pc, collapse = ", "), "\n")
  }
  if (!is.null(out$top_correlations)) {
    cat(sprintf("  (c) correlation: %d parameter pairs with |r| >= %.2f (trade-offs)\n",
                nrow(out$top_correlations), cor_threshold))
    if (nrow(out$top_correlations))
      for (k in seq_len(min(10, nrow(out$top_correlations))))
        cat(sprintf("      %s <-> %s : %+.2f\n",
                    out$top_correlations$p1[k], out$top_correlations$p2[k],
                    out$top_correlations$cor[k]))
  }
  
  invisible(out)
}

.ident_cor_heatmap <- function(Cp) {
  nm <- colnames(Cp)
  df <- expand.grid(p1 = nm, p2 = nm, stringsAsFactors = FALSE)
  df$cor <- as.vector(Cp)
  ggplot(df, aes(p1, p2, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_blank(), panel.grid = element_blank()) +
    labs(fill = "cor", title = "Posterior correlation")
}