# ============================================================
# Prior-vs-posterior overlay plots, one panel per parameter.
#
#   plot_prior_posterior(fit)                 # grid of all parameters (ggplot)
#   plot_prior_posterior(fit, "g_HA")         # a single parameter
#
# Posterior comes from the draws (natural scale). Prior density is
# reconstructed from prior_spec (the same bounds/distributions the model used).
# Where the posterior is much narrower than the prior, the data was informative.
# ============================================================

# Return a data frame of (x, density) for a parameter's PRIOR on natural scale.
.prior_density_df <- function(spec, n = 400) {
  d <- spec$dist; a <- as.numeric(spec$args)
  lo <- spec$lower; hi <- spec$upper
  
  # Choose an x-range to evaluate the prior over.
  rng <- switch(d,
                Uniform   = c(a[1], a[2]),
                Normal    = c(a[1] - 4*a[2], a[1] + 4*a[2]),
                LogNormal = c(0, exp(a[1] + 4*a[2])),
                Beta      = c(0, 1),
                Gamma     = c(0, qgamma(0.999, shape = a[1], scale = a[2])),
                c(lo, hi))
  # Respect finite truncation bounds if present.
  if (is.finite(lo)) rng[1] <- max(rng[1], lo)
  if (is.finite(hi)) rng[2] <- min(rng[2], hi)
  
  x <- seq(rng[1], rng[2], length.out = n)
  dens <- switch(d,
                 Uniform   = dunif(x, a[1], a[2]),
                 Normal    = dnorm(x, a[1], a[2]),
                 LogNormal = dlnorm(x, a[1], a[2]),
                 Beta      = dbeta(x, a[1], a[2]),
                 Gamma     = dgamma(x, shape = a[1], scale = a[2]),
                 rep(NA_real_, length(x)))
  
  # Renormalize density over a truncated support so it integrates to 1 there.
  if (d != "Uniform" && (is.finite(lo) || is.finite(hi))) {
    keep <- rep(TRUE, length(x))
    if (is.finite(lo)) keep <- keep & x >= lo
    if (is.finite(hi)) keep <- keep & x <= hi
    x <- x[keep]; dens <- dens[keep]
    area <- sum(0.5 * (head(dens,-1) + tail(dens,-1)) * diff(x))
    if (is.finite(area) && area > 0) dens <- dens / area
  }
  data.frame(x = x, density = dens)
}

# Look up the prior spec for a parameter name (searches params then states).
.lookup_spec <- function(prior_spec, name) {
  if (!is.null(prior_spec$params[[name]])) return(prior_spec$params[[name]])
  if (!is.null(prior_spec$states[[name]])) return(prior_spec$states[[name]])
  NULL
}

#' Plot priors against posteriors
#'
#' One panel per fitted parameter overlaying the prior density and the posterior
#' histogram, with a credible interval marked.
#'
#' @param fit A Bayesian `"compartmentalFit"` object.
#' @param which Optional subset of parameter names to plot.
#' @param ncol Number of facet columns.
#' @param bins Histogram bin count.
#' @param cri Central credible mass to mark (default `0.95`).
#' @return A list with per-parameter `plots` and (if available) a combined `grid`.
#' @examples
#' \dontrun{
#' plot_prior_posterior(fit_b)   # fit_b a method = "bayes" fit
#' }
#' @export
plot_prior_posterior <- function(fit, which = NULL, ncol = 4, bins = 30,
                                 cri = 0.95) {
  .check_bayes(fit)
  draws <- posterior_draws(fit)          # natural scale, _n suffix stripped
  ps    <- fit$samples$prior_spec
  
  # Which parameters to plot. Default: all that have a prior spec (skip e.g.
  # sigma/phi which have no entry in prior_spec$params/states).
  candidates <- names(draws)
  if (!is.null(which)) candidates <- intersect(candidates, which)
  params <- candidates[vapply(candidates,
                              function(nm) !is.null(.lookup_spec(ps, nm)),
                              logical(1))]
  if (length(params) == 0)
    stop("None of the requested parameters have a prior in prior_spec.")
  
  alpha <- (1 - cri) / 2
  
  # Assemble long data frames for posteriors, priors, and per-parameter stats.
  post_list  <- list()
  prior_list <- list()
  stat_list  <- list()
  for (nm in params) {
    x <- draws[[nm]]
    post_list[[nm]] <- data.frame(parameter = nm, value = x)
    
    spec <- .lookup_spec(ps, nm)
    pd <- .prior_density_df(spec)
    prior_list[[nm]] <- data.frame(parameter = nm, x = pd$x, density = pd$density)
    
    qs <- quantile(x, c(alpha, 0.5, 1 - alpha), na.rm = TRUE)
    stat_list[[nm]] <- data.frame(
      parameter = nm,
      mean   = mean(x, na.rm = TRUE),
      median = qs[2],
      lower  = qs[1],
      upper  = qs[3]
    )
  }
  post_df  <- do.call(rbind, post_list)
  prior_df <- do.call(rbind, prior_list)
  stat_df  <- do.call(rbind, stat_list)
  
  ggplot() +
    # 95% credible interval as a shaded band
    geom_rect(data = stat_df,
              aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
              fill = "grey70", alpha = 0.25, inherit.aes = FALSE) +
    # posterior as a density (area-normalized so it's comparable to the prior)
    geom_histogram(data = post_df,
                   aes(x = value, y = after_stat(density)),
                   bins = bins, fill = "steelblue", alpha = 0.55,
                   colour = "white") +
    # prior as a line
    geom_line(data = prior_df, aes(x = x, y = density),
              colour = "darkred", linewidth = 0.7) +
    # median (solid) and mean (dashed) verticals
    geom_vline(data = stat_df, aes(xintercept = median),
               colour = "black", linewidth = 0.6) +
    geom_vline(data = stat_df, aes(xintercept = mean),
               colour = "black", linewidth = 0.6, linetype = "dashed") +
    facet_wrap(~ parameter, scales = "free", ncol = ncol) +
    theme_minimal() +
    labs(title = "Prior (red) vs posterior (blue)",
         subtitle = sprintf("median = solid, mean = dashed, %d%% CrI = grey band",
                            round(cri * 100)),
         x = NULL, y = "density") +
    theme(strip.text = element_text(size = 8),
          axis.text  = element_text(size = 6),
          plot.subtitle = element_text(size = 8))
}