# ============================================================
# plots.R
# One entry point -- plot_fit() -- that turns a compartmentalFit into a set of
# nicely designed, colour-blind-safe plots, and ALSO returns the R code that
# produced them so you can copy, edit the design, and re-run.
#
#   res <- plot_fit(fit)
#   res$plots      # named list of ggplot objects (one per data stream + dummies)
#   res$grid       # all panels arranged in a grid (one object you can ggsave)
#   res$code       # character string: the editable source of the plot builder
#
# Design: Okabe-Ito colour-blind-safe palette + a minimal theme.
#   data points    : dark blue-grey
#   model line     : vermillion
#   predictive band: translucent vermillion (Bayes only)
#   censored limit : orange downward triangle (a 2nd marker; the 0 is plotted too)
#   dummy series   : muted grey
#
# Bands: for a Bayesian fit, plot_fit() draws a posterior-PREDICTIVE band by
# default (parameter uncertainty + observation noise), computed from a thinned
# set of draws. Turn off with bands = FALSE; switch to the mean-trajectory band
# (parameter uncertainty only) with band_type = "mean".
# ============================================================

# ---- Palette / theme (edit here to restyle everything) ---------------------
#' compfit colour palette
#'
#' Named list of hex colours (Okabe-Ito based) used across compfit plots: edit
#' this to restyle data, model, band, censor and dummy elements.
#'
#' @format A named list of hex colour strings.
#' @examples
#' cfit_palette$data
#' @export
cfit_palette <- list(
  data    = "#0072B2",  # blue        (Okabe-Ito)
  model   = "#D55E00",  # vermillion  (Okabe-Ito)
  band    = "#D55E00",  # vermillion, used translucently
  censor  = "#666666",  # mid grey -- censor/interval/asym markers (a neutral tone,
                        # distinguished from the data/model by shape and colour)
  dummy   = "#999999"   # neutral grey
)

#' compfit grayscale palette
#'
#' A grayscale counterpart to [cfit_palette], for print/monochrome output. Same
#' named slots (`data`, `model`, `band`, `censor`, `dummy`); elements are kept
#' apart by shade (and by shape/linetype in the plots). Select it with
#' `plot_fit(..., palette = "grey")` or globally via
#' `options(compfit.palette = "grey")`.
#'
#' @format A named list of hex colour strings.
#' @examples
#' cfit_palette_grey$data
#' @export
cfit_palette_grey <- list(
  data    = "#000000",  # black       (observed points)
  model   = "#737373",  # mid grey    (fitted curve)
  band    = "#737373",  # mid grey, used translucently
  censor  = "#3D3D3D",  # dark grey   (censoring markers)
  dummy   = "#B3B3B3"   # light grey  (dummy overlay)
)

# TRUE when a scheme NAME (string) selects the grayscale palette.
.cfit_is_grey_name <- function(x)
  tolower(as.character(x)[1]) %in% c("grey", "gray", "greyscale", "grayscale", "bw", "mono")

# Resolve a `palette` argument to a palette list. A list is used as-is (custom
# palette); a string picks the built-in scheme; NULL falls back to the global
# option `compfit.palette` (default "okabe"). So one option switches every plot,
# while a per-call `palette=` overrides it locally.
.cfit_resolve_palette <- function(palette = NULL) {
  if (is.list(palette)) return(palette)
  name <- if (is.null(palette)) getOption("compfit.palette", "okabe") else palette
  if (.cfit_is_grey_name(name)) cfit_palette_grey else cfit_palette
}

# The active palette per the global option (for plots without a palette= arg).
.cfit_active_palette <- function() .cfit_resolve_palette(NULL)

#' compfit ggplot2 theme
#'
#' A minimal ggplot2 theme used across compfit plots.
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + cfit_theme()
#' @export
cfit_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = base_size, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = "grey90"),
      axis.line        = ggplot2::element_line(linewidth = 0.3, colour = "grey60"),
      axis.text        = ggplot2::element_text(size = base_size - 3, colour = "grey30"),
      plot.margin      = ggplot2::margin(6, 8, 6, 6),
      legend.position  = "none"
    )
}

# ---- Single-stream plot builder (this source is returned as $code) ---------
# Kept deliberately self-contained and readable: it is the "recipe" users edit.
.cfit_stream_plot <- function(stream, label, data_points, model_eval,
                              band_df = NULL, mean_band_df = NULL,
                              spaghetti_df = NULL, cens_df = NULL, lcens_df = NULL,
                              interval_df = NULL, asym_df = NULL,
                              pal = cfit_palette, base_size = 11) {
  p <- ggplot2::ggplot()

  # Outer band: predictive (data-point) uncertainty -- lightest, drawn first.
  if (!is.null(band_df)) {
    p <- p + ggplot2::geom_ribbon(
      data = band_df,
      ggplot2::aes(x = date, ymin = lo, ymax = hi),
      fill = pal$band, alpha = 0.15)
  }

  # Inner band: fitted-curve (parameter) uncertainty -- darker, drawn on top.
  # Narrower than the predictive band; shows how well-determined the curve is.
  if (!is.null(mean_band_df)) {
    p <- p + ggplot2::geom_ribbon(
      data = mean_band_df,
      ggplot2::aes(x = date, ymin = lo, ymax = hi),
      fill = pal$band, alpha = 0.38)
  }

  # Spaghetti: one thin, transparent line per posterior draw (drawn beneath the
  # central trajectory). Preserves cross-time structure that bands collapse.
  if (!is.null(spaghetti_df) && nrow(spaghetti_df) > 0) {
    p <- p + ggplot2::geom_line(
      data = spaghetti_df,
      ggplot2::aes(x = date, y = value, group = draw),
      colour = pal$model, linewidth = 0.18, alpha = 0.12)
  }
  
  # Model trajectory.
  p <- p + ggplot2::geom_line(
    data = model_eval,
    ggplot2::aes(x = date, y = .data[[stream]]),
    colour = pal$model, linewidth = 0.7)
  
  # Observed data. Non-observed years (censored / interval / asym / missing) are
  # NA here so no dot is drawn -- the limit/interval markers carry that info.
  p <- p + ggplot2::geom_point(
    data = data_points,
    ggplot2::aes(x = date, y = .data[[stream]]),
    colour = pal$data, size = 1.4, na.rm = TRUE)
  
  # Censoring markers: a horizontal cap at the limit L plus an arrow in the OPEN
  # direction -- down for "<=L / <L" (at most L), up for ">=L / >L" (at least L).
  # The arrow length is a fixed fraction of the panel's y-range so it scales.
  cens_arrow <- grid::arrow(length = grid::unit(0.05, "inches"), type = "open")
  y_ref  <- c(model_eval[[stream]], data_points[[stream]],
              cens_df$limit, lcens_df$limit, interval_df$low, interval_df$high)
  y_span <- diff(range(y_ref, na.rm = TRUE))
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  arr_len <- 0.10 * y_span

  # Upper bound (<L / <=L): cap at L + arrow pointing DOWN.
  if (!is.null(cens_df) && any(!is.na(cens_df$limit))) {
    p <- p +
      ggplot2::geom_point(data = cens_df, ggplot2::aes(x = date, y = limit),
        colour = pal$censor, shape = 95, size = 4, na.rm = TRUE) +          # cap "-"
      ggplot2::geom_segment(data = cens_df,
        ggplot2::aes(x = date, xend = date, y = limit, yend = limit - arr_len),
        colour = pal$censor, linewidth = 0.6, arrow = cens_arrow, na.rm = TRUE)
  }

  # Lower bound (>L / >=L): cap at L + arrow pointing UP.
  if (!is.null(lcens_df) && any(!is.na(lcens_df$limit))) {
    p <- p +
      ggplot2::geom_point(data = lcens_df, ggplot2::aes(x = date, y = limit),
        colour = pal$censor, shape = 95, size = 4, na.rm = TRUE) +          # cap "-"
      ggplot2::geom_segment(data = lcens_df,
        ggplot2::aes(x = date, xend = date, y = limit, yend = limit + arr_len),
        colour = pal$censor, linewidth = 0.6, arrow = cens_arrow, na.rm = TRUE)
  }

  # Interval [A,B]: a vertical bracket from A to B with tick ends.
  if (!is.null(interval_df) && any(!is.na(interval_df$low))) {
    p <- p +
      ggplot2::geom_linerange(data = interval_df,
        ggplot2::aes(x = date, ymin = low, ymax = high),
        colour = pal$censor, linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = interval_df, ggplot2::aes(x = date, y = low),
        colour = pal$censor, shape = 95, size = 3, na.rm = TRUE) +   # "-" tick
      ggplot2::geom_point(data = interval_df, ggplot2::aes(x = date, y = high),
        colour = pal$censor, shape = 95, size = 3, na.rm = TRUE)
  }

  # Asymmetric A +/- dev: a capped interval line over [A, A + dir*dev] (like the
  # [A,B] marker) plus a data-sized grey dot at the recorded value A.
  if (!is.null(asym_df) && any(!is.na(asym_df$value))) {
    p <- p +
      ggplot2::geom_linerange(data = asym_df,
        ggplot2::aes(x = date,
                     ymin = pmin(value, value + dir * dev),
                     ymax = pmax(value, value + dir * dev)),
        colour = pal$censor, linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = asym_df,
        ggplot2::aes(x = date, y = pmin(value, value + dir * dev)),
        colour = pal$censor, shape = 95, size = 3, na.rm = TRUE) +   # bottom "-" cap
      ggplot2::geom_point(data = asym_df,
        ggplot2::aes(x = date, y = pmax(value, value + dir * dev)),
        colour = pal$censor, shape = 95, size = 3, na.rm = TRUE) +   # top "-" cap
      ggplot2::geom_point(data = asym_df, ggplot2::aes(x = date, y = value),
        colour = pal$censor, size = 1.4, na.rm = TRUE)               # grey dot at A
  }

  # Adaptive x-axis: pick a year step so there are ~6-8 date labels regardless
  # of horizon (avoids a wall of labels on long spans).
  yrs <- as.numeric(format(range(data_points$date, na.rm = TRUE), "%Y"))
  span <- diff(yrs)
  step <- if (span <= 8) 1 else if (span <= 16) 2 else if (span <= 30) 5 else 10
  brks <- seq(as.Date(sprintf("%d-12-31", yrs[1])),
              as.Date(sprintf("%d-12-31", yrs[2])), by = paste(step, "years"))
  
  p +
    ggplot2::scale_x_date(breaks = brks, date_labels = "%Y") +
    ggplot2::scale_y_continuous(
      n.breaks = 5,
      labels = function(v) format(v, big.mark = ",", scientific = FALSE,
                                  trim = TRUE),
      # Generous headroom/footroom so points and the curve don't touch the panel
      # edge: pad ~12% of the data range at the top, ~8% at the bottom.
      expand = ggplot2::expansion(mult = c(0.08, 0.12))) +
    ggplot2::labs(title = label, x = NULL, y = NULL) +
    cfit_theme(base_size)
}

# ---- Family-correct observation simulation ---------------------------------
# Given a stream's family, a mean trajectory mu (vector over the grid), the
# draw's noise parameter, and the stream scale, simulate ONE predictive
# observation vector from that family's actual distribution. n_rep>1 draws
# several per call to smooth the quantiles.
.cfit_simulate_obs <- function(family, mu, noise, scale, n_rep = 1) {
  mu <- pmax(mu, 0)                      # guard tiny negatives from the solver
  out <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    out[[r]] <- switch(family,
                       gaussian  = rnorm(length(mu),  mean = mu, sd = noise * scale),
                       lognormal = rlnorm(length(mu), meanlog = log(pmax(mu, 1e-9)), sdlog = noise),
                       poisson   = rpois(length(mu),  lambda = mu),
                       # R's rnbinom uses (size = dispersion phi, mu = mean): Var = mu + mu^2/phi.
                       negbin    = rnbinom(length(mu), size = noise, mu = mu),
                       rnorm(length(mu), mean = mu, sd = noise * scale)   # fallback: gaussian
    )
  }
  do.call(rbind, out)                    # n_rep x length(mu)
}

# ---- Posterior-predictive bands (Bayes) ------------------------------------
# Two modes:
#   band_type = "mean"        -> spread of the model trajectory mu across draws
#                                (parameter uncertainty only; family-independent)
#   band_type = "predictive"  -> spread of SIMULATED observations from each
#                                stream's actual likelihood (parameter +
#                                observation noise; correct per family).
# Returns a named list of data.frame(date, lo, hi) on the trajectory grid.
# Dispersion phi for stream i from a single posterior draw, tolerating BOTH the
# Julia naming (phi<i>, e.g. phi1) and the Stan naming (phi_<i>, e.g. phi_1);
# without this the Stan negbin band found no phi -> NA noise -> no data-stream
# CrI. NA when absent (a stream that is not negbin).
.cfit_phi_of <- function(row, i) {
  for (pn in c(paste0("phi_", i), paste0("phi", i)))
    if (pn %in% names(row)) return(row[[pn]])
  NA
}

.cfit_predictive_bands <- function(fit, n_draws = 200, band_type = "predictive",
                                   probs = c(0.025, 0.975), n_rep = 5,
                                   data_dummy = NULL) {
  if (is.null(fit$samples)) return(NULL)
  if (!exists("solve_and_evaluate"))
    stop("plot_fit() needs solve_and_evaluate() in scope to compute Bayes bands.")
  
  draws   <- posterior_draws(fit)
  S       <- nrow(draws)
  idx     <- if (S > n_draws) round(seq(1, S, length.out = n_draws)) else seq_len(S)
  streams <- fit$data$names_data_points
  dummies <- if (!is.null(data_dummy)) data_dummy$Formula else character(0)
  
  # Per-stream family + scale, for the predictive simulation.
  like_specs  <- fit$samples$like_specs
  fam_of      <- function(i) if (is.null(like_specs)) "gaussian" else like_specs[[i]]$family
  scale_cols  <- fit$samples$scale_cols
  if (is.null(scale_cols)) {
    # Fallback for fits made before scale_cols was attached: recompute it with
    # the SAME rule .fit_bayes uses -- each stream's mean absolute OBSERVED value
    # (obs_mask-gated), defaulting to 1 for empty/degenerate streams.
    ym   <- fit$data$matrix_data_points
    obsm <- fit$data$obs_mask
    if (is.null(obsm)) obsm <- matrix(1L, nrow(ym), ncol(ym))
    ym[is.na(ym)] <- 0
    scale_cols <- vapply(seq_len(ncol(ym)), function(j) {
      vals <- ym[obsm[, j] == 1, j]
      m <- mean(abs(vals))
      if (!is.finite(m) || m == 0) 1 else m
    }, numeric(1))
  }
  
  acc <- setNames(vector("list", length(c(streams, dummies))), c(streams, dummies))
  grid_dates  <- NULL
  n_failed    <- 0L

  for (s in idx) {
    row <- draws[s, , drop = FALSE]
    pe  <- .cfit_point_from_draw(fit, row)
    out <- tryCatch(solve_and_evaluate(fit, pe$initial_state, pe$parms, data_dummy),
                    error = function(e) NULL)
    if (is.null(out)) { n_failed <- n_failed + 1L; next }
    if (is.null(grid_dates)) grid_dates <- out$evaluation$date
    
    for (i in seq_along(streams)) {
      st  <- streams[i]
      mu  <- out$evaluation[[st]]
      if (band_type == "mean") {
        acc[[st]] <- rbind(acc[[st]], mu)               # trajectory only
      } else {
        fam   <- fam_of(i)
        # noise parameter per family: sigma for gaussian/lognormal; the stream's
        # phi for negbin; none for poisson.
        noise <- switch(fam,
                        gaussian  = row[["sigma"]],
                        lognormal = row[["sigma"]],
                        negbin    = .cfit_phi_of(row, i),
                        poisson   = NA, NA)
        sims <- .cfit_simulate_obs(fam, mu, noise, scale_cols[i], n_rep = n_rep)
        acc[[st]] <- rbind(acc[[st]], sims)              # n_rep rows per draw
      }
    }
    # Dummy series have NO likelihood, so they get TRAJECTORY (mean) bands only,
    # regardless of band_type -- the spread of the dummy value across draws.
    for (st in dummies) {
      v <- out$evaluation[[st]]
      if (!is.null(v)) acc[[st]] <- rbind(acc[[st]], v)
    }
  }
  if (n_failed > 0)
    warning(sprintf(
      "%d of %d posterior draws failed to solve (ODE error); excluded from band.",
      n_failed, length(idx)), call. = FALSE)
  if (is.null(grid_dates)) return(NULL)

  bands <- list()
  for (st in c(streams, dummies)) {
    M <- acc[[st]]
    if (is.null(M)) next
    lo <- apply(M, 2, quantile, probs[1], na.rm = TRUE)
    hi <- apply(M, 2, quantile, probs[2], na.rm = TRUE)
    df <- data.frame(date = grid_dates, lo = lo, hi = hi)
    df <- df[is.finite(df$lo) & is.finite(df$hi), ]
    if (nrow(df) > 0) bands[[st]] <- df
  }
  bands
}

# ---- Spaghetti: individual joint-draw trajectories -------------------------
# Returns a named list (per stream) of long data frames (date, value, draw) with
# one trajectory per posterior draw. This preserves the joint, cross-time
# structure that pointwise bands collapse away -- a more honest uncertainty view.
# Uses the MEAN trajectory mu per draw (no observation-noise layer); set
# predictive = TRUE to instead overlay one simulated observation series per draw.
.cfit_spaghetti_lines <- function(fit, n_draws = 100, predictive = FALSE,
                                  n_rep = 1, data_dummy = NULL) {
  if (is.null(fit$samples)) return(NULL)
  if (!exists("solve_and_evaluate"))
    stop("plot_fit() needs solve_and_evaluate() in scope for spaghetti plots.")
  
  draws   <- posterior_draws(fit)
  S       <- nrow(draws)
  idx     <- if (S > n_draws) round(seq(1, S, length.out = n_draws)) else seq_len(S)
  streams <- fit$data$names_data_points
  dummies <- if (!is.null(data_dummy)) data_dummy$Formula else character(0)
  
  like_specs <- fit$samples$like_specs
  fam_of     <- function(i) if (is.null(like_specs)) "gaussian" else like_specs[[i]]$family
  scale_cols <- fit$samples$scale_cols
  if (is.null(scale_cols)) scale_cols <- rep(1, length(streams))
  
  acc <- setNames(lapply(c(streams, dummies), function(.) list()),
                  c(streams, dummies))
  grid_dates <- NULL
  n_failed   <- 0L

  for (s in idx) {
    row <- draws[s, , drop = FALSE]
    pe  <- .cfit_point_from_draw(fit, row)
    out <- tryCatch(solve_and_evaluate(fit, pe$initial_state, pe$parms, data_dummy),
                    error = function(e) NULL)
    if (is.null(out)) { n_failed <- n_failed + 1L; next }
    if (is.null(grid_dates)) grid_dates <- out$evaluation$date
    
    for (i in seq_along(streams)) {
      st <- streams[i]
      mu <- out$evaluation[[st]]
      if (!predictive) {
        acc[[st]][[length(acc[[st]]) + 1]] <-
          data.frame(date = out$evaluation$date, value = mu, draw = s)
      } else {
        fam <- fam_of(i)
        noise <- switch(fam,
                        gaussian = row[["sigma"]], lognormal = row[["sigma"]],
                        negbin = .cfit_phi_of(row, i),
                        poisson = NA, NA)
        sims <- .cfit_simulate_obs(fam, mu, noise, scale_cols[i], n_rep = n_rep)
        for (r in seq_len(nrow(sims)))
          acc[[st]][[length(acc[[st]]) + 1]] <-
          data.frame(date = out$evaluation$date, value = sims[r, ],
                     draw = s * 1000L + r)
      }
    }
    
    # Dummy formulas: trajectory only (no likelihood -> no predictive noise).
    for (st in dummies) {
      v <- out$evaluation[[st]]
      if (!is.null(v))
        acc[[st]][[length(acc[[st]]) + 1]] <-
          data.frame(date = out$evaluation$date, value = v, draw = s)
    }
  }
  if (n_failed > 0)
    warning(sprintf(
      "%d of %d posterior draws failed to solve (ODE error); excluded from spaghetti.",
      n_failed, length(idx)), call. = FALSE)
  if (is.null(grid_dates)) return(NULL)

  lapply(acc, function(lst) if (length(lst)) do.call(rbind, lst) else NULL)
}

# Reconstruct a natural-scale (initial_state, parms) from one posterior draw row,
# reusing the model's recovery machinery so ordering matches the ODE.
.cfit_point_from_draw <- function(fit, row) {
  sap <- fit$model$sap
  tg  <- fit$time_grid
  nm  <- names(row)
  vals <- as.numeric(row[1, ])
  names(vals) <- nm
  
  # Drop likelihood hyperparameters (sigma, phi1/phi_1, ...) — observation-model
  # quantities the ODE does not take. Dropping anything else would silently corrupt
  # the parameter vector, so we validate the remainder against sap. The `_` covers
  # the Stan backend's phi_<i> naming as well as Julia's phi<i>.
  hyper <- grepl("^sigma$|^phi_?[0-9]*$", nm)
  best_par <- vals[!hyper]

  expected <- c(names(sap$params_fitted), names(sap$states_fitted))
  if (!setequal(names(best_par), expected))
    warning(sprintf(
      ".cfit_point_from_draw: unexpected columns after dropping hyperparameters.\n  expected: [%s]\n  got:      [%s]",
      paste(sort(expected), collapse = ", "),
      paste(sort(names(best_par)), collapse = ", ")
    ))

  .recover_solution(best_par, sap, tg)
}

# ---- Top-level entry point -------------------------------------------------
#' Plot a compartmental fit
#'
#' One panel per data stream: observed points, the fitted trajectory, and
#' (optionally) uncertainty bands or spaghetti draws for Bayesian fits.
#'
#' @param fit A `"compartmentalFit"` object.
#' @param bands Whether to draw uncertainty bands (Bayesian fits).
#' @param band_type Band style: predictive, mean, both, or spaghetti.
#' @param n_draws Posterior draws used for band construction.
#' @param n_rep Replicates per draw for predictive bands.
#' @param spaghetti_draws Number of spaghetti trajectories.
#' @param spaghetti_predictive If `TRUE`, spaghetti lines include observation
#'   noise.
#' @param base_size Base font size.
#' @param ncol Number of facet columns (default: automatic).
#' @param data_dummy Optional dummy-data data frame to overlay.
#' @param palette Colour scheme: `"okabe"` (default, colour-blind-safe) or
#'   `"grey"`/`"grayscale"` for monochrome output; also accepts a custom palette
#'   list (see [cfit_palette]). `NULL` honours `options(compfit.palette=)`.
#' @return A list with per-stream `plots` and (if patchwork is available) a
#'   combined `grid`, plus the plotting `code`.
#' @examples
#' \dontrun{
#' res <- plot_fit(fit, band_type = "both")       # fit from fitCompartmentalModel()
#' res$grid
#' res_bw <- plot_fit(fit, palette = "grey")      # monochrome
#' options(compfit.palette = "grey")              # ... or switch every plot
#' }
#' @export
plot_fit <- function(fit, bands = TRUE,
                     band_type = c("predictive", "mean", "both", "spaghetti"),
                     n_draws = 200, n_rep = 5, spaghetti_draws = 100,
                     spaghetti_predictive = FALSE,
                     base_size = 11, ncol = NULL, data_dummy = NULL,
                     palette = NULL) {
  .is_fit(fit)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_fit() requires ggplot2.")
  band_type <- match.arg(band_type)
  # Colour scheme: "okabe" (default) or "grey"/"grayscale", or a custom palette
  # list. NULL honours options(compfit.palette=). Threaded into every panel.
  pal <- .cfit_resolve_palette(palette)
  
  if (!exists("solve_and_evaluate"))
    stop("plot_fit() expects solve_and_evaluate() in scope (sourced from main.R).")
  
  # Model trajectory at the central estimate (MLE point, or posterior means).
  if (fit$method == "bayes") {
    pe <- .cfit_point_from_draw(fit, as.data.frame(t(posterior_means(fit))))
  } else {
    pe <- fit$point
  }
  model_eval <- solve_and_evaluate(fit, pe$initial_state, pe$parms, data_dummy)$evaluation
  
  data_points <- fit$data$data_points
  labels <- setNames(fit$data$data_combined$Label,
                     fit$data$data_combined$Formula)
  streams <- fit$data$names_data_points

  # Only draw a data dot on genuinely OBSERVED years. The plotting frame renders
  # censored / interval / asym / missing cells as 0; blank those to NA (per
  # stream, via obs_mask) so no dot appears -- the limit/interval markers show
  # what is known there instead.
  obsm <- fit$data$obs_mask
  if (!is.null(obsm) && nrow(data_points) == nrow(obsm)) {
    for (i in seq_along(streams)) {
      st <- streams[i]
      if (st %in% names(data_points) && i <= ncol(obsm))
        data_points[[st]][obsm[, i] != 1] <- NA
    }
  }
  
  # Uncertainty overlay (Bayes only): quantile band(s) or spaghetti lines.
  #   "mean"       -> fitted-curve (parameter) uncertainty only
  #   "predictive" -> data-point uncertainty (param + observation noise)
  #   "both"       -> BOTH: inner mean band + outer predictive band
  #   "spaghetti"  -> per-draw trajectories
  band_list      <- NULL   # the primary band (predictive, or mean)
  mean_band_list <- NULL   # the inner curve-uncertainty band (for "both")
  spaghetti_list <- NULL
  if (isTRUE(bands) && fit$method == "bayes") {
    if (band_type == "spaghetti") {
      spaghetti_list <- .cfit_spaghetti_lines(
        fit, n_draws = spaghetti_draws,
        predictive = spaghetti_predictive, n_rep = n_rep,
        data_dummy = data_dummy)
    } else if (band_type == "both") {
      band_list      <- .cfit_predictive_bands(
        fit, n_draws = n_draws, band_type = "predictive", n_rep = n_rep,
        data_dummy = data_dummy)
      mean_band_list <- .cfit_predictive_bands(
        fit, n_draws = n_draws, band_type = "mean", n_rep = n_rep,
        data_dummy = data_dummy)
    } else {
      band_list <- .cfit_predictive_bands(
        fit, n_draws = n_draws, band_type = band_type, n_rep = n_rep,
        data_dummy = data_dummy)
    }
  }
  
  # Censored-limit frames per stream.
  # For Bayesian fits: apply the same discrete-family limit adjustment (.fit_bayes
  # uses L-1 for strict "<L" on poisson/negbin) so the plotted marker matches
  # what the sampler actually saw as the censoring threshold.
  # Censored-limit frames per stream. Upper bound (<L) and lower bound (>L) each
  # get their own marker (downward vs upward triangle).
  cens_mask  <- fit$data$cens_mask
  limit_mat  <- fit$data$limit_mat
  lcens_mask <- fit$data$lcens_mask
  llimit_mat <- fit$data$llimit_mat
  # For Bayesian fits with discrete families, strict "<L" is shown at L-1 to
  # match what the sampler used (P(Y<L) = P(Y<=L-1) on integer support).
  limit_mat_plot <- limit_mat
  if (fit$method == "bayes" && !is.null(fit$samples$like_specs)) {
    inc_mask_plot <- fit$data$inc_mask
    for (j in seq_along(streams)) {
      if (fit$samples$like_specs[[j]]$family %in% .DISCRETE_FAMILIES) {
        strict_j <- which(cens_mask[, j] == 1 &
                            (is.na(inc_mask_plot[, j]) | inc_mask_plot[, j] == 0))
        limit_mat_plot[strict_j, j] <- limit_mat_plot[strict_j, j] - 1
      }
    }
  }
  dates <- data_points$date
  cens_frame <- function(j) {
    L <- limit_mat_plot[, j]; L[cens_mask[, j] != 1] <- NA
    data.frame(date = dates, limit = L)
  }
  lcens_frame <- function(j) {
    if (is.null(lcens_mask)) return(data.frame(date = dates, limit = NA_real_))
    L <- llimit_mat[, j]; L[lcens_mask[, j] != 1] <- NA
    data.frame(date = dates, limit = L)
  }
  # Interval / asymmetric frames (NULL for fits made before the grammar existed).
  interval_mask <- fit$data$interval_mask
  asym_mask     <- fit$data$asym_mask
  interval_frame <- function(j) {
    if (is.null(interval_mask))
      return(data.frame(date = dates, low = NA_real_, high = NA_real_))
    lo <- fit$data$ilow_mat[, j]; hi <- fit$data$iupp_mat[, j]
    keep <- interval_mask[, j] == 1; lo[!keep] <- NA; hi[!keep] <- NA
    data.frame(date = dates, low = lo, high = hi)
  }
  asym_frame <- function(j) {
    if (is.null(asym_mask))
      return(data.frame(date = dates, value = NA_real_, dev = NA_real_, dir = NA_real_))
    v <- fit$data$asym_val_mat[, j]; dv <- fit$data$asym_dev_mat[, j]; dr <- fit$data$asym_dir_mat[, j]
    keep <- asym_mask[, j] == 1; v[!keep] <- NA; dv[!keep] <- NA; dr[!keep] <- NA
    data.frame(date = dates, value = v, dev = dv, dir = dr)
  }

  plots <- list()
  for (i in seq_along(streams)) {
    st  <- streams[i]
    lab <- labels[[st]]; if (is.null(lab) || is.na(lab)) lab <- st
    plots[[st]] <- .cfit_stream_plot(
      stream = st, label = lab,
      data_points = data_points, model_eval = model_eval,
      band_df      = if (!is.null(band_list))      band_list[[st]]      else NULL,
      mean_band_df = if (!is.null(mean_band_list)) mean_band_list[[st]] else NULL,
      spaghetti_df = if (!is.null(spaghetti_list)) spaghetti_list[[st]] else NULL,
      cens_df  = cens_frame(i),
      lcens_df = lcens_frame(i),
      interval_df = interval_frame(i),
      asym_df     = asym_frame(i),
      pal = pal, base_size = base_size)
  }
  
  # Dummy (display-only) series: model trajectory + its uncertainty. Dummies have
  # no likelihood, so they get the TRAJECTORY (mean) band -- prefer the inner mean
  # band ("both"), else the primary band list -- plus spaghetti when chosen.
  if (!is.null(data_dummy)) {
    data_dummy   <- .default_label(data_dummy)   # missing Label -> Formula
    dummy_labels <- setNames(data_dummy$Label, data_dummy$Formula)
    for (st in data_dummy$Formula) {
      lab <- dummy_labels[[st]]; if (is.null(lab) || is.na(lab)) lab <- st
      p <- ggplot2::ggplot()
      # Band: prefer the inner mean band ("both"), else the primary band list.
      bdf <- if (!is.null(mean_band_list)) mean_band_list[[st]]
             else if (!is.null(band_list)) band_list[[st]] else NULL
      if (!is.null(bdf)) {
        p <- p + ggplot2::geom_ribbon(
          data = bdf, ggplot2::aes(x = date, ymin = lo, ymax = hi),
          fill = pal$dummy, alpha = 0.22)
      }
      sp <- if (!is.null(spaghetti_list)) spaghetti_list[[st]] else NULL
      if (!is.null(sp) && nrow(sp) > 0) {
        p <- p + ggplot2::geom_line(
          data = sp, ggplot2::aes(x = date, y = value, group = draw),
          colour = pal$dummy, linewidth = 0.18, alpha = 0.12)
      }
      plots[[st]] <- p +
        ggplot2::geom_line(data = model_eval,
                           ggplot2::aes(x = date, y = .data[[st]]),
                           colour = pal$dummy, linewidth = 0.7) +
        ggplot2::labs(title = lab, x = NULL, y = NULL) +
        cfit_theme(base_size)
    }
  }
  
  # Arrange into a grid. Page GROWS with panel count so each panel keeps a
  # constant, readable size (option A) instead of being squashed onto a fixed
  # page. Per-panel target ~3.4in wide x 2.5in tall; long horizons get a touch
  # more width.
  n_panels <- length(plots)
  if (is.null(ncol)) ncol <- min(5, max(1, ceiling(sqrt(n_panels))))
  nrow <- ceiling(n_panels / ncol)
  
  horizon <- fit$time_grid$endpoint - fit$time_grid$startpoint + 1
  panel_w <- 3.4 + (horizon > 15) * 0.6        # wider panels for long spans
  panel_h <- 2.5
  width  <- ncol * panel_w
  height <- nrow * panel_h
  
  grid <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    grid <- patchwork::wrap_plots(plots, ncol = ncol)
  }
  
  # Editable code: the builder's own source + the theme/palette, so the user can
  # paste, tweak colours/geoms, and rebuild any panel.
  code <- paste(
    "# --- Editable plot recipe (requires `fit` and ggplot2 in scope) ---",
    "# Palette and theme (the palette used for this plot):",
    paste("cfit_palette <-", paste(deparse(pal), collapse = "\n")),
    paste(deparse(cfit_theme), collapse = "\n"),
    "# Per-stream builder (edit colours, geoms, theme here):",
    paste(deparse(.cfit_stream_plot), collapse = "\n"),
    "# Example: rebuild one panel (pass pal = cfit_palette to keep this scheme)",
    'p <- .cfit_stream_plot("X3", "My label", fit$data$data_points,',
    '                       solve_and_evaluate(fit, fit$point$initial_state, fit$point$parms)$evaluation,',
    '                       pal = cfit_palette)',
    "print(p)",
    sep = "\n\n")
  
  list(plots = plots, grid = grid, code = code,
       width = width, height = height, ncol = ncol, nrow = nrow)
}