# Plot a compartmental fit

One panel per data stream: observed points, the fitted trajectory, and
(optionally) uncertainty bands or spaghetti draws for Bayesian fits.

## Usage

``` r
plot_fit(
  fit,
  bands = TRUE,
  band_type = c("predictive", "mean", "both", "spaghetti"),
  n_draws = 200,
  n_rep = 5,
  spaghetti_draws = 100,
  spaghetti_predictive = FALSE,
  base_size = 11,
  ncol = NULL,
  data_dummy = NULL,
  palette = NULL,
  endpoint = NULL
)
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

- bands:

  Whether to draw uncertainty bands (Bayesian fits).

- band_type:

  Band style: predictive, mean, both, or spaghetti.

- n_draws:

  Posterior draws used for band construction.

- n_rep:

  Replicates per draw for predictive bands.

- spaghetti_draws:

  Number of spaghetti trajectories.

- spaghetti_predictive:

  If \`TRUE\`, spaghetti lines include observation noise.

- base_size:

  Base font size.

- ncol:

  Number of facet columns (default: automatic).

- data_dummy:

  Optional dummy-data data frame to overlay.

- palette:

  Colour scheme: \`"okabe"\` (default, colour-blind-safe) or
  \`"grey"\`/\`"grayscale"\` for monochrome output; also accepts a
  custom palette list (see \[cfit_palette\]). \`NULL\` honours
  \`options(compfit.palette=)\`.

- endpoint:

  Optional end year, to draw the fit projected beyond the fitted
  horizon. The central trajectory and, for a Bayesian fit, the bands or
  spaghetti lines all run to \`endpoint\`, while the observed data stays
  where it is – so the curve continues past the last data point with its
  uncertainty fanning out. \`NULL\` (default) stops at the fitted
  endpoint. See \[solve_and_evaluate()\] for the caveat about
  time-varying \`Functions\` beyond the data.

## Value

A list with per-stream \`plots\` and (if patchwork is available) a
combined \`grid\`, plus the plotting \`code\`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- plot_fit(fit, band_type = "both")       # fit from fitCompartmentalModel()
res$grid
res_bw <- plot_fit(fit, palette = "grey")      # monochrome
options(compfit.palette = "grey")              # ... or switch every plot

# project to 2040, with predictive bands widening past the data
plot_fit(fit, endpoint = 2040)$grid
} # }
```
