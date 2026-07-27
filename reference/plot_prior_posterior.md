# Plot priors against posteriors

One panel per fitted parameter overlaying the prior density and the
posterior histogram, with a credible interval marked.

## Usage

``` r
plot_prior_posterior(fit, which = NULL, ncol = 4, bins = 30, cri = 0.95)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

- which:

  Optional subset of parameter names to plot.

- ncol:

  Number of facet columns.

- bins:

  Histogram bin count.

- cri:

  Central credible mass to mark (default \`0.95\`).

## Value

A list with per-parameter \`plots\` and (if available) a combined
\`grid\`.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_prior_posterior(fit_b)   # fit_b a method = "bayes" fit
} # }
```
