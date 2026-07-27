# Posterior credible intervals (natural scale)

Posterior credible intervals (natural scale)

## Usage

``` r
posterior_intervals(fit, p = 0.95)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

- p:

  Central credible mass (default \`0.95\`).

## Value

A data frame with lower/upper interval bounds per parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
posterior_intervals(fit_b, p = 0.9)   # fit_b a method = "bayes" fit
} # }
```
