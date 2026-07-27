# Posterior summary (natural scale)

Per-parameter posterior summary (mean, sd, rhat, ess) with means/sds
mapped to the natural scale.

## Usage

``` r
posterior_summary(fit)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

## Value

A data frame, one row per parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
posterior_summary(fit_b)   # fit_b a method = "bayes" fit
} # }
```
