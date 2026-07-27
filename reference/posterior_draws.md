# Posterior draws (natural scale)

Returns the posterior draws with each normalised \`name_n\` column
mapped back to its natural scale \`name\` on \`\[lo, hi\]\`;
non-normalised columns (e.g. \`sigma\`, \`phi\`) pass through unchanged.

## Usage

``` r
posterior_draws(fit)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

## Value

A data frame of natural-scale draws.

## Examples

``` r
if (FALSE) { # \dontrun{
d <- posterior_draws(fit_b)   # fit_b a method = "bayes" fit
head(d)
} # }
```
