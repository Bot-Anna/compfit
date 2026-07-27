# Compare a Bayesian fit to an MLE fit

Tabulates posterior mean and median against the MLE point estimate, with
relative differences for each.

## Usage

``` r
compare_to_mle(fit, fit_mle)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

- fit_mle:

  A maximum-likelihood \`"compartmentalFit"\` object.

## Value

A data frame with \`parameter\`, \`bayes_mean\`, \`bayes_median\`,
\`mle\`, and relative-difference columns.

## Examples

``` r
if (FALSE) { # \dontrun{
compare_to_mle(fit_b, fit_mle)   # Bayesian vs maximum-likelihood fit
} # }
```
