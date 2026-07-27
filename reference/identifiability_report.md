# Practical identifiability report

Diagnoses practical identifiability of a Bayesian fit:
prior-to-posterior narrowing, posterior coverage of the prior range, and
posterior parameter correlations, flagging weakly-identified or
correlated parameters.

## Usage

``` r
identifiability_report(
  fit,
  components = c("narrowing", "coverage", "correlation"),
  cor_threshold = 0.8,
  coverage_threshold = 0.7,
  cri = 0.95,
  plots = TRUE
)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

- components:

  Which diagnostics to run: narrowing, coverage, correlation.

- cor_threshold:

  Absolute correlation above which a pair is flagged.

- coverage_threshold:

  Posterior/prior coverage above which a parameter is flagged as poorly
  narrowed.

- cri:

  Central credible mass used for the diagnostics (default \`0.95\`).

- plots:

  If \`TRUE\`, include diagnostic plots.

## Value

Invisibly, a list of the diagnostic tables (and plots).

## Examples

``` r
if (FALSE) { # \dontrun{
identifiability_report(fit_b)   # fit_b a method = "bayes" fit
} # }
```
