# Report on a Bayesian posterior

Prints a posterior summary and convergence diagnostics (R-hat / ESS)
with threshold flags, and returns the summary invisibly.

## Usage

``` r
posterior_report(fit, p = 0.95, rhat_thresh = 1.01, ess_thresh = 400)
```

## Arguments

- fit:

  A Bayesian \`"compartmentalFit"\` object.

- p:

  Central credible mass (default \`0.95\`).

- rhat_thresh:

  R-hat threshold above which a parameter is flagged.

- ess_thresh:

  Effective-sample-size threshold below which a parameter is flagged.

## Value

The posterior summary data frame, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
posterior_report(fit_b)   # fit_b a method = "bayes" fit
} # }
```
