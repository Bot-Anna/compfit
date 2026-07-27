# Local (derivative) sensitivity at a fit

Finite-difference sensitivity of each model output to the fitted
parameters, evaluated at the fitted point. \`type\` selects absolute
derivatives, relative elasticities, or semi-relative (parameter-scaled)
sensitivities.

## Usage

``` r
local_sensitivity(
  fit,
  type = c("relative", "absolute", "semi"),
  output = c("both", "summary", "trajectory"),
  summary_fun = mean,
  eps = 1e-04,
  data_dummy = NULL,
  progress = TRUE
)
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

- type:

  Sensitivity type: relative (elasticity), absolute, or semi.

- output:

  What to return: both, summary, or trajectory.

- summary_fun:

  Function summarising each output trajectory to a scalar.

- eps:

  Finite-difference step on the normalised scale.

- data_dummy:

  Optional dummy-data data frame.

- progress:

  Show a progress bar.

## Value

A list with \`summary\` and/or \`trajectory\` tidy data frames.

## Examples

``` r
if (FALSE) { # \dontrun{
ls <- local_sensitivity(fit, type = "relative")   # fit from fitCompartmentalModel()
} # }
```
