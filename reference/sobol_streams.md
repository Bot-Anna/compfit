# Sobol sensitivity of model outputs

Variance-based (Sobol) global sensitivity of a summary of each model
output stream to the fitted parameters.

## Usage

``` r
sobol_streams(
  x,
  n = 500,
  summary_fun = mean,
  streams = NULL,
  dataCombined = NULL,
  solver = NULL,
  progress = TRUE,
  ...
)
```

## Arguments

- x:

  A \`"compartmentalFit"\` or \`"compartmentalModel"\` object.

- n:

  Base sample size for the Sobol design.

- summary_fun:

  Function summarising each output trajectory to a scalar.

- streams:

  Optional subset of output streams to analyse.

- dataCombined:

  Optional override data (defaults to the fit's data).

- solver:

  Optional solver settings (defaults to the fit's solver).

- progress:

  Show a progress bar.

- ...:

  Passed to \[sensitivity::soboljansen()\].

## Value

A tidy data frame of first-order/total indices per output and parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
sob <- sobol_streams(fit, n = 500)
} # }
```
