# Sobol sensitivity of the loss

Variance-based (Sobol) global sensitivity of the loss function to the
fitted parameters, via \[sensitivity::soboljansen()\].

## Usage

``` r
sobol_loss(
  x,
  n = 1000,
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

- dataCombined:

  Optional override data (defaults to the fit's data).

- solver:

  Optional solver settings (defaults to the fit's solver).

- progress:

  Show a progress bar.

- ...:

  Passed to \[sensitivity::soboljansen()\].

## Value

A tidy data frame of first-order and total indices per parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
sobol_loss(fit, n = 500)   # fit or model from build_compartmental_model()
} # }
```
