# Write fitted estimates back into a model sheet

Returns a copy of the model parameter sheet with fitted
States/Parameters replaced by their point estimates (fixing them as
\`\*name=value\`); structural columns are left intact.

## Usage

``` r
fill_params(fit, modelParams, digits = NULL, summary = c("median", "mean"))
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

- modelParams:

  The model parameter sheet (data frame) to fill.

- digits:

  Optional rounding for written values.

- summary:

  Which point estimate to write (\`"median"\` or \`"mean"\`).

## Value

The filled model parameter data frame.

## Examples

``` r
if (FALSE) { # \dontrun{
# fit + modelParams from a fitted scenario
filled <- fill_params(fit, sc$modelParams, digits = 4)
} # }
```
