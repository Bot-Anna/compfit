# Write a filled model sheet to a workbook

Fills the sheet via \[fill_params()\] and writes it to an \`.xlsx\` in
the scenario directory.

## Usage

``` r
write_filled_params(
  fit,
  scenario_dir,
  modelParams = fit$model$modelParams,
  out_file = "modelParams_filled.xlsx",
  digits = NULL,
  summary = c("median", "mean")
)
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

- scenario_dir:

  Directory to write the workbook into.

- modelParams:

  The model parameter sheet (defaults to the fit's stored sheet).

- out_file:

  Output workbook name.

- digits:

  Optional rounding for written values.

- summary:

  Which point estimate to write (\`"median"\` or \`"mean"\`).

## Value

The output path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
write_filled_params(fit, tempdir(), modelParams = sc$modelParams)
} # }
```
