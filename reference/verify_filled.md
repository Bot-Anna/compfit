# Verify a filled model sheet

Re-reads a filled workbook and checks it rebuilds the same model
structure as the original sheet (catches accidental structural edits).

## Usage

``` r
verify_filled(filled_path, modelParams)
```

## Arguments

- filled_path:

  Path to the filled \`.xlsx\`.

- modelParams:

  The original model parameter sheet (data frame).

## Value

\`TRUE\`, invisibly, if all checks pass; otherwise an error.

## Examples

``` r
if (FALSE) { # \dontrun{
verify_filled(file.path(tempdir(), "modelParams_filled.xlsx"), sc$modelParams)
} # }
```
