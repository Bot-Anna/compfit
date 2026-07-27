# Save a fit to disk

Serialises a \`"compartmentalFit"\` to an \`.rds\` file (the
non-portable loss closure is dropped; \[load_fit()\] rebuilds it).

## Usage

``` r
save_fit(fit, path)
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

- path:

  Output path (\`.rds\` appended if missing).

## Value

The output path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
save_fit(fit, file.path(tempdir(), "fit.rds"))   # fit from fitCompartmentalModel()
} # }
```
