# Load a fit from disk

Reads a \`"compartmentalFit"\` saved by \[save_fit()\] and rebuilds the
loss closure, re-registering the Julia ODE function in the running
session.

## Usage

``` r
load_fit(path)
```

## Arguments

- path:

  Path to the \`.rds\` file (\`.rds\` appended if missing).

## Value

A \`"compartmentalFit"\` object.

## Examples

``` r
if (FALSE) { # \dontrun{
setup_julia()
fit <- load_fit(file.path(tempdir(), "fit.rds"))
} # }
```
