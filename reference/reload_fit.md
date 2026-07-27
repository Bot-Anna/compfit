# Reload a fit (optionally re-registering Julia)

Entry point used by the counterfactual pipeline and interactive reloads.
With \`register = TRUE\` returns a fully re-executable fit (re-registers
the Julia ODE and rebuilds the loss via \[load_fit()\]); with \`register
= FALSE\` does a results-only read (object + draws, no Julia session).

## Usage

``` r
reload_fit(path, register = TRUE)
```

## Arguments

- path:

  Path to the \`.rds\` file (\`.rds\` appended if missing).

- register:

  If \`TRUE\`, re-register Julia and rebuild the loss.

## Value

A \`"compartmentalFit"\` object.

## Examples

``` r
if (FALSE) { # \dontrun{
fit  <- reload_fit(file.path(tempdir(), "fit.rds"))                  # with Julia
draws <- reload_fit(file.path(tempdir(), "fit.rds"), register = FALSE) # results only
} # }
```
