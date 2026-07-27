# Summarise a Sobol sensitivity analysis

Reports validity, per-stream top parameters, interaction structure, and
a global parameter ranking from a Sobol analysis. Accepts a precomputed
tidy table or a fit/model (computed via \[sobol_streams()\]).

## Usage

``` r
sobol_report(x, n = 500, cor_to = NULL, plots = TRUE, top_n = 5, ...)
```

## Arguments

- x:

  A tidy Sobol data frame (output, parameter, first_order, total) or a
  \`"compartmentalFit"\`/\`"compartmentalModel"\` object.

- n:

  Base sample size when computing from a fit/model.

- cor_to:

  Optional parameter name to correlate others against.

- plots:

  If \`TRUE\`, build a total-index heatmap.

- top_n:

  Number of top parameters to report per stream.

- ...:

  Passed to \[sobol_streams()\] when computing.

## Value

Invisibly, a list with the table, stream sums, parameter ranking,
not-converged flags, range/NA diagnostics, and (if \`plots\`) a heatmap.

## Examples

``` r
if (FALSE) { # \dontrun{
sobol_report(fit, n = 500)            # compute from a fit, or:
sob <- sobol_streams(fit, n = 500)
r   <- sobol_report(sob)              # report a precomputed table
} # }
```
