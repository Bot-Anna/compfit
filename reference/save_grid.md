# Save a plot grid

Writes the assembled \`\$grid\` of a plot result to disk. No-op (with a
message) when \`\$grid\` is \`NULL\` (e.g. patchwork unavailable).

## Usage

``` r
save_grid(res, path, limitsize = FALSE, ...)
```

## Arguments

- res:

  A plot result list containing a \`\$grid\` element.

- path:

  Output file path.

- limitsize:

  Passed to \[ggplot2::ggsave()\].

- ...:

  Passed to \[ggplot2::ggsave()\].

## Value

The output path, invisibly (or \`NULL\` if there is no grid).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- plot_fit(fit)            # fit from fitCompartmentalModel()
save_grid(res, file.path(tempdir(), "Plot.pdf"))
} # }
```
