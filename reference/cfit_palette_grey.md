# compfit grayscale palette

A grayscale counterpart to \[cfit_palette\], for print/monochrome
output. Same named slots (\`data\`, \`model\`, \`band\`, \`censor\`,
\`dummy\`); elements are kept apart by shade (and by shape/linetype in
the plots). Select it with \`plot_fit(..., palette = "grey")\` or
globally via \`options(compfit.palette = "grey")\`.

## Usage

``` r
cfit_palette_grey
```

## Format

A named list of hex colour strings.

## Examples

``` r
cfit_palette_grey$data
#> [1] "#000000"
```
