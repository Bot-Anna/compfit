# Next numbered plot path

Given a plot path, returns the next available numbered variant (e.g.
\`Plot_1.pdf\`, \`Plot_2.pdf\`, ...) so existing plots are never
overwritten.

## Usage

``` r
next_plot_path(path)
```

## Arguments

- path:

  A plot file path.

## Value

A character path with the next free numeric suffix.

## Examples

``` r
next_plot_path(file.path(tempdir(), "Plot.pdf"))
#> [1] "/tmp/RtmpF6wWtu/Plot_1.pdf"
```
