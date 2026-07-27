# compfit ggplot2 theme

A minimal ggplot2 theme used across compfit plots.

## Usage

``` r
cfit_theme(base_size = 11)
```

## Arguments

- base_size:

  Base font size.

## Value

A ggplot2 theme object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + cfit_theme()
```
