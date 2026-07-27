# Read a data file (csv or xlsx)

Reads a \`.csv\` or \`.xlsx\` data file, optionally forcing every column
to be read as text so numbers stored as text round-trip unchanged.

## Usage

``` r
read_data_file(path, text_cols = FALSE, ...)
```

## Arguments

- path:

  Path to a \`.csv\` or \`.xlsx\` file.

- text_cols:

  If \`TRUE\`, read all columns as character.

- ...:

  Passed to the underlying reader.

## Value

A data frame.

## Examples

``` r
f <- system.file("extdata", "minimal", "modelParams.csv", package = "compfit")
read_data_file(f, text_cols = TRUE)
#> # A tibble: 4 × 10
#>   `_Level1` Others     States Functions Parameters Conditions Linear1 Quadratic1
#>   <chr>     <chr>      <chr>  <chr>     <chr>      <chr>      <chr>   <chr>     
#> 1 1         startpoin… X1=[5… NA        k=[0,1]    NA         -k      0         
#> 2 2         endpoint=… *X2=0  NA        m=[0,1]    NA         0       0         
#> 3 NA        partition… NA     NA        NA         NA         NA      0         
#> 4 NA        cutoff=20… NA     NA        NA         NA         NA      0         
#> # ℹ 2 more variables: Linear2 <chr>, Quadratic2 <chr>
```
