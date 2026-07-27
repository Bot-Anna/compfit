# Plot local sensitivity

Visualises a \[local_sensitivity()\] result as a tornado plot (ranked
summary sensitivities) or a trajectory plot (sensitivity over time).

## Usage

``` r
plot_local_sensitivity(
  ls,
  kind = c("tornado", "trajectory"),
  top = NULL,
  base_size = 11,
  ncol = NULL
)
```

## Arguments

- ls:

  A \[local_sensitivity()\] result.

- kind:

  Plot kind: tornado or trajectory.

- top:

  Optional number of top parameters to show.

- base_size:

  Base font size.

- ncol:

  Number of facet columns (default: automatic).

## Value

A ggplot object (or a list of them).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_local_sensitivity(ls, kind = "tornado")   # ls from local_sensitivity()
} # }
```
