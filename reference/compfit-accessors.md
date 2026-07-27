# Accessors for a compartmental fit

Extract components of a \`"compartmentalFit"\` object.

## Usage

``` r
get_loss(fit)

get_model(fit)

get_data(fit)

get_time_grid(fit)

get_bounds(fit)

get_solver(fit)

get_compartmental_function(fit)

get_julia_code(fit)

get_expressions(fit)

get_states_and_params(fit)

get_point(fit)

get_samples(fit)
```

## Arguments

- fit:

  A \`"compartmentalFit"\` object.

## Value

The requested component of the fit.

## Examples

``` r
if (FALSE) { # \dontrun{
# fit from fitCompartmentalModel()
get_point(fit)
get_julia_code(fit)
} # }
```
