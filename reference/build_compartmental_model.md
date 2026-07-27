# Build a compartmental model (no fitting)

Runs everything \[fitCompartmentalModel()\] does except the fit
dispatch: builds the model, time grid, prepared data, bounds and loss
closure. Useful for tools that need the model and bounds without a prior
fit (e.g. Sobol sensitivity).

## Usage

``` r
build_compartmental_model(
  modelParams,
  dataCombined = NULL,
  solver = solver_control(),
  init = NULL,
  checkpoint_file = tempfile(fileext = ".rds")
)
```

## Arguments

- modelParams:

  Model parameter sheet (data frame).

- dataCombined:

  Combined observation data (data frame). Optional: pass \`NULL\` (the
  default) to build in simulation mode with no observed streams, in
  which case the loss closure is \`NULL\` (nothing to fit against). See
  \[simulate_model()\].

- solver:

  Solver settings, see \[solver_control()\].

- init:

  Optional initial point (named numeric) for the loss closure.

- checkpoint_file:

  Path for the best-solution checkpoint \`.rds\`.

## Value

A list of class \`"compartmentalModel"\` (model, sap, time_grid, data,
bounds, loss, best_state, best_start, solver). \`loss\` is \`NULL\` when
\`dataCombined\` is absent.

## Examples

``` r
mini <- system.file("extdata", "minimal", package = "compfit")
sc   <- load_scenario(mini, combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
# Build the machinery (no fitting). The R backend needs no Julia:
m <- build_compartmental_model(sc$modelParams, sc$dataCombined,
                               solver = solver_control(backend = "r"))
m$bounds
#> $lower
#> X1  k  m 
#> 50  0  0 
#> 
#> $upper
#>  X1   k   m 
#> 150   1   1 
#> 
#> $init_norm
#>  X1   k   m 
#> 0.5 0.5 0.5 
#> 
#> $random_init
#>    X1     k     m 
#> FALSE FALSE FALSE 
#> 
```
