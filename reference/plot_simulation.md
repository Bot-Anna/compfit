# Plot a simulated model

Plot the trajectories from a \[simulate_model()\] result: one panel per
compartment and/or per evaluated formula (dummy / overlay series), each
a line over the date grid. Mirrors \[plot_fit()\]'s return shape.

## Usage

``` r
plot_simulation(
  sim,
  which = c("both", "states", "series"),
  palette = NULL,
  ncol = 2,
  base_size = 11
)
```

## Arguments

- sim:

  A \`"compartmentalSim"\` object from \[simulate_model()\].

- which:

  Which panels to draw: \`"states"\` (compartment trajectories),
  \`"series"\` (evaluated dummy/overlay formulas), or \`"both"\`
  (default).

- palette:

  Colour palette: a scheme name (\`"okabe"\`/\`"grey"\`), a custom list
  like \[cfit_palette\], or \`NULL\` to use the \`compfit.palette\`
  option.

- ncol:

  Number of columns in the arranged grid.

- base_size:

  Base font size for the theme.

## Value

A list with \`plots\` (named list of ggplot objects) and \`grid\` (a
patchwork arrangement, or \`NULL\` if patchwork is unavailable).

## See also

\[simulate_model()\], \[plot_fit()\].

## Examples

``` r
sim_dir <- system.file("extdata", "SIR_sim", package = "compfit")
mp <- read_data_file(file.path(sim_dir, "modelParams.csv"))
dd <- read_data_file(file.path(sim_dir, "dataDummy.csv"))
sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
if (FALSE)  plot_simulation(sim)$grid  # \dontrun{}
```
