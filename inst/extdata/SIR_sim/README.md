# SIR_sim — simulate-only (no fitting, no data)

A fully-specified SIR (susceptible-infected-recovered) model with **every** state and parameter fixed
(`*name=value`). There is **no `dataCombined.csv`**: this scenario
demonstrates `simulate_model()`, which solves the model forward from the known
values and evaluates the `dataDummy` formulas on the trajectory — no observed
data, no estimation. It is the fixed twin of the `SIR` fitting scenario.

```r
mp <- read_data_file(system.file("extdata/SIR_sim/modelParams.csv", package = "compfit"))
dd <- read_data_file(system.file("extdata/SIR_sim/dataDummy.csv",   package = "compfit"))
sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
sim                       # a compact summary
plot_simulation(sim)$grid
```
