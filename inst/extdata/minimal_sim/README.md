# minimal_sim — simulate-only (no fitting, no data)

A fully-specified two-compartment linear model (`dX1=-k*X1`, `dX2=k*X1-m*X2`)
with **every** state and parameter fixed. There is **no `dataCombined.csv`**:
this scenario exists to demonstrate `simulate_model()`, which solves the model
forward and evaluates the `dataDummy` formulas on the trajectory — no data, no
estimation.

```r
mp <- read_data_file(system.file("extdata/minimal_sim/modelParams.csv", package = "compfit"))
dd <- read_data_file(system.file("extdata/minimal_sim/dataDummy.csv",   package = "compfit"))
sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
plot_simulation(sim)$grid
```
