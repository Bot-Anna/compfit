# Simulate a fully-specified compartmental model (no fitting)

Solve a compartmental ODE model whose states and parameters are ALL
fixed (\`\*name=value\`), with no observed data and no estimation. The
model is built, its fixed point recovered, solved over the time grid,
and any \`data_dummy\` (and optional \`dataCombined\`) formulas are
evaluated on the trajectory.

## Usage

``` r
simulate_model(
  modelParams,
  data_dummy = NULL,
  dataCombined = NULL,
  solver = solver_control()
)
```

## Arguments

- modelParams:

  Model parameter sheet (data frame). Every state and parameter must be
  fixed.

- data_dummy:

  Optional dummy-data data frame (a \`Formula\` column of expressions to
  evaluate on the trajectory; no likelihood, never fit).

- dataCombined:

  Optional observed data to also evaluate as overlay series (never fit).
  \`NULL\` (the default) simulates with dummy series only.

- solver:

  Solver settings, see \[solver_control()\]. Use
  \`solver_control(backend = "r")\` for a Julia-free solve.

## Value

An object of class \`"compartmentalSim"\`: a list with
\`initial_state\`, \`parms\`, \`sir_out\` (solved trajectory),
\`evaluation\` (formulas over the grid), \`time_grid\`, \`model\`,
\`solver\`, \`data_dummy\` and \`dataCombined\`.

## Details

Use this to run a scenario forward from known values – e.g. a
hand-specified SIR/SEIR model – without supplying \`dataCombined\`. If
any quantity is still fittable (a \`\[lo,hi\]\` box or a distributional
prior), the sheet is not fully specified and an error is raised naming
the offending quantities; fit those with \[fitCompartmentalModel()\]
instead.

## See also

\[plot_simulation()\], \[build_compartmental_model()\],
\[fitCompartmentalModel()\].

## Examples

``` r
sim_dir <- system.file("extdata", "SIR_sim", package = "compfit")
mp <- read_data_file(file.path(sim_dir, "modelParams.csv"))
dd <- read_data_file(file.path(sim_dir, "dataDummy.csv"))
sim <- simulate_model(mp, data_dummy = dd, solver = solver_control(backend = "r"))
head(sim$sir_out)
#>         X1       X2        X3       date   N1  N_1 total_pop gamma_t
#> 1 4975.000 25.00000  0.000000 2009-12-31 5000 5000      5000   0.400
#> 2 4966.755 30.46296  2.782406 2010-04-01 5000 5000      5000   0.405
#> 3 4956.735 37.05306  6.211768 2010-07-02 5000 5000      5000   0.410
#> 4 4944.587 44.98267 10.430049 2010-10-01 5000 5000      5000   0.415
#> 5 4929.895 54.49727 15.607698 2010-12-31 5000 5000      5000   0.420
#> 6 4912.174 65.87797 21.948435 2011-04-02 5000 5000      5000   0.425
sim$evaluation
#>          date       X1         X2          X3 X2/(X1+X2+X3)
#> 1  2009-12-31 4975.000   25.00000    0.000000   0.005000000
#> 2  2010-04-01 4966.755   30.46296    2.782406   0.006092592
#> 3  2010-07-02 4956.735   37.05306    6.211768   0.007410612
#> 4  2010-10-01 4944.587   44.98267   10.430049   0.008996533
#> 5  2010-12-31 4929.895   54.49727   15.607698   0.010899454
#> 6  2011-04-02 4912.174   65.87797   21.948435   0.013175594
#> 7  2011-07-02 4890.863   79.44299   29.694459   0.015888598
#> 8  2011-10-01 4865.321   95.54745   39.131974   0.019109491
#> 9  2011-12-31 4834.823  114.58066   50.596781   0.022916132
#> 10 2012-04-01 4798.561  136.95967   64.479609   0.027391934
#> 11 2012-07-01 4755.651  163.11793   81.230654   0.032623585
#> 12 2012-09-30 4705.150  193.48749  101.362627   0.038697497
#> 13 2012-12-31 4646.075  228.47343  125.451381   0.045694686
#> 14 2013-04-01 4577.448  268.41942  154.132989   0.053683885
#> 15 2013-07-01 4498.340  313.56415  188.095995   0.062712830
#> 16 2013-10-01 4407.943  363.98973  228.067499   0.072797946
#> 17 2013-12-31 4305.643  419.56517  274.791942   0.083913034
#> 18 2014-04-01 4191.108  479.89045  329.001881   0.095978090
#> 19 2014-07-02 4064.370  544.24931  391.380870   0.108849862
#> 20 2014-10-01 3925.899  611.58090  462.519762   0.122316180
#> 21 2014-12-31 3776.650  680.48069  542.869107   0.136096138
#> 22 2015-04-01 3618.069  749.23922  632.691827   0.149847843
#> 23 2015-07-02 3452.057  815.92207  732.021358   0.163184414
#> 24 2015-10-01 3280.882  878.48751  840.630881   0.175697502
#> 25 2015-12-31 3107.052  934.92984  958.018531   0.186985969
#> 26 2016-04-01 2933.158  983.43022 1083.411831   0.196686043
#> 27 2016-07-01 2761.714 1022.49359 1215.792018   0.204498718
#> 28 2016-09-30 2595.011 1051.05264 1353.936192   0.210210528
#> 29 2016-12-31 2435.002 1068.52568 1496.472808   0.213705136
#> 30 2017-04-01 2283.231 1074.82454 1641.944498   0.214964908
#> 31 2017-07-01 2140.811 1070.31697 1788.871836   0.214063394
#> 32 2017-09-30 2008.433 1055.75468 1935.812381   0.211150936
#> 33 2017-12-31 1886.408 1032.18103 2081.410733   0.206436205
#> 34 2018-04-01 1774.731 1000.83228 2224.437121   0.200166457
#> 35 2018-07-01 1673.142  963.04383 2363.813740   0.192608766
#> 36 2018-10-01 1581.202  920.16899 2498.629400   0.184033798
#> 37 2018-12-31 1498.341  873.51455 2628.143993   0.174702910
```
