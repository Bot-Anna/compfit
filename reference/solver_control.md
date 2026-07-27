# ODE solver settings

Bundle the Julia ODE solver and tolerances used for fitting and
plotting.

## Usage

``` r
solver_control(
  solver = "AutoTsit5(Rosenbrock23())",
  abstol = 1e-08,
  reltol = 1e-08,
  backend = c("julia", "r", "stan")
)
```

## Arguments

- solver:

  ODE solver. Its meaning is backend-specific, and the default
  \`"AutoTsit5(Rosenbrock23())"\` maps to each backend's auto-switching
  / non-stiff default, which suits almost every compartmental model.

  - \*\*Julia\*\*: the \`OrdinaryDiffEq\` constructor expression, passed
    verbatim. Common choices: \`"Tsit5()"\` (non-stiff),
    \`"Vern7()"\`/\`"Vern9()"\` (high-accuracy),
    \`"Rosenbrock23()"\`/\`"Rodas5()"\`/\`"KenCarp4()"\`/\`"TRBDF2()"\`
    (stiff), \`"AutoVern7(Rodas5())"\` (auto-switching). Any
    \`OrdinaryDiffEq\` solver works; see the [DifferentialEquations.jl
    docs](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).

  - \*\*R\*\* (deSolve): a \`deSolve::ode\` \`method\` name, e.g.
    \`"lsoda"\` (default, auto-switching), \`"radau"\`/\`"bdf"\`
    (stiff), \`"rk4"\`/\`"ode45"\` (non-stiff). Unrecognised strings
    fall back to \`"lsoda"\`.

  - \*\*Stan\*\*: \`"rk45"\` (default, non-stiff), \`"bdf"\` (stiff), or
    \`"adams"\` (non-stiff multistep). Set via \`solver = "bdf"\` etc.

- abstol:

  Absolute tolerance (all backends).

- reltol:

  Relative tolerance (all backends).

- backend:

  Solver backend: \`"julia"\` (default; fast, needs a Julia session via
  \[setup_julia()\]), \`"r"\` (pure R via deSolve; no Julia required,
  but slower), or \`"stan"\` (Bayesian sampling via Stan/rstan NUTS – no
  Julia; all R-side solves, e.g. plotting and the optional MAP init, use
  deSolve).

## Value

A named list of solver settings.

## Examples

``` r
solver_control(abstol = 1e-6, reltol = 1e-6)
#> $solver
#> [1] "AutoTsit5(Rosenbrock23())"
#> 
#> $abstol
#> [1] 1e-06
#> 
#> $reltol
#> [1] 1e-06
#> 
#> $backend
#> [1] "julia"
#> 
solver_control(backend = "r")   # Julia-free
#> $solver
#> [1] "AutoTsit5(Rosenbrock23())"
#> 
#> $abstol
#> [1] 1e-08
#> 
#> $reltol
#> [1] 1e-08
#> 
#> $backend
#> [1] "r"
#> 
```
