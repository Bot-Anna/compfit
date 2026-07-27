# Solve the ODE in Julia

Low-level solver bridge: solves the registered Julia ODE at a given
initial state, parameter vector and save times. Mirrors \[solveWithR()\]
(the deSolve backend) and returns the same shape. Exported because the
scripts emitted by \[extract_code()\] call it.

## Usage

``` r
solveWithJulia(
  X,
  t,
  p,
  time_integer,
  solver = "BS3()",
  abstol = 1e-08,
  reltol = 1e-08
)
```

## Arguments

- X:

  Named numeric initial-state vector.

- t:

  Length-2 numeric integration span \`c(t0, t1)\`.

- p:

  Numeric parameter vector (in ODE order).

- time_integer:

  Numeric vector of save times.

- solver:

  Julia solver expression.

- abstol, reltol:

  Solver tolerances.

## Value

A list with \`matrix\` (states by time) and \`t\` (the save times).
