# Medium fixture scenario

A 3-compartment **SIR** model that builds on the minimal fixture by adding every
feature the minimal one lacks.

```
S = X1, I = X2, R = X3   (one country)
dX1/dt = -beta*X1*X2/N                       (frequency-dependent infection)
dX2/dt =  beta*X1*X2/N - gamma_t*X2
dX3/dt =  gamma_t*X2
gamma_t = gamma*(1 + ramp*time)              (time-varying recovery; a Function)
N = X1 + X2 + X3

Initial states (parameter-dependent, function-defined):
  X1(0) = N0*(1 - init_inf)
  X2(0) = init_inf*N0
  X3(0) = 0

Fitted: beta, gamma.   Fixed: N0=1000, init_inf=0.01, ramp=0.05.
Condition: beta > gamma  (adds a penalty when violated).
```

Features exercised that `minimal/` does not:
- **Quadratic / second-order term** — the infection `beta*S*I/N`, encoded in
  `Quadratic1` as `*2*-beta` (defined in the susceptible equation, routed to the
  infected with a negative coefficient).
- **Functions column** — `gamma_t <- gamma*(1+ramp*time)`, a time-varying rate
  spliced into the ODE.
- **Conditions column** — `beta > gamma`, turned into a fitting penalty.
- **Parameter-dependent function-defined initial states** — `*X1=N0_0*(1-init_inf_0)`
  etc. (parameters referenced by their `_0`-suffixed names).

First-order (linear) coefficients use one column per compartment — `Linear1`,
`Linear2`, `Linear3` (column j of the n×n coefficient matrix) — interleaved with
the `Quadratic<j>` columns: `Linear1, Quadratic1, Linear2, Quadratic2, Linear3,
Quadratic3`. Column order is cosmetic (the model reads columns by name).

`test-unit-modelchain-medium.R` builds it without Julia; `test-integration-medium.R`
fits it end to end and checks conservation (S+I+R constant).
