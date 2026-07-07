# SIR example scenario — **named compartments**

The **same 3-compartment SIR** epidemic as the `SIR` scenario, but written with
**named compartments** (`S`, `I`, `R`) instead of the positional `X1`, `X2`,
`X3`. It demonstrates that the machinery is decoupled from the `X<i>` naming
convention.

```
dS/dt = -beta*S*I/N               (was X1)
dI/dt =  beta*S*I/N - gamma_t*I    (was X2)
dR/dt =  gamma_t*I                 (was X3)
gamma_t = gamma*(1 + ramp*time)
S(0) = N0*(1-init_inf),  I(0) = init_inf*N0,  R(0) = 0
```

What the naming feature looks like on the sheet:

- **`States`** lists the compartments as `*S=...`, `*I=...`, `*R=...`. The
  **order of this column is canonical**: compartment 1 is `S`, 2 is `I`, 3 is
  `R`. Everything downstream (`X[i]` in the ODE, solution columns, plots) uses
  that order.
- **`_Level1`** references each compartment by its **name** (`S`, `I`, `R`). It
  would accept the numbers `1`, `2`, `3` equally — names and 1-based indices are
  interchangeable everywhere a compartment is referenced.
- The rate columns are **`Linear<name>` / `Quadratic<name>`** (`LinearS`,
  `LinearI`, ..., `QuadraticS`, ...) instead of `Linear1` / `Quadratic1`. The
  numeric forms still work; you may even mix them.
- **`dataCombined` formulas** reference compartments by name (`I`, `R`) — e.g.
  the infected stream is just `I`.

Note that the second-order `*goto*coeff` **target** (the `2` in `*2*-beta`)
stays a numeric compartment index — numbers are always valid compartment
references.

Everything else — priors, families, censoring, the time-varying `Function`, the
`Condition` — is identical to the `SIR` scenario. Fitting this sheet gives the
same result as its `X1..Xn` twin; the names are purely for readability.

Truth: `beta=1.2`, `gamma=0.35`, `ramp=0.05`, `init_inf=0.002`, `N0=5000`, over
2010–2018 (partition 4). `beta`, `gamma`, `init_inf` fitted; `ramp`, `N0` fixed.
