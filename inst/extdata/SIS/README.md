# SIS example scenario

A **2-compartment SIS** epidemic (susceptible ↔ infected): mass-action infection
plus recovery back to susceptible at a **time-varying** rate. Data were simulated
at known "true" parameters and noised.

```
dX1/dt = -beta*X1*X2/N + gamma_t*X2     (X1 = S)
dX2/dt =  beta*X1*X2/N - gamma_t*X2     (X2 = I)
gamma_t = gamma*(1 + ramp*time)          (Functions column, time-varying rate)
X1(0) = N0*(1-init_inf),  X2(0) = init_inf*N0
```

Features demonstrated:

- a **time-varying rate `Function`** used in a linear position (`gamma_t`);
- a **`Condition`** turned into a fitting penalty (`beta>gamma`);
- the **gaussian** likelihood family;
- a **left-censored** cell (`<=L`, e.g. a year below a detection limit);
- an **asymmetric** cell — the `2003` reading is `89->109`.

### The asymmetric cell, worked

`89->109` means: the reported value is **89**, and `89` is a **hard lower bound**
(the truth is at least 89), but you believe the true value is *higher* — the `->`
target `109` sets the soft direction (up) and the deviation scale
`dev = |109 - 89| = 20`. The likelihood is a floor at 89 plus a one-sided
soft penalty: each `20` units above 89 costs one log-unit. So the plausible mass
sits in `[89, ~109]`, hugging 89 and trailing off past 109 (values beyond 109 are
allowed, just increasingly costly). The mode stays near the recorded value; `109`
is the *width* of the upward uncertainty, not a second most-likely point.

Equivalently you could declare a stream-wide scale in the `Likelihood` cell
(`gaussian; asym=20`) and write the cell as `89+`. Use `A-` / `A->B` with `B<A`
for a downward belief. Note this is *not* the tool for "reported 89 but I think
the truth is a band that excludes 89" — for a plateau over a band use the
interval cell `[A,B]` (see the `SIR`/`SEIR` examples).

Truth: `beta=0.6`, `gamma=0.2`, `ramp=0.03`, `init_inf=0.02`, `N0=1000`, over
2000–2010 (partition 4). `beta`, `gamma`, `init_inf` are fitted (box priors);
`ramp`, `N0` fixed. `Average` is set to `1/mean(stream)` so the residuals are
scale-normalised (see the note in the top-level example index).

Fits by both maximum likelihood and Bayesian sampling.
