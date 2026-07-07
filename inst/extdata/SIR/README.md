# SIR example scenario

A **3-compartment SIR** epidemic (susceptible → infected → recovered): mass-action
infection plus recovery to `R` at a **time-varying** rate. Two observed streams
with two different likelihood families. Data were simulated at known "true"
parameters and noised.

```
dX1/dt = -beta*X1*X2/N              (X1 = S)
dX2/dt =  beta*X1*X2/N - gamma_t*X2 (X2 = I)
dX3/dt =  gamma_t*X2                (X3 = R)
gamma_t = gamma*(1 + ramp*time)      (Functions column, time-varying recovery)
X1(0) = N0*(1-init_inf),  X2(0) = init_inf*N0,  X3(0) = 0
```

Features demonstrated:

- **two data streams with two families**: `X2` (infected) as overdispersed
  **negbin** counts, `X3` (recovered) as **gaussian**;
- **interval-censored** (`[A,B]`) and **right-censored** (`>=L`) cells on the
  gaussian `R` stream;
- a time-varying rate **`Function`** and a **`Condition`** (`beta>gamma`).

Note the CDF-based cells (`[A,B]`, `>=L`) sit on the **continuous** (`gaussian`)
stream, not the discrete `negbin` one: under the Julia/NUTS path `logcdf` is not
autodifferentiable for discrete families, so discrete-family streams carry plain
observed / missing cells only. (The R backend has no such restriction.)

Truth: `beta=1.2`, `gamma=0.35`, `ramp=0.05`, `init_inf=0.002`, `N0=5000`, over
2010–2018 (partition 4). `beta`, `gamma`, `init_inf` fitted (box priors); `ramp`,
`N0` fixed.

Fits by both maximum likelihood and Bayesian sampling.
