# SIR_priors example scenario

An **SIR with under-reporting**, built to exercise the **parameter-prior grammar**:
each fitted parameter uses a *different* prior distribution. Reported cases are a
fraction `rho` of the true infected (a coverage model), and the fully-observed
recovered stream anchors the scale so `rho` is identifiable.

```
dX1/dt = -beta*X1*X2/N              (X1 = S)
dX2/dt =  beta*X1*X2/N - gamma_t*X2 (X2 = I)
dX3/dt =  gamma_t*X2                (X3 = R)
gamma_t = gamma*(1 + ramp*time)      (time-varying recovery)
X1(0) = N0*(1-init_inf),  X2(0) = init_inf*N0,  X3(0) = 0
```

Prior distributions demonstrated (one per fitted parameter):

| Parameter | Prior | Kind |
|----|----|----|
| `beta`     | `StudentT(4,1.2,0.4)[0,3]` | Student-t (heavy-tailed), truncated |
| `gamma`    | `Normal(0.4,0.08)[0,1]`    | Normal, truncated |
| `init_inf` | `Beta(2,400)`              | Beta (a small proportion) |
| `rho`      | `[0.2,0.9]`                | Uniform (box) |

`ramp` and `N0` are fixed. **Likelihood families**: reported infections
(`rho*X2`) as overdispersed **negbin** counts; recovered (`X3`) as **gaussian**.
Data cells: the negbin stream carries plain counts plus a missing (`x`) year; the
gaussian stream carries a right-censored (`>=L`) year. (CDF-based cells stay on
the continuous stream — see the `SIR` example's note.)

This is primarily a **Bayesian** example: the informative priors are what
constrain `rho` vs the epidemic scale. It also fits by maximum likelihood, though
MLE (which ignores the priors) identifies `rho` only loosely.

Truth: `beta=1.2`, `gamma=0.4`, `ramp=0.05`, `N0=6000`, `init_inf=0.002`,
`rho=0.5`, over 2000–2009 (partition 4). Condition `beta>gamma`.

Together with `SEIR_priors` this covers every prior distribution
(`Uniform`/`Normal`/`LogNormal`/`Beta`/`Gamma`/`StudentT`) and every likelihood
family (`gaussian`/`poisson`/`negbin`/`lognormal`).
