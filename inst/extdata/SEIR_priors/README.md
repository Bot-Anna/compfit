# SEIR_priors example scenario

An **SEIR** built to exercise the remaining **parameter priors** and the full set
of **likelihood families**. Three observed streams, three families; four
different prior kinds across the fitted parameters.

```
dX1/dt = -beta*X1*X3/N               (X1 = S)
dX2/dt =  beta*X1*X3/N - sigma_ei*X2  (X2 = E)
dX3/dt =  sigma_ei*X2  - gamma_t*X3   (X3 = I)
dX4/dt =  gamma_t*X3                  (X4 = R)
gamma_t = gamma*(1 + ramp*time)        (time-varying recovery)
X1(0) = N0*(1-init_inf),  X3(0) = init_inf*N0,  X2(0) = X4(0) = 0
```

(`sigma_ei` is the `E→I` incubation rate — not named `sigma`, which a Bayesian
fit reserves for the observation-noise scale.)

Prior distributions demonstrated (one per fitted parameter):

| Parameter | Prior | Kind |
|----|----|----|
| `beta`     | `[0,3]`               | Uniform (box) |
| `sigma_ei` | `Gamma(5,0.1)`        | Gamma |
| `gamma`    | `LogNormal(-1.2,0.4)` | LogNormal |
| `init_inf` | `Beta(2,1000)`        | Beta |

`ramp` and `N0` are fixed. **Likelihood families**: exposed (`X2`) **gaussian**,
infected (`X3`) **poisson**, recovered (`X4`) **lognormal**. Data cells: a
right-censored (`>=L`) year on `E`; plain counts plus a missing (`x`) year on the
poisson `I`; a **soft interval** (`[A,B]~s`) on the lognormal `R`.

Fits by both maximum likelihood and Bayesian sampling. Truth: `beta=1.5`,
`sigma_ei=0.5`, `gamma=0.3`, `ramp=0.04`, `N0=10000`, `init_inf=0.001`, over
2005–2014 (partition 4). Condition `beta>gamma`.

Together with `SIR_priors` this covers every prior distribution
(`Uniform`/`Normal`/`LogNormal`/`Beta`/`Gamma`/`StudentT`) and every likelihood
family (`gaussian`/`poisson`/`negbin`/`lognormal`).
