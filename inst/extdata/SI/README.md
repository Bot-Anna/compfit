# SI example scenario

A **2-compartment SI** epidemic (susceptible → infected, no recovery), driven by
a mass-action infection term. Data were simulated at known "true" parameters and
lightly noised, so the fit should recover them.

```
dX1/dt = -beta*X1*X2/N      (X1 = S)
dX2/dt =  beta*X1*X2/N      (X2 = I)
X1(0) = N0*(1-init_inf),  X2(0) = init_inf*N0      (parameter-dependent init)
```

Features demonstrated:

- a **quadratic (mass-action) infection term** (`Quadratic1 = *2*-beta`);
- **parameter-dependent initial states** (`*X1=N0_0*(1-init_inf_0)`, ...);
- **fitted vs fixed parameters** — `beta`, `init_inf` fitted with box priors
  `[lo,hi]`; `N0` fixed;
- the **poisson** likelihood family (integer counts on `X2`);
- a **missing** data cell (`x`) for an unobserved year.

Truth: `beta=0.9`, `init_inf=0.005`, `N0=2000`, over 2000–2008 (partition 4).

Fits by both maximum likelihood (`method="lbfgsb"`, works on the R backend with
no Julia) and Bayesian sampling (`method="bayes"`).
