# SEIR example scenario

A **4-compartment SEIR** epidemic (susceptible → exposed → infected → recovered):
mass-action infection into a latent `E` class, progression `E→I` at rate `sigma_ei`,
and recovery `I→R` at a **time-varying** rate. Three observed streams, three
likelihood families. Data were simulated at known "true" parameters and noised.

```
dX1/dt = -beta*X1*X3/N               (X1 = S)
dX2/dt =  beta*X1*X3/N - sigma_ei*X2     (X2 = E)
dX3/dt =  sigma_ei*X2     - gamma_t*X3   (X3 = I)
dX4/dt =  gamma_t*X3                  (X4 = R)
gamma_t = gamma*(1 + ramp*time)        (Functions column, time-varying recovery)
X1(0) = N0*(1-init_inf),  X3(0) = init_inf*N0,  X2(0) = X4(0) = 0
```

Note the infection is `S + I -> E` (the quadratic term in `Quadratic1` routes the
new infection into `E`, compartment 2), and `I` is `X3`.

Features demonstrated:

- a **4-compartment chain** with a **latent class** (the `E→I` incubation rate,
  named `sigma_ei` rather than `sigma` because a Bayesian fit reserves `sigma`
  for the Gaussian/log-normal observation-noise scale — see the note below);
- **three data streams with three families**: `E` **gaussian**, `I`
  **negbin**, `R` **lognormal**;
- **right-censored** (`>=L`) and **missing** (`x`) cells on `E`, an
  **interval** (`[A,B]`) cell on `R` — all on continuous families (see the note
  in the SIR README on why CDF-based cells avoid the discrete `negbin` stream);
- a time-varying rate **`Function`**, a **`Condition`** (`beta>gamma`), and
  **parameter-dependent initial states**.

Observing `E` directly pins the incubation rate `sigma_ei` (SEIR is weakly
identifiable from `I`/`R` alone).

Truth: `beta=1.5`, `sigma_ei=0.5`, `gamma=0.3`, `ramp=0.04`, `init_inf=0.001`,
`N0=10000`, over 2005–2013 (partition 4). `beta`, `sigma_ei`, `gamma`, `init_inf`
fitted (box priors); `ramp`, `N0` fixed.

Fits by both maximum likelihood and Bayesian sampling.

**Reserved names.** A Bayesian fit uses `sigma` (Gaussian/log-normal noise scale)
and `phi`/`phi<i>` (negbin dispersion) as its own parameters, so a model
parameter or state may not share those names — hence `sigma_ei` here. Using a
reserved name raises a clear error on the `method = "bayes"` path (MLE is
unaffected).
