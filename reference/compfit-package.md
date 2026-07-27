# compfit: Fit, Analyse and Compare Compartmental ODE Models (R + Julia)

Build compartmental (SIR-type) ODE models from a spreadsheet
specification, fit them by maximum likelihood (L-BFGS-B / DEoptim /
hypercube) or Bayesian sampling, and analyse the results: posterior
summaries, prior-vs-posterior and predictive plots, global (Sobol) and
local (derivative) sensitivity, practical identifiability,
counterfactual scenarios, and self-contained code export. ODE solves and
Bayesian sampling run in Julia via the JuliaCall bridge; call
\[setup_julia()\] once per session before fitting.

## Author

**Maintainer**: Anna Bot <anna.bot@outlook.com>
([ORCID](https://orcid.org/0000-0001-7803-1518))
