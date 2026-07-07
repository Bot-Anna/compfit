# ============================================================
# setup.R
# Run-once environment setup: libraries, Julia initialisation,
# and sourcing of all building-block subscripts.
#
# This file is about WHERE THE CODE LIVES, not about data.
# Scenario data locations are supplied at run time in main.R via
# load_scenario().
# ============================================================
# ---- Code location ----
# The folder containing this project's R code. R_DIR holds the
# building-block scripts (utils.R, the model builders, etc.).
if (!exists("PROJECT_ROOT")) {
  stop("PROJECT_ROOT must be set before sourcing setup.R (set it in main.R).")
}
R_DIR <- file.path(PROJECT_ROOT, "R")
setwd(PROJECT_ROOT)
# ---- Libraries ----
c_packages <- c('dplyr',
                'deSolve',
                'tidyverse',
                'reshape2',
                'stringr',
                'stats',
                'Matrix',
                'rlang',
                'ggplot2',
                'readxl',
                'patchwork',
                'zoo',
                'parallel',
                'optimParallel',
                'FME',
                'diffeqr',
                'JuliaCall',
                'qrng',
                'spacefillr',
                'DEoptim',
                'tictoc',
                'pracma',
                'beepr',
                'sensitivity')
lapply(c_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
})
lapply(c_packages, function(pkg) {
  library(pkg, character.only = TRUE)
})
# ---- Source building blocks ----
# Order matters: utils first (everyone depends on it), then the model
# builders, then scenario loader + master function.
source(file.path(R_DIR, "julia_setup.R"))  # setup_julia() + lazy .compfit_ensure_julia()
source(file.path(R_DIR, "utils.R"))
source(file.path(R_DIR, "scenario.R"))
source(file.path(R_DIR, "numberOfComps.R"))
source(file.path(R_DIR, "validate.R"))               # upfront modelParams entry checks
source(file.path(R_DIR, "statesAndParams.R"))
source(file.path(R_DIR, "generateExpressions.R"))
source(file.path(R_DIR, "compartmentalFunction.R"))  # ODE function: build + Julia register/solve
source(file.path(R_DIR, "lossFunction.R"))
source(file.path(R_DIR, "hypercubeSampling.R"))
source(file.path(R_DIR, "priorSpec.R"))    # Bayesian: parallel prior-spec builder
source(file.path(R_DIR, "bayesJulia.R"))   # Bayesian: Turing model + sampler bridge (julia backend)
source(file.path(R_DIR, "bayesR.R"))       # Bayesian: gradient-free pure-R sampler (r backend)
source(file.path(R_DIR, "bayesPost.R"))    # Bayesian: posterior post-processing
source(file.path(R_DIR, "plotPriorPosterior.R"))  # Bayesian: prior vs posterior plots
source(file.path(R_DIR, "accessors.R"))    # accessor layer + summary method
source(file.path(R_DIR, "evaluate.R"))     # solve_and_evaluate (shared by main.R + plots)
source(file.path(R_DIR, "plots.R"))        # consolidated plot_fit()
source(file.path(R_DIR, "identifiability.R"))  # practical identifiability (posterior)
source(file.path(R_DIR, "sensitivity.R"))  # sensitivity: global (Sobol) + local (derivative)
source(file.path(R_DIR, "fitCompartmentalModel.R"))
source(file.path(R_DIR, "io.R"))                # save_fit / load_fit / reload_fit
source(file.path(R_DIR, "extract_code.R"))      # model/loss/fit as portable source
source(file.path(R_DIR, "fillParams.R"))        # write fitted estimates back to the sheet
source(file.path(R_DIR, "counterfactual.R"))    # coupled-pair counterfactual pipeline

# ---- Julia initialisation ----
# Bridge + one-time Julia package install (guarded). Same logic the package
# exposes as setup_julia(); doing it here keeps the source(setup.R) workflow
# fully self-contained.
setup_julia()