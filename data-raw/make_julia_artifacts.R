# =============================================================================
# data-raw/make_julia_artifacts.R
#
# Regenerates the Julia/NUTS artifacts used by the SIR tutorial article
# (vignettes/articles/sir-tutorial.Rmd):
#
#   vignettes/articles/julia_fit.rds        posterior: report + prior-vs-posterior
#   vignettes/articles/julia_spaghetti.svg  pre-rendered trajectory spaghetti
#
# WHY THIS EXISTS
#   The tutorial does NOT run Julia at knit time, so CI / pkgdown stay Julia-free.
#   These artifacts are produced ONCE here and committed; the article loads them
#   with reload_fit(register = FALSE) and include_graphics().
#
# WHEN TO RE-RUN
#   Whenever the SIR scenario, the model sheet, or the Bayesian settings change
#   and you want the tutorial's Julia results to match. Re-run, then commit the
#   two regenerated files.
#
# REQUIREMENTS
#   Julia (>= 1.6). The first setup_julia() installs the Julia packages
#   (OrdinaryDiffEq, Turing, Distributions) and takes a few minutes.
#
# RUN FROM THE PROJECT ROOT
#   Rscript data-raw/make_julia_artifacts.R      # or source() it in an R session
# =============================================================================

## ---- config ---------------------------------------------------------------
art_dir    <- "vignettes/articles"
fit_path   <- file.path(art_dir, "julia_fit.rds")
bands_path <- file.path(art_dir, "julia_bands.svg")       # predictive + mean CrI bands
spag_path  <- file.path(art_dir, "julia_spaghetti.svg")   # per-draw trajectories
seed       <- 1L          # fixed -> the committed artifact is reproducible

## ---- 1. package + Julia ---------------------------------------------------
if (!requireNamespace("devtools", quietly = TRUE)) stop("install.packages('devtools')")
devtools::load_all(".")            # current source (or: library(compfit))
# Run the chains in parallel: give Julia >= chains threads BEFORE it starts
# (the thread count is fixed at Julia startup, i.e. the setup_julia() call below).
Sys.setenv(JULIA_NUM_THREADS = 4)  # matches bayes_control(chains = 4) in step 3
setup_julia()                      # first call installs Julia pkgs (slow)

## ---- 2. the SAME scenario the tutorial uses -------------------------------
sc <- load_scenario(
  system.file("extdata", "SIR", package = "compfit"),
  combined_file = "dataCombined.csv",
  dummy_file    = "dataDummy.csv",
  params_file   = "modelParams.csv"
)

## ---- 3. fit with the Julia/NUTS backend (the default) ---------------------
fit_julia <- fitCompartmentalModel(
  sc$modelParams, sc$dataCombined,
  method = "bayes",
  bayes  = bayes_control(chains = 4, iter = 2000, seed = seed)
)

## ---- 4. sanity check: diagnostics should be HEALTHY -----------------------
message("\n--- posterior_report (expect R-hat ~1.0, ESS in the hundreds+) ---")
print(posterior_report(fit_julia))

## ---- 5. save the fit (the article loads this, results-only) ---------------
save_fit(fit_julia, fit_path)
message("wrote ", fit_path)

## ---- 6. pre-render the trajectory plots (need the live Julia session) ------
# plot_fit() re-solves the ODE, so these must run here where Julia is registered;
# the article cannot produce them from a register = FALSE reload.
bands <- plot_fit(fit_julia, data_dummy = sc$dataDummy, band_type = "both")
ggplot2::ggsave(bands_path, bands$grid, width = 7, height = 4.5)
message("wrote ", bands_path)

sp <- plot_fit(fit_julia, data_dummy = sc$dataDummy,
               band_type = "spaghetti", spaghetti_draws = 80)
ggplot2::ggsave(spag_path, sp$grid, width = 7, height = 4.5)
message("wrote ", spag_path)

## ---- 7. verify the results-only reload the article relies on --------------
chk <- reload_fit(fit_path, register = FALSE)   # must work WITHOUT Julia
stopifnot(!is.null(chk$samples))

message("\nOK: artifacts regenerated. Commit them with:\n  git add ",
        fit_path, " ", bands_path, " ", spag_path)
