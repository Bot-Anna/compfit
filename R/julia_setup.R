# ============================================================
# julia_setup.R
# Julia bridge initialisation for compfit.
#
# Julia is a SYSTEM REQUIREMENT (not an R package): the model is solved and
# sampled in Julia via JuliaCall/diffeqr. We deliberately do NOT start Julia at
# package load (it is slow and can fail) -- call setup_julia() once per session
# before fitting, or rely on the lazy .compfit_ensure_julia() that the bridge
# functions call. The heavy, truly one-time-per-machine part (registry update +
# Pkg.add) is guarded by an option so it runs at most once per session.
# ============================================================

# Package-internal session state (Julia readiness, the diffeqr handle).
.compfit_state <- new.env(parent = emptyenv())

#' Initialise the Julia bridge
#'
#' Starts the Julia session (via diffeqr + JuliaCall) and, on first use,
#' provisions the required Julia packages (SciMLBase, OrdinaryDiffEq, Turing,
#' Distributions) from a \emph{pinned} environment shipped with the package
#' (\code{inst/julia/Project.toml}), copied to a writable per-user cache and
#' instantiated so every user gets the same tested versions. If that fails it
#' falls back to installing the packages into the default Julia environment.
#' Loads OrdinaryDiffEq. Safe to call repeatedly: the session bind and the
#' one-time provisioning are each guarded.
#'
#' @return TRUE, invisibly.
#' @examples
#' \dontrun{
#' setup_julia()  # starts Julia; installs the Julia packages on first call
#' }
#' @export
setup_julia <- function() {
  # The Julia bridge packages are Suggests (the Julia backend is optional). Fail
  # with an actionable message rather than a cryptic one if they are absent.
  need <- c("JuliaCall", "diffeqr")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("The Julia backend needs the package(s): ", paste(miss, collapse = ", "),
         ".\nInstall them (install.packages(c(", paste(sprintf('\"%s\"', miss), collapse = ", "),
         "))) and a Julia runtime, or fit without Julia via solver_control(backend = \"r\").",
         call. = FALSE)
  # Bridge: bind the diffeqr handle + JuliaCall session once per process.
  if (!isTRUE(.compfit_state$bridge_ready) || is.null(.compfit_state$de)) {
    .compfit_state$de <- diffeqr::diffeq_setup()
    JuliaCall::julia_setup()
    .compfit_state$bridge_ready <- TRUE
  }
  # Heavy package setup: once per R process. Prefer the PINNED environment
  # (inst/julia/Project.toml) so every user gets the same tested Julia stack;
  # fall back to installing into the default environment if that fails.
  if (!isTRUE(getOption("julia_pkgs_ready"))) {
    JuliaCall::julia_command("using Pkg")
    if (!.compfit_activate_pinned_env()) {
      # Fallback: previous behaviour -- add into the active/default environment.
      JuliaCall::julia_command('Pkg.Registry.update()')
      JuliaCall::julia_command('Pkg.add(PackageSpec(name="SciMLBase", version="2.3.1"))')
      JuliaCall::julia_command('Pkg.add("OrdinaryDiffEq")')
      # Bayesian stack (only needed for method = "bayes"); Turing pulls Distributions.
      JuliaCall::julia_command('Pkg.add("Turing")')
      JuliaCall::julia_command('Pkg.add("Distributions")')
    }
    options(julia_pkgs_ready = TRUE)
  }
  # Needed once per session so solve()/ODEProblem are at top level; cheap to repeat.
  JuliaCall::julia_command("using OrdinaryDiffEq")
  .compfit_state$ready <- TRUE
  invisible(TRUE)
}

# Copy the shipped pinned Project.toml (inst/julia) to a writable per-user cache,
# activate it and Pkg.instantiate() so the whole Julia stack matches the tested
# versions. Returns TRUE on success; on any error warns and returns FALSE so
# setup_julia() can fall back to installing into the default environment. A
# writable copy is used because the installed package location may be read-only
# (instantiate writes a Manifest).
.compfit_activate_pinned_env <- function() {
  proj_src <- system.file("julia", "Project.toml", package = "compfit")
  if (!nzchar(proj_src) || !file.exists(proj_src)) return(FALSE)
  tryCatch({
    cache <- file.path(tools::R_user_dir("compfit", "cache"), "julia")
    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
    file.copy(proj_src, file.path(cache, "Project.toml"), overwrite = TRUE)
    JuliaCall::julia_command(sprintf('Pkg.activate(raw"%s")', cache))
    JuliaCall::julia_command("Pkg.instantiate()")
    TRUE
  }, error = function(e) {
    warning("compfit: could not use the pinned Julia environment (",
            conditionMessage(e), "); falling back to the default environment.",
            call. = FALSE)
    FALSE
  })
}

# Lazy, memoised guard called by the bridge functions (solve/register/sample) so
# package users need not remember to call setup_julia() first. No-op once ready.
.compfit_ensure_julia <- function() {
  if (!isTRUE(.compfit_state$ready)) setup_julia()
  invisible(TRUE)
}
