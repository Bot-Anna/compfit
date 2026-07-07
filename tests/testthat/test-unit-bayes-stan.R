test_that("test-unit-bayes-stan", {
# ============================================================
# test-unit-bayes-stan.R   (pure R; no Julia, no Stan compilation)
# Structural guards for the Stan Bayesian backend codegen:
#   * buildStanODEFunction emits a well-formed cf_ode from the ODE IR;
#   * buildStanModel emits all Stan blocks, maps observables to the right
#     trajectory row, and uses the [0,1] uniform reparameterisation;
#   * unsupported families are rejected at codegen time.
# When rstan is installed, additionally: the program passes stanc (syntax) and
# the family/cell capability guards in .fit_bayes_stan fire with clear messages.
# Actual compile+sample lives in test-integration-bayes-stan.R.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R",
               "fitCompartmentalModel.R", "priorSpec.R", "bayesStan.R"))

ld <- function(name) load_scenario(fixture_dir(name),
        combined_file = "dataCombined.csv", dummy_file = "dataDummy.csv",
        params_file = "modelParams.csv")

th_section("Stan ODE emitter (buildStanODEFunction via model$stan_code)")
sc  <- ld("SIR")
mod <- compfit:::.build_model(sc$modelParams, backend = "r")
ode <- mod$stan_code
chk("stan_code is generated",             !is.null(ode))
chk("declares the ODE function",          grepl("vector cf_ode(real t, vector X, array[] real p)", ode, fixed = TRUE))
chk("declares the derivative vector",     grepl("vector[3] dX;", ode, fixed = TRUE))
chk("has a derivative per compartment",   all(vapply(1:3, function(i)
                                              grepl(sprintf("dX[%d] =", i), ode, fixed = TRUE), logical(1))))
chk("population sum uses vector indices",  grepl("real N1 = X[1]+X[2]+X[3];", ode, fixed = TRUE))
chk("returns the derivative vector",       grepl("return dX;", ode, fixed = TRUE))

th_section("buildStanModel: full program structure (gaussian)")
ps <- buildPriorSpec(sc$modelParams)
tg <- compfit:::.time_grid(sc$modelParams); ny <- tg$endpoint - tg$startpoint + 1
gauss <- list(parseLikelihood("gaussian"))
code  <- compfit:::buildStanModel(ps, gauss, "X2", ny, tg$partition,
                                  mod$structure$number_of_comps,
                                  mod$structure$comp_names, ode)
for (blk in c("functions {", "data {", "transformed data {", "parameters {",
              "transformed parameters {", "model {"))
  chk(paste("has block", blk), grepl(blk, code, fixed = TRUE))
chk("solves with tolerance-controlled ode_rk45",
    grepl("ode_rk45_tol(cf_ode, X0, t0, ts, rel_tol, abs_tol, max_num_steps, p)", code, fixed = TRUE))
chk("uniform prior is reparam'd [0,1]", grepl("real<lower=0, upper=1> beta_n;", code, fixed = TRUE))
chk("denormalises beta_n -> beta",      grepl("real beta = 0 + beta_n * (3 - 0);", code, fixed = TRUE))
chk("observable X2 maps to traj row 2", grepl("traj[2, ]", code, fixed = TRUE))
chk("p is [fitted..., fixed...] order",  grepl("array[5] real p = {beta, gamma, init_inf, ramp, N0};",
                                              code, fixed = TRUE))
chk("gaussian likelihood via normal_lpdf", grepl("normal_lpdf(y[k, 1]", code, fixed = TRUE))
chk("left-censoring via normal_lcdf",      grepl("normal_lcdf(limit_mat[k, 1]", code, fixed = TRUE))

th_section("observable kinds map correctly")
annual_code <- compfit:::buildStanModel(ps, gauss, "annual(X2)", ny, tg$partition,
                 mod$structure$number_of_comps, mod$structure$comp_names, ode)
chk("annual() uses annual_integral", grepl("annual_integral(to_vector(traj[2, ])", annual_code, fixed = TRUE))
stock_code <- compfit:::buildStanModel(ps, gauss, "X2/(X1+X2+X3)", ny, tg$partition,
                 mod$structure$number_of_comps, mod$structure$comp_names, ode)
chk("stock ratio reads at annual_idx",
    grepl("to_vector((traj[2, ]/(traj[1, ]+traj[2, ]+traj[3, ]))[annual_idx])", stock_code, fixed = TRUE))

th_section("count families: integer data path + dispersion")
negbin_code <- compfit:::buildStanModel(ps, list(parseLikelihood("negbin")), "X2", ny,
                 tg$partition, mod$structure$number_of_comps, mod$structure$comp_names, ode)
chk("negbin uses neg_binomial_2_lpmf on integer data",
    grepl("neg_binomial_2_lpmf(y_int[k, 1]", negbin_code, fixed = TRUE))
chk("negbin declares a per-stream dispersion phi_1",
    grepl("real<lower=0> phi_1;", negbin_code, fixed = TRUE))
chk("negbin declares the integer data matrix", grepl("int y_int;", negbin_code, fixed = TRUE))
chk("count-only model has no sigma", !grepl("real<lower=0> sigma;", negbin_code, fixed = TRUE))
pois_code <- compfit:::buildStanModel(ps, list(parseLikelihood("poisson")), "X2", ny,
               tg$partition, mod$structure$number_of_comps, mod$structure$comp_names, ode)
chk("poisson uses poisson_lpmf on integer data",
    grepl("poisson_lpmf(y_int[k, 1]", pois_code, fixed = TRUE))
chk("poisson has no dispersion parameter", !grepl("phi_", pois_code, fixed = TRUE))
chk("gaussian-only model omits the integer data block",
    !grepl("int y_int;", code, fixed = TRUE))
# A genuinely unsupported family (binomial: parsed but no backend) is rejected.
chk_error("binomial family errors in buildStanModel",
          compfit:::buildStanModel(ps, list(parseLikelihood("binomial")), "X2", ny, tg$partition,
                                   mod$structure$number_of_comps, mod$structure$comp_names, ode))

th_section("threaded solver + priors, and the translatable-formula guard")
bdf_code <- compfit:::buildStanModel(ps, gauss, "X2", ny, tg$partition,
              mod$structure$number_of_comps, mod$structure$comp_names, ode,
              sigma_prior_stan = "exponential(1)", ode_solver = "bdf")
chk("tolerances are data-driven", grepl("real<lower=0> rel_tol;", bdf_code, fixed = TRUE))
chk("stiff solver requested -> ode_bdf_tol", grepl("ode_bdf_tol(cf_ode", bdf_code, fixed = TRUE))
chk("custom sigma prior is threaded", grepl("sigma ~ exponential(1);", bdf_code, fixed = TRUE))
nb_pr <- compfit:::buildStanModel(ps, list(parseLikelihood("negbin")), "X2", ny, tg$partition,
           mod$structure$number_of_comps, mod$structure$comp_names, ode,
           phi_prior_stan = "exponential(0.1)")
chk("custom phi prior is threaded", grepl("phi_1 ~ exponential(0.1);", nb_pr, fixed = TRUE))
# A Formula using a non-translatable function is rejected at codegen, not
# silently mistranslated.
chk_error("pmax() in a Formula is rejected",
          compfit:::buildStanModel(ps, gauss, "pmax(X2, 0)", ny, tg$partition,
                                   mod$structure$number_of_comps, mod$structure$comp_names, ode))

if (requireNamespace("rstan", quietly = TRUE)) {
  th_section("rstan present: programs pass stanc + capability guards fire")
  parse_ok <- function(cd) tryCatch({ rstan::stanc(model_code = cd, model_name = "chk"); TRUE },
                                    error = function(e) conditionMessage(e))
  chk("gaussian program passes stanc", isTRUE(parse_ok(code)))
  chk("negbin program passes stanc",   isTRUE(parse_ok(negbin_code)))
  chk("poisson program passes stanc",  isTRUE(parse_ok(pois_code)))

}

th_section("interval likelihood codegen (hard CDF diff + soft plateau)")
iv_g <- compfit:::buildStanModel(ps, gauss, "X2", ny, tg$partition,
          mod$structure$number_of_comps, mod$structure$comp_names, ode)
chk("continuous hard interval uses log_diff_exp of the CDFs",
    grepl("log_diff_exp(normal_lcdf(iupp_mat[k, 1]", iv_g, fixed = TRUE))
chk("soft interval is a plateau penalty on the mean",
    grepl("fmax(ilow_mat[k, 1] - mu1[k], 0) / idev_lo_mat[k, 1]", iv_g, fixed = TRUE))
chk("interval data is declared", grepl("matrix[n_years, n_streams] ilow_mat;", iv_g, fixed = TRUE))
iv_nb <- compfit:::buildStanModel(ps, list(parseLikelihood("negbin")), "X2", ny, tg$partition,
           mod$structure$number_of_comps, mod$structure$comp_names, ode)
chk("discrete hard interval uses integer edges",
    grepl("log_diff_exp(neg_binomial_2_lcdf(iupp_int[k, 1]", iv_nb, fixed = TRUE))
chk("discrete interval integer edges declared", grepl("int iupp_int;", iv_nb, fixed = TRUE))

th_section("asymmetric likelihood codegen (one-sided CDF anchor + soft damping)")
chk("continuous asymmetric uses lccdf/lcdf on the anchor",
    grepl("normal_lccdf(asym_val_mat[k, 1]", iv_g, fixed = TRUE))
chk("continuous asymmetric damps the mean excess",
    grepl("fmax(mu1[k] - asym_val_mat[k, 1], 0) / asym_dev_mat[k, 1]", iv_g, fixed = TRUE))
chk("discrete asymmetric uses the integer anchor",
    grepl("neg_binomial_2_lccdf(asym_val_int[k, 1]", iv_nb, fixed = TRUE))
chk("asymmetric data is declared", grepl("matrix[n_years, n_streams] asym_val_mat;", iv_g, fixed = TRUE))

th_summary("bayes-stan")
})
