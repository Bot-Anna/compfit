test_that("test-unit-r-backend", {
  # The pure-R (deSolve) backend: a full MLE fit that needs NO Julia, so it runs
  # on CI / CRAN. Exercises load_scenario -> fit -> get_point -> solve_and_evaluate
  # entirely in R, and checks the bayes+r guard.
  skip_if_not_installed("deSolve")

  dir <- fixture_dir("minimal")
  chk("minimal fixture exists", dir.exists(dir))

  sc <- load_scenario(dir, combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

  fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                               solver = solver_control(backend = "r"),
                               checkpoint_file = tempfile(fileext = ".rds"))

  chk("fit succeeds without Julia", isTRUE(fit$success))
  chk("backend recorded as r",      identical(fit$solver$backend, "r"))

  p <- get_point(fit)
  chk("X1 recovered in [50,150]", p$initial_state[["X1"]] >= 50 && p$initial_state[["X1"]] <= 150)
  chk("X2 fixed at 0",            p$initial_state[["X2"]] == 0)
  chk("k in [0,1]",               p$parms[["k"]] >= 0 && p$parms[["k"]] <= 1)
  chk("m in [0,1]",               p$parms[["m"]] >= 0 && p$parms[["m"]] <= 1)

  ev <- solve_and_evaluate(fit, p$initial_state, p$parms)
  chk("solve_and_evaluate spans the full grid",
      nrow(ev$evaluation) == length(fit$model$date))
})
