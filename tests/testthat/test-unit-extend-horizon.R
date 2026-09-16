test_that("test-unit-extend-horizon", {
  # Projecting past the fitted horizon: solve_and_evaluate(endpoint=) /
  # plot_fit(endpoint=) rebuild the solve grid for a later end year. Runs on the
  # pure-R backend, so no Julia is needed.
  #
  # The contract has two halves:
  #   (1) the extra years are ADDED -- the overlapping years are bit-identical,
  #       so extending never silently perturbs the fitted window;
  #   (2) the grid stays exactly what compartmentalFunction() would have built
  #       for that endpoint, so the projected part is on the same footing as a
  #       refit with a longer sheet.
  skip_if_not_installed("deSolve")

  dir <- fixture_dir("minimal")
  sc  <- load_scenario(dir, combined_file = "dataCombined.csv",
                       dummy_file = "dataDummy.csv", params_file = "modelParams.csv")

  fit <- fitCompartmentalModel(sc$modelParams, sc$dataCombined, method = "lbfgsb",
                               solver = solver_control(backend = "r"),
                               checkpoint_file = tempfile(fileext = ".rds"))
  chk("fit succeeds", isTRUE(fit$success))

  tg <- fit$time_grid
  E  <- tg$endpoint
  P  <- tg$partition

  # ---- .extend_time_grid: grid arithmetic -----------------------------------
  chk("NULL endpoint is a no-op",  identical(.extend_time_grid(fit, NULL), fit))
  chk("same endpoint is a no-op",  identical(.extend_time_grid(fit, E),    fit))

  ext <- .extend_time_grid(fit, E + 10)
  n   <- (E + 10) - tg$startpoint + 1
  chk("endpoint updated",      ext$time_grid$endpoint == E + 10)
  chk("grid length is P*n+1", length(ext$time_grid$time) == P * n + 1)
  chk("grid still starts at 0", ext$time_grid$time[1] == 0)
  chk("grid ends at n",         ext$time_grid$time[length(ext$time_grid$time)] == n)
  chk("dates match the grid",   length(ext$model$date) == P * n + 1)
  chk("date opens a year early (carries the IC)",
      format(ext$model$date[1]) == paste0(tg$startpoint - 1, "-12-31"))
  chk("date closes on the new endpoint",
      format(ext$model$date[length(ext$model$date)]) == paste0(E + 10, "-12-31"))
  chk("startpoint untouched", ext$time_grid$startpoint == tg$startpoint)
  chk("cutoff untouched",     identical(ext$time_grid$cutoff, tg$cutoff))

  # The extended grid must equal what a sheet with that endpoint would produce --
  # otherwise the projected years sit on a different footing than fitted ones.
  mp2 <- sc$modelParams
  mp2$Others[grepl("^ *endpoint *=", mp2$Others)] <- paste0("endpoint=", E + 10)
  chk("grid matches a sheet written with that endpoint",
      isTRUE(all.equal(.time_grid(mp2)$time, ext$time_grid$time)))

  # ---- rejected endpoints ---------------------------------------------------
  bad <- function(e) inherits(tryCatch(.extend_time_grid(fit, e),
                                       error = function(x) x), "error")
  chk("truncation rejected",   bad(E - 1))
  chk("fractional year rejected", bad(E + 0.5))
  chk("non-numeric rejected",  bad("2040"))
  chk("NA rejected",           bad(NA_real_))
  chk("vector rejected",       bad(c(E + 1, E + 2)))

  # ---- solve_and_evaluate: extra years are added, not substituted -----------
  p    <- get_point(fit)
  base <- solve_and_evaluate(fit, p$initial_state, p$parms, sc$dataDummy)
  long <- solve_and_evaluate(fit, p$initial_state, p$parms, sc$dataDummy,
                             endpoint = E + 10)

  k <- nrow(base$evaluation)
  chk("extended run is longer", nrow(long$evaluation) == P * n + 1)
  chk("extended run covers the base grid", nrow(long$evaluation) > k)

  num <- vapply(base$evaluation, is.numeric, logical(1))
  chk("overlapping years bit-identical",
      max(abs(as.matrix(base$evaluation[, num, drop = FALSE]) -
              as.matrix(long$evaluation[seq_len(k), num, drop = FALSE])),
          na.rm = TRUE) == 0)
  # Dates over the overlap agree to within a DAY, not exactly. The date column is
  # built by spreading the whole span linearly over days and rounding
  # (compartmentalFunction.R), so a longer span puts leap days in different
  # places and a snapshot can land a day earlier. Pre-existing behaviour of the
  # date model, not an artefact of extending: the solve `time` (in years) IS
  # identical over the overlap, which is why the values above match exactly.
  # Only visible when comparing two separate runs -- a single plot has one grid.
  chk("overlapping dates agree within a day",
      max(abs(as.numeric(long$evaluation$date[seq_len(k)] - base$evaluation$date))) <= 1)
  chk("overlapping solve times identical",
      identical(fit$time_grid$time,
                .extend_time_grid(fit, E + 10)$time_grid$time[seq_len(length(fit$time_grid$time))]))
  chk("trajectory extended too", nrow(long$sir_out) == P * n + 1)
  chk("projected values finite",
      all(is.finite(as.matrix(long$evaluation[, num, drop = FALSE]))))

  # The caller's fit must not be mutated -- .extend_time_grid works on a copy.
  chk("caller's fit untouched", fit$time_grid$endpoint == E)

  # ---- get_central_point: one accessor for both fit kinds -------------------
  cp <- get_central_point(fit)
  chk("central point == MLE point for an optimiser fit", identical(cp, p))
})
