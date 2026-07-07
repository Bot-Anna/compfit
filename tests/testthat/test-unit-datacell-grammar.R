test_that("test-unit-datacell-grammar", {
  # End-to-end coverage for the interval [A,B] and asymmetric A->B / A+/A- data
  # cells: .prepare_data masks, the MLE fit (pure-R backend, no Julia), the
  # A+/A- "no asym= declared" error, the no-op equivalences with right-censoring,
  # and the generated Turing @model branches (string-only, no Julia). Fast and
  # Julia-free, so it runs on CI/CRAN.
  skip_if_not_installed("deSolve")

  dir <- fixture_dir("minimal")
  sc <- load_scenario(dir, combined_file = "dataCombined.csv",
                      dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
  tg <- compfit:::.time_grid(sc$modelParams)

  ## ---- .prepare_data builds the interval / asymmetric masks ----
  dc <- sc$dataCombined
  dc$Likelihood[1] <- "gaussian; asym=5"          # stream X1 global asym deviation
  dc[["2016"]][1] <- "[75,90]"                     # interval (hard edges)
  dc[["2017"]][1] <- "67->72"                      # asym explicit (dev 5, up)
  dc[["2018"]][1] <- "55+"                         # asym global (dev 5, up)
  dc[["2019"]][1] <- ">=30"                        # right-censored (cap + up arrow)
  dc[["2015"]][2] <- "[80,95]~10"                  # soft interval on stream X2
  d <- compfit:::.prepare_data(dc, tg)
  chk("interval mask + edges", d$ilow_mat[2, 1] == 75 && d$iupp_mat[2, 1] == 90)
  chk("hard interval has zero shoulders", d$idev_lo_mat[2, 1] == 0 && d$idev_hi_mat[2, 1] == 0)
  chk("soft interval shoulders", d$interval_mask[1, 2] == 1 &&
        d$ilow_mat[1, 2] == 80 && d$iupp_mat[1, 2] == 95 &&
        d$idev_lo_mat[1, 2] == 10 && d$idev_hi_mat[1, 2] == 10)
  chk("asym explicit A->B",   d$asym_val_mat[3, 1] == 67 && d$asym_dev_mat[3, 1] == 5 && d$asym_dir_mat[3, 1] == 1)
  chk("asym global A+",       d$asym_val_mat[4, 1] == 55 && d$asym_dev_mat[4, 1] == 5 && d$asym_dir_mat[4, 1] == 1)

  ## ---- MLE fit (R backend) runs with the new cells ----
  fit <- fitCompartmentalModel(sc$modelParams, dc, method = "lbfgsb",
                               solver = solver_control(backend = "r"),
                               checkpoint_file = tempfile(fileext = ".rds"))
  chk("MLE fit with interval/asym cells succeeds", isTRUE(fit$success))
  chk("point recovered", !is.null(tryCatch(get_point(fit), error = function(e) NULL)))
  res <- chk_ok("plot_fit builds with interval/asym markers", plot_fit(fit))
  # No data dot on non-observed years; special-cell markers are grey.
  p1 <- res$plots[[1]]
  ptL <- Filter(function(L) inherits(L$geom, "GeomPoint") &&
                  identical(L$aes_params$colour, cfit_palette$data), p1$layers)[[1]]
  li  <- which(vapply(p1$layers, function(L) identical(L, ptL), logical(1)))[1]
  n_obs <- sum(fit$data$obs_mask[, 1] == 1)
  chk("data dots only on observed years", sum(!is.na(ggplot2::layer_data(p1, li)$y)) == n_obs)
  chk("special-cell markers are grey (not orange)",
      cfit_palette$censor == "#666666" &&
      any(vapply(p1$layers, function(L) identical(L$aes_params$colour, "#666666"), logical(1))))
  # Right-censored (>=) marker is an arrow (GeomSegment with an arrow grob), not a triangle.
  seg <- Filter(function(L) inherits(L$geom, "GeomSegment"), p1$layers)
  chk("censoring drawn as an arrow, not a triangle",
      length(seg) >= 1 && !is.null(seg[[1]]$geom_params$arrow) &&
      !any(vapply(p1$layers, function(L) identical(L$aes_params$shape, 2L) ||
                                          identical(L$aes_params$shape, 6L), logical(1))))

  ## ---- A+/A- with no stream asym= is an error naming the stream ----
  dc_bad <- sc$dataCombined; dc_bad[["2018"]][1] <- "55+"   # Likelihood plain "gaussian"
  chk_error("A+ without asym= errors", compfit:::.prepare_data(dc_bad, tg))

  ## ---- no-op equivalences (MLE): [A, huge] and A->huge behave like >=A ----
  mkfit <- function(cell) {
    d3 <- sc$dataCombined; d3[["2016"]][1] <- cell
    fitCompartmentalModel(sc$modelParams, d3, method = "lbfgsb",
                          solver = solver_control(backend = "r"),
                          checkpoint_file = tempfile(fileext = ".rds"))
  }
  x1 <- function(f) get_point(f)$initial_state[["X1"]]
  ref <- x1(mkfit(">=50"))
  chk("interval [A,huge] ~ right-censor >=A", abs(x1(mkfit("[50,1e9]")) - ref) < 1)
  chk("asym A->huge (weak damping) ~ right-censor >=A", abs(x1(mkfit("50->1e9")) - ref) < 5)

  ## ---- generated Turing @model carries the new branches (no Julia) ----
  ps <- buildPriorSpec(sc$modelParams)
  ls <- lapply(c("gaussian", "gaussian"), parseLikelihood)
  mc <- buildJuliaBayesModel(ps, ls, c("X1", "X2"), 6, 4, 2)
  chk("model_code interval branch", grepl("cf_logsubexp(logcdf", mc, fixed = TRUE))
  chk("model_code soft-interval branch (no CDF)", grepl("idev_lo_mat[k", mc, fixed = TRUE))
  chk("model_code asym branch",     grepl("asym_dir_mat[k",     mc, fixed = TRUE))
})
