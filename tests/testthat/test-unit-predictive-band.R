test_that("test-unit-predictive-band", {
# ============================================================
# test-unit-predictive-band.R   (pure R; no Julia)
# Regression: the negbin predictive band must find the per-stream dispersion
# under BOTH backends' naming -- Julia emits phi<i> (phi1), Stan emits phi_<i>
# (phi_1). The band lookup previously used phi<i> only, so a Stan negbin fit got
# NA dispersion -> all-NA simulated observations -> no data-stream CrI.
# ============================================================

th_load_pure(c("plots.R"))

th_section(".cfit_phi_of tolerates Julia and Stan phi naming")
stan_row  <- data.frame(beta = 0.4, gamma = 0.2, phi_1 = 3.5)   # Stan naming
julia_row <- data.frame(beta = 0.4, gamma = 0.2, phi1  = 3.5)   # Julia naming
stan_2    <- data.frame(beta = 0.4, phi_1 = 2.0, phi_2 = 7.0)   # two negbin streams
chk_equal("Stan  phi_1 found", .cfit_phi_of(stan_row, 1), 3.5)
chk_equal("Julia phi1  found", .cfit_phi_of(julia_row, 1), 3.5)
chk_equal("Stan  phi_2 found for stream 2", .cfit_phi_of(stan_2, 2), 7.0)
chk("absent phi -> NA", is.na(.cfit_phi_of(stan_row, 9)))

th_section("negbin simulation is finite once phi is found (the actual symptom)")
# With the Stan-named draw, the negbin observation simulation must produce
# finite draws (size = phi). Before the fix phi was NA -> rnbinom(size = NA) ->
# all NA -> the band collapsed and no CrI was drawn.
sims <- .cfit_simulate_obs("negbin", mu = c(10, 20, 30),
                           noise = .cfit_phi_of(stan_row, 1), scale = 1, n_rep = 3)
chk("Stan negbin simulation yields finite observations", all(is.finite(sims)))
# Contrast: a genuinely missing phi (NA) still yields NA, as expected
# (rnbinom(size = NA) legitimately warns "NAs produced" -- that is the point).
sims_na <- suppressWarnings(.cfit_simulate_obs("negbin", mu = c(10, 20, 30),
                              noise = NA_real_, scale = 1, n_rep = 3))
chk("missing phi still gives NA sims (no false finite)", all(is.na(sims_na)))

th_summary("predictive-band")
})
