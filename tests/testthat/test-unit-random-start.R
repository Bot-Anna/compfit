test_that("test-unit-random-start", {
# ============================================================
# test-unit-random-start.R   (pure R; no Julia)
# The '|random' (and '|rand') initial-start keyword on a box prior: the fitted
# quantity draws its L-BFGS-B start uniformly from its box at fit time, seeded
# by optim_control(seed=) for reproducibility. Covers parsing -> random mask ->
# .bounds alignment -> the seed-tied draw in .fit_optim, plus the box-only guard.
# ============================================================

th_load_pure(c("utils.R", "numberOfComps.R", "statesAndParams.R",
               "fitCompartmentalModel.R"))

mp <- data.frame(
  States     = c("S=[0,10]", "I=[0,10]"),
  Parameters = c("beta=[0,3]|random", "gamma=[0,1]|1.5"),
  check.names = FALSE, stringsAsFactors = FALSE)

sap <- statesAndParams(mp)

th_section("parse: |random flagged, explicit |init not")
chk("beta flagged random",        isTRUE(unname(sap$random_params["beta"])))
chk("gamma (explicit |1.5) not random", isFALSE(unname(sap$random_params["gamma"])))

th_section(".bounds: random_init aligned to fitted order")
b <- .bounds(sap)
chk("random_init present + named", !is.null(b$random_init) && "beta" %in% names(b$random_init))
chk("beta random, gamma not", isTRUE(b$random_init[["beta"]]) && isFALSE(b$random_init[["gamma"]]))
chk("random entry init is the box midpoint placeholder (0.5 normalised)",
    b$init_norm[["beta"]] == 0.5)

th_section(".fit_optim: seed-tied uniform draw, non-random untouched")
start <- c(beta = 0.5, gamma = 0.5); lo <- c(beta = 0, gamma = 0); hi <- c(beta = 1, gamma = 1)
mask  <- c(beta = TRUE, gamma = FALSE)
grab  <- function(seed) {
  seen <- new.env(); seen$p <- NULL
  lf <- function(par, bs, flag) { if (is.null(seen$p)) seen$p <<- par; 0 }
  .fit_optim(lf, new.env(), start, lo, hi,
             list(opt_method = "L-BFGS-B", maxit = 1, factr = 1e7, seed = seed), mask)
  seen$p
}
a <- grab(123); b2 <- grab(123)
chk("same seed reproduces the drawn start", isTRUE(all.equal(a, b2)))
chk("drawn start lies in [0,1]", a[["beta"]] >= 0 && a[["beta"]] <= 1)
chk("non-random entry left at its start", a[["gamma"]] == 0.5)
c1 <- grab(NULL); c2 <- grab(NULL)
chk("seed=NULL gives a fresh draw each time", !isTRUE(all.equal(c1[["beta"]], c2[["beta"]])))
# no mask -> start passed through unchanged
seen <- new.env(); seen$p <- NULL
lf <- function(par, bs, flag) { if (is.null(seen$p)) seen$p <<- par; 0 }
.fit_optim(lf, new.env(), start, lo, hi,
           list(opt_method = "L-BFGS-B", maxit = 1, factr = 1e7, seed = 123), NULL)
chk("no random mask leaves the start untouched", seen$p[["beta"]] == 0.5)

th_section("|rand alias + box-only guard")
mp2 <- mp; mp2$Parameters[1] <- "beta=[0,3]|rand"
chk("'|rand' is accepted as an alias",
    isTRUE(unname(statesAndParams(mp2)$random_params["beta"])))
mp3 <- mp; mp3$Parameters[1] <- "beta=Normal(0,1)|random"   # not a box
chk("|random on a non-box prior errors",
    isTRUE(tryCatch({ statesAndParams(mp3); FALSE },
                    error = function(e) grepl("box prior", conditionMessage(e)))))

th_summary("random-start")
})
