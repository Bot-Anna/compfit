test_that("test-unit-quadratic-target", {
# ============================================================
# test-unit-quadratic-target.R   (pure R; no Julia)
# The Quadratic<i> cross-compartment routing prefix *goto* accepts a compartment
# INDEX (*2*) or a State NAME (*I*), resolved through .comp_index -- the same
# name<->index registry the _Level and Linear/Quadratic columns already use.
# The generated ODE is identical either way; the validator accepts a valid name
# and rejects an unknown name or an out-of-range index.
# ============================================================

th_load_pure(c("utils.R", "scenario.R", "numberOfComps.R", "statesAndParams.R",
               "generateExpressions.R", "compartmentalFunction.R", "validate.R"))

base <- read.csv(system.file("extdata/SIR/modelParams.csv", package = "compfit"),
                 stringsAsFactors = FALSE, check.names = FALSE)
# Rename the compartments X1/X2/X3 -> S/I/R (the routing target must then be
# reachable by either the index or the new name).
rn <- function(mp) { mp$States <- sub("X1", "S", mp$States)
                     mp$States <- sub("X2", "I", mp$States)
                     mp$States <- sub("X3", "R", mp$States); mp }
dX <- function(mp) { m <- suppressWarnings(compfit:::.build_model(mp, backend = "r"))
                     grep("dX", strsplit(m$stan_code, "\n")[[1]], value = TRUE, fixed = TRUE) }

th_section("name target *I* builds the same ODE as index target *2*")
a   <- dX(base)                                   # X-names, *2* (the shipped form)
mpN <- rn(base); mpN$Quadratic1[2] <- "*I*-beta"  # S/I/R names, target by NAME
mpI <- rn(base); mpI$Quadratic1[2] <- "*2*-beta"  # S/I/R names, target by INDEX
chk("name *I* routes identically to index *2*", identical(a, dX(mpN)))
chk("index *2* still works alongside named compartments", identical(a, dX(mpI)))

th_section("validator: valid name accepted, bad targets rejected")
ok  <- function(mp) tryCatch({ suppressWarnings(validate_modelParams(mp)); TRUE },
                             error = function(e) FALSE)
chk("*I* (a State name) validates", ok(mpN))
mpZ <- rn(base); mpZ$Quadratic1[2] <- "*Z*-beta"  # unknown compartment name
mp9 <- rn(base); mp9$Quadratic1[2] <- "*9*-beta"  # index out of range (n = 3)
chk("*Z* (unknown compartment) is rejected", !ok(mpZ))
chk("*9* (index out of range) is rejected",  !ok(mp9))

th_summary("quadratic-target")
})
