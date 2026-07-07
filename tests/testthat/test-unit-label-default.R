test_that("test-unit-label-default", {
# ============================================================
# test-unit-label-default.R   (pure R; no Julia)
# A missing or blank Label defaults to the row's Formula, for both dataCombined
# (via .prepare_data) and dataDummy (via .default_label), so plots/reports always
# have a name to show.
# ============================================================
th_load_pure(c("utils.R", "scenario.R", "fitCompartmentalModel.R"))

sc <- load_scenario(fixture_dir("SIS"), combined_file = "dataCombined.csv",
                    dummy_file = "dataDummy.csv", params_file = "modelParams.csv")
tg <- compfit:::.time_grid(sc$modelParams)

th_section(".default_label helper")
d1 <- compfit:::.default_label(data.frame(Formula = c("X1", "X2/X3"),
                                          stringsAsFactors = FALSE))          # no Label column
chk("absent Label column filled from Formula", identical(d1$Label, c("X1", "X2/X3")))
d2 <- compfit:::.default_label(data.frame(Label = c("S", NA, ""), Formula = c("X1", "X2", "X3"),
                                          stringsAsFactors = FALSE))
chk("blank/NA cells filled, present kept", identical(d2$Label, c("S", "X2", "X3")))
d3 <- compfit:::.default_label(data.frame(x = 1))                            # no Formula
chk("no Formula column -> unchanged", identical(names(d3), "x"))

th_section("dataCombined Label defaults inside .prepare_data")
dc <- sc$dataCombined; dc$Label <- NULL
dat <- compfit:::.prepare_data(dc, tg)
chk("stored data_combined has Label = Formula",
    all(dat$data_combined$Label == dat$data_combined$Formula))

th_summary("label-default")
})
