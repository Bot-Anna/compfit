test_that("test-unit-palette", {
# ============================================================
# test-unit-palette.R   (pure R; no Julia)
# The colour-scheme switch: the Okabe-Ito default vs the grayscale palette,
# resolved from a per-call argument or the global options(compfit.palette=).
# ============================================================

th_load_pure(c("plots.R"))

th_section("palette lists")
chk("grey palette has the same slots as the default",
    identical(sort(names(cfit_palette)), sort(names(cfit_palette_grey))))
chk("default and grey differ", !identical(cfit_palette, cfit_palette_grey))
chk("grey slots are all hex", all(grepl("^#[0-9A-Fa-f]{6}$", unlist(cfit_palette_grey))))

th_section(".cfit_is_grey_name aliases")
chk("grey names recognised",
    all(vapply(c("grey", "gray", "greyscale", "grayscale", "bw", "mono"),
               compfit:::.cfit_is_grey_name, logical(1))))
chk("colour names not grey",
    !any(vapply(c("okabe", "colour", "color", "default"),
                compfit:::.cfit_is_grey_name, logical(1))))

th_section(".cfit_resolve_palette")
chk("NULL -> Okabe default",  identical(compfit:::.cfit_resolve_palette(NULL), cfit_palette))
chk("'okabe' -> default",     identical(compfit:::.cfit_resolve_palette("okabe"), cfit_palette))
chk("'grey' -> grayscale",    identical(compfit:::.cfit_resolve_palette("grey"), cfit_palette_grey))
chk("'grayscale' -> grayscale", identical(compfit:::.cfit_resolve_palette("grayscale"), cfit_palette_grey))
custom <- list(data = "#111111", model = "#222222", band = "#222222",
               censor = "#333333", dummy = "#444444")
chk("custom list passes through", identical(compfit:::.cfit_resolve_palette(custom), custom))

th_section("global option flips the active palette")
old <- options(compfit.palette = "grey")
on.exit(options(old), add = TRUE)
chk("active palette follows the option", identical(compfit:::.cfit_active_palette(), cfit_palette_grey))
options(compfit.palette = "okabe")
chk("active palette back to Okabe", identical(compfit:::.cfit_active_palette(), cfit_palette))

th_summary("unit-palette")
})
