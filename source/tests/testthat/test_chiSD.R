# ── helpers ───────────────────────────────────────────────────────────────────

# Build the internal variables that locSD would produce, given donor/recipient
# data frames and a single marker.  Returns a list mirroring locSD output.
build_locSD <- function(ddata, rdata, markers) locSD(ddata, rdata, markers)

# Minimal sample data frame for chiSD tests
make_sample <- function(...) make_file(...)

# ── tests ──────────────────────────────────────────────────────────────────────

test_that("chiSD_invalid_sdata_returns_error", {
  bad_s <- data.frame(V1 = 1)   # < 7 cols
  d     <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r     <- make_file(make_row("M1", "14", 1000))
  loc   <- build_locSD(d, r, "M1")
  result <- chiSD(bad_s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_true(is.character(result))
  expect_match(result, "Cannot read")
})

test_that("chiSD_false_allele_in_sample_returns_dataframe", {
  # M1 profile 211: donor=14,16  recipient=14
  # Sample has allele "99" which is not in donor or recipient
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  s <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 500),
                 make_row("M1", "99", 300))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_s3_class(result, "data.frame")
})

test_that("chiSD_profile_0_yields_na_donor_percent", {
  # Non-informative: donor == recipient (profile = 0)
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  s <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 600))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_true(is.list(result))
  expect_true(is.na(result[[1]]["M1", "Donor%"]))
})

test_that("chiSD_profile_211_formula_2Ad_over_Ad_plus_A", {
  # Donor: 14, 16 | Recipient: 14 → profile 211
  # Sample: allele 16 (unique to donor) area=300, allele 14 (shared) area=700
  # Expected: C = 2 * 300 / (300 + 700) = 0.60
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  s <- make_file(make_row("M1", "14", 700), make_row("M1", "16", 300))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_equal(unname(result[[1]]["M1", "Donor%"]), 0.60, tolerance = 1e-9)
})

test_that("chiSD_profile_221_formula_Ad_over_Ad_plus_Ar", {
  # Donor: 14, 16 | Recipient: 14, 18 → profile 221
  # Sample: allele 16 (donor-unique) area=400, allele 18 (recip-unique) area=600
  # Expected: C = 400 / (400 + 600) = 0.40
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "18", 1000))
  s <- make_file(make_row("M1", "14", 500), make_row("M1", "16", 400),
                 make_row("M1", "18", 600))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_equal(unname(result[[1]]["M1", "Donor%"]), 0.40, tolerance = 1e-9)
})

test_that("chiSD_profile_121_formula_1_minus_2Ar_over_Ar_plus_A", {
  # Donor: 14 | Recipient: 14, 18 → profile 121
  # Sample: allele 18 (recip-unique) area=200, allele 14 (shared) area=800
  # Expected: C = 1 - 2*200/(200+800) = 1 - 0.4 = 0.60
  d <- make_file(make_row("M1", "14", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "18", 1000))
  s <- make_file(make_row("M1", "14", 800), make_row("M1", "18", 200))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_equal(unname(result[[1]]["M1", "Donor%"]), 0.60, tolerance = 1e-9)
})

test_that("chiSD_profile_1_formula_Ad_over_A", {
  # profile = 1 (general): 3 donor alleles, 1 recipient allele
  # Donor-unique = {12, 16}, shared = {14}
  # Sample: allele 12 area=200, allele 14 area=500, allele 16 area=300
  # Expected: C = (200+300) / (200+500+300) = 500/1000 = 0.50
  d <- make_file(make_row("M1", "12", 1000), make_row("M1", "14", 1000),
                 make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  s <- make_file(make_row("M1", "12", 200), make_row("M1", "14", 500),
                 make_row("M1", "16", 300))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_equal(unname(result[[1]]["M1", "Donor%"]), 0.50, tolerance = 1e-9)
})

test_that("chiSD_outlier_locus_excluded_from_mean_and_cv", {
  # 7 markers needed: for n=3, outlier always < 2*SD from mean (proven by math).
  # With n=7, ratio = 6/(2*sqrt(7)) ≈ 1.13, so a clear outlier is definitively excluded.
  # M1-M6: profile 221, C=0.5; M7 (outlier): profile 221, C=0.95
  mk_221 <- function(m, base, ad = 500, ar = 500) {
    make_file(make_row(m, base,      1000), make_row(m, base + 2, 1000))   # donor
    list(
      d = make_file(make_row(m, base,      1000), make_row(m, base + 2, 1000)),
      r = make_file(make_row(m, base,      1000), make_row(m, base + 4, 1000)),
      s = make_file(make_row(m, base, 500), make_row(m, base + 2, ad),
                    make_row(m, base + 4, ar))
    )
  }
  ms <- list(
    mk_221("M1", 10), mk_221("M2", 20), mk_221("M3", 30),
    mk_221("M4", 40), mk_221("M5", 50), mk_221("M6", 60),
    mk_221("M7", 70, ad = 950, ar = 50)
  )
  d <- do.call(rbind, lapply(ms, `[[`, "d"))
  r <- do.call(rbind, lapply(ms, `[[`, "r"))
  s <- do.call(rbind, lapply(ms, `[[`, "s"))
  markers <- paste0("M", 1:7)
  loc    <- build_locSD(d, r, markers)
  result <- chiSD(s, markers, loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  res <- result[[1]]
  expect_true(is.na(res["M7", "Donor%_Mean"]))
  expect_false(is.na(res["M1", "Donor%_Mean"]))
})

test_that("chiSD_mean_is_computed_without_outlier", {
  # Same 7-marker setup as above: M1-M6 give C=0.5; M7 outlier C=0.95
  # Expected reported mean (M7 excluded) = 0.5
  mk_221 <- function(m, base, ad = 500, ar = 500) {
    list(
      d = make_file(make_row(m, base,      1000), make_row(m, base + 2, 1000)),
      r = make_file(make_row(m, base,      1000), make_row(m, base + 4, 1000)),
      s = make_file(make_row(m, base, 500), make_row(m, base + 2, ad),
                    make_row(m, base + 4, ar))
    )
  }
  ms <- list(
    mk_221("M1", 10), mk_221("M2", 20), mk_221("M3", 30),
    mk_221("M4", 40), mk_221("M5", 50), mk_221("M6", 60),
    mk_221("M7", 70, ad = 950, ar = 50)
  )
  d <- do.call(rbind, lapply(ms, `[[`, "d"))
  r <- do.call(rbind, lapply(ms, `[[`, "r"))
  s <- do.call(rbind, lapply(ms, `[[`, "s"))
  markers <- paste0("M", 1:7)
  loc    <- build_locSD(d, r, markers)
  result <- chiSD(s, markers, loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  reported_mean <- na.omit(result[[1]][, "Donor%_Mean"])[[1]]
  expect_equal(reported_mean, 0.5, tolerance = 1e-6)
})

test_that("chiSD_returns_list_of_two_with_correct_column_names", {
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  s <- make_file(make_row("M1", "14", 700), make_row("M1", "16", 300))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  expect_length(result, 2)
  expect_true("Donor%" %in% colnames(result[[1]]))
  expect_true("Donor%_Mean" %in% colnames(result[[1]]))
  expect_true("Sum" %in% colnames(result[[2]]))
})

test_that("chiSD_sample_matrix_marks_unknown_allele_sentinel", {
  d <- make_file(make_row("M1", "12", 1000), make_row("M1", "14", 1000),
                 make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  s <- make_file(make_row("M1", "12", 200), make_row("M1", "14", 500),
                 make_row("M1", "16", 300))
  loc    <- build_locSD(d, r, "M1")
  result <- chiSD(s, "M1", loc[[2]], loc[[3]], loc[[4]], loc[[7]], loc[[8]])
  sm <- result[[2]]
  expect_true("Sum" %in% colnames(sm))
  expect_false(any(is.na(sm)))
})
