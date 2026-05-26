test_that("locSD_invalid_ddata_too_few_cols_returns_error", {
  bad  <- data.frame(V1 = 1, V2 = 2)   # only 2 columns
  good <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  result <- locSD(bad, good, "M1")
  expect_true(is.character(result))
  expect_match(result, "Cannot read")
})

test_that("locSD_invalid_rdata_empty_returns_error", {
  good <- make_file(make_row("M1", "14", 1000))
  empty <- good[0, ]   # zero rows
  result <- locSD(good, empty, "M1")
  expect_true(is.character(result))
  expect_match(result, "Cannot read")
})

test_that("locSD_unknown_marker_returns_error", {
  d <- make_file(make_row("M1", "14", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  result <- locSD(d, r, c("M1", "MISSING"))
  expect_true(is.character(result))
  expect_match(result, "MISSING")
})

test_that("locSD_returns_list_of_8_elements", {
  # M1: donor=14,16 / recipient=14,16 (non-informative, profile=0)
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  result <- locSD(d, r, "M1")
  expect_type(result, "list")
  expect_length(result, 8)
})

test_that("locSD_profile_0_non_informative_locus", {
  # Donor and recipient share the same alleles → profile = 0
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  result <- locSD(d, r, "M1")
  profile <- result[[2]]
  expect_equal(unname(profile["M1"]), 0)
})

test_that("locSD_profile_1_general_informative_locus", {
  # Donor has 3 distinct alleles vs recipient 1 → profile = 1
  # ssum > 4 or sdiff > 1: falls into the "general" case
  d <- make_file(make_row("M1", "12", 1000), make_row("M1", "14", 1000),
                 make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  result <- locSD(d, r, "M1")
  profile <- result[[2]]
  expect_equal(unname(profile["M1"]), 1)
})

test_that("locSD_profile_211_two_donor_one_recipient_one_shared", {
  # Donor: 14, 16 (2 alleles) | Recipient: 14 (1 allele, shared with donor)
  # ssum = 3, sdiff = 1, and a column with sum=2 → 211
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000))
  result <- locSD(d, r, "M1")
  profile <- result[[2]]
  expect_equal(unname(profile["M1"]), 211)
})

test_that("locSD_profile_221_two_donor_two_recipient_one_shared", {
  # Donor: 14, 16 | Recipient: 14, 18 → one shared (14), each has one unique
  # ssum = 4, sdiff = 0 → 221
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "18", 1000))
  result <- locSD(d, r, "M1")
  profile <- result[[2]]
  expect_equal(unname(profile["M1"]), 221)
})

test_that("locSD_profile_121_one_donor_two_recipient_one_shared", {
  # Donor: 14 (1 allele) | Recipient: 14, 18 (2 alleles, one shared)
  # ssum = 3, sdiff = -1 → 121
  d <- make_file(make_row("M1", "14", 1000))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "18", 1000))
  result <- locSD(d, r, "M1")
  profile <- result[[2]]
  expect_equal(unname(profile["M1"]), 121)
})

test_that("locSD_noise_filtering_removes_low_area_peaks", {
  # Allele "99" has area 100 (< 1000/2 = 500) → should be filtered out
  # Only "14" (area=1000) and "16" (area=1000) survive
  d <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000),
                 make_row("M1", "99", 100))
  r <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  result <- locSD(d, r, "M1")
  # After filtering, donor and recipient have same alleles → non-informative
  expect_equal(unname(result[[2]]["M1"]), 0)
})
