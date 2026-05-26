test_that("locDD_invalid_d1data_returns_error", {
  bad  <- data.frame(V1 = 1)   # only 1 column
  good <- make_file(make_row("M1", "14", 1000))
  result <- locDD(bad, good, good, "M1")
  expect_true(is.character(result))
  expect_match(result, "Cannot read")
})

test_that("locDD_all_three_inputs_invalid_returns_first_error", {
  bad  <- data.frame(V1 = 1)
  result <- locDD(bad, bad, bad, "M1")
  expect_true(is.character(result))
  expect_length(result, 1)   # returns only the first error
})

test_that("locDD_unknown_marker_returns_error", {
  good <- make_file(make_row("M1", "14", 1000))
  result <- locDD(good, good, good, c("M1", "ABSENT"))
  expect_true(is.character(result))
  expect_match(result, "ABSENT")
})

test_that("locDD_returns_list_of_15_elements", {
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "18", 1000), make_row("M1", "20", 1000))
  r  <- make_file(make_row("M1", "14", 1000), make_row("M1", "22", 1000))
  result <- locDD(d1, d2, r, "M1")
  expect_type(result, "list")
  expect_length(result, 15)
})

test_that("locDD_profile_is_named_integer_vector", {
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  result <- locDD(d1, d2, r, "M1")
  profile <- result[[2]]
  expect_type(profile, "integer")
  expect_named(profile)
  expect_true("M1" %in% names(profile))
})

test_that("locDD_profile_informative_when_d1_has_unique_allele", {
  # d1 has allele 16 not present in d2 or r → locus is informative
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "14", 1000))
  r  <- make_file(make_row("M1", "14", 1000))
  result <- locDD(d1, d2, r, "M1")
  expect_equal(unname(result[[2]]["M1"]), 1L)
})

test_that("locDD_profile_noninformative_when_all_share_same_alleles", {
  # All three share exactly the same allele → no unique alleles anywhere
  d1 <- make_file(make_row("M1", "14", 1000))
  d2 <- make_file(make_row("M1", "14", 1000))
  r  <- make_file(make_row("M1", "14", 1000))
  result <- locDD(d1, d2, r, "M1")
  expect_equal(unname(result[[2]]["M1"]), 0L)
})

test_that("locDD_uses_area_col7_not_height_col6_for_noise_filtering", {
  # col6 (height) = 1 for allele 99 (would survive if col6 used)
  # col7 (area)   = 10 for allele 99, while allele 14 has area=1000
  # → allele 99 is noise and should be filtered (area 10 < 1000/2 = 500)
  d1 <- make_file(
    make_row("M1", "14", area = 1000, height = 1),
    make_row("M1", "99", area = 10,   height = 9999)  # high height, low area
  )
  d2 <- make_file(make_row("M1", "14", 1000))
  r  <- make_file(make_row("M1", "14", 1000))
  result <- locDD(d1, d2, r, "M1")
  # After filtering, d1 and r share allele 14 only → d1 has no unique allele
  expect_equal(unname(result[[2]]["M1"]), 0L)
})
