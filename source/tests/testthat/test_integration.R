EXAMPLE_DIR <- normalizePath(file.path(getwd(), "../example"), mustWork = FALSE)
if (!dir.exists(EXAMPLE_DIR)) {
  EXAMPLE_DIR <- file.path(getwd(), "example")
}
if (!dir.exists(EXAMPLE_DIR)) {
  EXAMPLE_DIR <- system.file("example", package = "Rchimerism")
}

load_example <- function(...) read.delim(file.path(EXAMPLE_DIR, ...))
load_markers <- function() {
  strsplit(readLines(file.path(EXAMPLE_DIR, "markers.csv")), ",")[[1]]
}

test_that("integration_single_donor_end_to_end", {
  skip_if_not(dir.exists(EXAMPLE_DIR), "example/ directory not found")
  markers <- load_markers()
  ddata   <- load_example("singleDonorExample", "ddata.txt")
  rdata   <- load_example("singleDonorExample", "rdata.txt")
  sdata   <- load_example("singleDonorExample", "sdata.txt")

  loc_out <- locSD(ddata, rdata, markers)
  expect_type(loc_out, "list")

  chi_out <- chiSD(sdata, markers,
                   loc_out[[2]], loc_out[[3]], loc_out[[4]],
                   loc_out[[7]], loc_out[[8]])
  expect_type(chi_out, "list")

  res         <- chi_out[[1]]
  donor_mean  <- na.omit(res[, "Donor%_Mean"])[[1]] * 100
  recip_mean  <- na.omit(res[, "Recipient%_Mean"])[[1]] * 100

  # Known expected values from example dataset
  expect_equal(donor_mean, 99.897, tolerance = 0.01)
  expect_equal(recip_mean,  0.103, tolerance = 0.01)
})

test_that("integration_double_donor_end_to_end", {
  skip_if_not(dir.exists(EXAMPLE_DIR), "example/ directory not found")
  markers <- load_markers()
  d1data  <- load_example("doubleDonorExample", "d1data.txt")
  d2data  <- load_example("doubleDonorExample", "d2data.txt")
  rdata   <- load_example("doubleDonorExample", "rdata.txt")
  sdata   <- load_example("doubleDonorExample", "sdata.txt")

  loc_out <- locDD(d1data, d2data, rdata, markers)
  expect_type(loc_out, "list")
  expect_length(loc_out, 15)

  chi_out <- chiDD(sdata, markers,
                   loc_out[[2]],  loc_out[[3]],  loc_out[[4]],
                   loc_out[[5]],  loc_out[[6]],  loc_out[[7]],
                   loc_out[[8]],  loc_out[[9]],  loc_out[[10]],
                   loc_out[[11]], loc_out[[12]])
  expect_type(chi_out, "list")

  res     <- chi_out[[1]]
  data_rows <- res[-nrow(res), ]   # exclude Info# row

  d1_mean <- na.omit(data_rows[, "Donor_1%_Mean"])[[1]] * 100
  d2_mean <- na.omit(data_rows[, "Donor_2%_Mean"])[[1]] * 100
  r_mean  <- na.omit(data_rows[, "Recipient%_Mean"])[[1]] * 100

  # Known expected values from example dataset
  expect_equal(d1_mean, 70.071, tolerance = 0.5)
  expect_equal(d2_mean, 18.708, tolerance = 0.5)
  expect_equal(r_mean,   0.0,   tolerance = 1.0)
})
