# ── helpers ───────────────────────────────────────────────────────────────────

build_locDD <- function(d1, d2, r, markers) locDD(d1, d2, r, markers)

# Unpack locDD output into named list for readability in tests
unpack_locDD <- function(loc) {
  list(markers = loc[[1]], profile = loc[[2]], ru = loc[[3]], rt = loc[[4]],
       rnn = loc[[5]], d1nn = loc[[6]], d2nn = loc[[7]],
       d1u = loc[[8]], d2u = loc[[9]], d1t = loc[[10]], d2t = loc[[11]],
       r   = loc[[12]])
}

# ── tests ──────────────────────────────────────────────────────────────────────

test_that("chiDD_invalid_sdata_returns_error", {
  bad_s <- data.frame(V1 = 1)
  d1 <- make_file(make_row("M1", "14", 1000))
  d2 <- make_file(make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(bad_s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  expect_true(is.character(result))
  expect_match(result, "Cannot read")
})

test_that("chiDD_false_allele_in_sample_returns_dataframe", {
  # Allele "99" not in d1, d2, or r → triggers false-call return
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  s  <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 400),
                  make_row("M1", "99", 300))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  expect_s3_class(result, "data.frame")
})

test_that("chiDD_dys391_locus_excluded_from_processing", {
  # DYS391 rows removed before building matrices; result should exclude it
  d1 <- make_file(make_row("M1", "14", 1000), make_row("DYS391", "11", 1000))
  d2 <- make_file(make_row("M1", "18", 1000), make_row("DYS391", "11", 1000))
  r  <- make_file(make_row("M1", "22", 1000), make_row("DYS391", "11", 1000))
  s  <- make_file(make_row("M1", "14", 800),  make_row("DYS391", "11", 500))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  # DYS391 should not appear as a row in the sample matrix
  expect_false("DYS391" %in% rownames(result[[2]]))
})

test_that("chiDD_donor1_percentage_formula_correct", {
  # d1={14,16}, d2={14,18}: d1 shares allele 14 with d2 → d1 unique=16 only, d1nn=2
  # ss = 800+300+200+100 = 1400, sad1u = area of allele 16 = 300
  # sd1 = 300 * 2 / 1400 ≈ 0.4286
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "14", 1000), make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  s  <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 300),
                  make_row("M1", "18", 200), make_row("M1", "22", 100))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  res <- result[[1]]
  expect_equal(unname(res["M1", "Donor_1%"]), 300 * 2 / 1400, tolerance = 1e-6)
})

test_that("chiDD_donor2_percentage_formula_correct", {
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  s  <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 300),
                  make_row("M1", "18", 200), make_row("M1", "22", 100))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  res <- result[[1]]
  # d2 unique allele is 18, area=200; d2 has 1 allele → d2nn[M1]=1
  # sd2 = 200 * 1 / 1400 ≈ 0.1429
  expect_equal(unname(res["M1", "Donor_2%"]), 200 * 1 / 1400, tolerance = 1e-6)
})

test_that("chiDD_recipient_percentage_formula_correct", {
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  s  <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 300),
                  make_row("M1", "18", 200), make_row("M1", "22", 100))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  res <- result[[1]]
  # r unique allele is 22, area=100; r has 1 allele → rnn[M1]=1
  # sr = 100 * 1 / 1400 ≈ 0.0714
  expect_equal(unname(res["M1", "Recipient%"]), 100 * 1 / 1400, tolerance = 1e-6)
})

test_that("chiDD_sd1_mean_uses_only_informative_loci_not_outliers", {
  # 7 markers needed: for n=3, outlier is never >2*SD from mean (proven by math).
  # With n=7, ratio ≈ 1.32, so the clear outlier is definitively excluded.
  # M1-M6: d1={Mk1,Mk2}, d2={Mk3}, r={Mk4}; sd1 = (area_Mk1+area_Mk2)/ss = 0.3
  # M7 (outlier): same structure but d1-unique areas sum to 900/1000 → sd1 = 0.9
  mk_dd <- function(m, b, d1a_tot = 300) {
    sa1 <- d1a_tot %/% 2; sa2 <- d1a_tot - sa1
    list(
      d1 = make_file(make_row(m, b,     1000), make_row(m, b + 1, 1000)),
      d2 = make_file(make_row(m, b + 2, 1000)),
      r  = make_file(make_row(m, b + 3, 1000)),
      s  = make_file(make_row(m, b,     sa1),  make_row(m, b + 1, sa2),
                     make_row(m, b + 2, 400),  make_row(m, b + 3, 300))
    )
  }
  ms <- list(
    mk_dd("M1", 100), mk_dd("M2", 110), mk_dd("M3", 120),
    mk_dd("M4", 130), mk_dd("M5", 140), mk_dd("M6", 150),
    mk_dd("M7", 160, d1a_tot = 900)
  )
  d1 <- do.call(rbind, lapply(ms, `[[`, "d1"))
  d2 <- do.call(rbind, lapply(ms, `[[`, "d2"))
  r  <- do.call(rbind, lapply(ms, `[[`, "r"))
  s  <- do.call(rbind, lapply(ms, `[[`, "s"))
  markers <- paste0("M", 1:7)
  loc <- build_locDD(d1, d2, r, markers); p <- unpack_locDD(loc)
  result <- chiDD(s, markers, p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  res <- result[[1]]
  expect_true(is.na(res["M7", "Donor_1%_Mean"]))
  expect_false(is.na(res["M1", "Donor_1%_Mean"]))
  m1_mean  <- res["M1", "Donor_1%_Mean"]
  m1_d1pct <- res["M1", "Donor_1%"]
  m2_d1pct <- res["M2", "Donor_1%"]
  expect_equal(m1_mean, mean(c(m1_d1pct, m2_d1pct, res["M3","Donor_1%"],
                               res["M4","Donor_1%"], res["M5","Donor_1%"],
                               res["M6","Donor_1%"])), tolerance = 1e-6)
})

test_that("chiDD_sd_computation_uses_only_informative_loci", {
  # Same 7-marker setup: SD of Donor_1% uses only non-outlier loci (M1-M6),
  # not M7 — validates the 4 corrected SD/mean lines in chiDD.
  mk_dd <- function(m, b, d1a_tot = 300) {
    sa1 <- d1a_tot %/% 2; sa2 <- d1a_tot - sa1
    list(
      d1 = make_file(make_row(m, b,     1000), make_row(m, b + 1, 1000)),
      d2 = make_file(make_row(m, b + 2, 1000)),
      r  = make_file(make_row(m, b + 3, 1000)),
      s  = make_file(make_row(m, b,     sa1),  make_row(m, b + 1, sa2 + 10),
                     make_row(m, b + 2, 400),  make_row(m, b + 3, 290))
    )
  }
  ms <- list(
    mk_dd("M1", 100), mk_dd("M2", 110), mk_dd("M3", 120),
    mk_dd("M4", 130), mk_dd("M5", 140), mk_dd("M6", 150),
    mk_dd("M7", 160, d1a_tot = 900)
  )
  d1 <- do.call(rbind, lapply(ms, `[[`, "d1"))
  d2 <- do.call(rbind, lapply(ms, `[[`, "d2"))
  r  <- do.call(rbind, lapply(ms, `[[`, "r"))
  s  <- do.call(rbind, lapply(ms, `[[`, "s"))
  markers <- paste0("M", 1:7)
  loc <- build_locDD(d1, d2, r, markers); p <- unpack_locDD(loc)
  result <- chiDD(s, markers, p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  res <- result[[1]]
  m1_sd <- res["M1", "Donor_1%_SD"]
  expect_false(is.na(m1_sd))
  pcts <- res[paste0("M", 1:6), "Donor_1%"]
  expect_equal(m1_sd, sd(pcts), tolerance = 1e-6)
})

test_that("chiDD_info_row_counts_nonempty_loci", {
  d1 <- make_file(make_row("M1", "14", 1000), make_row("M1", "16", 1000))
  d2 <- make_file(make_row("M1", "18", 1000))
  r  <- make_file(make_row("M1", "22", 1000))
  s  <- make_file(make_row("M1", "14", 800), make_row("M1", "16", 300),
                  make_row("M1", "18", 200), make_row("M1", "22", 100))
  loc <- build_locDD(d1, d2, r, "M1"); p <- unpack_locDD(loc)
  result <- chiDD(s, "M1", p$profile, p$ru, p$rt, p$rnn, p$d1nn, p$d2nn,
                  p$d1u, p$d2u, p$d1t, p$d2t, p$r)
  res <- result[[1]]
  expect_equal(rownames(res)[nrow(res)], "Info#")
  # Info# for Donor_1% should count 1 (M1 is informative)
  expect_equal(unname(res["Info#", "Donor_1%"]), 1)
})
