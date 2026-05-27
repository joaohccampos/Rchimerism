parse_genemapper_tsv      <- Rchimerism:::parse_genemapper_tsv
list_genemapper_samples   <- Rchimerism:::list_genemapper_samples
detect_genemapper_role    <- Rchimerism:::detect_genemapper_role
find_auto_role            <- Rchimerism:::find_auto_role
extract_genemapper_sample <- Rchimerism:::extract_genemapper_sample

GM_HEADER <- paste(
  "Dye/Sample Peak", "Sample File Name", "Marker", "Allele",
  "Size", "Height", "Area", "Data Point",
  sep = "\t"
)

make_gm_row <- function(sample, marker, allele, area = 10000, height = 5000) {
  paste("B,1", sample, marker, allele, "120.0", height, area, "2000", sep = "\t")
}

write_gm_tsv <- function(...) {
  path <- tempfile(fileext = ".txt")
  writeLines(c(GM_HEADER, ...), path)
  path
}

# ── detect_genemapper_role ─────────────────────────────────────────────────

test_that("detect_genemapper_role_donor_DOADOR", {
  expect_equal(detect_genemapper_role("DOADOR_001.fsa"), "donor")
})

test_that("detect_genemapper_role_donor_DONOR", {
  expect_equal(detect_genemapper_role("DONOR_A.fsa"), "donor")
})

test_that("detect_genemapper_role_donor_DON", {
  expect_equal(detect_genemapper_role("DON_X"), "donor")
})

test_that("detect_genemapper_role_is_case_insensitive", {
  expect_equal(detect_genemapper_role("doador_001.fsa"), "donor")
})

test_that("detect_genemapper_role_recipient_PRE", {
  expect_equal(detect_genemapper_role("PRE_TX_001.fsa"), "recipient")
})

test_that("detect_genemapper_role_recipient_RECEP", {
  expect_equal(detect_genemapper_role("RECEP_B.fsa"), "recipient")
})

test_that("detect_genemapper_role_recipient_RECEPTOR", {
  expect_equal(detect_genemapper_role("RECEPTOR_C.fsa"), "recipient")
})

test_that("detect_genemapper_role_recipient_RECIPIENT", {
  expect_equal(detect_genemapper_role("RECIPIENT_D.fsa"), "recipient")
})

test_that("detect_genemapper_role_sample_QUI", {
  expect_equal(detect_genemapper_role("QUI_001.fsa"), "sample")
})

test_that("detect_genemapper_role_sample_CHIM", {
  expect_equal(detect_genemapper_role("CHIM_POS.fsa"), "sample")
})

test_that("detect_genemapper_role_sample_POS", {
  expect_equal(detect_genemapper_role("POS_030.fsa"), "sample")
})

test_that("detect_genemapper_role_sample_POST", {
  expect_equal(detect_genemapper_role("POST_TX.fsa"), "sample")
})

test_that("detect_genemapper_role_returns_null_for_unknown", {
  expect_null(detect_genemapper_role("UNRELATED_SAMPLE.fsa"))
})

test_that("detect_genemapper_role_returns_null_for_empty_string", {
  expect_null(detect_genemapper_role(""))
})

# ── find_auto_role ─────────────────────────────────────────────────────────

test_that("find_auto_role_returns_matching_name", {
  names <- c("DOADOR_001.fsa", "PRE_TX.fsa", "QUI_001.fsa")
  expect_equal(find_auto_role(names, "donor"), "DOADOR_001.fsa")
})

test_that("find_auto_role_returns_empty_string_when_no_match", {
  names <- c("DOADOR_001.fsa", "PRE_TX.fsa")
  expect_equal(find_auto_role(names, "sample"), "")
})

test_that("find_auto_role_returns_first_match_when_multiple_candidates", {
  names <- c("DOADOR_001.fsa", "DOADOR_002.fsa", "PRE_TX.fsa")
  expect_equal(find_auto_role(names, "donor"), "DOADOR_001.fsa")
})

test_that("find_auto_role_empty_vector_returns_empty_string", {
  expect_equal(find_auto_role(character(0), "donor"), "")
})

# ── parse_genemapper_tsv ───────────────────────────────────────────────────

test_that("parse_genemapper_tsv_contains_expected_columns", {
  path <- write_gm_tsv(make_gm_row("DOADOR.fsa", "M1", "14"))
  df   <- parse_genemapper_tsv(path)
  expect_true(all(c("Sample.File.Name", "Marker", "Allele", "Area") %in% colnames(df)))
})

test_that("parse_genemapper_tsv_drops_extra_columns", {
  header <- paste(GM_HEADER, "ExtraCol", sep = "\t")
  path   <- tempfile(fileext = ".txt")
  writeLines(c(header, paste(make_gm_row("DOADOR.fsa", "M1", "14"), "extra", sep = "\t")), path)
  df <- parse_genemapper_tsv(path)
  expect_false("ExtraCol" %in% colnames(df))
})

test_that("parse_genemapper_tsv_reads_allele_as_character", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("DOADOR.fsa", "M1", "OL")
  )
  df <- parse_genemapper_tsv(path)
  expect_type(df[["Allele"]], "character")
})

test_that("parse_genemapper_tsv_handles_trailing_tab_like_existing_files", {
  path <- tempfile(fileext = ".txt")
  writeLines(
    c(paste0(GM_HEADER, "\t"), paste0(make_gm_row("DOADOR.fsa", "M1", "14"), "\t")),
    path
  )
  df <- parse_genemapper_tsv(path)
  expect_equal(nrow(df), 1)
  expect_true("Marker" %in% colnames(df))
})

test_that("parse_genemapper_tsv_returns_subset_when_some_columns_absent", {
  path <- tempfile(fileext = ".txt")
  writeLines(c("Sample File Name\tMarker\tAllele\tArea", "DOADOR.fsa\tM1\t14\t10000"), path)
  df <- parse_genemapper_tsv(path)
  expect_true("Marker" %in% colnames(df))
  expect_true("Sample.File.Name" %in% colnames(df))
  expect_false("Dye.Sample.Peak" %in% colnames(df))
})

test_that("parse_genemapper_tsv_preserves_all_rows", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("PRE_TX.fsa", "M1", "16"),
    make_gm_row("QUI.fsa",    "M1", "14")
  )
  df <- parse_genemapper_tsv(path)
  expect_equal(nrow(df), 3)
})

# ── list_genemapper_samples ────────────────────────────────────────────────

test_that("list_genemapper_samples_deduplicates_names", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("DOADOR.fsa", "M1", "16"),
    make_gm_row("PRE_TX.fsa", "M1", "14")
  )
  df <- parse_genemapper_tsv(path)
  expect_equal(list_genemapper_samples(df), c("DOADOR.fsa", "PRE_TX.fsa"))
})

test_that("list_genemapper_samples_preserves_insertion_order", {
  path <- write_gm_tsv(
    make_gm_row("PRE_TX.fsa", "M1", "14"),
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("QUI.fsa",    "M1", "14")
  )
  df <- parse_genemapper_tsv(path)
  expect_equal(list_genemapper_samples(df), c("PRE_TX.fsa", "DOADOR.fsa", "QUI.fsa"))
})

# ── extract_genemapper_sample ──────────────────────────────────────────────

test_that("extract_genemapper_sample_returns_only_matching_rows", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("DOADOR.fsa", "M1", "16"),
    make_gm_row("PRE_TX.fsa", "M1", "14")
  )
  df     <- parse_genemapper_tsv(path)
  result <- extract_genemapper_sample(df, "DOADOR.fsa")
  expect_equal(nrow(result), 2)
  expect_true(all(result[["Sample.File.Name"]] == "DOADOR.fsa"))
})

test_that("extract_genemapper_sample_returns_empty_df_for_unknown_name", {
  path   <- write_gm_tsv(make_gm_row("DOADOR.fsa", "M1", "14"))
  df     <- parse_genemapper_tsv(path)
  result <- extract_genemapper_sample(df, "NONEXISTENT.fsa")
  expect_equal(nrow(result), 0)
})

test_that("extract_genemapper_sample_has_enough_columns_for_locSD", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("DOADOR.fsa", "M1", "16")
  )
  df     <- parse_genemapper_tsv(path)
  result <- extract_genemapper_sample(df, "DOADOR.fsa")
  expect_gte(ncol(result), 7L)
})

test_that("extract_genemapper_sample_preserves_allele_as_character", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14"),
    make_gm_row("PRE_TX.fsa", "M1", "OL")
  )
  df     <- parse_genemapper_tsv(path)
  result <- extract_genemapper_sample(df, "DOADOR.fsa")
  expect_type(result[["Allele"]], "character")
})

test_that("extract_genemapper_sample_result_feeds_locSD_without_error", {
  path <- write_gm_tsv(
    make_gm_row("DOADOR.fsa", "M1", "14", area = 8000),
    make_gm_row("DOADOR.fsa", "M1", "16", area = 6000),
    make_gm_row("PRE_TX.fsa", "M1", "14", area = 7000),
    make_gm_row("QUI.fsa",    "M1", "14", area = 5000)
  )
  df    <- parse_genemapper_tsv(path)
  ddata <- extract_genemapper_sample(df, "DOADOR.fsa")
  rdata <- extract_genemapper_sample(df, "PRE_TX.fsa")
  result <- Rchimerism:::locSD(ddata, rdata, "M1")
  expect_true(is.list(result))
  expect_false(is.character(result))
})
