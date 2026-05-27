GENEMAPPER_EXPECTED_COLUMNS <- c(
  "Dye.Sample.Peak",
  "Sample.File.Name",
  "Marker",
  "Allele",
  "Size",
  "Height",
  "Area",
  "Data.Point"
)

GENEMAPPER_AUTODETECT_PATTERNS <- list(
  donor     = c("DOADOR", "DONOR", "DON"),
  recipient = c("PRE", "RECEP", "RECIPIENT", "RECEPTOR"),
  sample    = c("QUI", "CHIM", "POS", "POST")
)

parse_genemapper_tsv <- function(path) {
  raw <- read.delim(
    path,
    header           = TRUE,
    sep              = "\t",
    na.strings       = c("", " "),
    colClasses       = c(Allele = "character"),
    fill             = TRUE,
    check.names      = TRUE,
    stringsAsFactors = FALSE
  )
  available <- intersect(GENEMAPPER_EXPECTED_COLUMNS, colnames(raw))
  raw[, available, drop = FALSE]
}

list_genemapper_samples <- function(df) {
  unique(as.character(df[["Sample.File.Name"]]))
}

detect_genemapper_role <- function(sample_name) {
  upper <- toupper(sample_name)
  for (role in names(GENEMAPPER_AUTODETECT_PATTERNS)) {
    patterns <- GENEMAPPER_AUTODETECT_PATTERNS[[role]]
    if (any(sapply(patterns, function(p) grepl(p, upper, fixed = TRUE)))) {
      return(role)
    }
  }
  NULL
}

find_auto_role <- function(sample_names, role) {
  for (name in sample_names) {
    detected <- detect_genemapper_role(name)
    if (!is.null(detected) && detected == role) return(name)
  }
  ""
}

extract_genemapper_sample <- function(df, sample_name) {
  df[df[["Sample.File.Name"]] == sample_name, , drop = FALSE]
}
