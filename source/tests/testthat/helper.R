# Helpers to build minimal GeneMapper-format data frames for testing.
# Column layout matches ABI GeneMapper export:
#   col1=Dye/Peak, col2=SampleFile, col3=Marker, col4=Allele,
#   col5=Size, col6=Height, col7=Area, col8=DataPoint, col9=extra
#
# locSD/locDD grep-filter keeps rows where col4 is NOT purely alphabetic,
# so numeric allele strings (e.g. "14") pass through.

make_row <- function(marker, allele, area, height = 500) {
  data.frame(
    Dye.Sample.Peak  = "B,1",
    Sample.File.Name = "sample",
    Marker           = marker,
    Allele           = as.character(allele),
    Size             = 100.0,
    Height           = height,
    Area             = area,
    Data.Point       = 1000L,
    Extra            = "",
    stringsAsFactors = FALSE
  )
}

make_file <- function(...) {
  rows <- list(...)
  do.call(rbind, rows)
}

# A minimal valid markers vector used across tests
TEST_MARKERS <- c("M1", "M2")

# Load the package's internal functions (not exported)
locSD  <- Rchimerism:::locSD
locDD  <- Rchimerism:::locDD
chiSD  <- Rchimerism:::chiSD
chiDD  <- Rchimerism:::chiDD
