
#' Rchimerism chiSD
#'
#' An internal function to determine donor percentages for a single donor
#'
#'
#' @param sdata Sample data input text file
#' @param markers List of locus markers
#' @param profile,rt,dt,d,r Internal variables for matrix operations
#'
#' @return A data frame with donor percentage results, and sample allele matrix
#'
#'

chiSD <- function(sdata, markers, profile, rt, dt, d, r,
                  ignore_unknown_alleles = FALSE) {

sample_data <- sdata

ci <- c(validate_input_format(sample_data, "sample data"))
if (length(ci) > 0) {
  return(ci[[1]])
}

sample_peaks <- sample_data[grep("[^[:alpha:]]", sample_data[, 4]), c(3:4, 7)]
sample_peaks <- droplevels(sample_peaks)
sample_peaks$Allele <- as.factor(sample_peaks$Allele)
sample_presence <- rt
sample_presence[, ] <- 0

sample_peaks <- sample_peaks[(sample_peaks[, 1] %in% rownames(sample_presence)), ]

if (nrow(sample_peaks[!(sample_peaks[, 2] %in% colnames(sample_presence)), ]) != 0) {
  false_calls <- sample_peaks[!(sample_peaks[, 2] %in% colnames(sample_presence)), ]
  if (!ignore_unknown_alleles) return(false_calls)
  sample_peaks <- sample_peaks[sample_peaks[, 2] %in% colnames(sample_presence), ]
}

sample_presence[as.matrix(sample_peaks[, 1:2])] <- 1
sample_presence[(dt + rt) == 0 & sample_presence != 0] <- UNKNOWN_ALLELE_SENTINEL

sample_area <- rt
sample_area[, ] <- 0
sample_area[as.matrix(sample_peaks[, 1:2])] <- sample_peaks[, 3]

chimerism <- profile

for (m in markers) {

  if (profile[m] == 211) {
    area_donor  <- sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] == setdiff(d[d[, 1] == m, 2], r[r[, 1] == m, 2]), 3]
    area_shared <- sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] == intersect(d[d[, 1] == m, 2], r[r[, 1] == m, 2]), 3]
    if (length(area_donor) == 0) { area_donor <- 0 }
    chimerism[m] <- 2 * area_donor / (area_donor + area_shared)
  }

  if (profile[m] == 221) {
    area_donor  <- sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] == setdiff(d[d[, 1] == m, 2], r[r[, 1] == m, 2]), 3]
    area_recipient <- sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] == setdiff(r[r[, 1] == m, 2], d[d[, 1] == m, 2]), 3]
    if (length(area_recipient) == 0) { area_recipient <- 0 }
    if (length(area_donor) == 0) { area_donor <- 0 }
    chimerism[m] <- area_donor / (area_donor + area_recipient)
  }

  if (profile[m] == 121) {
    area_recipient <- sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] == setdiff(r[r[, 1] == m, 2], d[d[, 1] == m, 2]), 3]
    area_shared    <- sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] == intersect(r[r[, 1] == m, 2], d[d[, 1] == m, 2]), 3]
    if (length(area_recipient) == 0) { area_recipient <- 0 }
    chimerism[m] <- 1 - (2 * area_recipient / (area_recipient + area_shared))
  }

  if (profile[m] == 1) {
    area_donor  <- sum(sample_peaks[sample_peaks[, 1] == m & sample_peaks[, 2] %in% setdiff(d[d[, 1] == m, 2], r[r[, 1] == m, 2]), 3])
    area_total  <- sum(sample_peaks[sample_peaks[, 1] == m, 3])
    if (length(area_donor) == 0) { area_donor <- 0 }
    chimerism[m] <- area_donor / area_total
  }

  if (profile[m] == 0) {
    chimerism[m] <- NA
  }
}

results <- cbind(profile, chimerism, NA, NA, NA, NA)
locus_mean <- mean(results[, 2], na.rm = TRUE)
sd_val     <- sd(results[, 2], na.rm = TRUE)

informative <- !((abs(results[, 2] - locus_mean) > 2 * sd_val) | is.na(results[, 2]))
results[informative, 3] <- mean(results[informative, 2])
results[informative, 4] <- sd(results[informative, 2])
results[informative, 5] <- results[informative, 4] / results[informative, 3]
results[informative, 6] <- 1 - results[informative, 3]

colnames(results)[1] <- 'Profile'
colnames(results)[2] <- 'Donor%'
colnames(results)[3] <- 'Donor%_Mean'
colnames(results)[4] <- 'Donor%_SD'
colnames(results)[5] <- 'Donor%_CV'
colnames(results)[6] <- 'Recipient%_Mean'

sample_matrix <- cbind(sample_presence, apply(sample_presence, 1, sum))
colnames(sample_matrix)[length(colnames(sample_matrix))] <- 'Sum'

return(list(results, sample_matrix))
}
