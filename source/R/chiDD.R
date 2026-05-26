
#' Rchimerism chiDD
#'
#' An internal function to determine donor percentages for a double donor
#'
#'
#' @param sdata Sample data input text file
#' @param markers List of locus markers
#' @param profile,ru,rt,rnn,d1nn,d2nn,d1u,d2u,d1t,d2t,r Internal variables for matrix operations
#'
#' @return A data frame with donor percentage results, and sample allele matrix
#'
#'

chiDD <- function(sdata, markers, profile, ru, rt, rnn, d1nn, d2nn, d1u, d2u, d1t, d2t, r) {

sample_data <- sdata

ci <- c(validate_input_format(sample_data, "sample data"))
if (length(ci) > 0) {
  return(ci[[1]])
}

sample_peaks <- sample_data[grep("[^[:alpha:]]", sample_data[, 4]), c(3:4, 7)]
sample_peaks <- sample_peaks[sample_peaks[, 1] != EXCLUDED_MARKERS, ]
sample_peaks <- droplevels(sample_peaks)
sample_peaks$Allele <- as.factor(sample_peaks$Allele)

sample_presence <- ru
sample_presence[, ] <- 0

sample_peaks <- sample_peaks[(sample_peaks[, 1] %in% rownames(sample_presence)), ]

if (nrow(sample_peaks[!(sample_peaks[, 2] %in% colnames(sample_presence)), ]) != 0) {
  false_calls <- sample_peaks[!(sample_peaks[, 2] %in% colnames(sample_presence)), ]
  return(false_calls)
}

sample_presence[as.matrix(sample_peaks[, 1:2])] <- 1
sample_presence[(d1t + d2t + rt) == 0 & sample_presence != 0] <- UNKNOWN_ALLELE_SENTINEL
allele_count <- apply(sample_presence, 1, sum)

sample_area <- ru
sample_area[, ] <- 0
sample_area[as.matrix(sample_peaks[, 1:2])] <- sample_peaks[, 3]

area_recipient_unique <- sample_area
area_recipient_unique[ru == 0] <- 0
area_donor1_unique <- sample_area
area_donor1_unique[d1u == 0] <- 0
area_donor2_unique <- sample_area
area_donor2_unique[d2u == 0] <- 0

total_area <- apply(sample_area, 1, sum)
area_recipient_sum <- apply(area_recipient_unique, 1, sum)
area_recipient_sum[apply(ru, 1, sum) == 0] <- NA
area_donor1_sum <- apply(area_donor1_unique, 1, sum)
area_donor1_sum[apply(d1u, 1, sum) == 0] <- NA
area_donor2_sum <- apply(area_donor2_unique, 1, sum)
area_donor2_sum[apply(d2u, 1, sum) == 0] <- NA

recipient_pct <- area_recipient_sum * rnn / total_area
recipient_mean <- mean(recipient_pct, na.rm = TRUE)
sd_val <- sd(recipient_pct, na.rm = TRUE)
keep <- !((abs(recipient_pct - recipient_mean) > 2 * sd_val) | is.na(recipient_pct))
recipient_pct_mean <- recipient_pct
recipient_pct_mean[!keep] <- NA
recipient_pct_mean[keep] <- mean(recipient_pct_mean, na.rm = TRUE)
recipient_pct_sd <- recipient_pct
recipient_pct_sd[!keep] <- NA
recipient_pct_sd[keep] <- sd(recipient_pct[keep])
recipient_pct_cv <- recipient_pct_sd / recipient_pct_mean

donor1_pct <- area_donor1_sum * d1nn / total_area
donor1_mean <- mean(donor1_pct, na.rm = TRUE)
sd_val <- sd(donor1_pct, na.rm = TRUE)
keep <- !((abs(donor1_pct - donor1_mean) > 2 * sd_val) | is.na(donor1_pct))
donor1_pct_mean <- donor1_pct
donor1_pct_mean[!keep] <- NA
donor1_pct_mean[keep] <- mean(donor1_pct_mean, na.rm = TRUE)
donor1_pct_sd <- donor1_pct
donor1_pct_sd[!keep] <- NA
donor1_pct_sd[keep] <- sd(donor1_pct[keep])
donor1_pct_cv <- donor1_pct_sd / donor1_pct_mean

donor2_pct <- area_donor2_sum * d2nn / total_area
donor2_mean <- mean(donor2_pct, na.rm = TRUE)
sd_val <- sd(donor2_pct, na.rm = TRUE)
keep <- !((abs(donor2_pct - donor2_mean) > 2 * sd_val) | is.na(donor2_pct))
donor2_pct_mean <- donor2_pct
donor2_pct_mean[!keep] <- NA
donor2_pct_mean[keep] <- mean(donor2_pct_mean, na.rm = TRUE)
donor2_pct_sd <- donor2_pct
donor2_pct_sd[!keep] <- NA
donor2_pct_sd[keep] <- sd(donor2_pct[keep])
donor2_pct_cv <- donor2_pct_sd / donor2_pct_mean

results <- cbind(
  donor1_pct, donor1_pct_mean, donor1_pct_sd, donor1_pct_cv,
  donor2_pct, donor2_pct_mean, donor2_pct_sd, donor2_pct_cv,
  recipient_pct, recipient_pct_mean, recipient_pct_sd, recipient_pct_cv,
  donor1_pct + donor2_pct + recipient_pct
)
results <- rbind(results, apply(!is.na(results), 2, sum))
rownames(results)[nrow(results)] <- "Info#"
colnames(results) <- c(
  'Donor_1%', 'Donor_1%_Mean', 'Donor_1%_SD', 'Donor_1%_CV',
  'Donor_2%', 'Donor_2%_Mean', 'Donor_2%_SD', 'Donor_2%_CV',
  'Recipient%', 'Recipient%_Mean', 'Recipient%_SD', 'Recipient%_CV',
  'Sum'
)

sample_matrix <- cbind(sample_presence, allele_count)
colnames(sample_matrix)[length(colnames(sample_matrix))] <- 'Sum'

return(list(results, sample_matrix))
}
