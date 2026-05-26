
#' Rchimerism locSD
#'
#' An internal function to determine informative loci for a single donor
#'
#'
#' @param ddata Donor data input text file
#' @param rdata Recipient data input text file
#' @param markers List of locus markers
#'
#' @return Internal variables used by chiSD.R

locSD <- function(ddata, rdata, markers) {

donor_data <- ddata
recipient_data <- rdata

ci <- c(
  validate_input_format(donor_data, "donor data"),
  validate_input_format(recipient_data, "recipient data")
)
if (length(ci) > 0) {
  return(ci[[1]])
}

donor_peaks <- donor_data[grep("[^[:alpha:]]", donor_data[, 4]), c(3:4, 7)]
recipient_peaks <- recipient_data[grep("[^[:alpha:]]", recipient_data[, 4]), c(3:4, 7)]

donor_clean <- droplevels(donor_peaks)
recipient_clean <- droplevels(recipient_peaks)

max_donor <- tapply(donor_clean[, 3], donor_clean[, 1], max) / 2
donor_peaks <- donor_clean[donor_clean[, 3] > max_donor[donor_clean[, 1]], ]
max_recipient <- tapply(recipient_clean[, 3], recipient_clean[, 1], max) / 2
recipient_peaks <- recipient_clean[recipient_clean[, 3] > max_recipient[recipient_clean[, 1]], ]

recipient_peaks[, 4] <- 'r'
donor_peaks[, 4] <- 'd'
combined <- rbind(recipient_peaks, donor_peaks)

xtra_in_markers <- setdiff(markers, combined[, 1])
if (length(xtra_in_markers) != 0) {
  return(paste("'", xtra_in_markers, "' from markers not found in input data", sep = ""))
}

allele_table <- table(combined[, c(1, 2, 4)])
allele_names <- dimnames(allele_table)
rt <- matrix(allele_table[, , 'r'], nrow = dim(allele_table)[1],
  ncol = dim(allele_table)[2], dimnames = allele_names[1:2])
dt <- matrix(allele_table[, , 'd'], nrow = dim(allele_table)[1],
  ncol = dim(allele_table)[2], dimnames = allele_names[1:2])
rt <- rt[markers, , drop = FALSE]
rt <- rt[, sort(colnames(rt)), drop = FALSE]
dt <- dt[markers, , drop = FALSE]
dt <- dt[, sort(colnames(dt)), drop = FALSE]

allele_sum <- dt + rt
allele_diff <- dt - rt

profile <- allele_diff[, 1]
names(profile) <- rownames(allele_diff)
profile[apply(allele_diff, 1, any)] <- 1
profile[!apply(allele_diff, 1, any)] <- 0

allele_total <- apply(allele_sum, 1, sum)
net_diff <- apply(allele_diff, 1, sum)

locus_221 <- apply(allele_sum == 1, 1, any) & apply(allele_sum == 2, 1, any) &
  allele_total == 4 & net_diff == 0
profile[locus_221] <- 221

locus_121 <- allele_total == 3 & net_diff == -1 & apply(allele_sum == 2, 1, any)
profile[locus_121] <- 121

locus_211 <- allele_total == 3 & net_diff == 1 & apply(allele_sum == 2, 1, any)
profile[locus_211] <- 211

profile <- profile[markers]

dm <- cbind(dt, apply(dt, 1, sum), profile)
colnames(dm)[length(colnames(dm)) - 1] <- 'Sum'
colnames(dm)[length(colnames(dm))] <- 'Profile'
rm <- cbind(rt, apply(rt, 1, sum))
colnames(rm)[length(colnames(rm))] <- 'Sum'

return(list(markers, profile, rt, dt, dm, rm, donor_peaks, recipient_peaks))
}
