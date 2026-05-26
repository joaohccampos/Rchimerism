
#' Rchimerism locDD
#'
#' An internal function to determine informative loci for double donor cases
#'
#'
#' @param donor1_data Donor 1 data input text file
#' @param donor2_data Donor 2 data input text file
#' @param recipient_data Recipient data input text file
#' @param markers List of locus markers
#'
#'
#' @return Internal variables used by chiDD.R

locDD <- function(donor1_data, donor2_data, recipient_data, markers) {

d1_data <- donor1_data
d2_data <- donor2_data
r_data <- recipient_data

ci <- c(
  validate_input_format(d1_data, "donor 1 data"),
  validate_input_format(d2_data, "donor 2 data"),
  validate_input_format(r_data, "recipient data")
)
if (length(ci) > 0) {
  return(ci[[1]])
}

recipient_peaks <- r_data[grep("[^[:alpha:]]", r_data[, 4]), c(3:4, 7)]
donor1_peaks <- d1_data[grep("[^[:alpha:]]", d1_data[, 4]), c(3:4, 7)]
donor2_peaks <- d2_data[grep("[^[:alpha:]]", d2_data[, 4]), c(3:4, 7)]

recipient_clean <- droplevels(recipient_peaks)
donor1_clean <- droplevels(donor1_peaks)
donor2_clean <- droplevels(donor2_peaks)

max_recipient <- tapply(recipient_clean[, 3], recipient_clean[, 1], max) / 2
recipient_peaks <- recipient_clean[recipient_clean[, 3] > max_recipient[recipient_clean[, 1]], ]
max_donor1 <- tapply(donor1_clean[, 3], donor1_clean[, 1], max) / 2
donor1_peaks <- donor1_clean[donor1_clean[, 3] > max_donor1[donor1_clean[, 1]], ]
max_donor2 <- tapply(donor2_clean[, 3], donor2_clean[, 1], max) / 2
donor2_peaks <- donor2_clean[donor2_clean[, 3] > max_donor2[donor2_clean[, 1]], ]

recipient_peaks[, 4] <- 'r'
donor1_peaks[, 4] <- 'd1'
donor2_peaks[, 4] <- 'd2'
combined <- rbind(recipient_peaks, donor1_peaks, donor2_peaks)

xtra_in_markers <- setdiff(markers, combined[, 1])
if (length(xtra_in_markers) != 0) {
  return(paste("'", xtra_in_markers, "' from markers not found in input data", sep = ""))
}

allele_table <- table(combined[, c(1, 2, 4)])
allele_names <- dimnames(allele_table)
rt <- matrix(allele_table[, , 'r'], nrow = dim(allele_table)[1],
  ncol = dim(allele_table)[2], dimnames = allele_names[1:2])
d1t <- matrix(allele_table[, , 'd1'], nrow = dim(allele_table)[1],
  ncol = dim(allele_table)[2], dimnames = allele_names[1:2])
d2t <- matrix(allele_table[, , 'd2'], nrow = dim(allele_table)[1],
  ncol = dim(allele_table)[2], dimnames = allele_names[1:2])
rt <- rt[markers, , drop = FALSE]
rt <- rt[, sort(colnames(rt)), drop = FALSE]
d1t <- d1t[markers, , drop = FALSE]
d1t <- d1t[, sort(colnames(d1t)), drop = FALSE]
d2t <- d2t[markers, , drop = FALSE]
d2t <- d2t[, sort(colnames(d2t)), drop = FALSE]

ru <- rt - d1t - d2t; ru[ru != 1] <- 0; run <- apply(ru, 1, sum)
d1u <- d1t - d2t - rt; d1u[d1u != 1] <- 0; d1un <- apply(d1u, 1, sum)
d2u <- d2t - d1t - rt; d2u[d2u != 1] <- 0; d2un <- apply(d2u, 1, sum)

rn <- apply(rt, 1, sum)
d1n <- apply(d1t, 1, sum)
d2n <- apply(d2t, 1, sum)

rur <- sum(run) / sum(rn)
d1ur <- sum(d1un) / sum(d1n)
d2ur <- sum(d2un) / sum(d2n)

rnn <- rn; rnn[!((rn == 2) & (run == 1))] <- 1
d1nn <- d1n; d1nn[!((d1n == 2) & (d1un == 1))] <- 1
d2nn <- d2n; d2nn[!((d2n == 2) & (d2un == 1))] <- 1

profile <- as.integer(apply(d1u != 0, 1, any) | apply(d2u != 0, 1, any))
names(profile) <- markers

d1m <- cbind(d1t, d1n, d1un, d1nn)
colnames(d1m)[length(colnames(d1m)) - 2] <- 'Sum'
colnames(d1m)[length(colnames(d1m)) - 1] <- 'Unique'
colnames(d1m)[length(colnames(d1m))] <- 'Factor'
d2m <- cbind(d2t, d2n, d2un, d2nn); colnames(d2m) <- colnames(d1m)
rm <- cbind(rt, rn, run, rnn); colnames(rm) <- colnames(d1m)

return(list(markers, profile, ru, rt, rnn, d1nn, d2nn,
  d1u, d2u, d1t, d2t, recipient_peaks, d1m, d2m, rm))
}
