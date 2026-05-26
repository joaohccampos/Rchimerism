validate_input_format <- function(data, label) {
  if (ncol(data) < 7L || nrow(data) < 1L) {
    return(paste("Cannot read", label))
  }
  NULL
}
