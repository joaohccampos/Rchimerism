test_file_path <- function(path, file_name) {
  if (!file.exists(path)) {
    t1 <- paste("Missing ", file_name, " input file at path:")
    t2 <- paste("<i>", trim_path(path), "</i>")
    shiny::showModal(shiny::modalDialog(
      title = paste("Missing ", file_name, " input file"),
      shiny::HTML(paste(t1, t2, sep = '<br/>'))
    ))
    return(FALSE)
  }
}

test_sdata_path <- function(path, file_name) {
  if (!file.exists(path)) {
    t1 <- paste("Missing ", file_name, " input file")
    t2 <- paste("<i>", path, "</i>")
    shiny::showModal(shiny::modalDialog(
      title = paste("Missing ", file_name, " input file"),
      shiny::HTML(paste(t1, t2, sep = '<br/>'))
    ))
    return(FALSE)
  }
}

trim_path <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  path_parts <- unlist(strsplit(path, .Platform$file.sep))
  path_parts <- path_parts[path_parts != ""]
  n <- length(path_parts)
  if (n >= 2) {
    path_parts <- path_parts[-(n - 1)]
    path_parts <- path_parts[-(length(path_parts))]
  }
  paste(path_parts, collapse = "/")
}

validate_path <- function(path, file_name) {
  shiny::validate(test_file_path(path, file_name))
}

bad_input <- function(input_file, ext) {
  if (is.null(input_file)) {
    shiny::showModal(shiny::modalDialog(title = "Missing input file"))
    return(FALSE)
  } else if (ext != tools::file_ext(input_file$datapath)) {
    shiny::showModal(shiny::modalDialog(
      title = paste("Wrong File Extension for '", input_file[1],
                    "', requires '.", ext, "'", sep = "")
    ))
    return(FALSE)
  } else {
    return(NULL)
  }
}

check_input <- function(input_file, ext) {
  shiny::validate(bad_input(input_file, ext))
}

safe_read_delim <- function(path, label) {
  tryCatch(
    read.delim(path),
    error = function(error) {
      shiny::showModal(shiny::modalDialog(
        title = paste("Error reading", label),
        paste("Could not read file:", conditionMessage(error))
      ))
      shiny::validate(shiny::need(FALSE, paste("Could not read", label)))
    }
  )
}
