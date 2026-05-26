#' Rchimerism
#'
#' Use Rchimerism through the shiny interface
#'
#' This function allows the user to input data files
#' and view the chimerism percentages along with the input data matrices
#' through the browser interface, with the options to download the results
#'
#' @name Rchimerism
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny br
#' @importFrom shiny fluidRow
#' @importFrom shiny column
#' @importFrom shiny fileInput
#' @importFrom shiny actionButton
#' @importFrom shiny radioButtons
#' @importFrom shiny conditionalPanel
#' @importFrom DT dataTableOutput
#' @importFrom shiny h4
#' @importFrom shiny verbatimTextOutput
#' @importFrom shiny downloadButton
#' @importFrom shiny h5
#' @importFrom DT renderDataTable
#' @importFrom DT datatable
#' @importFrom DT formatStyle
#' @importFrom DT styleEqual
#' @importFrom shinyFiles getVolumes
#' @importFrom shinyFiles shinyDirChoose
#' @importFrom shinyFiles shinyDirButton
#' @importFrom shinyFiles parseDirPath
#' @importFrom tools file_ext
#' @importFrom shiny reactive
#' @importFrom shiny validate
#' @importFrom shiny need
#' @importFrom shiny downloadHandler
#' @importFrom shiny runApp
#' @importFrom shiny stopApp
#' @importFrom shiny shinyApp
#' @importFrom shiny shinyUI
#'
#' @return Returns nothing

#' @export

Rchimerism <- function() {
  shiny::runApp(
    shiny::shinyApp(build_ui(), build_server()),
    quiet = TRUE,
    launch.browser = TRUE
  )
}
