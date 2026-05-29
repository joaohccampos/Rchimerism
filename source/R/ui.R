build_ui <- function() {
  shiny::shinyUI({
    shiny::fluidPage(

    shiny::titlePanel("R Chimerism"),
    shiny::br(),
    shiny::br(),

    shiny::fluidRow(
      shiny::column(3,
         shiny::radioButtons("directory_mode", label = NULL,
                             choices = list("Normal Mode" = 1, "Directory Mode" = 2,
                                            "GeneMapper Mode" = 3),
                             selected = 1),

         shiny::fileInput("markers", label = "Marker File (default: GlobalFiler 23-loci)",
                          accept = c("csv", ".csv")),

         shiny::radioButtons("donor_type", label = NULL,
                             choices = list("Single Donor" = 1, "Double Donor" = 2),
                             selected = 1),

        shiny::conditionalPanel("input.directory_mode==1",

        shiny::fileInput("ddata", label = "Choose Donor Data Input",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv")),

        shiny::conditionalPanel(condition = "input.donor_type == 2",
          shiny::fileInput("d2data", label = "Choose Second Donor Data Input",
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv"))
        ),

        shiny::fileInput("rdata", label = "Choose Recipient Data Input",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv")),

        shiny::fileInput("sdata", label = "Choose Sample Data Input",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".tsv"))
        ),

        shiny::conditionalPanel("input.directory_mode == 2",
          shinyFiles::shinyDirButton('directory_select', 'Folder select', 'Please select a folder'),
          shiny::br(),
          shiny::br()
        ),

        shiny::conditionalPanel("input.directory_mode == 3",
          shiny::fileInput("gm_file", label = "GeneMapper Export File (.txt / .tsv)",
                           accept = c(".txt", ".tsv")),
          shiny::uiOutput("gm_role_ui")
        ),

        shiny::checkboxInput("ignore_unknown_alleles",
                             "Ignore unrecognized sample alleles",
                             value = FALSE),

        shiny::actionButton("run_locSD_button", "Read input files"),
        shiny::br(),
        shiny::br(),
        shiny::actionButton("new_analysis_button", "New Analysis")
      ),

      shiny::conditionalPanel("input.donor_type == 1",
        shiny::column(3, offset = 1, shiny::h4("Final Results"),
          shiny::br(),
          shiny::fluidRow(DT::dataTableOutput("results_table"))
        ),
        shiny::column(2, offset = 1,
          shiny::h4(shiny::strong("Donor% Mean")),
          shiny::br(),
          shiny::fluidRow(shiny::strong(shiny::verbatimTextOutput("donor_p_mean"))),
          shiny::h5("Recipient% Mean"),
          shiny::fluidRow(shiny::verbatimTextOutput("recip_p_mean")),
          shiny::h5("Donor% SD"),
          shiny::fluidRow(shiny::verbatimTextOutput("donor_p_SD")),
          shiny::h5("Donor% CV"),
          shiny::fluidRow(shiny::verbatimTextOutput("donor_p_CV")),
          shiny::br(),
          shiny::fluidRow(shiny::downloadButton("check_file_txt", "Download Check File")),
          shiny::fluidRow(shiny::downloadButton("results_xls", "Download Results Excel File"))
        )
      ),

      shiny::conditionalPanel("input.donor_type == 2",
        shiny::column(6, shiny::h4("Final Results"),
          shiny::br(),
          shiny::fluidRow(DT::dataTableOutput("dd_results_table")),
          shiny::fluidRow(shiny::h4("% Mean"), shiny::br(),
            shiny::column(4,
              shiny::fluidRow(
                shiny::column(5, shiny::strong("Donor_1 % Mean")),
                shiny::column(6, shiny::verbatimTextOutput("donor_1_p_mean"))
              ),
              shiny::fluidRow(
                shiny::column(5, shiny::strong("Donor_2 % Mean")),
                shiny::column(6, shiny::verbatimTextOutput("donor_2_p_mean"))
              ),
              shiny::fluidRow(
                shiny::column(5, shiny::strong("Recipient % Mean")),
                shiny::column(6, shiny::verbatimTextOutput("recipient_p_mean"))
              )
            ),
            shiny::column(4,
              shiny::fluidRow(
                shiny::column(4, shiny::h5("Donor_1 % SD")),
                shiny::column(6, shiny::verbatimTextOutput("donor_1_SD"))
              ),
              shiny::fluidRow(
                shiny::column(4, shiny::h5("Donor_2 % SD")),
                shiny::column(6, shiny::verbatimTextOutput("donor_2_SD"))
              ),
              shiny::fluidRow(
                shiny::column(4, shiny::h5("Recipient % SD")),
                shiny::column(6, shiny::verbatimTextOutput("recipient_SD"))
              )
            ),
            shiny::column(4,
              shiny::fluidRow(
                shiny::column(4, shiny::h5("Donor_1 % CV")),
                shiny::column(6, shiny::verbatimTextOutput("donor_1_CV"))
              ),
              shiny::fluidRow(
                shiny::column(4, shiny::h5("Donor_2 % CV")),
                shiny::column(6, shiny::verbatimTextOutput("donor_2_CV"))
              ),
              shiny::fluidRow(
                shiny::column(4, shiny::h5("Recipient % CV")),
                shiny::column(6, shiny::verbatimTextOutput("recipient_CV"))
              )
            )
          )
        ),
        shiny::column(2, offset = 1,
          shiny::br(), shiny::br(), shiny::br(), shiny::br(),
          shiny::fluidRow(shiny::downloadButton("dd_check_file_txt", "Download Check File")),
          shiny::fluidRow(shiny::downloadButton("dd_results_xls", "Download Results Excel File"))
        )
      )
    ),

    shiny::conditionalPanel("input.donor_type ==1",
      shiny::fluidRow(
        shiny::column(4, shiny::h4("Donor Matrix"),    DT::dataTableOutput("d_matrix")),
        shiny::column(4, shiny::h4("Recipient Matrix"), DT::dataTableOutput("r_matrix")),
        shiny::column(4, shiny::h4("Sample Matrix"),   DT::dataTableOutput("s_matrix"))
      )
    ),

    shiny::conditionalPanel("input.donor_type ==2",
      shiny::fluidRow(
        shiny::column(6, shiny::h4("Donor_1 Matrix"),   DT::dataTableOutput("d1_matrix")),
        shiny::column(6, shiny::h4("Donor_2 Matrix"),   DT::dataTableOutput("d2_matrix"))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::h4("Recipient Matrix"), DT::dataTableOutput("r_dd_matrix")),
        shiny::column(6, shiny::h4("Sample Matrix"),    DT::dataTableOutput("s_dd_matrix"))
      )
    )

  )
  })
}
