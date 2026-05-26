build_server <- function() {
  function(input, output, session) {

    volumes <- c(
      Home = "~",
      DirectoryExample = system.file("extdata/Directory_example", package = "Rchimerism"),
      shinyFiles::getVolumes()()
    )
    shinyFiles::shinyDirChoose(input, 'directory_select', roots = volumes, session = session)
    dirname <- shiny::reactive({ shinyFiles::parseDirPath(volumes, input$directory_select) })

    initialize_dir_mode_vars <- function() {
      ddata_path  <- paste0(dirname(), "/../../ddata.txt")
      d1data_path <- paste0(dirname(), "/../../d1data.txt")
      d2data_path <- paste0(dirname(), "/../../d2data.txt")
      rdata_path  <- paste0(dirname(), "/../../rdata.txt")
      sdata_path  <- paste0(dirname(), "/sdata.txt")

      if (input$donor_type == 1) {
        validate_path(ddata_path, "ddata.txt")
      } else {
        validate_path(d1data_path, "d1data.txt")
        validate_path(d2data_path, "d2data.txt")
      }

      validate_path(rdata_path, "rdata.txt")
      shiny::validate(test_sdata_path(sdata_path, "sdata.txt"))

      output_vars <- c(ddata_path, d1data_path, d2data_path, rdata_path, sdata_path)
      names(output_vars) <- c("ddata", "d1data", "d2data", "rdata", "sdata")
      return(output_vars)
    }

    shiny::observeEvent(input$run_locSD_button, {

      if (!is.null(input$markers)) {
        check_input(input$markers, "csv")
      }

      if (input$directory_mode == 1) {
        check_input(input$ddata, "txt")
        check_input(input$rdata, "txt")
        check_input(input$sdata, "txt")

        ddata <- safe_read_delim(input$ddata$datapath, "donor data")
        if (input$donor_type == 2) {
          check_input(input$d2data, "txt")
          d2data <- safe_read_delim(input$d2data$datapath, "donor 2 data")
        }
        rdata <- safe_read_delim(input$rdata$datapath, "recipient data")
        sdata <- safe_read_delim(input$sdata$datapath, "sample data")
      } else {
        dir_mode_vars <- initialize_dir_mode_vars()
        if (input$donor_type == 1) {
          ddata <- safe_read_delim(dir_mode_vars[["ddata"]], "donor data")
        } else {
          ddata  <- safe_read_delim(dir_mode_vars[["d1data"]], "donor 1 data")
          d2data <- safe_read_delim(dir_mode_vars[["d2data"]], "donor 2 data")
        }
        rdata <- safe_read_delim(dir_mode_vars[["rdata"]], "recipient data")
        sdata <- safe_read_delim(dir_mode_vars[["sdata"]], "sample data")
      }

      markers_path <- if (!is.null(input$markers)) {
        input$markers$datapath
      } else {
        system.file("extdata/globalfiler_markers.csv", package = "Rchimerism")
      }
      markers <- as.character(read.csv(
        markers_path, header = FALSE, nrows = 1,
        stringsAsFactors = FALSE
      ))

      incoherent_input <- function(any_input) {
        if (is.character(any_input)) {
          shiny::showModal(shiny::modalDialog(title = any_input))
          return(FALSE)
        } else {
          return(NULL)
        }
      }

      is_coherent_input <- function(any_output) {
        shiny::validate(incoherent_input(any_output))
      }

      check_sample_data <- function(any_chi_output) {
        if (length(any_chi_output) == 3) {
          shiny::showModal(shiny::modalDialog(
            title = "Sample Data Format Error",
            DT::renderDataTable(any_chi_output)
          ))
        }
        shiny::validate(shiny::need(length(any_chi_output) != 3, "Sample Data Format Error"))
      }

      if (input$donor_type == 2) {
        loc_dd_output <- locDD(ddata, d2data, rdata, markers)
        is_coherent_input(loc_dd_output)

        profile <- loc_dd_output[[2]]
        ru      <- loc_dd_output[[3]]
        rt      <- loc_dd_output[[4]]
        rnn     <- loc_dd_output[[5]]
        d1nn    <- loc_dd_output[[6]]
        d2nn    <- loc_dd_output[[7]]
        d1u     <- loc_dd_output[[8]]
        d2u     <- loc_dd_output[[9]]
        d1t     <- loc_dd_output[[10]]
        d2t     <- loc_dd_output[[11]]
        r       <- loc_dd_output[[12]]
        d1m     <- loc_dd_output[[13]]
        d2m     <- loc_dd_output[[14]]
        rm      <- loc_dd_output[[15]]

        chi_dd_output <- chiDD(sdata, markers, profile,
                               ru, rt, rnn, d1nn, d2nn, d1u, d2u, d1t, d2t, r)
        is_coherent_input(chi_dd_output)
        check_sample_data(chi_dd_output)

        results <- chi_dd_output[[1]]
        sm      <- chi_dd_output[[2]]
      } else {
        loc_sd_output <- locSD(ddata, rdata, markers)
        is_coherent_input(loc_sd_output)

        profile <- loc_sd_output[[2]]
        rt      <- loc_sd_output[[3]]
        dt      <- loc_sd_output[[4]]
        dm      <- loc_sd_output[[5]]
        rm      <- loc_sd_output[[6]]
        d       <- loc_sd_output[[7]]
        r       <- loc_sd_output[[8]]

        chi_sd_output <- chiSD(sdata, markers, profile, rt, dt, d, r)
        is_coherent_input(chi_sd_output)
        check_sample_data(chi_sd_output)

        results <- chi_sd_output[[1]]
        sm      <- chi_sd_output[[2]]
      }

      if (input$donor_type == 1) {
        printed_result      <- results[, 1:3]
        printed_result[, 1] <- sapply(printed_result[, 1], as.character)
        printed_result[, 2] <- results[, 2] * 100

        rdt    <- DT::datatable(printed_result, rownames = TRUE,
                                options = list(pageLength = 5, searching = FALSE,
                                               lengthChange = FALSE,
                                               columnDefs = list(list(visible = FALSE, targets = c(3)))))
        frdt   <- DT::formatStyle(rdt, columns = 2, valueColumns = 3,
                                  target = "row", color = DT::styleEqual(NA, "gray"))
        ifrdt  <- DT::formatStyle(frdt, columns = 0, fontStyle = 'italic')
        output$results_table <- DT::renderDataTable({ ifrdt })
      } else {
        dd_printed_result <- results[-nrow(results), c(1, 2, 5, 6, 9, 10, 13)]
        dd_printed_result <- dd_printed_result[, c(1:7)] * 100

        dd_rdt   <- DT::datatable(dd_printed_result, rownames = TRUE,
                                  options = list(pageLength = 5, searching = FALSE,
                                                 lengthChange = FALSE,
                                                 columnDefs = list(list(visible = FALSE, targets = c(2, 4, 6)))))
        dd_frdt  <- DT::formatStyle(dd_rdt, columns = c(1, 3, 5), valueColumns = c(2, 4, 6),
                                    color = DT::styleEqual(c(NA), c("gray")))
        dd_ifrdt <- DT::formatStyle(dd_frdt, columns = 0, fontStyle = 'italic')
        output$dd_results_table <- DT::renderDataTable(dd_ifrdt)
      }

      stat_pct <- function(col_index) {
        val <- na.omit(results[, col_index])
        if (length(val) == 0) return(NaN)
        round(val[[1]] * 100, 3)
      }

      if (input$donor_type == 1) {
        output$donor_p_mean <- shiny::renderText(stat_pct(3))
        output$recip_p_mean <- shiny::renderText(stat_pct(6))
        output$donor_p_SD   <- shiny::renderText(stat_pct(4))
        output$donor_p_CV   <- shiny::renderText(stat_pct(5))
      } else {
        p_output <- function(num) {
          shiny::renderText(format(round(
            max(results[-nrow(results), num], na.rm = TRUE) * 100, 3
          )))
        }
        output$donor_1_p_mean  <- p_output(2)
        output$donor_2_p_mean  <- p_output(6)
        output$recipient_p_mean <- p_output(10)
        output$donor_1_SD      <- p_output(3)
        output$donor_2_SD      <- p_output(7)
        output$recipient_SD    <- p_output(11)
        output$donor_1_CV      <- p_output(4)
        output$donor_2_CV      <- p_output(8)

        if (length(na.omit(results[-length(results[, 12]), 12])) == 0) {
          output$recipient_CV <- shiny::renderText(NaN)
        } else {
          output$recipient_CV <- p_output(12)
        }
      }

      matrix_output <- function(internal_matrix, rcol, srch) {
        if (missing(srch)) { srch <- "" }
        DT::renderDataTable({
          DT::formatStyle(
            DT::datatable(internal_matrix,
                          rownames = TRUE, extensions = "FixedColumns",
                          options = list(pageLength = 5, scrollX = TRUE,
                                         fixedColumns = list(leftColumns = 1, rightColumns = rcol),
                                         autoWidth = TRUE, search = list(search = srch))),
            columns = c(0), fontStyle = 'italic')
        })
      }

      if (input$donor_type == 1) {
        dm_minus_profile_col  <- dm[, -ncol(dm)]
        output$d_matrix       <- matrix_output(dm_minus_profile_col, 1)
        output$r_matrix       <- matrix_output(rm, 1)

        if (UNKNOWN_ALLELE_SENTINEL %in% sm) {
          output$s_matrix <- matrix_output(sm, 1, as.character(UNKNOWN_ALLELE_SENTINEL))
          shiny::showModal(shiny::modalDialog(
            paste("Clear '", UNKNOWN_ALLELE_SENTINEL, "' from Sample Matrix search box to remove filter"),
            title = paste("Found '", UNKNOWN_ALLELE_SENTINEL, "' value(s) in Sample Allele Matrix, review sample data")
          ))
        } else {
          output$s_matrix <- matrix_output(sm, 1)
        }
      } else {
        output$d1_matrix    <- matrix_output(d1m, 3)
        output$d2_matrix    <- matrix_output(d2m, 3)
        output$r_dd_matrix  <- matrix_output(rm,  3)

        if (UNKNOWN_ALLELE_SENTINEL %in% sm) {
          output$s_dd_matrix <- matrix_output(sm, 1, as.character(UNKNOWN_ALLELE_SENTINEL))
          shiny::showModal(shiny::modalDialog(
            paste("Clear '", UNKNOWN_ALLELE_SENTINEL, "' from Sample Matrix search box to remove filter"),
            title = paste("Found '", UNKNOWN_ALLELE_SENTINEL, "' value(s) in Sample Allele Matrix, review sample data")
          ))
        } else {
          output$s_dd_matrix <- matrix_output(sm, 1)
        }
      }

      return_xls <- function() {
        shiny::downloadHandler(
          filename = function() { "results.xls" },
          content  = function(fname) write.table(results, fname, sep = "\t", col.names = FALSE)
        )
      }

      return_txt <- function() {
        shiny::downloadHandler(
          filename = function() { "check.txt" },
          content  = function(file) {
            sink(file)
            print(getwd(), quote = FALSE)
            if (input$donor_type == 1) {
              print("Donor Allele Matrix", quote = FALSE); print(dm)
            } else {
              print("Donor 1 Allele Matrix", quote = FALSE); print(d1m)
              print("Donor 2 Allele Matrix", quote = FALSE); print(d2m)
            }
            print("Recipient Allele Matrix", quote = FALSE); print(rm)
            print("Sample Allele Matrix",    quote = FALSE); print(sm)
            print("Final Results",           quote = FALSE); print(results)
            sink()
          }
        )
      }

      if (input$donor_type == 1) {
        output$results_xls    <- return_xls()
        output$check_file_txt <- return_txt()
      } else {
        output$dd_results_xls    <- return_xls()
        output$dd_check_file_txt <- return_txt()
      }

    })

    session$onSessionEnded(shiny::stopApp)
  }
}
