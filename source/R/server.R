build_server <- function() {
  function(input, output, session) {

    volumes <- c(
      Home = "~",
      DirectoryExample = system.file("extdata/Directory_example", package = "Rchimerism"),
      shinyFiles::getVolumes()()
    )
    shinyFiles::shinyDirChoose(input, 'directory_select', roots = volumes, session = session)
    dirname <- shiny::reactive({ shinyFiles::parseDirPath(volumes, input$directory_select) })

    gm_parsed <- shiny::reactive({
      if (input$directory_mode != 3) return(NULL)
      if (is.null(input$gm_file)) return(NULL)
      parse_genemapper_tsv(input$gm_file$datapath)
    })

    output$gm_role_ui <- shiny::renderUI({
      df <- gm_parsed()
      if (is.null(df)) return(NULL)

      if (!"Sample.File.Name" %in% colnames(df)) {
        return(shiny::p("Error: file does not contain a 'Sample File Name' column."))
      }

      sample_names <- list_genemapper_samples(df)
      if (length(sample_names) == 0) {
        return(shiny::p("No samples found in file."))
      }

      build_label <- function(name) {
        role <- detect_genemapper_role(name)
        if (!is.null(role)) paste0(name, "  [auto: ", role, "]") else name
      }
      choices     <- stats::setNames(sample_names, sapply(sample_names, build_label))
      all_choices <- c("--- select ---" = "", choices)

      if (input$donor_type == 1) {
        shiny::tagList(
          shiny::selectInput("gm_donor",     "Donor",
                             choices = all_choices,
                             selected = find_auto_role(sample_names, "donor")),
          shiny::selectInput("gm_recipient", "Recipient",
                             choices = all_choices,
                             selected = find_auto_role(sample_names, "recipient")),
          shiny::selectInput("gm_sample",    "Sample (chimera)",
                             choices = all_choices,
                             selected = find_auto_role(sample_names, "sample"))
        )
      } else {
        shiny::tagList(
          shiny::selectInput("gm_donor1",    "Donor 1",
                             choices = all_choices,
                             selected = find_auto_role(sample_names, "donor")),
          shiny::selectInput("gm_donor2",    "Donor 2",
                             choices = all_choices,
                             selected = ""),
          shiny::selectInput("gm_recipient", "Recipient",
                             choices = all_choices,
                             selected = find_auto_role(sample_names, "recipient")),
          shiny::selectInput("gm_sample",    "Sample (chimera)",
                             choices = all_choices,
                             selected = find_auto_role(sample_names, "sample"))
        )
      }
    })

    shiny::observeEvent(input$new_analysis_button, {
      session$reload()
    })

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
      } else if (input$directory_mode == 2) {
        dir_mode_vars <- initialize_dir_mode_vars()
        if (input$donor_type == 1) {
          ddata <- safe_read_delim(dir_mode_vars[["ddata"]], "donor data")
        } else {
          ddata  <- safe_read_delim(dir_mode_vars[["d1data"]], "donor 1 data")
          d2data <- safe_read_delim(dir_mode_vars[["d2data"]], "donor 2 data")
        }
        rdata <- safe_read_delim(dir_mode_vars[["rdata"]], "recipient data")
        sdata <- safe_read_delim(dir_mode_vars[["sdata"]], "sample data")
      } else {
        df <- gm_parsed()
        shiny::validate(shiny::need(!is.null(df), "Upload a GeneMapper export file first."))

        if (input$donor_type == 1) {
          shiny::validate(
            shiny::need(!is.null(input$gm_donor) && nchar(input$gm_donor) > 0,
                        "Select a donor sample."),
            shiny::need(!is.null(input$gm_recipient) && nchar(input$gm_recipient) > 0,
                        "Select a recipient sample."),
            shiny::need(!is.null(input$gm_sample) && nchar(input$gm_sample) > 0,
                        "Select a chimera sample.")
          )
          ddata <- extract_genemapper_sample(df, input$gm_donor)
          rdata <- extract_genemapper_sample(df, input$gm_recipient)
          sdata <- extract_genemapper_sample(df, input$gm_sample)
        } else {
          shiny::validate(
            shiny::need(!is.null(input$gm_donor1) && nchar(input$gm_donor1) > 0,
                        "Select donor 1."),
            shiny::need(!is.null(input$gm_donor2) && nchar(input$gm_donor2) > 0,
                        "Select donor 2."),
            shiny::need(!is.null(input$gm_recipient) && nchar(input$gm_recipient) > 0,
                        "Select a recipient sample."),
            shiny::need(!is.null(input$gm_sample) && nchar(input$gm_sample) > 0,
                        "Select a chimera sample.")
          )
          ddata  <- extract_genemapper_sample(df, input$gm_donor1)
          d2data <- extract_genemapper_sample(df, input$gm_donor2)
          rdata  <- extract_genemapper_sample(df, input$gm_recipient)
          sdata  <- extract_genemapper_sample(df, input$gm_sample)
        }
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
                               ru, rt, rnn, d1nn, d2nn, d1u, d2u, d1t, d2t, r,
                               ignore_unknown_alleles = input$ignore_unknown_alleles)
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

        chi_sd_output <- chiSD(sdata, markers, profile, rt, dt, d, r,
                               ignore_unknown_alleles = input$ignore_unknown_alleles)
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
            title = "Unclassified allele(s) in sample data",
            shiny::tags$p(
              "One or more alleles in the sample were found at positions where",
              shiny::tags$b("neither donor nor recipient"),
              "carries that allele (marked as",
              shiny::tags$b(as.character(UNKNOWN_ALLELE_SENTINEL)),
              "in the Sample Matrix)."
            ),
            shiny::tags$p(
              "Possible causes: stutter peaks, pull-up artifacts, or a genuine allele",
              "not present in the reference samples."
            ),
            shiny::tags$p(
              "The analysis has proceeded normally.",
              "Review the highlighted rows in the Sample Matrix.",
              paste0("To remove the filter, clear '", UNKNOWN_ALLELE_SENTINEL,
                     "' from the search box above the table.")
            )
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
            title = "Unclassified allele(s) in sample data",
            shiny::tags$p(
              "One or more alleles in the sample were found at positions where",
              shiny::tags$b("neither donor nor recipient"),
              "carries that allele (marked as",
              shiny::tags$b(as.character(UNKNOWN_ALLELE_SENTINEL)),
              "in the Sample Matrix)."
            ),
            shiny::tags$p(
              "Possible causes: stutter peaks, pull-up artifacts, or a genuine allele",
              "not present in the reference samples."
            ),
            shiny::tags$p(
              "The analysis has proceeded normally.",
              "Review the highlighted rows in the Sample Matrix.",
              paste0("To remove the filter, clear '", UNKNOWN_ALLELE_SENTINEL,
                     "' from the search box above the table.")
            )
          ))
        } else {
          output$s_dd_matrix <- matrix_output(sm, 1)
        }
      }

      sample_name <- if (input$directory_mode == 1) {
        tools::file_path_sans_ext(input$sdata$name)
      } else if (input$directory_mode == 2) {
        basename(dirname())
      } else {
        input$gm_sample
      }

      names_list <- if (input$directory_mode == 1) {
        if (input$donor_type == 1) {
          list(donor     = tools::file_path_sans_ext(input$ddata$name),
               recipient = tools::file_path_sans_ext(input$rdata$name),
               sample    = tools::file_path_sans_ext(input$sdata$name))
        } else {
          list(donor1    = tools::file_path_sans_ext(input$ddata$name),
               donor2    = tools::file_path_sans_ext(input$d2data$name),
               recipient = tools::file_path_sans_ext(input$rdata$name),
               sample    = tools::file_path_sans_ext(input$sdata$name))
        }
      } else if (input$directory_mode == 2) {
        if (input$donor_type == 1) {
          list(donor     = "ddata",
               recipient = "rdata",
               sample    = basename(dirname()))
        } else {
          list(donor1    = "d1data",
               donor2    = "d2data",
               recipient = "rdata",
               sample    = basename(dirname()))
        }
      } else {
        if (input$donor_type == 1) {
          list(donor     = input$gm_donor,
               recipient = input$gm_recipient,
               sample    = input$gm_sample)
        } else {
          list(donor1    = input$gm_donor1,
               donor2    = input$gm_donor2,
               recipient = input$gm_recipient,
               sample    = input$gm_sample)
        }
      }

      return_xls <- function() {
        shiny::downloadHandler(
          filename = function() { paste0(sample_name, "_results.xlsx") },
          content  = function(fname) {
            export <- as.data.frame(results)
            export$Marker <- rownames(results)
            export <- export[, c(ncol(export), 1:(ncol(export) - 1))]
            rownames(export) <- NULL

            base_rows <- nrow(export)

            if (input$donor_type == 1) {
              donor_mean    <- na.omit(results[, "Donor%_Mean"])[1] * 100
              recip_mean    <- na.omit(results[, "Recipient%_Mean"])[1] * 100
              donor_sd      <- na.omit(results[, "Donor%_SD"])[1] * 100
              donor_cv      <- na.omit(results[, "Donor%_CV"])[1] * 100
              n_informative <- sum(!is.na(results[, "Donor%_Mean"]))

              labels <- c("", "SUMMARY",
                          "Donor", "Recipient", "Sample (chimera)",
                          "Donor%_Mean (%)", "Recipient%_Mean (%)",
                          "Donor%_SD (%)", "Donor%_CV (%)", "Informative Loci (N)")
              export[base_rows + seq_along(labels), 1] <- labels
              export[base_rows + 3,  2] <- names_list$donor
              export[base_rows + 4,  2] <- names_list$recipient
              export[base_rows + 5,  2] <- names_list$sample
              export[base_rows + 6,  4] <- donor_mean
              export[base_rows + 7,  7] <- recip_mean
              export[base_rows + 8,  5] <- donor_sd
              export[base_rows + 9,  6] <- donor_cv
              export[base_rows + 10, 2] <- n_informative
            } else {
              data_only  <- results[-nrow(results), ]
              d1_mean    <- na.omit(data_only[, "Donor_1%_Mean"])[1] * 100
              d2_mean    <- na.omit(data_only[, "Donor_2%_Mean"])[1] * 100
              r_mean     <- na.omit(data_only[, "Recipient%_Mean"])[1] * 100
              d1_sd      <- na.omit(data_only[, "Donor_1%_SD"])[1] * 100
              d2_sd      <- na.omit(data_only[, "Donor_2%_SD"])[1] * 100
              r_sd       <- na.omit(data_only[, "Recipient%_SD"])[1] * 100
              d1_cv      <- na.omit(data_only[, "Donor_1%_CV"])[1] * 100
              d2_cv      <- na.omit(data_only[, "Donor_2%_CV"])[1] * 100
              r_cv       <- na.omit(data_only[, "Recipient%_CV"])[1] * 100
              n_inf_d1   <- as.integer(results["Info#", "Donor_1%_Mean"])
              n_inf_d2   <- as.integer(results["Info#", "Donor_2%_Mean"])
              n_inf_r    <- as.integer(results["Info#", "Recipient%_Mean"])

              labels <- c(
                "", "SUMMARY",
                "Donor 1", "Donor 2", "Recipient", "Sample (chimera)",
                "Donor_1%_Mean (%)", "Donor_2%_Mean (%)", "Recipient%_Mean (%)",
                "Donor_1%_SD (%)", "Donor_2%_SD (%)", "Recipient%_SD (%)",
                "Donor_1%_CV (%)", "Donor_2%_CV (%)", "Recipient%_CV (%)",
                "Inf. Loci D1 (N)", "Inf. Loci D2 (N)", "Inf. Loci R (N)"
              )
              export[base_rows + seq_along(labels), 1] <- labels
              export[base_rows + 3,  2] <- names_list$donor1
              export[base_rows + 4,  2] <- names_list$donor2
              export[base_rows + 5,  2] <- names_list$recipient
              export[base_rows + 6,  2] <- names_list$sample
              export[base_rows + 7,  3]  <- d1_mean
              export[base_rows + 8,  7]  <- d2_mean
              export[base_rows + 9,  11] <- r_mean
              export[base_rows + 10, 4]  <- d1_sd
              export[base_rows + 11, 8]  <- d2_sd
              export[base_rows + 12, 12] <- r_sd
              export[base_rows + 13, 5]  <- d1_cv
              export[base_rows + 14, 9]  <- d2_cv
              export[base_rows + 15, 13] <- r_cv
              export[base_rows + 16, 2]  <- n_inf_d1
              export[base_rows + 17, 6]  <- n_inf_d2
              export[base_rows + 18, 10] <- n_inf_r
            }

            writexl::write_xlsx(export, fname)
          }
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

  }
}
