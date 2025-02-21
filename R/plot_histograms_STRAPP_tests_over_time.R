
#' @title Wrapper function to plot histogram of STRAPP test statistics over time-steps
#'
#' @description Plot an histogram of the distribution of the test statistics
#'   obtained from a STRAPP test carried out for each focal time in `$time_steps`.
#'   (See [deepSTRAPP::run_STRAPP_tests_over_time()]).
#'
#'   Returns one histogram for overall tests for each focal time in `$time_steps`.
#'   If `plot_posthoc_tests = TRUE`, it will return one faceted plot with an histogram
#'   per post hoc tests for each focal time in `$time_steps`.
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plots will be saved directly in a PDF file,
#'   with one page per focal time in `$time_steps`.
#'
#' @param STRAPP_tests_over_time List of elements generated with [deepSTRAPP::run_STRAPP_tests_over_time()],
#'   that summarize the results of multiple STRAPP tests across `$time_steps`. It needs to include the `$STRAPP_results_over_time`
#'   element with `$perm_data_df` obtained when setting both `return_STRAPP_results = TRUE` and `return_perm_data = TRUE`.
#' @param display_plots Logical. Whether to display the histograms generated in the R console. Default is `TRUE`.
#' @param plot_posthoc_tests Logical. For multinominal data only. Whether to plot the histograms for the overall Kruskal-Wallis test across all states (`plot_posthoc_tests = FALSE`),
#'   or plot the histograms for all the pairwise post hoc Dunn's tests across pairs of states (`plot_posthoc_tests = TRUE`). Default is `FALSE`.
#' @param PDF_file_path Character string. If provided, the plots will be saved in a unique PDF file following the path provided here. The path must end with '.pdf'.
#'   Each page of the PDF corresponds to a focal time in `$time_steps`.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes geom_vline labs ggtitle theme element_line element_rect element_text unit margin annotate annotation_custom
#' @importFrom grid gpar
#' @importFrom cowplot plot_grid save_plot
#' @importFrom grDevices pdf dev.off
#'
#' @details Histograms are build based on the distribution of the test statistics.
#'   Such distributions are recorded in the outputs of STRAPP tests carried out with [deepSTRAPP::run_STRAPP_tests_over_time()]
#'   when `return_STRAPP_results = TRUE` AND `return_perm_data = TRUE`. The `$STRAPP_results_over_time` objects provided within the input are lists that must contain
#'   a `$perm_data_df` element that summarizes test statistics computed across posterior samples.
#'
#'   For multinominal data (categorical or biogeographic data with more than 2 states), it is possible to plot the histograms of post hoc pairwise tests.
#'   Set `plot_posthoc_tests = TRUE` to generate histograms for all the pairwise post hoc Dunn's test across pairs of states.
#'   To achieve this, the `$STRAPP_results_over_time` objects must contain a `$posthoc_pairwise_tests$perm_data_array` element that summarizes test statistics
#'   computed across posterior samples for all pairwise post hoc tests. This is obtained from [deepSTRAPP::run_STRAPP_tests_over_time()] when setting
#'   `return_STRAPP_results = TRUE` to return the STRAPP results, `posthoc_pairwise_tests = TRUE` to carry out post hoc tests,
#'   and `return_perm_data = TRUE` to record distributions of test statistics.
#'
#' @return By default, the function returns a list of sub-lists of classes `gg` and `ggplot` ordered as in `$time_steps`.
#'   Each sub-list corresponds to a ggplot for a given `focal_time` i that can be displayed on the console with `print(output[[i]])`.
#'   If `display_plots = TRUE`, the histograms are being displayed on the console one by one while generated.
#'
#'   If using multinominal data and set `plot_posthoc_tests = TRUE`, the function will return a list of sub-lists of objects ordered as in `$time_steps`.
#'   Each sub-list is a list of the ggplots associated with pairwise post hoc tests carried out for this a given `focal_time`.
#'   For a given `focal_time` i, to plot each histogram j individually, use `print(output_list[[i]][[j]])`.
#'   To plot all histograms of a given `focal_time` i at once in a multifaceted plot, as displayed sequentially on the console if `display_plots = TRUE`,
#'   use `cowplot::plot_grid(plotlist = output_list[[i]])`.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plots with one page per `$time_steps`.
#'   For post hoc tests, this will save the multifaceted plots.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::run_STRAPP_tests_over_time()] [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()]
#'
#' @examples
#' ## Load results of run_STRAPP_tests_over_time()
#' data(STRAPP_tests_over_time_temp_example, package = "deepSTRAPP")
#'
#' ## Plot histograms of STRAPP overall test results
#' histogram_ggplots <- plot_histograms_STRAPP_tests_over_time(
#'   STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example,
#'   display_plots = TRUE,
#'   # PDF_file_path = "./plot_STRAPP_histograms_overall_tests_over_time.pdf",
#'   plot_posthoc_tests = FALSE)
#' # Print histogram for time-step 1
#' print(histogram_ggplots[[1]])
#' # Adjust aesthetics of plot for time-step 1 a posteriori
#' histogram_ggplot_adj <- histogram_ggplots[[1]] +
#'     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#' print(histogram_ggplot_adj)
#'
#' ## Plot histograms of STRAPP post hoc test results
#' histograms_ggplots_list <- plot_histograms_STRAPP_tests_over_time(
#'   STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example,
#'   display_plots = TRUE,
#'   # PDF_file_path = "./test_multi_histo.pdf"
#'   plot_posthoc_tests = TRUE)
#' # Print all histograms for time-step 1 one by one
#' print(histograms_ggplots_list[[1]])
#' # Plot all histograms on one faceted plot
#' cowplot::plot_grid(plotlist = histograms_ggplots_list[[1]])
#'


plot_histograms_STRAPP_tests_over_time <- function (STRAPP_tests_over_time,
                                                    display_plots = TRUE,
                                                    plot_posthoc_tests = FALSE,
                                                    PDF_file_path = NULL)

{
  ### Check input validity

  # STRAPP_tests_over_time must have $STRAPP_results_over_time
  # Make a special warning if $STRAPP_results_over_time is missing to ask to select return_STRAPP_results = TRUE in compute_STRAPP_test_over_time() to save the STRAPP_results needed for the histogram plot

  # STRAPP_tests_over_time must have $STRAPP_results_over_time with the $perm_data_df
  # Make a special warning if $perm_data_df is missing to ask to select return_STRAPP_results = TRUE AND return_perm_data = TRUE in compute_STRAPP_test_over_time() to save the raw data needed for the histogram plot

  # plot_posthoc_tests = TRUE only for "multinominal" data
  # Check if $posthoc_pairwise_tests$perm_data_array is present in $STRAPP_results_over_time.
  # Make a special warning if $posthoc_pairwise_tests$perm_data_array is missing to select return_STRAPP_results = TRUE, posthoc_pairwise_tests = TRUE AND return_perm_data = TRUE in compute_STRAPP_test_over_time() to save the raw data needed for the histogram plot of post hoc tests.

  # If provided, PDF_file_path must end with ".pdf"

  if (!plot_posthoc_tests)
  {
    ## Case for overall test plot

    # Initiate list of ggplots
    ggplot_histo_list <- list()

    # Loop per time-steps
    for (i in seq_along(STRAPP_tests_over_time$time_steps))
    {
      # i <- 1

      # Extract STRAPP_results
      STRAPP_results_i <- STRAPP_tests_over_time$STRAPP_results_over_time[[i]]

      # Extract name of the statistic
      stat_name <- names(STRAPP_results_i$perm_data_df)[4]
      # Extract null distribution of the statistic
      stat_null_distri <- STRAPP_results_i$perm_data_df[,stat_name]
      # Extract quantile of the critical threshold
      estimate_quantile <- names(STRAPP_results_i$estimate)

      # Build ggplot object
      ggplot_histo_i <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

        # Generate histogram
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = stat_null_distri),
                                bins = 30, # Set to avoid displaying automatic message from ggplot
                                fill = "grey80", color = "grey50", alpha = 0.8) +

        # Add vline for x = 0 (evaluation threshold)
        ggplot2::geom_vline(xintercept = 0, color = "black",
                            linetype = "dashed", linewidth = 1.0) +

        # Add vline for x = estimate (significance threshold)
        ggplot2::geom_vline(xintercept = STRAPP_results_i$estimate, color = "red",
                            linetype = "dashed", linewidth = 1.0) +

        # Add test summary
        # alpha value : Estimate,  p-value
        annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                     label = paste0("Q", estimate_quantile, " = ", round(STRAPP_results_i$estimate, digits = 3), "\n",
                                    "P-value = ", round(STRAPP_results_i$p_value, digits = 3))) +

        # Adjust axis labels
        ggplot2::labs(x = paste0("Distribution of the test statistics"),
                      y = "Counts across posterior samples") +

        # Add title
        ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results_i$method," test\n",
                                "Focal-time = ", STRAPP_results_i$focal_time)) +

        # Adjust aesthetics
        ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
                       panel.background = ggplot2::element_rect(fill = NA, color = NA),
                       plot.margin = ggplot2::unit(c(0.3, 0.5, 0.5, 0.5), "inches"),
                       plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                                          margin = ggplot2::margin(b = 15, t = 5)),
                       axis.title = ggplot2::element_text(size = 20, color = "black"),
                       axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                       axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
                       axis.line = ggplot2::element_line(linewidth = 1.5),
                       axis.text = ggplot2::element_text(size = 18, color = "black"),
                       axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

      ## Display plot at each time-step if requested
      if (display_plots)
      {
        print(ggplot_histo_i)
      }

      ## Store in list of ggplots
      ggplot_histo_list[[i]] <- ggplot_histo_i
    }

    ## Export all plots in one PDF if requested
    if (!is.null(PDF_file_path))
    {
      grDevices::pdf(file = PDF_file_path, height = 8, width = 10, onefile = TRUE)
      for (i in seq_along(ggplot_histo_list))
      {
        print(ggplot_histo_list[[i]])
      }
      grDevices::dev.off()
    }

    ## Return list of ggplots
    return(invisible(ggplot_histo_list))

  } else {

    ## Case for post hoc tests plot

    # Initiate list of lists of ggplots
    ggplot_histo_list_of_lists <- list()

    # Loop per time-steps
    for (i in seq_along(STRAPP_tests_over_time$time_steps))
    {
      # i <- 1

      # Extract STRAPP_results
      STRAPP_results_i <- STRAPP_tests_over_time$STRAPP_results_over_time[[i]]

      # Create list of pairwise plots for time-step i
      ggplots_histo_list_i <- list()

      # Extract pairs
      all_pairs <- dimnames(STRAPP_results_i$posthoc_pairwise_tests$perm_data_array)$pairs
      # Compute size factor to adjust size of everything
      size_factor <- sqrt(length(all_pairs))

      for (j in seq_along(all_pairs))
      {
        # j <- 1

        # Extract perm data for pair j
        perm_data_df_j <- as.data.frame(STRAPP_results_i$posthoc_pairwise_tests$perm_data_array[j, , ])

        # Extract pair names
        pair_j <- all_pairs[j]
        # # Replace "!=" with unicode
        # pair_j <- gsub(pattern = "!=", replacement = "\u2260", x = pair_j )

        # Extract name of the statistic
        stat_name <- names(perm_data_df_j)[3]
        # Extract null distribution of the statistic
        stat_null_distri <- perm_data_df_j[,stat_name]
        # Extract quantile of the critical threshold (same the one used for the main test)
        estimate_quantile <- names(STRAPP_results_i$estimate)
        # Extract estimate of the critical threshold
        estimate_value <- STRAPP_results_i$posthoc_pairwise_tests$summary_df$estimates[j]
        # Extract p-value
        p_value <- STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values[j]

        # Detect if adjusted p-values differ from p-value
        p_value_adj_to_plot <- any(STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values != STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values_adjusted)
        # Extract p-value adjusted if different
        if (p_value_adj_to_plot)
        {
          p_value_adj <- STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values_adjusted[j]
        }

        # Build ggplot object for pair j
        ggplot_histo_j <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

          # Generate histogram
          ggplot2::geom_histogram(mapping = ggplot2::aes(x = stat_null_distri),
                                  bins = 30, # Set to avoid displaying automatic message from ggplot
                                  fill = "grey80", color = "grey50", alpha = 0.8) +

          # Add vline for x = 0 (evaluation threshold)
          ggplot2::geom_vline(xintercept = 0, color = "black",
                              linetype = "dashed", linewidth = 1.0/size_factor) +

          # Add vline for x = estimate (significance threshold)
          ggplot2::geom_vline(xintercept = estimate_value, color = "red",
                              linetype = "dashed", linewidth = 1.0/size_factor) +

          # Add test summary
          # alpha value : Estimate,  p-value
          annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18/size_factor),
                       label = paste0("Q", estimate_quantile, " = ", round(estimate_value, digits = 3), "\n",
                                      "P-value = ", round(p_value, digits = 3))) +

          # Adjust axis labels
          ggplot2::labs(x = paste0("Distribution of the test statistics"),
                        y = "Counts across posterior samples") +

          # Add title
          ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results_i$posthoc_pairwise_tests$method," test\n",
                                  "Focal-time = ", STRAPP_results_i$focal_time, "\n",
                                  "Hypothesis = ",pair_j, "\n")) +

          # Adjust aesthetics
          ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3/size_factor),
                         panel.background = ggplot2::element_rect(fill = NA, color = NA),
                         plot.margin = ggplot2::unit(c(0.3/(size_factor*2), 0.5/(size_factor*2), 0.5/(size_factor*2), 0.5/(size_factor*2)), "inches"),
                         plot.title = ggplot2::element_text(size = 20/size_factor, hjust = 0.5, color = "black",
                                                            margin = ggplot2::margin(b = 15/size_factor, t = 5/size_factor)),
                         axis.title = ggplot2::element_text(size = 20/size_factor, color = "black"),
                         axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10/size_factor)),
                         axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12/size_factor)),
                         axis.line = ggplot2::element_line(linewidth = 1.5/size_factor),
                         axis.text = ggplot2::element_text(size = 18/size_factor, color = "black"),
                         axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5/size_factor)),
                         axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5/size_factor)))

        # Add adjusted p-value if present
        if (p_value_adj_to_plot)
        {
          ggplot_histo_j <- ggplot_histo_j +
            # Add test summary including adjusted p-value
            annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18/size_factor),
                         label = paste0("Q", estimate_quantile, " = ", round(estimate_value, digits = 3), "\n",
                                        "P-value = ", round(p_value, digits = 3), "\n",
                                        "P-value adj = ", round(p_value_adj, digits = 3)))
        }

        # Store ggplot for pairs j
        ggplots_histo_list_i[[j]] <- ggplot_histo_j

      }

      ## Display plot at each time-step if requested
      if (display_plots)
      {
        cow_plot_to_print_i <- cowplot::plot_grid(plotlist = ggplots_histo_list_i)
        print(cow_plot_to_print_i)
      }

      ## Store in list of lists of ggplots
      ggplot_histo_list_of_lists[[i]] <- ggplots_histo_list_i
    }

    ## Export all plots in one PDF if requested
    if (!is.null(PDF_file_path))
    {
      # # Find the number of rows and columns
      # nb_plots <- length(cow_plot_to_print_i$layers)
      # nb_cols <- ceiling(sqrt(nb_plots))
      # nb_rows <- ceiling(nb_plots/nb_cols)

      # Export PDF witn one page per plot
      grDevices::pdf(file = PDF_file_path, height = 8, width = 10, onefile = TRUE)
      # pdf(file = PDF_file_path, height = 8*nb_rows, width = 10*nb_cols, onefile = TRUE)
      for (i in seq_along(ggplot_histo_list_of_lists))
      {
        cow_plot_to_print_i <- cowplot::plot_grid(plotlist = ggplot_histo_list_of_lists[[i]])
        print(cow_plot_to_print_i)
      }
      grDevices::dev.off()
    }

    ## Return list of lists of ggplots
    return(invisible(ggplot_histo_list_of_lists))
  }
}

