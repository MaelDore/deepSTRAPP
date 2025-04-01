
#' @title Plot histogram of STRAPP test statistics to assess results
#'
#' @description Plot an histogram of the distribution of the test statistics
#'   obtained from a STRAPP test carried out for a unique `focal_time`.
#'   (See [deepSTRAPP::compute_STRAPP_test_for_focal_time()]).
#'
#'   Returns a single histogram for overall tests.
#'   If `plot_posthoc_tests = TRUE`, it will return a faceted plot with an histogram per post hoc tests.
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plot will be saved directly in a PDF file.
#'
#' @param STRAPP_results List of elements generated with [deepSTRAPP::compute_STRAPP_test_for_focal_time()],
#'   that summarize the results of a STRAPP test for a specific time in the past (i.e. the `focal_time`).
#' @param display_plot Logical. Whether to display the histogram(s) generated in the R console. Default is `TRUE`.
#' @param plot_posthoc_tests Logical. For multinominal data only. Whether to plot the histogram for the overall Kruskal-Wallis test across all states (`plot_posthoc_tests = FALSE`),
#'   or plot the histograms for all the pairwise post hoc Dunn's test across pairs of states (`plot_posthoc_tests = TRUE`). Default is `FALSE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with '.pdf'.
#'
#' @export
#' @importFrom BAMMtools plot.bammdata
#' @importFrom ggplot2 ggplot geom_histogram aes geom_vline labs ggtitle theme element_line element_rect element_text unit margin annotate annotation_custom
#' @importFrom grid gpar textGrob
#' @importFrom cowplot plot_grid save_plot
#'
#' @details Histograms are build based on the distribution of the test statistics.
#'   Such distributions are recorded in the outputs of STRAPP tests carried out with [deepSTRAPP::compute_STRAPP_test_for_focal_time()]
#'   when `return_perm_data = TRUE`. The `STRAPP_results` object provided as input is a list that must contain
#'   a `$perm_data_df` element that summarizes test statistics computed across posterior samples.
#'
#'   For multinominal data (categorical or biogeographic data with more than 2 states), it is possible to plot the histograms of post hoc pairwise tests.
#'   Set `plot_posthoc_tests = TRUE` to generate histograms for all the pairwise post hoc Dunn's test across pairs of states.
#'   To achieve this, the `STRAPP_results` input object must contain a `$posthoc_pairwise_tests$perm_data_array` element that summarizes test statistics
#'   computed across posterior samples for all pairwise post hoc tests. This is obtained from [deepSTRAPP::compute_STRAPP_test_for_focal_time()] when setting both
#'   `posthoc_pairwise_tests = TRUE` to carry out post hoc tests, and `return_perm_data = TRUE` to record distributions of test statistics.
#'
#' @return By default, the function returns a list of classes `gg` and `ggplot`.
#'   This object is a ggplot that can be displayed on the console with `print(output)`.
#'   It corresponds to the histogram being displayed on the console when the function is run, if `display_plot = TRUE`, and can be further
#'   modify for aesthetics using the ggplot2 grammar.
#'
#'   If using multinominal data and set `plot_posthoc_tests = TRUE`, the function will return a list of objects.
#'   Each object is the ggplot associated with a pairwise post hoc test.
#'   To plot each histogram i individually, use `print(output_list[[i]])`.
#'   To plot all histograms at once in a multifaceted plot, as displayed on the console if `display_plot = TRUE`, use `cowplot::plot_grid(plotlist = output_list)`.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plot. For post hoc tests, this will save the multifaceted plot.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::compute_STRAPP_test_for_focal_time()] [deepSTRAPP::run_STRAPP_test_for_focal_time()]
#'
#' @examples
#' # ------ Prepare data ------ #
#'
#' ## Load the BAMM_object summarizing 1000 posterior samples of BAMM with diversification rates
#' # for ponerine ants extracted for 10My ago.
#' data(Ponerinae_BAMM_object_10My, package = "deepSTRAPP")
#'
#' # Plot the associated phylogeny with mapped rates
#' BAMMtools::plot.bammdata(Ponerinae_BAMM_object_10My)
#'
#' ## Load the object containing head width trait data for ponerine ants extracted for 10My ago.
#' data(Ponerinae_trait_data_10My, package = "deepSTRAPP")
#'
#' # Categorize continuous trait data into three states to create multinomial data
#' trait_data_continuous <- Ponerinae_trait_data_10My
#' trait_data_multinominal <- trait_data_continuous
#' trait_data_multinominal$trait_data[trait_data_continuous$trait_data < 0] <- "state_B"
#' trait_data_multinominal$trait_data[trait_data_continuous$trait_data < -1] <- "state_A"
#' trait_data_multinominal$trait_data[trait_data_continuous$trait_data >= 0] <- "state_C"
#' trait_data_multinominal$trait_data_type <- "categorical"
#'
#' table(trait_data_multinominal$trait_data)
#'
#' \dontrun{  (May take several minutes to run)
#' # ------ Compute STRAPP test ------ #
#'
#' # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
#' # between "net_diversification" rates and trait data
#' STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_multinominal,
#'    posthoc_pairwise_tests = TRUE,
#'    two_tailed = TRUE,
#'    return_perm_data = TRUE)
#' str(STRAPP_results, max.level = 2)
#' # Data from the posterior samples for the overall Kruskal-Wallis test is available
#' # in STRAPP_results$perm_data_df
#' head(STRAPP_results$perm_data_df)
#' # Data from the posterior samples for the post hoc Dunn's tests is available
#' # in STRAPP_results$posthoc_pairwise_tests$perm_data_array
#' head(STRAPP_results$posthoc_pairwise_tests$perm_data_array[1,,])
#'
#' # ------ Plot histogram of STRAPP overall test results ------ #
#'
#' histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
#'                         STRAPP_results = STRAPP_results,
#'                         display_plot = TRUE,
#'                         # PDF_file_path = "./plot_STRAPP_histogram_overall_test.pdf",
#'                         plot_posthoc_tests = FALSE)
#' # Adjust aesthetics a posteriori
#' histogram_ggplot_adj <- histogram_ggplot +
#'    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#' print(histogram_ggplot_adj)
#'
#' # ------ Plot histograms of STRAPP post hoc test results ------ #
#'
#' histograms_ggplot_list <- plot_histogram_STRAPP_test_for_focal_time(
#'                               STRAPP_results = STRAPP_results,
#'                               display_plot = TRUE,
#'                               # PDF_file_path = "./plot_STRAPP_histograms_posthoc_tests.pdf",
#'                               plot_posthoc_tests = TRUE)
#' # Plot all histograms one by one
#' print(histograms_ggplot_list)
#' # Plot all histograms on one faceted plot
#' cowplot::plot_grid(plotlist = histograms_ggplot_list)}
#'


plot_histogram_STRAPP_test_for_focal_time <- function (STRAPP_results,
                                                       display_plot = TRUE,
                                                       plot_posthoc_tests = FALSE,
                                                       PDF_file_path = NULL)

{
  ### Check input validity
  {
    ## STRAPP_results
    # STRAPP_results must have all the needed elements: $focal_time, $estimate, $p_value, $method, $trait_data_type_for_stats, $perm_data_df
    if (is.null(STRAPP_results$focal_time))
    {
      stop(paste0("'$focal_time' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                  "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects."))
    }
    if (is.null(STRAPP_results$estimate))
    {
      stop(paste0("'$estimate' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                  "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects."))
    }
    if (is.null(STRAPP_results$p_value))
    {
      stop(paste0("'$p_value' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                  "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects."))
    }
    if (is.null(STRAPP_results$method))
    {
      stop(paste0("'$method' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                  "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects."))
    }
    if (is.null(STRAPP_results$trait_data_type_for_stats))
    {
      stop(paste0("'$trait_data_type_for_stats' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                  "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects."))
    }
    if (is.null(STRAPP_results$perm_data_df))
    {
      stop(paste0("'$perm_data_df' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                  "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects.\n",
                  "Especially, check if you used 'return_perm_data = TRUE' to save the permutated data needed for the histogram plot."))
    }

    ## plot_posthoc_tests
    if (plot_posthoc_tests)
    {
      # plot_posthoc_tests = TRUE only for "multinominal" data
      if (STRAPP_results$trait_data_type_for_stats != "multinominal")
      {
        stop(paste0("'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                    "Set 'posthoc_pairwise_tests = FALSE', or provide 'STRAPP_results' for a trait with more than two states/ranges'.\n",
                    "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects."))
      }
      # Check if STRAPP_results$posthoc_pairwise_tests$summary_df is present.
      if (is.null(STRAPP_results$posthoc_pairwise_tests$summary_df))
      {
        stop(paste0("'$posthoc_pairwise_tests$summary_df' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                    "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects.\n",
                    "Especially, check if you used 'posthoc_pairwise_tests = TRUE' to compute pairwise tests.\n",
                    "In addition, you must also select 'return_perm_data = TRUE' to save the permutated data needed for the histogram plots."))
      }
      # Check if $posthoc_pairwise_tests$perm_data_array is present.
      if (is.null(STRAPP_results$posthoc_pairwise_tests$perm_data_array))
      {
        stop(paste0("'$posthoc_pairwise_tests$perm_data_array' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
                    "See ?deepSTRAPP::run_STRAPP_test_for_focal_time() to learn how to generate those objects.\n",
                    "Especially, check if you used 'posthoc_pairwise_tests = TRUE' AND 'return_perm_data = TRUE' to compute pairwise tests ",
                    "and save the permutated data needed for the histogram plots."))
      }
    }

    ## PDF_file_path
    # If provided, PDF_file_path must end with ".pdf"
    if (!is.null(PDF_file_path))
    {
      if (length(grep(pattern = "\\.pdf$", x = PDF_file_path)) != 1)
      {
        stop("'PDF_file_path' must end with '.pdf'")
      }
    }
  }

  if (!plot_posthoc_tests)
  {
    ## Case for overall test plot

    # Extract name of the statistic
    stat_name <- names(STRAPP_results$perm_data_df)[4]
    # Extract null distribution of the statistic
    stat_null_distri <- STRAPP_results$perm_data_df[,stat_name]
    # Extract quantile of the critical threshold
    estimate_quantile <- names(STRAPP_results$estimate)

    # Build ggplot object
    ggplot_histo <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

      # Generate histogram
      ggplot2::geom_histogram(mapping = ggplot2::aes(x = stat_null_distri),
                              bins = 30, # Set to avoid displaying automatic message from ggplot
                              fill = "grey80", color = "grey50", alpha = 0.8) +

      # Add vline for x = 0 (evaluation threshold)
      ggplot2::geom_vline(xintercept = 0, color = "black",
                          linetype = "dashed", linewidth = 1.0) +

      # Add vline for x = estimate (significance threshold)
      ggplot2::geom_vline(xintercept = STRAPP_results$estimate, color = "red",
                          linetype = "dashed", linewidth = 1.0) +

      # Add test summary
      # alpha value : Estimate,  p-value
      annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                   label = paste0("Q", estimate_quantile, " = ", round(STRAPP_results$estimate, digits = 3), "\n",
                                  "P-value = ", round(STRAPP_results$p_value, digits = 3))) +

      # Adjust axis labels
      ggplot2::labs(x = paste0("Distribution of the test statistics"),
                    y = "Counts across posterior samples") +

      # Add title
      ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results$method," test\n",
                              "Focal-time = ", STRAPP_results$focal_time)) +

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

    ## Display plot if requested
    if (display_plot)
    {
      print(ggplot_histo)
    }

    ## Export plot if requested
    if (!is.null(PDF_file_path))
    {
      cowplot::save_plot(plot = ggplot_histo,
                         filename = PDF_file_path,
                         base_height = 8, base_width = 10)
    }

    ## Return ggplot
    return(invisible(ggplot_histo))

   } else {

    ## Case for post hoc tests plot

    # Create list of pairwise plots
    ggplots_histo_list <- list()

    # Extract pairs
    all_pairs <- dimnames(STRAPP_results$posthoc_pairwise_tests$perm_data_array)$pairs
    # Compute size factor to adjust size of everything
    size_factor <- sqrt(length(all_pairs))

    for (i in seq_along(all_pairs))
    {
      # i <- 1

      # Extract perm data for pair i
      perm_data_df_i <- as.data.frame(STRAPP_results$posthoc_pairwise_tests$perm_data_array[i, , ])

      # Extract pair names
      pair_i <- all_pairs[i]
      # # Replace "!=" with unicode
      # pair_i <- gsub(pattern = "!=", replacement = "\u2260", x = pair_i )

      # Extract name of the statistic
      stat_name <- names(perm_data_df_i)[3]
      # Extract null distribution of the statistic
      stat_null_distri <- perm_data_df_i[,stat_name]
      # Extract quantile of the critical threshold (same the one used for the main test)
      estimate_quantile <- names(STRAPP_results$estimate)
      # Extract estimate of the critical threshold
      estimate_value <- STRAPP_results$posthoc_pairwise_tests$summary_df$estimates[i]
      # Extract p-value
      p_value <- STRAPP_results$posthoc_pairwise_tests$summary_df$p_values[i]

      # Detect if adjusted p-values differ from p-value
      p_value_adj_to_plot <- any(STRAPP_results$posthoc_pairwise_tests$summary_df$p_values != STRAPP_results$posthoc_pairwise_tests$summary_df$p_values_adjusted)
      # Extract p-value adjusted if different
      if (p_value_adj_to_plot)
      {
        p_value_adj <- STRAPP_results$posthoc_pairwise_tests$summary_df$p_values_adjusted[i]
      }

      # Build ggplot object
      ggplot_histo_i <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

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
        ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results$posthoc_pairwise_tests$method," test\n",
                                "Focal-time = ", STRAPP_results$focal_time, "\n",
                                "Hypothesis = ",pair_i, "\n")) +

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
        ggplot_histo_i <- ggplot_histo_i +
          # Add test summary including adjusted p-value
          annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18/size_factor),
                       label = paste0("Q", estimate_quantile, " = ", round(estimate_value, digits = 3), "\n",
                                      "P-value = ", round(p_value, digits = 3), "\n",
                                      "P-value adj = ", round(p_value_adj, digits = 3)))
      }

      # Store ggplot for pairs i
      ggplots_histo_list[[i]] <- ggplot_histo_i

    }

    ## Display plot if requested
    if (display_plot)
    {
      cow_plot_to_print <- cowplot::plot_grid(plotlist = ggplots_histo_list)
      print(cow_plot_to_print)
    }

    ## Export plot if requested
    if (!is.null(PDF_file_path))
    {
      cowplot::save_plot(plot = cowplot::plot_grid(plotlist = ggplots_histo_list),
                         filename = PDF_file_path,
                         base_height = 8, base_width = 10)
    }

    ## Return list of ggplots
    return(invisible(ggplots_histo_list))
  }
}



### Helper function to enable the use of "npc" units in ggplot2::annotate()

#' @noRd

annotate_npc <- function(label, x, y, ...)
{
  ggplot2::annotation_custom(
    grob = grid::textGrob(x = ggplot2::unit(x, "npc"),
                          y = ggplot2::unit(y, "npc"),
                          label = label, ...))
}
