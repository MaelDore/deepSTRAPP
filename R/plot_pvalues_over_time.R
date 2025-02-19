
#' @title Plot evolution of p-values of STRAPP tests over time
#'
#' @description Plot the evolution of the p-values yield STRAPP tests
#'   carried out for across multiple `time_steps`, obtained from
#'   [deepSTRAPP::run_STRAPP_test_over_time()].
#'
#'   By default, return a plot with a single line for p-values of overall tests.
#'   If `plot_posthoc_tests = TRUE`, it will return a plot with multiple lines, one per pair in post hoc tests
#'   (only for multinominal data, with more than two states).
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plot will be saved directly in a PDF file.
#'
#' @param STRAPP_tests_over_time List of elements generated with [deepSTRAPP::run_STRAPP_test_over_time()],
#'   that summarize the results of multiple STRAPP tests across `$time_steps`.
#' @param time_range Vector of two numerical values. Time boundaries used for the plot.
#'   If `NULL` (the default), the range of data provided in `STRAPP_tests_over_time` will be used.
#' @param alpha Numerical. Significance level to display as a red dashed line on the plot. If set to `NULL`, no line will be added. Default is `0.05`.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param plot_posthoc_tests Logical. For multinominal data only. Whether to plot the p-values for the overall Kruskal-Wallis test across all states (`plot_posthoc_tests = FALSE`),
#'   or plot the p-values for the pairwise post hoc Dunn's test across pairs of states (`plot_posthoc_tests = TRUE`). Default is `FALSE`.
#'   This is only possible if `STRAPP_tests_over_time` contains the `$pvalues_summary_df_for_posthoc_pairwise_tests` element returned by
#'   [deepSTRAPP::run_STRAPP_test_over_time()] when `posthoc_pairwise_tests = TRUE`.
#' @param select_posthoc_pairs Vector of character strings used to specify the pairs to include in the plot. Names of pairs must match the pairs found in
#'   `STRAPP_tests_over_time$pvalues_summary_df_for_posthoc_pairwise_tests$pair`. Default is "all" to include all pairs.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with '.pdf'.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_line aes geom_hline scale_y_continuous scale_x_continuous scale_color_discrete xlab ylab ggtitle theme element_line element_rect element_text unit margin
#' @importFrom cowplot save_plot
#'
#' @details Plots are build based on the p-values recorded in summary_df provided by [deepSTRAPP::run_STRAPP_test_over_time()].
#'
#'   For overall tests, those p-values are found in `$pvalues_summary_df`.
#'
#'   For multinominal data (categorical or biogeographic data with more than 2 states), it is possible to plot p-values of post hoc pairwise tests.
#'   Set `plot_posthoc_tests = TRUE` to generate plots for the pairwise post hoc Dunn's test across pairs of states.
#'   To achieve this, the `STRAPP_tests_over_time` input object must contain a `$pvalues_summary_df_for_posthoc_pairwise_tests` element that summarizes p-values
#'   computed across pairs of states for all post hoc tests. This is obtained from [deepSTRAPP::run_STRAPP_test_over_time()] when setting
#'   `posthoc_pairwise_tests = TRUE` to carry out post hoc tests.
#'
#' @return The function returns a list of classes `gg` and `ggplot`.
#'   This object is a ggplot that can be displayed on the console with `print(output)`.
#'   It corresponds to the plot being displayed on the console when the function is run, if `display_plot = TRUE`,
#'   and can be further modify for aesthetics using the ggplot2 grammar.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plot.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::run_STRAPP_test_over_time()]
#'
#' @examples
#' ## Load results of run_STRAPP_test_over_time()
#' data(STRAPP_tests_over_time_temp_example, package = "deepSTRAPP")
#'
#' ## Plot results of overall Kruskal-Wallis test across all tests
#' plot_STRAPP_pvalues_over_time(
#'    STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example,
#'    alpha = 0.1,
#'    time_range = c(20, 150))
#'
#' ## Plot results of post hoc pairwise Dunn's tests between pairs of tests
#' plot_STRAPP_pvalues_over_time(
#'    STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example,
#'    plot_posthoc_tests = TRUE,
#'    # PDF_file_path = "./pvalues_over_time.pdf",
#'    select_posthoc_pairs = c("state_A != state_B",
#'                             "state_A != state_C"))
#'


plot_STRAPP_pvalues_over_time <-  function (
    STRAPP_tests_over_time,
    time_range = NULL,
    alpha = 0.05,
    display_plot = TRUE,
    plot_posthoc_tests = FALSE,
    select_posthoc_pairs = "all",
    PDF_file_path = NULL
)
{
  ### Check input validity

  # STRAPP_tests_over_time must have $pvalues_summary_df, with $focal_time and $p_value columns

  # If provided, PDF_file_path must end with ".pdf"

  # plot_posthoc_tests = TRUE only for "multinominal" data
  # Check if $pvalues_summary_df_for_posthoc_pairwise_tests is present, with $focal_time, $pair, and $p_value columns
  # Make a special warning to select posthoc_pairwise_tests = TRUE in run_STRAPP_test_over_time().

  # Check that select_posthoc_pairs is "all" or one of the pair in STRAPP_tests_over_time$pvalues_summary_df_for_posthoc_pairwise_tests$pair

  # Check that time_range is strictly positive, ordered in increasing age, and encompass multiple data points (per pair, if post hoc)

  # Check that alpha is NULL or a numerical between 0 and 1.

  ## Remove "_" from $rate_type
  rate_type <- gsub(pattern = "_", replacement = " ", x = STRAPP_tests_over_time$rate_type)

  if (!plot_posthoc_tests) ## Case for overall test plot
  {
    # Extract summary df
    pvalues_summary_df <- STRAPP_tests_over_time$pvalues_summary_df

    # Extract data for the selected time range
    if (!is.null(time_range))
    {
      pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$focal_time <= time_range[2], ]
      pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$focal_time >= time_range[1], ]
    } else {
      # Extract time range from data
      time_range <- range(pvalues_summary_df$focal_time)
    }

    # Extract data to avoid 'binding warning'
    p_value <- pvalues_summary_df$p_value
    focal_time <- pvalues_summary_df$focal_time

    # Build ggplot
    pvalues_plot <- ggplot2::ggplot(data = pvalues_summary_df,
                                    mapping = ggplot2::aes(y = p_value, x = focal_time)) +

      # Plot p_values line
      ggplot2::geom_line(col = "black",
                         alpha = 1.0,
                         linewidth = 1.5) +

      # Set plot title +
      ggplot2::ggtitle(label = paste0("STRAPP tests\nDifferences in ",rate_type," rates through time\n")) +

      # Set axes labels
      ggplot2::xlab("Time") +
      ggplot2::ylab("P-value") +

      # Reverse time scale
      ggplot2::scale_x_continuous(
        transform = "reverse",
        expand = c(0, 0),
        limits = rev(time_range) # Set limits
      ) +

      # Reverse p-value scale
      ggplot2::scale_y_continuous(
        transform = "reverse",
        expand = c(0, 0),
        limits = c(1, 0) # Set limits
      ) +

      # Adjust aesthetics
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
        panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
        panel.background = ggplot2::element_rect(fill = NA, color = NA),
        plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                           margin = ggplot2::margin(b = 10, t = 5)),
        axis.title = ggplot2::element_text(size = 20, color = "black"),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
        axis.line = ggplot2::element_line(linewidth = 1.0),
        axis.ticks.length = ggplot2::unit(8, "pt"),
        axis.text = ggplot2::element_text(size = 18, color = "black"),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

    # Add significance line if requested
    if (!is.null(alpha))
    {
      pvalues_plot <- pvalues_plot +
        ggplot2::geom_hline(yintercept = alpha,
                            col = "red",
                            alpha = 1.0,
                            linetype = "dashed",
                            linewidth = 1.0)
    }

  } else { ## Case for post hoc tests plot

    # Extract summary df
    pvalues_summary_df <- STRAPP_tests_over_time$pvalues_summary_df_for_posthoc_pairwise_tests

    # Extract data for the selected time range
    if (!is.null(time_range))
    {
      pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$focal_time <= time_range[2], ]
      pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$focal_time >= time_range[1], ]
    } else {
      # Extract time range from data
      time_range <- range(pvalues_summary_df$focal_time)
    }

    # Extract data for the selected pairs
    if (!("all" %in% select_posthoc_pairs))
    {
      pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$pair %in% select_posthoc_pairs, ]
    }

    # Extract data to avoid 'binding warning'
    p_value <- pvalues_summary_df$p_value
    focal_time <- pvalues_summary_df$focal_time
    pair <- pvalues_summary_df$pair

    # Build ggplot
    pvalues_plot <- ggplot2::ggplot(data = pvalues_summary_df,
                                    mapping = ggplot2::aes(y = p_value, x = focal_time, color = pair)) +

      # Plot p_values line
      ggplot2::geom_line(alpha = 1.0,
                         linewidth = 1.5) +

      # Set legend title
      ggplot2::scale_color_discrete(name = "Pairs") +

      # Set plot title +
      ggplot2::ggtitle(label = paste0("STRAPP tests\nDifferences in ",rate_type," rates through time\n")) +

      # Set axes labels
      ggplot2::xlab("Time") +
      ggplot2::ylab("P-value") +

      # Reverse time scale
      ggplot2::scale_x_continuous(
        transform = "reverse",
        expand = c(0, 0),
        limits = rev(time_range) # Set limits
      ) +

      # Reverse p-value scale
      ggplot2::scale_y_continuous(
        transform = "reverse",
        expand = c(0, 0),
        limits = c(1, 0) # Set limits
      ) +

      # Adjust aesthetics
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
        panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
        panel.background = ggplot2::element_rect(fill = NA, color = NA),
        plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                           margin = ggplot2::margin(b = 10, t = 5)),
        legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.2),
        legend.text = ggplot2::element_text(size = 12),
        legend.key = ggplot2::element_rect(colour = NA, fill = NA, linewidth = 5),
        legend.key.size = ggplot2::unit(1.8, "line"),
        legend.spacing.y = ggplot2::unit(0.5, "line"),
        axis.title = ggplot2::element_text(size = 20, color = "black"),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
        axis.line = ggplot2::element_line(linewidth = 1.0),
        axis.ticks.length = ggplot2::unit(8, "pt"),
        axis.text = ggplot2::element_text(size = 18, color = "black"),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

    # Add significance line if requested
    if (!is.null(alpha))
    {
      pvalues_plot <- pvalues_plot +
        ggplot2::geom_hline(yintercept = alpha,
                            col = "red",
                            alpha = 1.0,
                            linetype = "dashed",
                            linewidth = 1.0)
    }
  }

  ## Display plot if requested
  if (display_plot)
  {
    print(pvalues_plot)
  }

  ## Export plot if requested
  if (!is.null(PDF_file_path))
  {
    cowplot::save_plot(plot = pvalues_plot,
                       filename = PDF_file_path,
                       base_height = 8, base_width = 10)
  }

  ## Return ggplot
  return(invisible(pvalues_plot))

}


