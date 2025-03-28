## Functions to plot rates through time in relation with trait values
# One master function to prepare data and select the proper test function according to data type
# Three sub-functions carrying out tests according to data type

#' @title Plot evolution of diversification rates in relation to trait values over time
#'
#' @description Plot the evolution of diversification rates in relation to trait values
#'   extracted for multiple `time_steps` with [deepSTRAPP::run_STRAPP_tests_over_time()].
#'
#'   Rates are averaged across branches at each time step (i.e., `focal_time`).
#'   * For continuous data, branches are grouped by ranges of trait values defined by `quantile_ranges`.
#'   * For categorical data, branches are grouped by trait states.
#'   * For biogeographic data, branches are grouped by ranges.
#'
#' @param STRAPP_tests_over_time List of elements generated with [deepSTRAPP::run_STRAPP_tests_over_time()],
#'   that summarize the results of multiple STRAPP tests across `$time_steps`. The list needs to include two data.frame:
#'   `$trait_data_df_over_time` and `$diversification_data_df_over_time` by setting `extract_trait_data_melted_df = TRUE`
#'   and `extract_diversification_data_melted_df = TRUE`.
#' @param rate_type A character string specifying the type of diversification rates to use.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#'   Even if the `STRAPP_tests_over_time` object was generated with [deepSTRAPP::run_STRAPP_tests_over_time()]
#'   for testing another type of rates, the `$trait_data_df_over_time` and `$diversification_data_df_over_time` data frames
#'   will contain data for all types of rates.
#' @param quantile_ranges Vector of numerical. Only for continuous trait data. Quantiles used as thresholds to group branches
#'  by trait values. It must start with 0 and finish with 1. Default is `c(0, 0.25, 0.5, 0.75, 1.0)`
#'  which produces four balanced quantile groups.
#' @param select_trait_states (Vector of) character string. Only for categorical and biogeographic trait data.
#'  To provide a list of a subset of states/ranges to plot. Names must match the ones found in
#'  `STRAPP_tests_over_time$trait_data_df_over_time$trait_value`. Default is `all` which means all states/ranges will be plotted.
#' @param time_range Vector of two numerical values. Time boundaries used for the plot.
#'   If `NULL` (the default), the range of data provided in `STRAPP_tests_over_time` will be used.
#' @param plot_CI Logical. Whether to plot a confidence interval (CI) based on the distribution of rates found in posterior samples.
#' @param CI_type Character string. To select the type of confidence interval (CI) to plot.
#'  * `fuzzy` (default): to overlay the evolution of rates found in all posterior samples with high transparency levels.
#'  * `quantiles_rect`: to add a polygon encompassing a proportion of the rate values found in posterior samples.
#'   This proportion is defined with `CI_quantiles`.
#' @param CI_quantiles Numerical. Proportion of rate values encompassed by the confidence interval. Default is `0.95`.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with '.pdf'.
#' @param return_mean_data_per_samples_df Logical. Whether to include in the output the data.frame of mean rates per trait values computed for
#'   each posterior sample at each time-step (aggregated across groups of branches based on trait data). This is used to draw the confidence interval.
#' @param return_median_data_across_samples_df Logical. Whether to include in the output the data.frame of median rates per trait values
#'  across posterior samples computed for at each time-step (aggregated across groups of branches based on trait data AND posterior samples).
#'  This is used to draw the lines on the plot.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_line aes geom_hline geom_polygon scale_y_continuous scale_x_continuous scale_color_discrete scale_color_brewer scale_fill_brewer xlab ylab ggtitle theme element_line element_rect element_text unit margin
#' @importFrom dplyr left_join join_by group_by reframe summarise ungroup mutate arrange select filter
#' @importFrom cowplot save_plot
#' @importFrom stringr str_to_title
#' @importFrom stats quantile
#'
#' @return The function returns a list with at least one element.
#'
#'   * `rates_TT_ggplot` An object of classes `gg` and `ggplot`. This is a ggplot that can be displayed
#'     on the console with `print(output$rates_TT_ggplot)`. It corresponds to the plot being displayed on the console
#'     when the function is run, if `display_plot = TRUE`, and can be further modify for aesthetics using the ggplot2 grammar.
#'
#'   Optional summary data frames:
#'   * `mean_data_per_samples_df` A data.frame with four columns providing the `$mean_rates` observed along branches
#'     with a similar `$trait_value` (if categorical or biogeographic) or falling into the same `$quantile_ranges`.
#'     Data are extracted for each posterior sample (`$BAMM_sample_ID`) at each time-step (i.e., `$focal_time`).
#'     This is used to draw the confidence interval.
#'   * `$median_data_across_samples_df` A data.frame with three columns providing the `$median_rates`
#'   observed across all posterior samples in `$mean_data_per_samples_df`. This is used to draw the lines on the plot.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plot.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::run_STRAPP_tests_over_time()]
#'
#' @examples
#' ## Load results of run_STRAPP_tests_over_time()
#' data(STRAPP_tests_over_time_temp_example_2, package = "deepSTRAPP")
#'
#' ## Plot rates through time for continuous data
#'
#' # Visualize trait data
#' hist(STRAPP_tests_over_time_temp_example_2$trait_data_df_over_time$trait_value)
#'
#' # Generate plot
#' plotTT_continuous <- plot_rates_through_time(
#'   STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example_2,
#'   quantile_ranges = c(0, 0.25, 0.5, 0.75, 1.0),
#'   time_range = c(0, 15),
#'   plot_CI = TRUE,
#'   CI_type = "quantiles_rect",
#'   CI_quantiles = 0.9,
#'   # PDF_file_path = "./plotTT_continuous.pdf",
#'   return_mean_data_per_samples_df = TRUE,
#'   return_median_data_across_samples_df = TRUE
#'   )
#'
#' # Explore output
#' str(plotTT_continuous, max.level = 1)
#'
#' # Plot again
#' print(plotTT_continuous$rates_TT_ggplot)
#' # Adjust aesthetics of plot a posteriori
#' plotTT_continuous_adj <- plotTT_continuous$rates_TT_ggplot +
#'     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#' print(plotTT_continuous_adj)
#'
#' ## Plot rates through time for categorical data
#'
#' # Turn trait data into multiple states
#' trait_data_continuous <- STRAPP_tests_over_time_temp_example_2$trait_data_df_over_time$trait_value
#' trait_data_multinominal <- trait_data_continuous
#' trait_data_multinominal[trait_data_continuous < 0] <- "state_B"
#' trait_data_multinominal[trait_data_continuous < -1] <- "state_A"
#' trait_data_multinominal[trait_data_continuous >= 0] <- "state_C"
#'
#' # Visualize trait data
#' table(trait_data_multinominal)
#'
#' # Change trait data for categorical
#' STRAPP_tests_over_time_temp_example_2$trait_data_df_over_time$trait_value <- trait_data_multinominal
#' STRAPP_tests_over_time_temp_example_2$trait_data_type <- "categorical"
#'
#' # Generate plot only for "state_A" and "state_C"
#' plotTT_categorical <- plot_rates_through_time(
#'   STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example_2,
#'   select_trait_states = c("state_A", "state_C"),
#'   time_range = c(0, 10),
#'   plot_CI = TRUE,
#'   CI_type = "quantiles_rect",
#'   CI_quantiles = 0.9,
#'   # PDF_file_path = "./plotTT_categorical.pdf",
#'   return_mean_data_per_samples_df = TRUE,
#'   return_median_data_across_samples_df = TRUE
#' )
#'
#' # Explore output
#' str(plotTT_categorical, max.level = 1)
#' # Plot again
#' print(plotTT_categorical$rates_TT_ggplot)
#'
#' ## Plot rates through time for biogeographic data
#'
#' # Turn trait data into multiple ranges
#' trait_data_biogeographic <- trait_data_multinominal
#' trait_data_biogeographic[trait_data_multinominal == "state_A"] <- "range_A"
#' trait_data_biogeographic[trait_data_multinominal == "state_B"] <- "range_B"
#' trait_data_biogeographic[trait_data_multinominal == "state_C"] <- "range_C"
#'
#' # Visualize trait data
#' table(trait_data_biogeographic)
#'
#' # Change trait data for biogeographic
#' trait_data_df_over_time <- STRAPP_tests_over_time_temp_example_2$trait_data_df_over_time
#' trait_data_df_over_time$trait_value <- trait_data_biogeographic
#' STRAPP_tests_over_time_temp_example_2$trait_data_df_over_time <- trait_data_df_over_time
#' STRAPP_tests_over_time_temp_example_2$trait_data_type <- "biogeographic"
#'
#' plotTT_biogeographic <- plot_rates_through_time(
#'   STRAPP_tests_over_time = STRAPP_tests_over_time_temp_example_2,
#'   select_trait_states = "all",
#'   time_range = c(0, 10),
#'   plot_CI = TRUE,
#'   CI_type = "quantiles_rect",
#'   CI_quantiles = 0.9,
#'   # PDF_file_path = "./plotTT_biogeographic.pdf",
#'   return_mean_data_per_samples_df = TRUE,
#'   return_median_data_across_samples_df = TRUE
#' )
#'
#' # Explore output
#' str(plotTT_biogeographic, max.level = 1)
#' # Plot again
#' print(plotTT_biogeographic$rates_TT_ggplot)
#'


### Master function to prepare data and select the proper test function according to data type ####

plot_rates_through_time <- function (
    STRAPP_tests_over_time,
    rate_type = "net_diversification",
    quantile_ranges = c(0, 0.25, 0.5, 0.75, 1.0),
    select_trait_states = "all",
    time_range = NULL,
    plot_CI = FALSE,
    CI_type = "fuzzy",
    CI_quantiles = 0.95,
    display_plot = TRUE,
    PDF_file_path = NULL,
    return_mean_data_per_samples_df = FALSE,
    return_median_data_across_samples_df = FALSE
)
{
  ### Check input validity

  # STRAPP_tests_over_time must contain $trait_data_df_over_time and $diversification_data_df_over_time
  # Special warning to use extract_trait_data_melted_df = TRUE and extract_diversification_data_melted_df = TRUE to get them

  # If provided, PDF_file_path must end with ".pdf"

  # Other checks are carried in sub-functions

  ## Detect the type of trait data
  trait_data_type <- STRAPP_tests_over_time$trait_data_type

  ## Compute the appropriate internal function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous =   { # Case for continuous data
           # Need quantile_ranges to define groups of branches per trait values
           plotTT_output <- plot_rates_through_time_for_continuous_data(
             STRAPP_tests_over_time = STRAPP_tests_over_time,
             rate_type = rate_type,
             quantile_ranges = quantile_ranges,
             time_range = time_range,
             plot_CI = plot_CI,
             CI_type = CI_type,
             CI_quantiles = CI_quantiles,
             display_plot = display_plot,
             PDF_file_path = PDF_file_path,
             return_mean_data_per_samples_df = return_mean_data_per_samples_df,
             return_median_data_across_samples_df = return_median_data_across_samples_df
           )
         },
         categorical =  { # Case for categorical data
           # Can select the states to plot
           plotTT_output <- plot_rates_through_time_for_categorical_data(
             STRAPP_tests_over_time = STRAPP_tests_over_time,
             rate_type = rate_type,
             select_trait_states = select_trait_states,
             time_range = time_range,
             plot_CI = plot_CI,
             CI_type = CI_type,
             CI_quantiles = CI_quantiles,
             display_plot = display_plot,
             PDF_file_path = PDF_file_path,
             return_mean_data_per_samples_df = return_mean_data_per_samples_df,
             return_median_data_across_samples_df = return_median_data_across_samples_df
           )
         },
         biogeographic = { # Case for biogeographic data
           # Can select the states/ranges to plot
           plotTT_output <- plot_rates_through_time_for_biogeographic_data(
             STRAPP_tests_over_time = STRAPP_tests_over_time,
             rate_type = rate_type,
             select_trait_states = select_trait_states,
             time_range = time_range,
             plot_CI = plot_CI,
             CI_type = CI_type,
             CI_quantiles = CI_quantiles,
             display_plot = display_plot,
             PDF_file_path = PDF_file_path,
             return_mean_data_per_samples_df = return_mean_data_per_samples_df,
             return_median_data_across_samples_df = return_median_data_across_samples_df
           )
         }
  )

  ## Export the output
  return(invisible(plotTT_output))
}



### Sub-function to handle continuous data ####

plot_rates_through_time_for_continuous_data <- function (
    STRAPP_tests_over_time,
    rate_type = "net_diversification",
    quantile_ranges = c(0, 0.25, 0.5, 0.75, 1.0),
    time_range = NULL,
    plot_CI = FALSE,
    CI_type = "fuzzy",
    CI_quantiles = 0.95,
    display_plot = TRUE,
    PDF_file_path = NULL,
    return_mean_data_per_samples_df = FALSE,
    return_median_data_across_samples_df = FALSE
)
{
  ### Check input validity

  # STRAPP_tests_over_time must contain $trait_data_df_over_time and $diversification_data_df_over_time
  # Special warning to use extract_trait_data_melted_df = TRUE and extract_diversification_data_melted_df = TRUE to get them

  # quantile_ranges contains numerical between 0 and 1, in increasing order.
  # Special warning that display the provided argument if error.
  # If range(quantile_ranges) is not c(0, 1), and the 0, 1 boundaries to the list of quantile_ranges and print a WARNING to say they have been added.

  # Check that time_range is strictly positive, ordered in increasing age, and encompass multiple data points
  # Can also order by doing time_range <- range(time_range) instead of stopping

  # CI_type is either "fuzzy" or "quantiles_rect"

  # CI_quantiles is a numerical between 0 and 1

  # If provided, PDF_file_path must end with ".pdf"

  ## Adjust rate_type for labels
  rate_type_label <- stringr::str_to_title(rate_type)
  rate_type_label <- gsub(pattern = "_", replacement = " ", x = rate_type_label)

  ## Create binding of new variables to avoid Notes
  tip_ID <- BAMM_sample_ID <- focal_time <- quant_traits <- NULL
  trait_value <- rates <- median_rates <- mean_rates <- NULL
  n_points <- points_ID <- quant_rates <- NULL

  ## Merge diversification and trait data
  # Trait data are copied across BAMM samples
  data_per_samples_df <- dplyr::left_join(
    x = STRAPP_tests_over_time$diversification_data_df_over_time,
    y = STRAPP_tests_over_time$trait_data_df_over_time,
    by = dplyr::join_by(focal_time, tip_ID))

  ## Filter data for selected rate_type
  data_per_samples_df <- data_per_samples_df[data_per_samples_df$rate_type == rate_type, ]

  # Filter data for the selected time range
  if (!is.null(time_range))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$focal_time <= time_range[2], ]
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$focal_time >= time_range[1], ]
  } else {
    # Extract time range from data
    time_range <- range(data_per_samples_df$focal_time)
  }

  if (nrow(data_per_samples_df) == 0)
  {
    stop("No data found in the time range c(",time_range[1],", ", time_range[2],").\n")
  }

  ## Compute quantile thresholds for each focal_time
  quantiles_data_df <- data_per_samples_df |>
    dplyr::group_by(focal_time) |>
    # Compute quantiles of trait value
    dplyr::reframe(quant_traits = stats::quantile(trait_value, probs = quantile_ranges, na.rm = T))

  ## Attribute a quantile range to trait values found across branches for each focal_time
  # (will be copy across BAMM samples as trait value do not change across BAMM samples)

  data_per_samples_df$quantile_ranges <- NA
  # Loop per focal_time
  focal_time_list <- unique(data_per_samples_df$focal_time)
  for (i in seq_along(focal_time_list))
  {
    # i <- 1

    # Extract focal_time
    focal_time_i <- focal_time_list[i]

    # Get thresholds for this focal_time
    quantile_thresholds <- quantiles_data_df$quant_traits[quantiles_data_df$focal_time == focal_time_i]

    # Loop per quantile ranges
    # From highest to lowest
    for (j in length(quantile_ranges):2)
    {
      # j <- 2

      # Get quantile range name
      quantile_range_name <- paste0("Q",quantile_ranges[j-1]*100,"% - Q",quantile_ranges[j]*100,"%")

      # Get max threshold for this quantile range
      threshold_j <- quantile_thresholds[j]

      # Inform $quantile_ranges by attributing all branches with rates <= than the max threshold
      data_per_samples_df$quantile_ranges[(data_per_samples_df$focal_time == focal_time_i) & (data_per_samples_df$trait_value <= threshold_j)] <- quantile_range_name
    }
  }
  # table(data_per_samples_df$quantile_ranges) # Should be roughly equally distributed

  ## Aggregate across tip_ID (branches), per quantile ranges
  mean_data_per_samples_df <- data_per_samples_df |>
    dplyr::group_by(focal_time, BAMM_sample_ID, quantile_ranges) |>
    dplyr::summarise(mean_rates = mean(rates), .groups = "keep") |>
    dplyr::ungroup()

  ## Aggregate across BAMM samples
  median_data_across_samples_df <- mean_data_per_samples_df |>
    dplyr::group_by(focal_time, quantile_ranges) |>
    dplyr::summarise(median_rates = median(mean_rates), .groups = "keep") |>
    dplyr::ungroup()

  ## Case for plot without CI
  if (!plot_CI)
  {
    rates_TT_ggplot <- ggplot2::ggplot(data = median_data_across_samples_df) +

      # Plot mean lines
      ggplot2::geom_line(mapping = aes(y = median_rates, x = focal_time,
                                       group = quantile_ranges, col = quantile_ranges),
                         alpha = 1.0,
                         linewidth = 1.5) +

      # Plot div = 0 line
      ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

      # Set plot title +
      ggplot2::ggtitle(label = paste0(rate_type_label, " rates per trait values through time")) +

      # Set axes labels
      ggplot2::xlab("Time") +
      ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

      # Prevent rate Y-scale to expand
      ggplot2::scale_y_continuous(expand = c(0, 0)) +

      # Reverse time scale
      ggplot2::scale_x_continuous(transform = "reverse",
                                  limits = rev(time_range)) +

      # Adjust color scheme and legend
      ggplot2::scale_color_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +

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

  } else { ## Case for plot with CI

    if (CI_type == "fuzzy")
    {
      ## Plot with fuzzy CI

      rates_TT_ggplot <- ggplot2::ggplot(data = mean_data_per_samples_df) +

        # Plot line replicates for all samples
        ggplot2::geom_line(data = mean_data_per_samples_df,
                           mapping = aes(y = mean_rates, x = focal_time,
                                         group = interaction(quantile_ranges, BAMM_sample_ID),
                                         col = quantile_ranges),
                           alpha = 0.01,
                           linewidth = 3.0) +

        # Plot mean lines
        ggplot2::geom_line(data = median_data_across_samples_df,
                           mapping = aes(y = median_rates, x = focal_time,
                                         group = quantile_ranges, col = quantile_ranges),
                           alpha = 1.0,
                           linewidth = 1.5) +

        # Plot div = 0 line
        ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

        # Set plot title +
        ggplot2::ggtitle(label = paste0(rate_type_label, " rates per trait values through time")) +

        # Set axes labels
        ggplot2::xlab("Time") +
        ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

        # Prevent rate Y-scale to expand
        ggplot2::scale_y_continuous(expand = c(0, 0)) +

        # Reverse time scale
        ggplot2::scale_x_continuous(transform = "reverse",
                                    limits = rev(time_range)) +

        # Adjust color scheme and legend
        ggplot2::scale_color_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +

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


    } else {

      ## Plot with quantiles_rect CI

      ## Convert CI quantiles to the proportion of data to NOT include
      CI_quantiles_inv <- (1 - CI_quantiles)

      ## Create data.frame for quantile polygons
      quantiles_mean_data_df <- mean_data_per_samples_df |>
        dplyr::group_by(focal_time, quantile_ranges) |>
        # Compute quantiles
        dplyr::reframe(quant_rates = stats::quantile(mean_rates, probs = c(CI_quantiles_inv/2, (1 - CI_quantiles_inv/2)), na.rm = T)) |>
        dplyr::group_by(focal_time, quantile_ranges) |>
        dplyr::mutate(quantile = c(CI_quantiles_inv/2, (1 - CI_quantiles_inv/2))) |>
        # Assign points ID (order for drawing the polygon)
        dplyr::group_by(quantile_ranges) |>
        dplyr::arrange(quantile_ranges, quantile) |>
        dplyr::mutate(n_points = dplyr::n()) |> # Count the number of points in a polygon
        dplyr::mutate(points_ID = c(1:(dplyr::first(n_points)/2), dplyr::first(n_points):((dplyr::first(n_points)/2) + 1))) |>
        dplyr::select(-n_points) |>
        # Reorder by points ID
        dplyr::arrange(quantile_ranges, points_ID) |>
        # Filter for NA
        dplyr::filter(!is.na(quant_rates)) |>
        # Reattribute points_ID after filtering
        dplyr::mutate(points_ID = dplyr::row_number()) |>
        dplyr::ungroup()

      rates_TT_ggplot <- ggplot2::ggplot(data = quantiles_mean_data_df) +

        # Plot quantile polygons
        ggplot2::geom_polygon(data = quantiles_mean_data_df,
                              mapping = aes(y = quant_rates, x = focal_time,
                                            group = quantile_ranges,
                                            fill = quantile_ranges),
                              alpha = 0.3,
                              linewidth = 1.0) +

        # Plot mean lines
        ggplot2::geom_line(data = median_data_across_samples_df,
                           mapping = aes(y = median_rates, x = focal_time,
                                         group = quantile_ranges, col = quantile_ranges),
                           alpha = 1.0,
                           linewidth = 1.5) +

        # Plot div = 0 line
        ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

        # Set plot title +
        ggplot2::ggtitle(label = paste0(rate_type_label, " rates per trait values through time")) +

        # Set axes labels
        ggplot2::xlab("Time") +
        ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

        # Prevent rate Y-scale to expand
        ggplot2::scale_y_continuous(expand = c(0, 0)) +

        # Reverse time scale
        ggplot2::scale_x_continuous(transform = "reverse",
                                    limits = rev(time_range)) +

        # Adjust fill scheme and legend
        ggplot2::scale_fill_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +

        # Adjust color scheme and legend
        ggplot2::scale_color_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +

        # Remove fill legend
        ggplot2::guides(fill = "none") +

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

    }
  }

  ## Display plot if requested
  if (display_plot)
  {
    print(rates_TT_ggplot)
  }

  ## Export plot if requested
  if (!is.null(PDF_file_path))
  {
    cowplot::save_plot(plot = rates_TT_ggplot,
                       filename = PDF_file_path,
                       base_height = 8, base_width = 10)
  }

  ## Build output
  output <- list()

  ## Store ggplot
  output$rates_TT_ggplot <- rates_TT_ggplot

  ## Store melted df if requested
  if (return_mean_data_per_samples_df)
  {
    output$mean_data_per_samples_df <- as.data.frame(mean_data_per_samples_df)
  }
  if (return_median_data_across_samples_df)
  {
    output$median_data_across_samples_df <- as.data.frame(median_data_across_samples_df)
  }

  ## Return output
  return(invisible(output))
}


### Sub-function to handle categorical data ####

plot_rates_through_time_for_categorical_data <- function (
    STRAPP_tests_over_time,
    rate_type = "net_diversification",
    select_trait_states = "all",
    time_range = NULL,
    plot_CI = FALSE,
    CI_type = "fuzzy",
    CI_quantiles = 0.95,
    display_plot = TRUE,
    PDF_file_path = NULL,
    return_mean_data_per_samples_df = FALSE,
    return_median_data_across_samples_df = FALSE
)
{
  ### Check input validity

  # STRAPP_tests_over_time must contain $trait_data_df_over_time and $diversification_data_df_over_time
  # Special warning to use extract_trait_data_melted_df = TRUE and extract_diversification_data_melted_df = TRUE to get them

  # select_trait_states must be "all" or states that are found in the $trait_data_df_over_time and $diversification_data_df_over_time data.frame
  # Special warning that display the provided argument and the observed states

  # Check that time_range is strictly positive, ordered in increasing age, and encompass multiple data points
  # Can also order by doing time_range <- range(time_range) instead of stopping

  # CI_type is either "fuzzy" or "quantiles_rect"

  # CI_quantiles is a numerical between 0 and 1.0

  # If provided, PDF_file_path must end with ".pdf"

  ## Create binding of new variables to avoid Notes
  tip_ID <- BAMM_sample_ID <- focal_time <- quant_traits <- NULL
  trait_value <- rates <- median_rates <- mean_rates <- NULL
  n_points <- points_ID <- quant_rates <- NULL

  ## Adjust rate_type for labels
  rate_type_label <- stringr::str_to_title(rate_type)
  rate_type_label <- gsub(pattern = "_", replacement = " ", x = rate_type_label)

  ## Merge diversification and trait data
  # Trait data are copied across BAMM samples
  data_per_samples_df <- dplyr::left_join(
    x = STRAPP_tests_over_time$diversification_data_df_over_time,
    y = STRAPP_tests_over_time$trait_data_df_over_time,
    by = dplyr::join_by(focal_time, tip_ID))

  ## Filter data for selected rate_type
  data_per_samples_df <- data_per_samples_df[data_per_samples_df$rate_type == rate_type, ]

  ## Filter data for selected states
  if (!("all" %in% select_trait_states))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$trait_value %in% select_trait_states, ]
  }

  # Filter data for the selected time range
  if (!is.null(time_range))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$focal_time <= time_range[2], ]
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$focal_time >= time_range[1], ]
  } else {
    # Extract time range from data
    time_range <- range(data_per_samples_df$focal_time)
  }

  if (nrow(data_per_samples_df) == 0)
  {
    stop("No data found in the time range c(",time_range[1],", ", time_range[2],") for ",paste(select_trait_states, collapse = ", ")," states.\n")
  }

  ## Aggregate across tip_ID (branches), per trait states
  mean_data_per_samples_df <- data_per_samples_df |>
    dplyr::group_by(focal_time, BAMM_sample_ID, trait_value) |>
    dplyr::summarise(mean_rates = mean(rates), .groups = "keep") |>
    dplyr::ungroup()

  ## Aggregate across BAMM samples
  median_data_across_samples_df <- mean_data_per_samples_df |>
    dplyr::group_by(focal_time, trait_value) |>
    dplyr::summarise(median_rates = median(mean_rates), .groups = "keep") |>
    dplyr::ungroup()

  ## Case for plot without CI
  if (!plot_CI)
  {
    rates_TT_ggplot <- ggplot2::ggplot(data = median_data_across_samples_df) +

      # Plot mean lines
      ggplot2::geom_line(mapping = aes(y = median_rates, x = focal_time,
                                       group = trait_value, col = trait_value),
                         alpha = 1.0,
                         linewidth = 1.5) +

      # Plot div = 0 line
      ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

      # Set plot title +
      ggplot2::ggtitle(label = paste0(rate_type_label, " rates per trait states through time")) +

      # Set axes labels
      ggplot2::xlab("Time") +
      ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

      # Prevent rate Y-scale to expand
      ggplot2::scale_y_continuous(expand = c(0, 0)) +

      # Reverse time scale
      ggplot2::scale_x_continuous(transform = "reverse",
                                  limits = rev(time_range)) +

      # Adjust color scheme and legend
      ggplot2::scale_color_discrete(name = "States") +

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

  } else { ## Case for plot with CI

    if (CI_type == "fuzzy")
    {
      ## Plot with fuzzy CI

      rates_TT_ggplot <- ggplot2::ggplot(data = mean_data_per_samples_df) +

        # Plot line replicates for all samples
        ggplot2::geom_line(data = mean_data_per_samples_df,
                           mapping = aes(y = mean_rates, x = focal_time,
                                         group = interaction(trait_value, BAMM_sample_ID),
                                         col = trait_value),
                           alpha = 0.01,
                           linewidth = 3.0) +

        # Plot mean lines
        ggplot2::geom_line(data = median_data_across_samples_df,
                           mapping = aes(y = median_rates, x = focal_time,
                                         group = trait_value, col = trait_value),
                           alpha = 1.0,
                           linewidth = 1.5) +

        # Plot div = 0 line
        ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

        # Set plot title +
        ggplot2::ggtitle(label = paste0(rate_type_label, " rates per trait states through time")) +

        # Set axes labels
        ggplot2::xlab("Time") +
        ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

        # Prevent rate Y-scale to expand
        ggplot2::scale_y_continuous(expand = c(0, 0)) +

        # Reverse time scale
        ggplot2::scale_x_continuous(transform = "reverse",
                                    limits = rev(time_range)) +

        # Adjust color scheme and legend
        ggplot2::scale_color_discrete(name = "States") +

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


    } else {

      ## Plot with quantiles_rect CI

      ## Convert CI quantiles to the proportion of data to not include
      CI_quantiles_inv <- (1 - CI_quantiles)

      ## Create data.frame for quantile polygons
      quantiles_mean_data_df <- mean_data_per_samples_df |>
        dplyr::group_by(focal_time, trait_value) |>
        # Compute quantiles
        dplyr::reframe(quant_rates = stats::quantile(mean_rates, probs = c(CI_quantiles_inv/2, (1 - CI_quantiles_inv/2)), na.rm = T)) |>
        dplyr::group_by(focal_time, trait_value) |>
        dplyr::mutate(quantile = c(CI_quantiles_inv/2, (1 - CI_quantiles_inv/2))) |>
        # Assign points ID (order for drawing the polygon)
        dplyr::group_by(trait_value) |>
        dplyr::arrange(trait_value, quantile) |>
        dplyr::mutate(n_points = dplyr::n()) |> # Count the number of points in a polygon
        dplyr::mutate(points_ID = c(1:(dplyr::first(n_points)/2), dplyr::first(n_points):((dplyr::first(n_points)/2) + 1))) |>
        dplyr::select(-n_points) |>
        # Reorder by points ID
        dplyr::arrange(trait_value, points_ID) |>
        # Filter for NA
        dplyr::filter(!is.na(quant_rates)) |>
        # Reattribute points_ID after filtering
        dplyr::mutate(points_ID = dplyr::row_number()) |>
        dplyr::ungroup()

      rates_TT_ggplot <- ggplot2::ggplot(data = quantiles_mean_data_df) +

        # Plot quantile polygons
        ggplot2::geom_polygon(data = quantiles_mean_data_df,
                              mapping = aes(y = quant_rates, x = focal_time,
                                            group = trait_value,
                                            fill = trait_value),
                              alpha = 0.3,
                              linewidth = 1.0) +

        # Plot mean lines
        ggplot2::geom_line(data = median_data_across_samples_df,
                           mapping = aes(y = median_rates, x = focal_time,
                                         group = trait_value, col = trait_value),
                           alpha = 1.0,
                           linewidth = 1.5) +

        # Plot div = 0 line
        ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

        # Set plot title +
        ggplot2::ggtitle(label = paste0(rate_type_label, " rates per trait states through time")) +

        # Set axes labels
        ggplot2::xlab("Time") +
        ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

        # Prevent rate Y-scale to expand
        ggplot2::scale_y_continuous(expand = c(0, 0)) +

        # Reverse time scale
        ggplot2::scale_x_continuous(transform = "reverse",
                                    limits = rev(time_range)) +

        # Adjust fill scheme and legend
        ggplot2::scale_fill_discrete("States") +

        # Adjust color scheme and legend
        ggplot2::scale_color_discrete(name = "States") +

        # Remove fill legend
        ggplot2::guides(fill = "none") +

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

    }
  }

  ## Display plot if requested
  if (display_plot)
  {
    print(rates_TT_ggplot)
  }

  ## Export plot if requested
  if (!is.null(PDF_file_path))
  {
    cowplot::save_plot(plot = rates_TT_ggplot,
                       filename = PDF_file_path,
                       base_height = 8, base_width = 10)
  }

  ## Build output
  output <- list()

  ## Store ggplot
  output$rates_TT_ggplot <- rates_TT_ggplot

  ## Store melted df if requested
  if (return_mean_data_per_samples_df)
  {
    output$mean_data_per_samples_df <- as.data.frame(mean_data_per_samples_df)
  }
  if (return_median_data_across_samples_df)
  {
    output$median_data_across_samples_df <- as.data.frame(median_data_across_samples_df)
  }

  ## Return output
  return(invisible(output))
}


### Sub-function to handle biogeographic data ####

plot_rates_through_time_for_biogeographic_data <- function (
    STRAPP_tests_over_time,
    rate_type = "net_diversification",
    select_trait_states = "all",
    time_range = NULL,
    plot_CI = FALSE,
    CI_type = "fuzzy",
    CI_quantiles = 0.95,
    display_plot = TRUE,
    PDF_file_path = NULL,
    return_mean_data_per_samples_df = FALSE,
    return_median_data_across_samples_df = FALSE
)
{
  ### Check input validity

  # STRAPP_tests_over_time must contain $trait_data_df_over_time and $diversification_data_df_over_time
  # Special warning to use extract_trait_data_melted_df = TRUE and extract_diversification_data_melted_df = TRUE to get them

  # select_trait_states must be "all" or states that are found in the $trait_data_df_over_time and $diversification_data_df_over_time data.frame
  # Special warning that display the provided argument and the observed states

  # Check that time_range is strictly positive, ordered in increasing age, and encompass multiple data points
  # Can also order by doing time_range <- range(time_range) instead of stopping

  # CI_type is either "fuzzy" or "quantiles_rect"

  # CI_quantiles is a numerical between 0 and 1.0

  # If provided, PDF_file_path must end with ".pdf"

  ## Create binding of new variables to avoid Notes
  tip_ID <- BAMM_sample_ID <- focal_time <- quant_traits <- NULL
  trait_value <- rates <-  median_rates <- mean_rates <- NULL
  n_points <- points_ID <- quant_rates <- NULL

  ## Adjust rate_type for labels
  rate_type_label <- stringr::str_to_title(rate_type)
  rate_type_label <- gsub(pattern = "_", replacement = " ", x = rate_type_label)

  ## Merge diversification and trait data
  # Trait data are copied across BAMM samples
  data_per_samples_df <- dplyr::left_join(
    x = STRAPP_tests_over_time$diversification_data_df_over_time,
    y = STRAPP_tests_over_time$trait_data_df_over_time,
    by = dplyr::join_by(focal_time, tip_ID))

  ## Filter data for selected rate_type
  data_per_samples_df <- data_per_samples_df[data_per_samples_df$rate_type == rate_type, ]

  ## Filter data for selected states/ranges
  if (!("all" %in% select_trait_states))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$trait_value %in% select_trait_states, ]
  }

  # Filter data for the selected time range
  if (!is.null(time_range))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$focal_time <= time_range[2], ]
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$focal_time >= time_range[1], ]
  } else {
    # Extract time range from data
    time_range <- range(data_per_samples_df$focal_time)
  }

  if (nrow(data_per_samples_df) == 0)
  {
    stop("No data found in the time range c(",time_range[1],", ", time_range[2],") in ",paste(select_trait_states, collapse = ", ")," ranges.\n")
  }

  ## Aggregate across tip_ID (branches), per trait ranges
  mean_data_per_samples_df <- data_per_samples_df |>
    dplyr::group_by(focal_time, BAMM_sample_ID, trait_value) |>
    dplyr::summarise(mean_rates = mean(rates), .groups = "keep") |>
    dplyr::ungroup()

  ## Aggregate across BAMM samples
  median_data_across_samples_df <- mean_data_per_samples_df |>
    dplyr::group_by(focal_time, trait_value) |>
    dplyr::summarise(median_rates = median(mean_rates), .groups = "keep") |>
    dplyr::ungroup()

  ## Case for plot without CI
  if (!plot_CI)
  {
    rates_TT_ggplot <- ggplot2::ggplot(data = median_data_across_samples_df) +

      # Plot mean lines
      ggplot2::geom_line(mapping = aes(y = median_rates, x = focal_time,
                                       group = trait_value, col = trait_value),
                         alpha = 1.0,
                         linewidth = 1.5) +

      # Plot div = 0 line
      ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

      # Set plot title +
      ggplot2::ggtitle(label = paste0(rate_type_label, " rates per ranges through time")) +

      # Set axes labels
      ggplot2::xlab("Time") +
      ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

      # Prevent rate Y-scale to expand
      ggplot2::scale_y_continuous(expand = c(0, 0)) +

      # Reverse time scale
      ggplot2::scale_x_continuous(transform = "reverse",
                                  limits = rev(time_range)) +

      # Adjust color scheme and legend
      ggplot2::scale_color_discrete(name = "Ranges") +

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

  } else { ## Case for plot with CI

    if (CI_type == "fuzzy")
    {
      ## Plot with fuzzy CI

      rates_TT_ggplot <- ggplot2::ggplot(data = mean_data_per_samples_df) +

        # Plot line replicates for all samples
        ggplot2::geom_line(data = mean_data_per_samples_df,
                           mapping = aes(y = mean_rates, x = focal_time,
                                         group = interaction(trait_value, BAMM_sample_ID),
                                         col = trait_value),
                           alpha = 0.01,
                           linewidth = 3.0) +

        # Plot mean lines
        ggplot2::geom_line(data = median_data_across_samples_df,
                           mapping = aes(y = median_rates, x = focal_time,
                                         group = trait_value, col = trait_value),
                           alpha = 1.0,
                           linewidth = 1.5) +

        # Plot div = 0 line
        ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

        # Set plot title +
        ggplot2::ggtitle(label = paste0(rate_type_label, " rates per ranges through time")) +

        # Set axes labels
        ggplot2::xlab("Time") +
        ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

        # Prevent rate Y-scale to expand
        ggplot2::scale_y_continuous(expand = c(0, 0)) +

        # Reverse time scale
        ggplot2::scale_x_continuous(transform = "reverse",
                                    limits = rev(time_range)) +

        # Adjust color scheme and legend
        ggplot2::scale_color_discrete(name = "Ranges") +

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


    } else {

      ## Plot with quantiles_rect CI

      ## Convert CI quantiles to the proportion of data to not include
      CI_quantiles_inv <- (1 - CI_quantiles)

      ## Create data.frame for quantile polygons
      quantiles_mean_data_df <- mean_data_per_samples_df |>
        dplyr::group_by(focal_time, trait_value) |>
        # Compute quantiles
        dplyr::reframe(quant_rates = stats::quantile(mean_rates, probs = c(CI_quantiles_inv/2, (1 - CI_quantiles_inv/2)), na.rm = T)) |>
        dplyr::group_by(focal_time, trait_value) |>
        dplyr::mutate(quantile = c(CI_quantiles_inv/2, (1 - CI_quantiles_inv/2))) |>
        # Assign points ID (order for drawing the polygon)
        dplyr::group_by(trait_value) |>
        dplyr::arrange(trait_value, quantile) |>
        dplyr::mutate(n_points = dplyr::n()) |> # Count the number of points in a polygon
        dplyr::mutate(points_ID = c(1:(dplyr::first(n_points)/2), dplyr::first(n_points):((dplyr::first(n_points)/2) + 1))) |>
        dplyr::select(-n_points) |>
        # Reorder by points ID
        dplyr::arrange(trait_value, points_ID) |>
        # Filter for NA
        dplyr::filter(!is.na(quant_rates)) |>
        # Reattribute points_ID after filtering
        dplyr::mutate(points_ID = dplyr::row_number()) |>
        dplyr::ungroup()

      rates_TT_ggplot <- ggplot2::ggplot(data = quantiles_mean_data_df) +

        # Plot quantile polygons
        ggplot2::geom_polygon(data = quantiles_mean_data_df,
                              mapping = aes(y = quant_rates, x = focal_time,
                                            group = trait_value,
                                            fill = trait_value),
                              alpha = 0.3,
                              linewidth = 1.0) +

        # Plot mean lines
        ggplot2::geom_line(data = median_data_across_samples_df,
                           mapping = aes(y = median_rates, x = focal_time,
                                         group = trait_value, col = trait_value),
                           alpha = 1.0,
                           linewidth = 1.5) +

        # Plot div = 0 line
        ggplot2::geom_hline(yintercept = 0, linewidth = 1.0, linetype = "dashed") +

        # Set plot title +
        ggplot2::ggtitle(label = paste0(rate_type_label, " rates per ranges through time")) +

        # Set axes labels
        ggplot2::xlab("Time") +
        ggplot2::ylab(paste0(rate_type_label, " rates\n[Events / lineage / My]")) +

        # Prevent rate Y-scale to expand
        ggplot2::scale_y_continuous(expand = c(0, 0)) +

        # Reverse time scale
        ggplot2::scale_x_continuous(transform = "reverse",
                                    limits = rev(time_range)) +

        # Adjust fill scheme and legend
        ggplot2::scale_fill_discrete("Ranges") +

        # Adjust color scheme and legend
        ggplot2::scale_color_discrete(name = "Ranges") +

        # Remove fill legend
        ggplot2::guides(fill = "none") +

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

    }
  }

  ## Display plot if requested
  if (display_plot)
  {
    print(rates_TT_ggplot)
  }

  ## Export plot if requested
  if (!is.null(PDF_file_path))
  {
    cowplot::save_plot(plot = rates_TT_ggplot,
                       filename = PDF_file_path,
                       base_height = 8, base_width = 10)
  }

  ## Build output
  output <- list()

  ## Store ggplot
  output$rates_TT_ggplot <- rates_TT_ggplot

  ## Store melted df if requested
  if (return_mean_data_per_samples_df)
  {
    output$mean_data_per_samples_df <- as.data.frame(mean_data_per_samples_df)
  }
  if (return_median_data_across_samples_df)
  {
    output$median_data_across_samples_df <- as.data.frame(median_data_across_samples_df)
  }

  ## Return output
  return(invisible(output))
}

