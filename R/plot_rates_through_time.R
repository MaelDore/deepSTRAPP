## Functions to plot rates through time in relation with trait values
# One master function to prepare data and select the proper test function according to data type
# Three sub-functions carrying out tests according to data type

#' @title Plot evolution of diversification rates in relation to trait values over time
#'
#' @description Plot the evolution of diversification rates in relation to trait values
#'   extracted for multiple `time_steps` with [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'   Rates are averaged across branches at each time step (i.e., `focal_time`).
#'   * For continuous data, branches are grouped by ranges of trait values defined by `quantile_ranges`.
#'   * For categorical data, branches are grouped by trait states.
#'   * For biogeographic data, branches are grouped by ranges.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   that summarize the results of multiple STRAPP tests across `$time_steps`. The list needs to include two data.frame:
#'   `$trait_data_df_over_time` and `$diversification_data_df_over_time` by setting `extract_trait_data_melted_df = TRUE`
#'   and `extract_diversification_data_melted_df = TRUE`.
#' @param rate_type A character string specifying the type of diversification rates to use.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#'   Even if the `deepSTRAPP_outputs` object was generated with [deepSTRAPP::run_deepSTRAPP_over_time()]
#'   for testing another type of rates, the `$trait_data_df_over_time` and `$diversification_data_df_over_time` data frames
#'   will contain data for all types of rates.
#' @param quantile_ranges Vector of numerical. Only for continuous trait data. Quantiles used as thresholds to group branches
#'  by trait values. It must start with 0 and finish with 1. Default is `c(0, 0.25, 0.5, 0.75, 1.0)`
#'  which produces four balanced quantile groups.
#' @param select_trait_levels (Vector of) character string. Only for categorical and biogeographic trait data.
#'  To provide a list of a subset of states/ranges to plot. Names must match the ones found in
#'  `deepSTRAPP_outputs$trait_data_df_over_time$trait_value`. Default is `all` which means all states/ranges will be plotted.
#' @param time_range Vector of two numerical values. Time boundaries used for the plot.
#'   If `NULL` (the default), the range of data provided in `deepSTRAPP_outputs` will be used.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   to display the quantile groups used to discretize the continuous trait data. From lowest values to highest values. Only for continuous data.
#'   Default = `NULL` will use the 'Spectral' color palette in [RColorBrewer::brewer.pal()].
#' @param colors_per_levels Named character string. To set the colors to use to plot rates of each state/range. Names = states/ranges; values = colors.
#'   If `NULL` (default), the default ggplot2 color palette ([scales::hue_pal()]) will be used. Only for categorical and biogeographic data.
#' @param plot_CI Logical. Whether to plot a confidence interval (CI) based on the distribution of rates found in posterior samples. Default is `FALSE`.
#' @param CI_type Character string. To select the type of confidence interval (CI) to plot.
#'  * `fuzzy` (default): to overlay the evolution of rates found in all posterior samples with high transparency levels.
#'  * `quantiles_rect`: to add a polygon encompassing a proportion of the rate values found in posterior samples.
#'   This proportion is defined with `CI_quantiles`.
#' @param CI_quantiles Numerical. Proportion of rate values across posterior samples encompassed by the confidence interval. Only if `CI_type = "quantiles_rect"`. Default is `0.95`.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#' @param return_mean_data_per_samples_df Logical. Whether to include in the output the data.frame of mean rates per trait values computed for
#'   each posterior sample at each time-step (aggregated across groups of branches based on trait data). This is used to draw the confidence interval. Default is `FALSE`.
#' @param return_median_data_across_samples_df Logical. Whether to include in the output the data.frame of median rates per trait values
#'  across posterior samples computed for at each time-step (aggregated across groups of branches based on trait data AND posterior samples).
#'  This is used to draw the lines on the plot. Default is `FALSE`.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_line aes geom_hline geom_polygon scale_y_continuous scale_x_continuous scale_color_discrete scale_color_brewer scale_fill_brewer xlab ylab ggtitle theme element_line element_rect element_text unit margin
#' @importFrom dplyr left_join join_by group_by reframe summarise ungroup mutate arrange select filter
#' @importFrom cowplot save_plot
#' @importFrom stringr str_to_title
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales hue_pal
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
#'     This is used to draw the confidence interval. Included if `return_mean_data_per_samples_df = TRUE`.
#'   * `$median_data_across_samples_df` A data.frame with three columns providing the `$median_rates`
#'     observed across all posterior samples in `$mean_data_per_samples_df`. This is used to draw the lines on the plot.
#'     Included if `return_median_data_across_samples_df = TRUE`.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plot.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::run_deepSTRAPP_over_time()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("plot_rates_through_time", package = "deepSTRAPP")}
#'
#' @examples
#'
#' # ------ Example 1: Plot rates through time for continuous data ------ #
#'
#' if (deepSTRAPP:::is_dev_version())
#' {
#'   ## Load results of run_deepSTRAPP_over_time()
#'   data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
#'   # This dataset is only available in development versions installed from GitHub.
#'   # It is not available in CRAN versions.
#'   # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'   # Visualize trait data
#'   hist(Ponerinae_deepSTRAPP_cont_old_calib_0_40$trait_data_df_over_time$trait_value,
#'      xlab = "Trait values", main = NULL)
#'
#'   # Generate plot
#'   plotTT_continuous <- plot_rates_through_time(
#'      deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'      quantile_ranges = c(0, 0.25, 0.5, 0.75, 1.0),
#'      time_range = c(0, 50), # Control range of the X-axis
#'      # color_scale = c("limegreen", "red"),
#'      plot_CI = TRUE,
#'      CI_type = "quantiles_rect",
#'      CI_quantiles = 0.9,
#'      display_plot = FALSE,
#'      # PDF_file_path = "./plotTT_continuous.pdf",
#'      return_mean_data_per_samples_df = TRUE,
#'      return_median_data_across_samples_df = TRUE)
#'
#'   # Explore output
#'   # str(plotTT_continuous, max.level = 1)
#'
#'   # Plot
#'   print(plotTT_continuous$rates_TT_ggplot)
#'   # Adjust aesthetics of plot a posteriori
#'   plotTT_continuous_adj <- plotTT_continuous$rates_TT_ggplot +
#'       ggplot2::theme(
#'          plot.title = ggplot2::element_text(color = "red", size = 15),
#'          axis.title = ggplot2::element_text(size = 14),
#'          axis.text = ggplot2::element_text(size = 12))
#'   # Plot again
#'   print(plotTT_continuous_adj)
#' }
#'
#'
#' # ------ Example 2: Plot rates through time for categorical data ------ #
#'
#' if (deepSTRAPP:::is_dev_version())
#' {
#'   ## Load results of run_deepSTRAPP_over_time()
#'   data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, package = "deepSTRAPP")
#'   # This dataset is only available in development versions installed from GitHub.
#'   # It is not available in CRAN versions.
#'   # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'   # Explore trait data
#'   table(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40$trait_data_df_over_time$trait_value)
#'
#'   # Set colors to use
#'   colors_per_states <- c("forestgreen", "sienna", "goldenrod")
#'   names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
#'
#'   # Generate plot only for "arboreal" and "terricolous"
#'   plotTT_categorical <- plot_rates_through_time(
#'       deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#'       select_trait_levels = c("arboreal", "terricolous"),
#'       time_range = c(0, 50),
#'       colors_per_levels = colors_per_states,
#'       plot_CI = TRUE,
#'       CI_type = "quantiles_rect",
#'       CI_quantiles = 0.9,
#'       display_plot = FALSE,
#'       # PDF_file_path = "./plotTT_categorical.pdf",
#'       return_mean_data_per_samples_df = TRUE,
#'       return_median_data_across_samples_df = TRUE)
#'
#'   # Explore output
#'   # str(plotTT_categorical, max.level = 1)
#'
#'   # Adjust aesthetics of plot a posteriori
#'   plotTT_categorical_adj <- plotTT_categorical$rates_TT_ggplot +
#'       ggplot2::theme(
#'          plot.title = ggplot2::element_text(size = 15),
#'          axis.title = ggplot2::element_text(size = 14),
#'          axis.text = ggplot2::element_text(size = 12))
#'   print(plotTT_categorical_adj)
#' }
#'
#' # ------ Example 3: Plot rates through time for biogeographic data ------ #
#'
#' if (deepSTRAPP:::is_dev_version())
#' {
#'   ## Load results of run_deepSTRAPP_over_time()
#'   data(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40, package = "deepSTRAPP")
#'   # This dataset is only available in development versions installed from GitHub.
#'   # It is not available in CRAN versions.
#'   # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'   # Explore range data
#'   table(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$trait_data_df_over_time$trait_value)
#'
#'   # Set colors to use
#'   colors_per_ranges <- c("mediumpurple2", "peachpuff2")
#'   names(colors_per_ranges) <- c("N", "O")
#'
#'   plotTT_biogeographic <- plot_rates_through_time(
#'       deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
#'       select_trait_levels = "all",
#'       time_range = c(0, 50),
#'       colors_per_levels = colors_per_ranges,
#'       plot_CI = TRUE,
#'       CI_type = "quantiles_rect",
#'       CI_quantiles = 0.9,
#'       display_plot = FALSE,
#'       # PDF_file_path = "./plotTT_biogeographic.pdf",
#'       return_mean_data_per_samples_df = TRUE,
#'       return_median_data_across_samples_df = TRUE)
#'
#'   # Explore output
#'   # str(plotTT_biogeographic, max.level = 1)
#'
#'   # Adjust aesthetics of plot a posteriori
#'   plotTT_biogeographic_adj <- plotTT_biogeographic$rates_TT_ggplot +
#'       ggplot2::theme(
#'          plot.title = ggplot2::element_text(size = 15),
#'          axis.title = ggplot2::element_text(size = 14),
#'          axis.text = ggplot2::element_text(size = 12))
#'   print(plotTT_biogeographic_adj)
#' }
#'


### Master function to prepare data and select the proper test function according to data type ####

plot_rates_through_time <- function (
    deepSTRAPP_outputs,
    rate_type = "net_diversification",
    quantile_ranges = c(0, 0.25, 0.5, 0.75, 1.0),
    select_trait_levels = "all",
    time_range = NULL,
    color_scale = NULL,
    colors_per_levels = NULL,
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
  {
    ## deepSTRAPP_outputs
    # deepSTRAPP_outputs must have element $trait_data_df_over_time
    if (is.null(deepSTRAPP_outputs$trait_data_df_over_time))
    {
      stop(paste0("'$trait_data_df_over_time' is missing from 'deepSTRAPP_outputs'. You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                  "Especially, check if you used 'extract_trait_data_melted_df = TRUE' to save the summary data.frame of trait values ",
                  "found along branches at each time-step, needed for the RTT plot."))
    }
    # deepSTRAPP_outputs must have element $diversification_data_df_over_time
    if (is.null(deepSTRAPP_outputs$diversification_data_df_over_time))
    {
      stop(paste0("'$diversification_data_df_over_time' is missing from 'deepSTRAPP_outputs'. You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                  "Especially, check if you used 'extract_diversification_data_melted_df = TRUE' to save the summary data.frame of diversification rates ",
                  "found along branches at each time-step, needed for the RTT plot."))
    }

    ## rate_type must be either "speciation", "extinction" or "net_diversification"
    if (!(rate_type %in% c("speciation", "extinction", "net_diversification")))
    {
      stop("'rate_type' can only be 'speciation', 'extinction', or 'net_diversification'.")
    }

    ## time_range
    if (!is.null(time_range))
    {
      # Check that two values are provided for time_range
      if (length(time_range) != 2)
      {
        stop(paste0("'time_range' must be a vector of two positive numerical values providing the time boundaries used for the plot."))
      }
      # Check that time_range is strictly positive
      if (!identical(time_range, abs(time_range)))
      {
        stop(paste0("'time_range' must be strictly positive numerical values providing the time boundaries used for the plot."))
      }
      # Ensure that time_range are properly ordered in increasing values
      time_range <- range(time_range)
      # Check that time_range encompass multiple focal times to be able to draw a line
      focal_times_in_trait_df <- unique(deepSTRAPP_outputs$trait_data_df_over_time$focal_time)
      focal_times_in_diversification_df <- unique(deepSTRAPP_outputs$diversification_data_df_over_time$focal_time)
      shared_focal_times <- intersect(focal_times_in_trait_df, focal_times_in_diversification_df)
      shared_focal_times_in_range <- (shared_focal_times >= time_range[1]) & (shared_focal_times <= time_range[2])
      if (sum(shared_focal_times_in_range) < 2)
      {
        stop(paste0("'time_range' must encompass at least two focal_time recorded in the summary data.frames ",
                    "for trait data ('deepSTRAPP_outputs$trait_data_df_over_time') and diversification rates ('deepSTRAPP_outputs$diversification_data_df_over_time').\n",
                    "Current values of 'time_range' = ", paste(time_range, collapse = ", "), ".\n",
                    "'focal_time' with data recorded for traits and rates are: ", paste(shared_focal_times, collapse = ", "),"."))
      }
    }

    ## CI_type
    # CI_type is either "fuzzy" or "quantiles_rect"
    if (!(CI_type %in% c("fuzzy", "quantiles_rect")))
    {
      stop("'CI_type' can only be 'fuzzy', or 'quantiles_rect'.")
    }

    ## CI_quantiles
    # CI_quantiles should be a numerical between 0 and 1
    if ((CI_quantiles < 0) | (CI_quantiles > 1))
    {
      stop(paste0("'CI_quantiles' reflects the proportion of rate values encompass by the confidence interval.\n",
                  "This is used to display the variance in rates observed across posterior samples. It must be between 0 and 1.\n",
                  "Current value of 'CI_quantiles' is ",CI_quantiles,"."))
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

    ## Other checks are carried in dedicated sub-functions
  }


  ## Detect the type of trait data
  trait_data_type <- deepSTRAPP_outputs$trait_data_type

  ## Compute the appropriate internal function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous =   { # Case for continuous data
           # Need quantile_ranges to define groups of branches per trait values
           plotTT_output <- plot_rates_through_time_for_continuous_data(
             deepSTRAPP_outputs = deepSTRAPP_outputs,
             rate_type = rate_type,
             quantile_ranges = quantile_ranges,
             time_range = time_range,
             color_scale = color_scale,
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
             deepSTRAPP_outputs = deepSTRAPP_outputs,
             rate_type = rate_type,
             select_trait_levels = select_trait_levels,
             time_range = time_range,
             colors_per_levels = colors_per_levels,
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
             deepSTRAPP_outputs = deepSTRAPP_outputs,
             rate_type = rate_type,
             select_trait_levels = select_trait_levels,
             time_range = time_range,
             colors_per_levels = colors_per_levels,
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
    deepSTRAPP_outputs,
    rate_type = "net_diversification",
    quantile_ranges = c(0, 0.25, 0.5, 0.75, 1.0),
    time_range = NULL,
    color_scale = NULL,
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
  {
    ## quantile_ranges
    # Check that quantile_ranges are between 0 and 1
    if (any(quantile_ranges < 0))
    {
      stop(paste0("'quantile_ranges' must be numerical values between [0, 1] providing trait data quantiles to use to aggregate rate values.\n",
                  "Current values are ", paste(quantile_ranges, collapse = ", ")),".")
    }
    if (any(quantile_ranges > 1))
    {
      stop(paste0("'quantile_ranges' must be numerical values between [0, 1] providing trait data quantiles to use to aggregate rate values.\n",
                  "Current values are ", paste(quantile_ranges, collapse = ", ")),".")
    }
    # Ensure that quantile_ranges are properly ordered in increasing values
    quantile_ranges <- quantile_ranges[order(quantile_ranges)]
    # Ensure that quantile_ranges include the c(0,1) boundaries
    initial_quantile_ranges <- quantile_ranges
    if (range(quantile_ranges)[1] != 0)
    {
      cat(paste0("WARNING: 'quantile_ranges' does not start with 0. The lower boundary was added.\n"))
      quantile_ranges <- c(0, quantile_ranges)
    }
    if (range(quantile_ranges)[2] != 1)
    {
      cat(paste0("WARNING: 'quantile_ranges' does not end with 1. The upper boundary was added.\n"))
      quantile_ranges <- c(quantile_ranges, 1)
    }
    if (!(identical(range(initial_quantile_ranges), c(0, 1))))
    {
      cat(paste0("WARNING: New 'quantile_ranges' are ", paste(quantile_ranges, collapse = ", "),".\n"))
    }

    ## color_scale
    # Check whether all colors are valid
    if (!is.null(color_scale))
    {
      if (!all(is_color(color_scale)))
      {
        invalid_colors <- color_scale[!is_color(color_scale)]
        stop(paste0("Some color names in 'color_scale' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }

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
    x = deepSTRAPP_outputs$diversification_data_df_over_time,
    y = deepSTRAPP_outputs$trait_data_df_over_time,
    by = dplyr::join_by(focal_time, tip_ID))

  ## Filter data for selected rate_type
  if (rate_type == "speciation") { rate_type <- "lambda" }
  if (rate_type == "extinction") { rate_type <- "mu" }
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

  ## Prepare colors to use for quantile groups
  nb_groups <- length(levels(as.factor(median_data_across_samples_df$quantile_ranges)))
  if (!is.null(color_scale))
  {
    # Use the provided color to build the color palette
    col_fn <- grDevices::colorRampPalette(colors = color_scale)
    colors_per_groups <- col_fn(n = nb_groups)
  } else {
    # Default: use the 'Spectral' palette from RColorBrewer
    colors_per_groups <- rev(RColorBrewer::brewer.pal(name = "Spectral", n = nb_groups))
  }
  names(colors_per_groups) <- levels(as.factor(median_data_across_samples_df$quantile_ranges))

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
      # ggplot2::scale_color_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +
      ggplot2::scale_color_manual(name = "Trait quantile groups", values = colors_per_groups) +

      # Adjust aesthetics
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
        panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
        panel.background = ggplot2::element_rect(fill = NA, color = NA),
        plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                           margin = ggplot2::margin(b = 10, t = 5)),
        legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
        legend.position = "right",
        # legend.position = "inside",
        # legend.position.inside = c(0.15, 0.2),
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
        # ggplot2::scale_color_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +
        ggplot2::scale_color_manual(name = "Trait quantile groups", values = colors_per_groups) +

        # Adjust aesthetics
        ggplot2::theme(
          plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
          panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
          panel.background = ggplot2::element_rect(fill = NA, color = NA),
          plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                             margin = ggplot2::margin(b = 10, t = 5)),
          legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
          legend.position = "right",
          # legend.position = "inside",
          # legend.position.inside = c(0.15, 0.2),
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
        # ggplot2::scale_fill_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +
        ggplot2::scale_fill_manual(name = "Trait quantile groups", values = colors_per_groups) +

        # Adjust color scheme and legend
        # ggplot2::scale_color_brewer(name = "Trait quantile groups", palette = "Spectral", direction = -1) +
        ggplot2::scale_color_manual(name = "Trait quantile groups", values = colors_per_groups) +

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
          legend.position = "right",
          # legend.position = "inside",
          # legend.position.inside = c(0.15, 0.2),
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
                       base_height = 8, base_width = 14)
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
    deepSTRAPP_outputs,
    rate_type = "net_diversification",
    select_trait_levels = "all",
    time_range = NULL,
    colors_per_levels = NULL,
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
  {
    ## Extract state levels
    states_in_trait_df <- unique(deepSTRAPP_outputs$trait_data_df_over_time$trait_value)
    states_in_trait_df <- states_in_trait_df[order(states_in_trait_df)]

    ## select_trait_levels
    if (!any(select_trait_levels == "all"))
    {
      # Check that select_trait_levels are all found in the summary data.frame $trait_data_df_over_time

      if (!all(select_trait_levels %in% states_in_trait_df))
      {
        stop(paste0("Some states listed in 'select_trait_levels' are not found in the summary data.frame for trait data ('deepSTRAPP_outputs$trait_data_df_over_time').\n",
                    "'select_trait_levels' = ",paste(select_trait_levels[order(select_trait_levels)], collapse = ", "),".\n",
                    "Observed states in trait data = ", paste(states_in_trait_df, collapse = ", ")),".")
      }
    }

    # Update list of states to keep only the selected ones
    if (!any(select_trait_levels == "all"))
    {
      states_in_trait_df <- select_trait_levels
    }

    ## colors_per_levels
    # Check whether all colors are valid
    if (!is.null(colors_per_levels))
    {
      # Check that the color match the selected states
      if (!all(states_in_trait_df %in% names(colors_per_levels)))
      {
        missing_states <- states_in_trait_df[!(states_in_trait_df %in% names(colors_per_levels))]
        stop(paste0("Not all selected states are found in 'colors_per_levels'.\n",
                    "Missing states: ", paste(missing_states, collapse = ", "), "."))
      }
      if (!all(is_color(colors_per_levels)))
      {
        invalid_colors <- colors_per_levels[!is_color(colors_per_levels)]
        stop(paste0("Some color names in 'colors_per_levels' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }

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
    x = deepSTRAPP_outputs$diversification_data_df_over_time,
    y = deepSTRAPP_outputs$trait_data_df_over_time,
    by = dplyr::join_by(focal_time, tip_ID))

  ## Filter data for selected rate_type
  if (rate_type == "speciation") { rate_type <- "lambda" }
  if (rate_type == "extinction") { rate_type <- "mu" }
  data_per_samples_df <- data_per_samples_df[data_per_samples_df$rate_type == rate_type, ]

  ## Filter data for selected states
  if (!("all" %in% select_trait_levels))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$trait_value %in% select_trait_levels, ]
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
    stop("No data found in the time range c(",time_range[1],", ", time_range[2],") for ",paste(select_trait_levels, collapse = ", ")," states.\n")
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

  ## Prepare colors_per_levels to use in plots
  if (is.null(colors_per_levels))
  {
    nb_groups <- length(levels(as.factor(median_data_across_samples_df$trait_value)))
    # Default: use the default ggplot palette from scales
    col_fn <- scales::hue_pal()
    colors_per_levels <- col_fn(n = nb_groups)
    names(colors_per_levels) <- levels(as.factor(median_data_across_samples_df$trait_value))
  }

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
      # ggplot2::scale_color_discrete(name = "States") +
      ggplot2::scale_color_manual(name = "States", values = colors_per_levels) +

      # Adjust aesthetics
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
        panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
        panel.background = ggplot2::element_rect(fill = NA, color = NA),
        plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                           margin = ggplot2::margin(b = 10, t = 5)),
        legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
        legend.position = "right",
        # legend.position = "inside",
        # legend.position.inside = c(0.15, 0.2),
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
        # ggplot2::scale_color_discrete(name = "States") +
        ggplot2::scale_color_manual(name = "States", values = colors_per_levels) +

        # Adjust aesthetics
        ggplot2::theme(
          plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
          panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
          panel.background = ggplot2::element_rect(fill = NA, color = NA),
          plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                             margin = ggplot2::margin(b = 10, t = 5)),
          legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
          legend.position = "right",
          # legend.position = "inside",
          # legend.position.inside = c(0.15, 0.2),
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
        # ggplot2::scale_fill_discrete("States") +
        ggplot2::scale_fill_manual(name = "States", values = colors_per_levels) +

        # Adjust color scheme and legend
        # ggplot2::scale_color_discrete(name = "States") +
        ggplot2::scale_color_manual(name = "States", values = colors_per_levels) +

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
          legend.position = "right",
          # legend.position = "inside",
          # legend.position.inside = c(0.15, 0.2),
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
                       base_height = 8, base_width = 14)
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
    deepSTRAPP_outputs,
    rate_type = "net_diversification",
    select_trait_levels = "all",
    time_range = NULL,
    colors_per_levels = NULL,
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
  {
    ## Extract range levels
    ranges_in_trait_df <- unique(deepSTRAPP_outputs$trait_data_df_over_time$trait_value)
    ranges_in_trait_df <- ranges_in_trait_df[order(ranges_in_trait_df)]

    ## select_trait_levels
    if (!any(select_trait_levels == "all"))
    {
      # Check that select_trait_levels are all found in the summary data.frame $trait_data_df_over_time
      if (!all(select_trait_levels %in% ranges_in_trait_df))
      {
        stop(paste0("Some ranges listed in 'select_trait_levels' are not found in the summary data.frame for trait data ('deepSTRAPP_outputs$trait_data_df_over_time').\n",
                    "'select_trait_levels' = ",paste(select_trait_levels[order(select_trait_levels)], collapse = ", "),".\n",
                    "Observed ranges in trait data = ", paste(ranges_in_trait_df, collapse = ", ")),".")
      }
    }

    # Update list of ranges to keep only the selected ones
    if (!any(select_trait_levels == "all"))
    {
      ranges_in_trait_df <- select_trait_levels
    }

    ## colors_per_levels
    # Check whether all colors are valid
    if (!is.null(colors_per_levels))
    {
      # Check that the color match the states
      if (!all(ranges_in_trait_df %in% names(colors_per_levels)))
      {
        missing_ranges <- ranges_in_trait_df[!(ranges_in_trait_df %in% names(colors_per_levels))]
        stop(paste0("Not all selected ranges are found in 'colors_per_levels'.\n",
                    "Missing ranges: ", paste(missing_ranges, collapse = ", "), "."))
      }
      if (!all(is_color(colors_per_levels)))
      {
        invalid_colors <- colors_per_levels[!is_color(colors_per_levels)]
        stop(paste0("Some color names in 'colors_per_levels' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }


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
    x = deepSTRAPP_outputs$diversification_data_df_over_time,
    y = deepSTRAPP_outputs$trait_data_df_over_time,
    by = dplyr::join_by(focal_time, tip_ID))

  ## Filter data for selected rate_type
  if (rate_type == "speciation") { rate_type <- "lambda" }
  if (rate_type == "extinction") { rate_type <- "mu" }
  data_per_samples_df <- data_per_samples_df[data_per_samples_df$rate_type == rate_type, ]

  ## Filter data for selected states/ranges
  if (!("all" %in% select_trait_levels))
  {
    data_per_samples_df <- data_per_samples_df[data_per_samples_df$trait_value %in% select_trait_levels, ]
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
    stop("No data found in the time range c(",time_range[1],", ", time_range[2],") in ",paste(select_trait_levels, collapse = ", ")," ranges.\n")
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

  ## Prepare colors_per_levels to use in plots
  if (is.null(colors_per_levels))
  {
    nb_groups <- length(levels(as.factor(median_data_across_samples_df$trait_value)))
    # Default: use the default ggplot palette from scales
    col_fn <- scales::hue_pal()
    colors_per_levels <- col_fn(n = nb_groups)
    names(colors_per_levels) <- levels(as.factor(median_data_across_samples_df$trait_value))
  }

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
      # ggplot2::scale_color_discrete(name = "Ranges") +
      ggplot2::scale_color_manual(name = "Ranges", values = colors_per_levels) +

      # Adjust aesthetics
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
        panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
        panel.background = ggplot2::element_rect(fill = NA, color = NA),
        plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                           margin = ggplot2::margin(b = 10, t = 5)),
        legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
        legend.position = "right",
        # legend.position = "inside",
        # legend.position.inside = c(0.15, 0.2),
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
        # ggplot2::scale_color_discrete(name = "Ranges") +
        ggplot2::scale_color_manual(name = "Ranges", values = colors_per_levels) +

        # Adjust aesthetics
        ggplot2::theme(
          plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
          panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
          panel.background = ggplot2::element_rect(fill = NA, color = NA),
          plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                             margin = ggplot2::margin(b = 10, t = 5)),
          legend.title = ggplot2::element_text(size  = 16, margin = ggplot2::margin(b = 5)),
          legend.position = "right",
          # legend.position = "inside",
          # legend.position.inside = c(0.15, 0.2),
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
        # ggplot2::scale_fill_discrete("Ranges") +
        ggplot2::scale_fill_manual(name = "Ranges", values = colors_per_levels) +

        # Adjust color scheme and legend
        # ggplot2::scale_color_discrete(name = "Ranges") +
        ggplot2::scale_color_manual(name = "Ranges", values = colors_per_levels) +

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
          legend.position = "right",
          # legend.position = "inside",
          # legend.position.inside = c(0.15, 0.2),
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
                       base_height = 8, base_width = 14)
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

