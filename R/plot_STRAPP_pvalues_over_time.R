
#' @title Plot evolution of p-values of STRAPP tests over time
#'
#' @description Plot the evolution of the p-values of STRAPP tests
#'   carried out for across multiple `time_steps`, obtained from
#'   [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'   By default, return a plot with a single line for p-values of overall tests.
#'   If `plot_posthoc_tests = TRUE`, it will return a plot with multiple lines, one per pair in post hoc tests
#'   (only for multinominal data, with more than two states).
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plot will be saved directly in a PDF file.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   that summarize the results of multiple deepSTRAPP across `$time_steps`.
#' @param time_range Vector of two numerical values. Time boundaries used for X-axis the plot.
#'   If `NULL` (the default), the range of data provided in `deepSTRAPP_outputs` will be used.
#' @param pvalues_max Numerical. Set the max boundary used for the Y-axis of the plot.
#'   If `NULL` (the default), the maximum p-value provided in `deepSTRAPP_outputs` will be used.
#' @param alpha Numerical. Significance level to display as a red dashed line on the plot. If set to `NULL`, no line will be added. Default is `0.05`.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param plot_significant_time_frame Logical. Whether to display a green band over the time frame that yields significant results according to the chosen alpha level. Default is `TRUE`.
#' @param plot_posthoc_tests Logical. For multinominal data only. Whether to plot the p-values for the overall Kruskal-Wallis test across all states (`plot_posthoc_tests = FALSE`),
#'   or plot the p-values for the pairwise post hoc Dunn's test across pairs of states (`plot_posthoc_tests = TRUE`). Default is `FALSE`.
#'   This is only possible if `deepSTRAPP_outputs` contains the `$pvalues_summary_df_for_posthoc_pairwise_tests` element returned by
#'   [deepSTRAPP::run_deepSTRAPP_over_time()] when `posthoc_pairwise_tests = TRUE`.
#' @param select_posthoc_pairs Vector of character strings used to specify the pairs to include in the plot. Names of pairs must match the pairs found in
#'   `deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests$pair`. Default is "all" to include all pairs.
#' @param plot_adjusted_pvalues Logical. Whether to display the p-values adjusted for multiple testing rather than the raw p-values.
#'   See argument 'p.adjust_method' in [deepSTRAPP::run_deepSTRAPP_for_focal_time()] or [deepSTRAPP::run_deepSTRAPP_over_time()]. Default is `FALSE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#'
#' @export
#' @importFrom ggplot2 ggplot geom_line aes geom_hline scale_y_continuous scale_x_continuous scale_color_discrete xlab ylab ggtitle theme element_line element_rect element_text unit margin
#' @importFrom cowplot save_plot
#' @importFrom stats complete.cases
#'
#' @details Plots are build based on the p-values recorded in summary_df provided by [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'   For overall tests, those p-values are found in `$pvalues_summary_df`.
#'
#'   For multinominal data (categorical or biogeographic data with more than 2 states), it is possible to plot p-values of post hoc pairwise tests.
#'   Set `plot_posthoc_tests = TRUE` to generate plots for the pairwise post hoc Dunn's test across pairs of states.
#'   To achieve this, the `deepSTRAPP_outputs` input object must contain a `$pvalues_summary_df_for_posthoc_pairwise_tests` element that summarizes p-values
#'   computed across pairs of states for all post hoc tests. This is obtained from [deepSTRAPP::run_deepSTRAPP_over_time()] when setting
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
#' @seealso [deepSTRAPP::run_deepSTRAPP_over_time()]
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'   ## Load results of run_deepSTRAPP_over_time() for categorical data with 3-levels
#'   Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- readRDS(system.file("extdata",
#'      "Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", package = "deepSTRAPP"))
#'   ## This dataset is only available in development versions installed from GitHub.
#'   # It is not available in CRAN versions.
#'   # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'   ## Plot results of overall Kruskal-Wallis / Mann-Whitney-Wilcoxon tests across all time-steps
#'   plot_overall <- plot_STRAPP_pvalues_over_time(
#'      deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#'      alpha = 0.10,
#'      time_range = c(0, 30), # Adjust time range if needed
#'      display_plot = FALSE)
#'
#'   # Adjust aesthetics a posteriori
#'   plot_overall <- plot_overall +
#'      ggplot2::theme(
#'         plot.title = ggplot2::element_text(size = 16))
#'
#'   print(plot_overall)
#'
#'   ## Plot results of post hoc pairwise Dunn's tests between selected pairs of states
#'   plot_posthoc <- plot_STRAPP_pvalues_over_time(
#'      deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#'      alpha = 0.10,
#'      plot_posthoc_tests = TRUE,
#'      # PDF_file_path = "./pvalues_over_time.pdf",
#'      select_posthoc_pairs = c("arboreal != subterranean",
#'                               "arboreal != terricolous"),
#'      display_plot = FALSE)
#'
#'   # Adjust aesthetics a posteriori
#'   plot_posthoc <- plot_posthoc +
#'      ggplot2::theme(
#'         plot.title = ggplot2::element_text(size = 16),
#'         legend.title = ggplot2::element_text(size  = 14),
#'         legend.position.inside = c(0.25, 0.25))
#'
#'   print(plot_posthoc)
#' }
#'


plot_STRAPP_pvalues_over_time <-  function (
    deepSTRAPP_outputs,
    time_range = NULL,
    pvalues_max = NULL,
    alpha = 0.05,
    display_plot = TRUE,
    plot_significant_time_frame = TRUE,
    plot_posthoc_tests = FALSE,
    select_posthoc_pairs = "all",
    plot_adjusted_pvalues = FALSE,
    PDF_file_path = NULL
)
{
  ### Check input validity
  {
    ## deepSTRAPP_outputs
    # deepSTRAPP_outputs must have element $pvalues_summary_df
    if (is.null(deepSTRAPP_outputs$pvalues_summary_df))
    {
      stop(paste0("'$pvalues_summary_df' is missing from 'deepSTRAPP_outputs'. You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects."))
    }

    ## time_range
    if (!is.null(time_range))
    {
      # Check that two values are provided for time_range
      if (length(time_range) != 2)
      {
        stop(paste0("'time_range' must be a vector of two positive numerical values providing the time boundaries of the X-axis used for the plot."))
      }
      # Check that time_range is strictly positive
      if (!identical(time_range, abs(time_range)))
      {
        stop(paste0("'time_range' must be strictly positive numerical values providing the time boundaries of the X-axis used for the plot."))
      }
      # Ensure that time_range are properly ordered in increasing values
      time_range <- range(time_range)
    } else {
      # Extract time range from data if not provided
      time_range <- range(deepSTRAPP_outputs$pvalues_summary_df$focal_time, na.rm = TRUE)
    }
    # Check that time_range encompass multiple focal-time with recorded p-values to be able to draw a line
    pvalues_summary_df_no_NA <- deepSTRAPP_outputs$pvalues_summary_df[stats::complete.cases(deepSTRAPP_outputs$pvalues_summary_df), ]
    focal_times_in_pvalues_df <- unique(pvalues_summary_df_no_NA$focal_time)
    focal_times_in_range <- (focal_times_in_pvalues_df >= time_range[1]) & (focal_times_in_pvalues_df <= time_range[2])
    if (sum(focal_times_in_range) < 2)
    {
      stop(paste0("'time_range' must encompass at least two focal_time with p-value recorded in 'deepSTRAPP_outputs$pvalues_summary_df'.\n",
                  "Current values of 'time_range' = ", paste(time_range, collapse = ", "), ".\n",
                  "'focal_time' with p-values recorded are: ", paste(focal_times_in_pvalues_df, collapse = ", "),"."))
    }

    ## pvalues_max
    if (!is.null(pvalues_max))
    {
      # Check that pvalues_max is comprised within [0,1]
      if ((pvalues_max < 0) | (pvalues_max > 1))
      {
        stop(paste0("'pvalues_max' must be strictly a positive numerical value between 0 and 1 providing the maximum boundary of the Y-axis used for the plot."))
      }
    }

    ## alpha
    # alpha must be set between 0 and 1, or NULL.
    if (!is.null(alpha))
    {
      if ((alpha < 0) | (alpha > 1))
      {
        stop(paste0("'alpha' is the proportion of type I error tolerated to assess significance of the test. It must be between 0 and 1.\n",
                    "Current value of 'alpha' is ",alpha,"."))
      }
    }

    ## plot_significant_time_frame
    if ((plot_significant_time_frame) & is.null(alpha))
    {
      stop(paste0("You need to provide a value for 'alpha' if you wish to plot the significant time frame ('plot_significant_time_frame' = TRUE)."))
    }

    ## plot_posthoc_tests
    if (plot_posthoc_tests)
    {
      # Extract $trait_data_type_for_stats from deepSTRAPP_outputs$STRAPP_results_over_time
      trait_data_type_for_stats <- deepSTRAPP_outputs$trait_data_type_for_stats
      # plot_posthoc_tests = TRUE only for "multinominal" data
      if (trait_data_type_for_stats != "multinominal")
      {
        stop(paste0("'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                    "Set 'posthoc_pairwise_tests = FALSE', or provide 'deepSTRAPP_outputs' for a trait with more than two states/ranges'.\n",
                    "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects."))
      }
      # Check if $pvalues_summary_df_for_posthoc_pairwise_tests is present in deepSTRAPP_outputs.
      if (is.null(deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests))
      {
        stop(paste0("'$pvalues_summary_df_for_posthoc_pairwise_tests' is missing from 'deepSTRAPP_outputs'. You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                    "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                    "Especially, check if you used 'posthoc_pairwise_tests = TRUE' to run post hoc tests."))
      }

      # Extract posthoc pairwise tests with recorded p-values
      pvalues_summary_df_for_posthoc_no_NA <- deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests[stats::complete.cases(deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests), ]

      ## select_posthoc_pairs
      # Check that select_posthoc_pairs is "all" or match the pairs in deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests$pair
      if (!("all" %in% select_posthoc_pairs))
      {
        available_pairs <- unique(pvalues_summary_df_for_posthoc_no_NA$pair)
        available_pairs <- available_pairs[order(available_pairs)]
        if (!all(select_posthoc_pairs %in% available_pairs))
        {
          stop(paste0("Some pairs of states/ranges listed in 'select_posthoc_pairs' are not found in the summary data.frame for posthoc tests ('deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests').\n",
                      "'select_posthoc_pairs' = ",paste(select_posthoc_pairs[order(select_posthoc_pairs)], collapse = ", "),".\n",
                      "Observed pairs of states/ranges with recorded p-values = ", paste(available_pairs, collapse = ", ")),".")
        }
      }

      ## time_range
      # Check that time_range encompass multiple focal-time with recorded p-values per pair (!) to be able to draw a line
      # Split pvalues_df per pairs
      pvalues_summary_df_per_pairs <- split(x = pvalues_summary_df_for_posthoc_no_NA, f = pvalues_summary_df_for_posthoc_no_NA$pair)
      # Extract focal time per pairs
      focal_times_per_pairs <- lapply(X = pvalues_summary_df_per_pairs, FUN = function (x) { unique(x$focal_time) } )
      # Detect focal time within range, per pairs
      focal_times_in_range_per_pairs <- lapply(X = focal_times_per_pairs, FUN = function (x) { (x >= time_range[1]) & (x <= time_range[2]) } )
      focal_times_in_range_per_pairs_counts <- unlist(lapply(X = focal_times_in_range_per_pairs, FUN = sum))
      if (any(focal_times_in_range_per_pairs_counts < 2))
      {
        pairs_with_issue <- names(focal_times_in_range_per_pairs_counts)[focal_times_in_range_per_pairs_counts < 2]
        focal_times_for_pair_with_issue <- focal_times_per_pairs[[which.max(focal_times_in_range_per_pairs_counts < 2)]]
        stop(paste0("'time_range' must encompass at least two focal_time with recorded p-value for each pair in 'deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests'.\n",
                    "Current values of 'time_range' = ", paste(time_range, collapse = ", "), ".\n",
                    "'focal_time' with p-values recorded for pair '",pairs_with_issue[1],"' are: ", paste(focal_times_for_pair_with_issue, collapse = ", "),"."))
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

  ## Save the graphical parameters on entry and, on exit, restore only those this function actually changed.

  # List graphical parameters describing the plot that now exists (Never to restore as they describe the current plot state)
  plot_state_par_names <- c("usr", "plt", "fig", "mfg", "new", "xaxp", "yaxp", "pin")
  # Save current graphic parameters
  entry_par <- par(no.readonly = TRUE)

  # On exit, compare changes in par and restore the ones that have change, but the plot_state_par
  on.exit({
    # List par on exit
    exit_par <- par(no.readonly = TRUE)
    # Identify differences with initial parameters
    changed_par_names <- names(entry_par)[!vapply(X = names(entry_par),
                                                  FUN = function (x) { isTRUE(all.equal(entry_par[[x]], exit_par[[x]])) },
                                                  FUN.VALUE = logical(1))]
    # Exclude plot state parameters from the restore
    changed_par_names <- setdiff(changed_par_names, plot_state_par_names)
    # Restore parameters if any has changed
    if (length(changed_par_names) > 0)
    {
      par(entry_par[changed_par_names])
    }

    # Force new to FALSE, so the next plot calling plot.new() start on a new clean window,
    # But only if we actually plot something
    if (isTRUE(display_plot) || (isTRUE(exit_par$new) && !isTRUE(entry_par$new)))
    {
      par(new = FALSE)
    }
  }, add = TRUE)


  ## Remove "_" from $rate_type
  rate_type <- gsub(pattern = "_", replacement = " ", x = deepSTRAPP_outputs$rate_type)

  if (!plot_posthoc_tests) ## Case for overall test plot
  {
    # Extract summary df
    pvalues_summary_df <- deepSTRAPP_outputs$pvalues_summary_df

    # Remove NA
    pvalues_summary_df <- pvalues_summary_df[stats::complete.cases(pvalues_summary_df), ]

    # Extract data for the selected time range
    pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$focal_time <= time_range[2], ]
    pvalues_summary_df <- pvalues_summary_df[pvalues_summary_df$focal_time >= time_range[1], ]

    # Extract data to avoid 'binding warning'
    p_value <- pvalues_summary_df$p_value
    focal_time <- time <- pvalues_summary_df$focal_time
    y_axis <- poly_ID <- NA

    # Set pvalues_max
    if (is.null(pvalues_max))
    {
      pvalues_max <- max(p_value)
    }

    # Prepare polygons for significant timeframes if needed
    if (plot_significant_time_frame)
    {
      # Filter time to avoid failed polygons when plotting
      pvalues_summary_df_filtered <- pvalues_summary_df |>
        dplyr::filter(focal_time >= time_range[1]) |>
        dplyr::filter(focal_time <= time_range[2])

      ## Interpolate p-values to thinner time scale
      x_increment <- (time_range[2] - time_range[1])/100
      pvalues_df_thin <- data.frame()
      for (i in 1:(nrow(pvalues_summary_df_filtered)-1))
      {
        # i <- 1

        # Extract coarse time boundaries
        start_time_i <- pvalues_summary_df_filtered$focal_time[i]
        end_time_i <- pvalues_summary_df_filtered$focal_time[i+1]

        # Extract coarse time p-values
        start_pval_i <- pvalues_summary_df_filtered$p_value[i]
        end_pval_i <- pvalues_summary_df_filtered$p_value[i+1]

        # Generate thinner time steps
        time_scale_thinner_i <- seq(from = start_time_i,
                                    to = end_time_i,
                                    by = x_increment)
        # If any NA, fill with NA
        if (any(is.na(c(start_pval_i, end_pval_i))))
        {
          pval_thinner_i <- rep(NA, length(time_scale_thinner_i))
        } else {
          # Interpolate p-values
          pval_increment_i <- (end_pval_i - start_pval_i) / (length(time_scale_thinner_i) - 1)
          pval_thinner_i <- seq(from = start_pval_i,
                                to = end_pval_i,
                                by = pval_increment_i)
        }

        # Store results
        pvalues_df_thin_i <- data.frame(time = time_scale_thinner_i,
                                        p_value = pval_thinner_i)

        pvalues_df_thin <- rbind(pvalues_df_thin, pvalues_df_thin_i)
      }
      # Remove duplicates
      pvalues_df_thin <- pvalues_df_thin |>
        dplyr::distinct(time, .keep_all = TRUE)

      # Add last value
      last_value_ID <- nrow(pvalues_summary_df_filtered)
      last_values <- c(pvalues_summary_df_filtered$focal_time[last_value_ID], pvalues_summary_df_filtered$p_value[last_value_ID])
      pvalues_df_thin <- rbind(pvalues_df_thin, last_values)

      ## Prepare data for significance gradient
      # One polygon per x-time
      nb_poly <- 0
      time_x_data <- c()

      # Loop per time
      for (i in 1:(nrow(pvalues_df_thin-1)))
      {
        time_x <- pvalues_df_thin$time[i]
        time_x2 <- pvalues_df_thin$time[i+1]

        time_x_data_i <- c(time_x, time_x2, time_x2, time_x)
        time_x_data <- c(time_x_data, time_x_data_i)

        nb_poly <- nb_poly + 1
      }

      significance_area_poly_df <- data.frame(time = time_x_data,
                                              y_axis = rep(x = c(pvalues_max, pvalues_max, 0, 0), times = nb_poly),
                                              poly_ID = rep(1:nb_poly, each = 4))
      # Join p-values
      significance_area_poly_df <- left_join(x = significance_area_poly_df,
                                             y = pvalues_df_thin,
                                             by = join_by(time))

      # Remove polygons with p-value > alpha
      significance_area_poly_df <- significance_area_poly_df |>
        filter(!(p_value > alpha))
    } else {
      significance_area_poly_df <- data.frame(time = numeric(), y_axis = numeric(), pvalue = numeric())
    }

    # Build ggplot
    pvalues_plot <- ggplot2::ggplot(data = pvalues_summary_df,
                                    mapping = ggplot2::aes(y = p_value, x = focal_time)) +

      # Plot significance area
      geom_polygon(data = significance_area_poly_df,
                   mapping = aes(y = y_axis, x = time,
                                 group = poly_ID),
                   fill = "limegreen", col = NA,
                   alpha = 0.7,
                   linewidth = 1.0, show.legend = FALSE) +

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
        limits = rev(c(0, pvalues_max)) # Set limits
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
    pvalues_summary_df <- deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests

    # Remove NA
    pvalues_summary_df <- pvalues_summary_df[stats::complete.cases(pvalues_summary_df), ]

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

    # Replace p-values with adjusted p-values if needed
    if (plot_adjusted_pvalues)
    {
      pvalues_summary_df$p_value <- pvalues_summary_df$p_value_adjusted
    }

    # Extract data to avoid 'binding warning'
    p_value <- pvalues_summary_df$p_value
    focal_time <- time <- pvalues_summary_df$focal_time
    pair <- pvalues_summary_df$pair
    y_axis <- poly_ID <- NA

    # Set pvalues_max
    if (is.null(pvalues_max))
    {
      pvalues_max <- max(p_value)
    }

    # Prepare polygons for significant timeframes if needed
    if (plot_significant_time_frame)
    {
      # Filter time to avoid failed polygons when plotting
      pvalues_summary_df_filtered <- pvalues_summary_df |>
        dplyr::filter(focal_time >= time_range[1]) |>
        dplyr::filter(focal_time <= time_range[2])

      ## Interpolate p-values to thinner time scale
      x_increment <- (time_range[2] - time_range[1])/100
      pvalues_df_thin <- data.frame()
      focal_time_list <- unique(pvalues_summary_df_filtered$focal_time)
      for (i in 1:(length(focal_time_list)-1))
      {
        # i <- 1

        # Extract coarse time boundaries
        start_time_i <- focal_time_list[i]
        end_time_i <- focal_time_list[i+1]

        # Generate thinner time steps
        time_scale_thinner_i <- seq(from = start_time_i,
                                    to = end_time_i,
                                    by = x_increment)

        ## Loop per pairs
        pvalues_df_thin_i <- data.frame()
        pairs_list <- unique(pvalues_summary_df_filtered$pair)
        for (j in seq_along(pairs_list))
        {
          # j <- 1

          # Extract pair
          pair_j <- pairs_list[j]
          # Extract p-values for pair j
          pvalues_summary_df_filtered_j <- pvalues_summary_df_filtered[pvalues_summary_df_filtered$pair == pair_j, ]

          # Extract coarse time p-values
          start_pval_ij <- pvalues_summary_df_filtered_j$p_value[i]
          end_pval_ij <- pvalues_summary_df_filtered_j$p_value[i+1]

          # If any NA, fill with NA
          if (any(is.na(c(start_pval_ij, end_pval_ij))))
          {
            pval_thinner_ij <- rep(NA, length(time_scale_thinner_i))
          } else {
            # Interpolate p-values
            pval_increment_ij <- (end_pval_ij - start_pval_ij) / (length(time_scale_thinner_i) - 1)
            pval_thinner_ij <- seq(from = start_pval_ij,
                                   to = end_pval_ij,
                                   by = pval_increment_ij)
          }

          # Store results
          pvalues_df_thin_ij <- data.frame(time = time_scale_thinner_i,
                                           p_value = pval_thinner_ij,
                                           pair = pair_j)

          pvalues_df_thin_i <- rbind(pvalues_df_thin_i, pvalues_df_thin_ij)
        }
        # Keep only the highest value at each time
        pvalues_df_thin_i <- pvalues_df_thin_i |>
          dplyr::group_by(time) |>
          dplyr::arrange(time, p_value) |>
          dplyr::mutate(rank = dplyr::row_number()) |>
          dplyr::ungroup() |>
          dplyr::filter(rank == 1) |>
          dplyr::select(time, p_value)

        # Store results
        pvalues_df_thin <- rbind(pvalues_df_thin, pvalues_df_thin_i)
      }
      # Remove duplicates
      pvalues_df_thin <- pvalues_df_thin |>
        dplyr::distinct(time, .keep_all = TRUE)

      # Add last value
      last_step_ID <- which(pvalues_summary_df_filtered$focal_time == focal_time_list[length(focal_time_list)])
      last_p_value <- min(pvalues_summary_df_filtered$p_value[last_step_ID])
      last_step_values <- c(focal_time_list[length(focal_time_list)], last_p_value)
      pvalues_df_thin <- rbind(pvalues_df_thin, last_step_values)

      ## Prepare data for significance gradient
      # One polygon per x-time
      nb_poly <- 0
      time_x_data <- c()

      # Loop per time
      for (i in 1:(nrow(pvalues_df_thin-1)))
      {
        time_x <- pvalues_df_thin$time[i]
        time_x2 <- pvalues_df_thin$time[i+1]

        time_x_data_i <- c(time_x, time_x2, time_x2, time_x)
        time_x_data <- c(time_x_data, time_x_data_i)

        nb_poly <- nb_poly + 1
      }

      significance_area_poly_df <- data.frame(time = time_x_data,
                                              y_axis = rep(x = c(pvalues_max, pvalues_max, 0, 0), times = nb_poly),
                                              poly_ID = rep(1:nb_poly, each = 4))
      # Join p-values
      significance_area_poly_df <- left_join(x = significance_area_poly_df,
                                             y = pvalues_df_thin,
                                             by = join_by(time))

      # Remove polygons with p-value > alpha
      significance_area_poly_df <- significance_area_poly_df |>
        filter(!(p_value > alpha))
    } else {
      significance_area_poly_df <- data.frame(time = numeric(), y_axis = numeric(), pvalue = numeric(), poly_ID = numeric())
    }

    # Build ggplot
    pvalues_plot <- ggplot2::ggplot(data = pvalues_summary_df,
                                    mapping = ggplot2::aes(y = p_value, x = focal_time, color = pair)) +

      # Plot significance area
      geom_polygon(data = significance_area_poly_df,
                   mapping = aes(y = y_axis, x = time,
                                 group = poly_ID),
                   fill = "limegreen", col = NA,
                   alpha = 0.7,
                   linewidth = 1.0, show.legend = FALSE) +

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
        limits = rev(c(0, pvalues_max)) # Set limits
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


