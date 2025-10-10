
#' @title Plot rates vs. trait data for a given focal time
#'
#' @description Plot rates vs. trait data as extracted for a given focal time.
#'   Data are extracted from the output of a deepSTRAPP run carried out with
#'   [deepSTRAPP::run_deepSTRAPP_for_focal_time()] or
#'   [deepSTRAPP::run_deepSTRAPP_over_time()]).
#'
#'   Returns a single plot showing rates vs. trait data for a given focal time.
#'   If the trait data are 'continuous', the plot is a scatter plot.
#'   If the trait data are 'categorical' or 'biogeographic', the plot is a boxplot.
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plot will be saved directly in a PDF file.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_for_focal_time()],
#'   that summarize the results of a STRAPP test for a specific time in the past (i.e. the `focal_time`).
#'   `deepSTRAPP_outputs` can also be extracted from the output of [deepSTRAPP::run_deepSTRAPP_over_time()] that
#'   runs the whole deepSTRAPP workflow over multiple time-steps.
#' @param focal_time Numerical. (Optional) If `deepSTRAPP_outputs` comprises results over multiple time-steps
#'   (i.e., output of [deepSTRAPP::run_deepSTRAPP_over_time()], this is the time of the STRAPP test targeted for plotting.
#' @param rate_type A character string specifying the type of diversification rates to plot.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#'   Even if the `deepSTRAPP_outputs` object was generated with [deepSTRAPP::run_deepSTRAPP_over_time()]
#'   for testing another type of rates, the object will contain data for all types of rates.
#' @param select_trait_levels (Vector of) character string. Only for categorical and biogeographic trait data.
#'  To provide a list of a subset of states/ranges to plot. Names must match the ones found in the `deepSTRAPP_outputs`.
#'  Default is `all` which means all states/ranges will be plotted.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   to display the points. Color scale from lowest values to highest rate values. Only for continuous data.
#'   Default = `NULL` will use the 'Spectral' color palette in [RColorBrewer::brewer.pal()].
#' @param colors_per_levels Named character string. To set the colors to use to plot data points and box for each state/range. Names = states/ranges; values = colors.
#'   If `NULL` (default), the default ggplot2 color palette ([scales::hue_pal()]) will be used. Only for categorical and biogeographic data.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#' @param return_mean_rates_vs_trait_data_df Logical. Whether to include in the output the data.frame of mean rates per trait values/states/ranges computed for
#'   each posterior sample at the focal time. Default is `FALSE`.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_jitter aes geom_point geom_boxplot scale_color_gradientn scale_color_manual ylab xlab guides ggtitle theme element_line element_rect element_text margin
#' @importFrom cowplot save_plot
#'
#' @details The main input `deepSTRAPP_outputs` is the typical output of [deepSTRAPP::run_deepSTRAPP_for_focal_time()].
#'   It provides information on results of a STRAPP test performed at a given `focal_time`.
#'
#'   Plots are built based on both trait data and diversification data as extracted for the given `focal_time`.
#'   Such data are recorded in the outputs of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_for_focal_time()]
#'   when `return_updated_trait_data_with_Map = TRUE` for trait data, and `extract_diversification_data_melted_df = TRUE` for diversification data.
#'   Please ensure to select those arguments when running deepSTRAPP.
#'
#'   Alternatively, the main input `deepSTRAPP_outputs` can be the output of [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   providing results of STRAPP tests over multiple time-steps. In this case, you must provide a `focal_time` to select the
#'   unique time-step used for plotting.
#'   * `return_updated_trait_data_with_Map` must be set to `TRUE` so that the trait data used to compute the tests are returned among the outputs
#'     under `$updated_trait_data_with_Map_over_time`. Alternatively, and more efficiently, `extract_trait_data_melted_df` can be set to `TRUE`
#'     so that trait data are already returned in a melted data.frame among the outputs under `$trait_data_df_over_time`.
#'   * `extract_diversification_data_melted_df` must be set to `TRUE` so that the diversification rates are returned
#'     among the outputs under `$diversification_data_df_over_time`.
#'
#'  For plotting all time-steps at once, see [deepSTRAPP::plot_rates_vs_trait_data_over_time()].
#'
#' @return The function returns a list with at least one element.
#'
#'   * `rates_vs_trait_ggplot` An object of classes `gg` and `ggplot`. This is a ggplot that can be displayed
#'     on the console with `print(output$rates_vs_trait_ggplot)`. It corresponds to the plot being displayed on the console
#'     when the function is run, if `display_plot = TRUE`, and can be further modify for aesthetics using the ggplot2 grammar.
#'
#'   If the trait data are 'continuous', the plot is a scatter plot showing how diversification rates varies with trait values.
#'   If the trait data are 'categorical' or 'biogeographic', the plot is a boxplot showing diversification rates per states/ranges.
#'
#'   Each plot also displays summary statistics for the STRAPP test associated with the data displayed:
#'   * An observed statistic computed across the mean traits/ranges and rates values shown on the plot. This is not the statistic of the STRAPP test itself,
#'     which is conducted across all BAMM posterior samples.
#'   * The quantile of null statistic distribution at the significant threshold used to define test significance. The test will be considered significant
#'     (i.e., the null hypothesis is rejected) if this value is higher than zero.
#'   * The p-value of the associated STRAPP test.
#'
#'   Optional summary data.frame:
#'   * `mean_rates_vs_trait_data_df` A data.frame with three columns providing the `$mean_rates` and `$trait_value`
#'     observed along branches at `focal_time`. Rates are averaged across all BAMM posterior samples.
#'     This is the raw data used to draw the plot. Included if `return_mean_rates_vs_trait_data_df = TRUE`.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plot.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::run_deepSTRAPP_for_focal_time()] [deepSTRAPP::plot_rates_vs_trait_data_over_time()]
#'
#' @examples
#'
#'
#' # ----- Example 1: Continuous trait ----- #
#'
#' # Load fake trait df
#' data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#' # Load phylogeny with old calibration
#' data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
#' # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#' data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'
#' ## Prepare trait data
#'
#' # Extract continuous trait data as a named vector
#' Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#'                                     nm = Ponerinae_trait_tip_data$Taxa)
#'
#' # Select a color scheme from lowest to highest values
#' color_scale = c("darkgreen", "limegreen", "orange", "red")
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)
#'
#' # Run a Stochastic Mapping based on a Brownian Motion model
#' # to interpolate values along branches and obtain a "contMap" object
#' Ponerinae_contMap <- phytools::contMap(Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data,
#'                                        res = 100, # Number of time steps
#'                                        plot = FALSE)
#' # Plot contMap = stochastic mapping of continuous trait
#' plot_contMap(contMap = Ponerinae_contMap,
#'              color_scale = color_scale)
#'
#' ## Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#' Ponerinae_deepSTRAPP_cont_old_calib_10My <- run_deepSTRAPP_for_focal_time(
#'    contMap = Ponerinae_contMap,
#'    ace = Ponerinae_ACE,
#'    tip_data = Ponerinae_cont_tip_data,
#'    trait_data_type = "continuous",
#'    BAMM_object = Ponerinae_BAMM_object_old_calib,
#'    focal_time = focal_time,
#'    rate_type = "net_diversification",
#'    return_perm_data = TRUE,
#'    # Need to be set to TRUE to save diversification data
#'    extract_diversification_data_melted_df = TRUE,
#'    # Need to be set to TRUE to save trait data
#'    return_updated_trait_data_with_Map = TRUE,
#'    return_updated_BAMM_object = TRUE)
#'
#' ## Explore output
#' str(Ponerinae_deepSTRAPP_cont_old_calib_10My, max.level = 1)
#'
#' # ------ Plot histogram of STRAPP overall test results from run_deepSTRAPP_for_focal_time() ------ #
#'
#' # Get plot
#' rates_vs_trait_output <- plot_rates_vs_trait_data_for_focal_time(
#'    deepSTRAPP_outputs = deepPonerinae_deepSTRAPP_cont_old_calib_10My,
#'    color_scale = c("grey80", "orange"),
#'    display_plot = TRUE,
#'    # PDF_file_path = "./plot_rates_vs_trait_10My.pdf"
#'    return_mean_rates_vs_trait_data_df = TRUE)
#' # Adjust aesthetics a posteriori
#' rates_vs_trait_ggplot_adj <- rates_vs_trait_output$rates_vs_trait_ggplot +
#'    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#' print(rates_vs_trait_ggplot_adj)
#'
#' # Explore melted data.frame of mean rates and trait values extracted for the given focal time.
#' head(rates_vs_trait_output$mean_rates_vs_trait_data_df) }
#'
#'
#' # ------ Plot histogram of STRAPP overall test results from run_deepSTRAPP_over_time() ------ #
#'
#' ## Load directly outputs from run_deepSTRAPP_over_time()
#' data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
#'
#' # Select focal_time = 10My
#' focal_time <- 10
#'
#' # Get plot
#' rates_vs_trait_output <- plot_rates_vs_trait_data_for_focal_time(
#'    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'    focal_time = focal_time,
#'    color_scale = c("grey80", "purple"),
#'    display_plot = TRUE)
#'    # PDF_file_path = "./plot_rates_vs_trait_10My.pdf"
#'
#' # Adjust aesthetics a posteriori
#' rates_vs_trait_ggplot_adj <- rates_vs_trait_output$rates_vs_trait_ggplot +
#'    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#' print(rates_vs_trait_ggplot_adj)
#'
#' # ----- Example 2: Categorical trait ----- #
#'
#' ## Load data
#'
#' # Load phylogeny
#' data(Ponerinae_tree, package = "deepSTRAPP")
#' # Load trait df
#' data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#' # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#' data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'
#' ## Prepare trait data
#'
#' # Extract categorical data with 3-levels
#' Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#'                                         nm = Ponerinae_trait_tip_data$Taxa)
#' table(Ponerinae_cat_3lvl_tip_data)
#'
# Select color scheme for states
#' colors_per_states <- c("forestgreen", "sienna", "goldenrod")
#' names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
#'
#' \dontrun{  (May take several minutes to run)
#' ## Produce densityMaps using stochastic character mapping based on an ARD Mk model
#' Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
#'    tip_data = Ponerinae_cat_3lvl_tip_data,
#'    phylo = Ponerinae_tree_old_calib,
#'    trait_data_type = "categorical",
#'    colors_per_states = colors_per_states,
#'    evolutionary_models = "ARD", # Use default ARD model
#'    nb_simulations = 100, # Reduce number of simulations to save time
#'    seed = 1234, # Seet seed for reproducibility
#'    return_best_model_fit = TRUE,
#'    return_model_selection_df = TRUE,
#'    plot_map = FALSE) }
#'
#' # Load directly output
#' data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
#'
#' ## Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#' Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My <- run_deepSTRAPP_for_focal_time(
#'     densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
#'     ace = Ponerinae_cat_3lvl_data_old_calib$ace,
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     trait_data_type = "categorical",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     focal_time = focal_time,
#'     rate_type = "net_diversification",
#'     posthoc_pairwise_tests = TRUE,
#'     return_perm_data = TRUE,
#'     # Need to be set to TRUE to save diversification data
#'     extract_diversification_data_melted_df = TRUE,
#'     # Need to be set to TRUE to save trait data
#'     return_updated_trait_data_with_Map = TRUE,
#'     return_updated_BAMM_object = TRUE)
#'
#' ## Explore output
#' str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My, max.level = 1)
#'
#' ## Plot rates vs. states
#' rates_vs_trait_output <- plot_rates_vs_trait_data_for_focal_time(
#'    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My,
#'    focal_time = 10,
#'    select_trait_levels = c("arboreal", "terricolous"), # Select only two levels
#'    colors_per_levels = colors_per_states[c("arboreal", "terricolous")], # Adjust colors
#'    display_plot = TRUE,
#'    # PDF_file_path = "./plot_rates_vs_trait_10My.pdf",
#'    return_mean_rates_vs_trait_data_df = TRUE)
#'
#' # Adjust aesthetics a posteriori
#' rates_vs_trait_ggplot_adj <- rates_vs_trait_output$rates_vs_trait_ggplot +
#'    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#' print(rates_vs_trait_ggplot_adj)
#'
#' # Explore melted data.frame of mean rates and states extracted for the given focal time.
#' head(rates_vs_trait_output$mean_rates_vs_trait_data_df) }
#'


plot_rates_vs_trait_data_for_focal_time <- function (deepSTRAPP_outputs,
                                                     focal_time = NULL,
                                                     rate_type = "net_diversification",
                                                     select_trait_levels = "all",
                                                     color_scale = NULL,
                                                     colors_per_levels = NULL,
                                                     display_plot = TRUE,
                                                     PDF_file_path = NULL,
                                                     return_mean_rates_vs_trait_data_df = FALSE)

{
  ### Check input validity
  {
    ## deepSTRAPP_outputs
    # Check presence of diversification_data_df
    if (is.null(deepSTRAPP_outputs$diversification_data_df) & is.null(deepSTRAPP_outputs$diversification_data_df_over_time))
    {
      stop(paste0("`deepSTRAPP_outputs` must have a `$diversification_data_df` or `$diversification_data_df_over_time` element.\n",
                  "Be sure to set `extract_diversification_data_melted_df = TRUE` in [deepSTRAPP::run_deepSTRAPP_for_focal_time] or [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                  "This element is needed to plot rates vs. trait data."))
    }

    ## Identify the type of inputs
    if (is.null(deepSTRAPP_outputs$diversification_data_df_over_time))
    {
      inputs_over_time <- FALSE
    } else {
      inputs_over_time <- TRUE
    }

    ## focal_time
    # Ensure a focal_time is provided if output is from [deepSTRAPP::run_deepSTRAPP_over_time()]
    if (inputs_over_time)
    {
      if (is.null(focal_time))
      {
        stop(paste0("You provided as input a `deepSTRAPP_outputs` object with multiple time-steps resulting from [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                    "You must provide a `focal_time` to select the appropriate time-step to be plotted.\n",
                    "For plotting all time-steps at once, please see [deepSTRAPP::plot_rates_vs_trait_data_over_time]."))
      }
      # Ensure focal_time match (any) time in deepSTRAPP_outputs
      if (!(focal_time %in% deepSTRAPP_outputs$time_steps))
      {
        stop(paste0("You provided as input a `deepSTRAPP_outputs` object with multiple time-steps resulting from [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                    "You must provide a `focal_time` that matches with the `$time_steps` recorded in the `deepSTRAPP_outputs` object."))
      }
    } else {
      # Ensure focal_time match with the focal_time recorded in deepSTRAPP_outputs
      if (!is.null(focal_time))
      {
        if (!(focal_time %in% deepSTRAPP_outputs$focal_time))
        {
          stop(paste0("You provided as input a `deepSTRAPP_outputs` object with a unique time-step resulting from [deepSTRAPP::run_deepSTRAPP_for_focal_time].\n",
                      "However, the `focal_time` you provided does not that match with the `$focal_time` recorded in the `deepSTRAPP_outputs` object.\n",
                      "focal_time provided: ",focal_time,".\n",
                      "focal_time recorded in `deepSTRAPP_outputs`: ",deepSTRAPP_outputs$focal_time,"."))
        }
      } else {
        # Extract focal time if not provided
        focal_time <- deepSTRAPP_outputs$focal_time
      }
    }

    ## rate_type must be either "speciation", "extinction" or "net_diversification"
    if (!(rate_type %in% c("speciation", "extinction", "net_diversification")))
    {
      stop("'rate_type' can only be 'speciation', 'extinction', or 'net_diversification'.")
    }

    ## Extract diversification_data_df
    if (!inputs_over_time)
    {
      # For outputs from run_deepSTRAPP_for_focal_time
      diversification_data_df <- deepSTRAPP_outputs$diversification_data_df
    } else {
      # For outputs from run_deepSTRAPP_over_time
      diversification_data_df <- deepSTRAPP_outputs$diversification_data_df_over_time
      diversification_data_df <- diversification_data_df[diversification_data_df$focal_time == focal_time, ]
    }

    ## Extract trait_data_df
    if (!inputs_over_time)
    {
      # Case for input from [deepSTRAPP::run_deepSTRAPP_for_focal_time]

      # Check presence of $updated_trait_data_with_Map
      if (is.null(deepSTRAPP_outputs$updated_trait_data_with_Map))
      {
        stop(paste0("`deepSTRAPP_outputs` must have a `$updated_trait_data_with_Map` element.\n",
                    "Be sure to set `return_updated_trait_data_with_Map = TRUE` in [deepSTRAPP::run_deepSTRAPP_for_focal_time].\n",
                    "This element is needed to plot rates vs. trait data."))
      } else {

        # Extract trait data from $updated_trait_data_with_Map
        trait_data_df <- as.data.frame(deepSTRAPP_outputs$updated_trait_data_with_Map$trait_data)
        trait_data_df$tip_ID <- names(deepSTRAPP_outputs$updated_trait_data_with_Map$trait_data)
        names(trait_data_df) <- c("trait_value", "tip_ID")
        trait_data_df <- trait_data_df[, c("tip_ID", "trait_value")]
      }

    } else {
      # Case for input from [deepSTRAPP::run_deepSTRAPP_over_time]

      # Check presence of $trait_data_df_over_time
      if (!is.null(deepSTRAPP_outputs$trait_data_df_over_time))
      {
        # Extract trait_df
        trait_data_df <- deepSTRAPP_outputs$trait_data_df_over_time
        trait_data_df <- trait_data_df[trait_data_df$focal_time == focal_time, ]
        trait_data_df <- trait_data_df[, c("tip_ID", "trait_value")]
      } else {
        # If absent, check presence of $updated_trait_data_with_Map_over_time instead
        if (!is.null(deepSTRAPP_outputs$updated_trait_data_with_Map_over_time))
        {
          stop(paste0("`deepSTRAPP_outputs` must have a `$trait_data_df_over_time` or `$updated_trait_data_with_Map_over_time` element.\n",
                      "Be sure to set `extract_trait_data_melted_df = TRUE` or at least `return_updated_trait_data_with_Map = TRUE`",
                      "in [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                      "One of these elements is needed to plot rates vs. trait data."))
        } else {
          # Extract trait_df from $updated_trait_data_with_Map_over_time
          focal_time_ID <- which(focal_time == deepSTRAPP_outputs$time_steps)

          updated_trait_data_with_Map <- deepSTRAPP_outputs$updated_trait_data_with_Map_over_time[[focal_time_ID]]

          trait_data_df <- as.data.frame(updated_trait_data_with_Map$trait_data)
          trait_data_df$tip_ID <- names(updated_trait_data_with_Map$trait_data)
          names(trait_data_df) <- c("trait_value", "tip_ID")
          trait_data_df <- trait_data_df[, c("tip_ID", "trait_value")]
        }
      }
    }

    ## Extract type of trait
    if (!inputs_over_time)
    {
      trait_data_type <- deepSTRAPP_outputs$updated_trait_data_with_Map$trait_data_type
    } else {
      trait_data_type <- deepSTRAPP_outputs$trait_data_type
    }

    if (trait_data_type == "continuous")
    {
      ## Case for "continuous" traits

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
    } else {
      ## Case for "categorical" and "biogeographic" traits

      ## Extract trait levels
      states_in_trait_data_df <- unique(trait_data_df$trait_value)
      states_in_trait_data_df <- states_in_trait_data_df[order(states_in_trait_data_df)]

      ## select_trait_levels
      if (!any(select_trait_levels == "all"))
      {
        # Check that select_trait_levels are all found in trait_data_df

        if (!all(select_trait_levels %in% states_in_trait_data_df))
        {
          stop(paste0("Some states/ranges listed in 'select_trait_levels' are not found in among the trait data.\n",
                      "'select_trait_levels' = ",paste(select_trait_levels[order(select_trait_levels)], collapse = ", "),".\n",
                      "Observed states/ranges in trait data = ", paste(states_in_trait_data_df, collapse = ", ")),".")
        }
      }

      # Update list of states/ranges to keep only the selected ones
      if (!any(select_trait_levels == "all"))
      {
        states_in_trait_data_df <- select_trait_levels
      }

      ## colors_per_levels
      # Check whether all colors are valid
      if (!is.null(colors_per_levels))
      {
        # Check that the color match the selected states/ranges
        if (!all(states_in_trait_data_df %in% names(colors_per_levels)))
        {
          missing_states <- states_in_trait_data_df[!(states_in_trait_data_df %in% names(colors_per_levels))]
          stop(paste0("Not all selected states/ranges are found in 'colors_per_levels'.\n",
                      "Missing states/ranges: ", paste(missing_states, collapse = ", "), "."))
        }
        if (!all(is_color(colors_per_levels)))
        {
          invalid_colors <- colors_per_levels[!is_color(colors_per_levels)]
          stop(paste0("Some color names in 'colors_per_levels' are not valid.\n",
                      "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
        }
      }
    }

    ## Extract STRAPP_results
    if (!inputs_over_time)
    {
      # For outputs from run_deepSTRAPP_for_focal_time
      STRAPP_results <- deepSTRAPP_outputs$STRAPP_results
    } else {
      # For outputs from run_deepSTRAPP_over_time
      focal_time_ID <- which(deepSTRAPP_outputs$time_steps == focal_time)
      STRAPP_results <- deepSTRAPP_outputs$STRAPP_results_over_time[[focal_time_ID]]
    }

    ## STRAPP_results

    # STRAPP_results must have recorded a STRAPP test to display its results.
    STRAPP_results_available <- TRUE
    if (STRAPP_results$trait_data_type_for_stats == "none")
    {
      if (!inputs_over_time)
      {
        warning(paste0("STRAPP test results are missing from 'deepSTRAPP_outputs$STRAPP_results'.\n",
                       "A unique ML state/range was inferred across branches for 'focal_time' = ",STRAPP_results$focal_time,".\n",
                       "No STRAPP test for differences in rates between states/ranges can be computed with a unique state/range.\n",
                       "Therefore, no STRAPP results can be associoted to the plot for this 'focal_time'."))
      } else {
        warning(paste0("STRAPP test results are missing from 'deepSTRAPP_outputs$STRAPP_results_over_time' for the given 'focal_time'.\n",
                       "A unique ML state/range was inferred across branches for 'focal_time' = ",STRAPP_results$focal_time,".\n",
                       "No STRAPP test for differences in rates between states/ranges can be computed with a unique state/range.\n",
                       "Therefore, no STRAPP results can be associoted to the plot for this 'focal_time'."))
      }
      STRAPP_results_available <- FALSE
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

  # Bind column names to prevent Notes
  tip_ID <- rates <- mean_rates <- trait_value <- NA

  ## Extract mean diversification data, only the selected rate_type
  diversification_data_df <- diversification_data_df |>
    dplyr::filter(diversification_data_df$rate_type == rate_type) |>
    dplyr::group_by(tip_ID) |>
    dplyr::summarize(mean_rates = mean(rates))

  ## Merge rates and traits in a unique data.frame
  data_melted_df <- dplyr::left_join(x = diversification_data_df, y = trait_data_df,
                                     by = dplyr::join_by(tip_ID))

  ## Extract test summary
  if (STRAPP_results_available)
  {
    # Extract quantile of the critical threshold
    estimate_quantile <- names(STRAPP_results$estimate)
    # Extract value at the critical threshold
    quantile_value <- round(STRAPP_results$estimate, digits = 3)
    # Extract p-value
    p_value <- round(STRAPP_results$p_value, digits = 3)
  }

  ## Compute observed stat across mean values (not used directly for the STRAPP test)
  if (STRAPP_results_available)
  {
    if (STRAPP_results$trait_data_type_for_stats == "continuous")
    {
      stat_name <- "\u03C1" # Rho (unicode)

      ## Wrapped-up function to extract rho stats from Spearman's correlation test
      spearman_test <- function(rates, trait_data)
      {
        if (stats::sd(rates, na.rm = TRUE) == 0)
        { # Case with no variance in rates. Rho = 0.
          return(0)
        } else { # Default case
          test_output <- stats::cor.test(rates, trait_data, method = "spearman", exact = FALSE)
          return(test_output$estimate)
        }
      }

      # Compute observed stat across mean data
      stat_estimate <- spearman_test(rates = data_melted_df$mean_rates, trait_data = data_melted_df$trait_value)
    }

    if (STRAPP_results$trait_data_type_for_stats == "binary")
    {
      stat_name <- "U-stat"

      # Check the type of test
      two_tailed <- STRAPP_results$two_tailed # Type of test: two-tailed or not
      one_tailed_hypothesis <- STRAPP_results$one_tailed_hypothesis # Type of hypothesis if one-tailed test

      # Parse one_tailed_hypothesis
      if (!is.null(one_tailed_hypothesis))
      {
        one_tailed_hypothesis_parsed <- gsub(pattern = " ", replacement = "", x = one_tailed_hypothesis)
        trait_states <- as.character(unlist(strsplit(x = one_tailed_hypothesis_parsed, split = ">")))
      } else {
        trait_states <- NULL
      }

      ## Wrapped-up function to extract U-stats from Mann-Whitney-Wilcoxon's rank-sum test
      mann_whitney_wilcoxon_test <- function(rates, trait_data, two_tailed, trait_states)
      {
        if (two_tailed)
        { # Case for two-tailed test
          test_output <- stats::wilcox.test(formula = rates ~ trait_data, exact = FALSE)
        } else { # Case for one-tailed test
          test_output <- stats::wilcox.test(x = rates[which(trait_data == trait_states[1])], # State with the higher ranked rates in Ha
                                            y = rates[which(trait_data == trait_states[2])], # State with the lower ranked rates in Ha
                                            exact = FALSE)
        }
        return(test_output$statistic)
      }

      # Compute observed stat across mean data
      stat_estimate <- mann_whitney_wilcoxon_test(rates = data_melted_df$mean_rates,
                                                  trait_data = data_melted_df$trait_value,
                                                  two_tailed = two_tailed,
                                                  trait_states = trait_states)
      # Center stats around location shift of the null hypothesis (mu)
      # Null hypothesis is that ranks of the values of the two groups are random
      # Compute location shift (mu) from state frequencies as average of the products of frequencies
      trait_data_counts <- table(data_melted_df$trait_value)
      trait_data_counts <- trait_data_counts[!is.na(names(trait_data_counts))] # Remove NA
      stat_mu <- prod(trait_data_counts)/2
      # Center U-stats to get an estimate of how greater/lower (far away) than the null hypothesis (mu) is the observed U-stats
      stat_estimate <- stat_estimate - stat_mu
      # In two-tailed tests, the absolute deviation to the null-expectation is used.
      # But for coherency with the correlation tests, better to provide a negative stat when 'state 1' has lower rates than 'state 2', and respectively.
      # if (two_tailed) { stat_estimate <- abs(stat_estimate) }
    }

    if (STRAPP_results$trait_data_type_for_stats == "multinominal")
    {
      stat_name <- "H-stat"

      ## Wrapped-up function to extract H-stats from Kruskal-Wallis's one-way ANOVA on ranks test
      kruskal_wallis_test <- function(rates, trait_data)
      {
        # Compute the Kruskal-Wallis test
        test_output <- stats::kruskal.test(rates ~ trait_data)

        # If the test failed to provide a statistic because the value is reaching the ceiling for computation,
        # use the Khi-squared approximation by setting an extremely high p-value
        if (is.na(test_output$statistic))
        {
          H_approximation <- stats::qchisq(p = 1 - 10^-9, df = test_output$parameter)
          return(H_approximation)
        } else { # Otherwise, provide the computed H-stats
          return(test_output$statistic)
        }
      }

      # Compute observed stat across mean data
      stat_estimate <- kruskal_wallis_test(rates = data_melted_df$mean_rates, trait_data = data_melted_df$trait_value)
    }
  }

  ## Set label for y-lab
  y_label <- stringr::str_to_sentence(paste0(sub(x = rate_type, pattern = "_", replacement = " "), " rates"))

  ## Case for continuous trait data
  if (trait_data_type == "continuous")
  {

    ## Prepare color scale
    if (!is.null(color_scale))
    {
      # Use the provided color to build the color palette
      col_fn <- grDevices::colorRampPalette(colors = color_scale)
      colors_per_values <- col_fn(n = 10)
    } else {
      # Default: use the 'Spectral' palette from RColorBrewer
      colors_per_values <- rev(RColorBrewer::brewer.pal(name = "Spectral", n = 10))
    }

    ## Plot rates vs. trait values across all BAMM samples
    ggplot_rates_vs_traits <- ggplot2::ggplot(data = data_melted_df) +

      # Plot points for all samples
      ggplot2::geom_point(mapping = ggplot2::aes(y = mean_rates,
                                                 x = trait_value,
                                                 color = trait_value),
                          alpha = 0.80, size = 3) +

      # Adjust point colors
      ggplot2::scale_color_gradientn(colours = colors_per_values) +
      ggplot2::guides(colour = "none") +

      # Adjust axis titles
      ggplot2::ylab(y_label) +
      ggplot2::xlab("Trait values") +

      # Add title
      ggplot2::ggtitle(paste0("Mean rates vs. trait values\n",
                              "Focal time = ", focal_time)) +

      # Adjust aesthetics
      ggplot2::theme(
         plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
         panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
         panel.background = ggplot2::element_rect(fill = NA, color = NA),
         plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                            margin = ggplot2::margin(b = 15, t = 5)),
         axis.title = ggplot2::element_text(size = 20, color = "black"),
         axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
         axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
         axis.line = ggplot2::element_line(linewidth = 1.0),
         axis.text = ggplot2::element_text(size = 18, color = "black"),
         axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
         axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

    # Add test summary if available
    if (STRAPP_results_available)
    {
      ggplot_rates_vs_traits <- ggplot_rates_vs_traits +

      # Observed stats, Q%: Estimate,  p-value
      annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                   label = paste0(stat_name," obs = ", round(stat_estimate, digits = 3), "\n",
                                  "Q", estimate_quantile, " = ", quantile_value, "\n",
                                  "P-value = ", p_value))
    }

  } else {
    ## Case for categorical/biogeographic trait data

    # Filter data to keep only the selected states/ranges
    data_melted_df <- data_melted_df[data_melted_df$trait_value %in% states_in_trait_data_df, ]

    ## Prepare colors_per_levels to use in plots
    if (is.null(colors_per_levels))
    {
      nb_groups <- length(levels(as.factor(data_melted_df$trait_value)))
      # Default: use the default ggplot palette from scales
      col_fn <- scales::hue_pal()
      colors_per_levels <- col_fn(n = nb_groups)
      names(colors_per_levels) <- levels(as.factor(data_melted_df$trait_value))
    }

    if (trait_data_type == "categorical")
    {
      ## Case for categorical trait data

      ## Plot rates vs. states across all BAMM samples
      ggplot_rates_vs_traits <- ggplot2::ggplot(data = data_melted_df) +

        # Plot boxplot per states
        ggplot2::geom_boxplot(mapping = ggplot2::aes(y = mean_rates, x = trait_value,
                                                     fill = trait_value)) +
        # Plot points for all samples
        ggplot2::geom_jitter(mapping = ggplot2::aes(y = mean_rates, x = trait_value,
                                                    fill = trait_value),
                             alpha = 0.50, size = 3, width = 0.25, shape = 21, color = "black") +

        # Adjust legend
        ggplot2::scale_fill_manual(name = "States",
                                   breaks = names(colors_per_levels),
                                   values = colors_per_levels) +
        ggplot2::guides(fill = "none") +

        # Adjust axis titles
        ggplot2::ylab(y_label) +
        ggplot2::xlab("States") +

        # Add title
        ggplot2::ggtitle(paste0("Mean rates vs. states\n",
                                "Focal time = ", focal_time)) +

        # Adjust aesthetics
        ggplot2::theme(
          plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
          panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
          panel.background = ggplot2::element_rect(fill = NA, color = NA),
          plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                             margin = ggplot2::margin(b = 15, t = 5)),
          axis.title = ggplot2::element_text(size = 20, color = "black"),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
          axis.line = ggplot2::element_line(linewidth = 1.0),
          axis.text = ggplot2::element_text(size = 18, color = "black"),
          axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
          axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

      # Add test summary if available
      if (STRAPP_results_available)
      {
        ggplot_rates_vs_traits <- ggplot_rates_vs_traits +

          # Observed stats, Q%: Estimate,  p-value
          annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                       label = paste0(stat_name," obs = ", round(stat_estimate, digits = 3), "\n",
                                      "Q", estimate_quantile, " = ", quantile_value, "\n",
                                      "P-value = ", p_value))
      }

    } else {

      ## Case for biogeographic data

      ## Plot rates vs. states across all BAMM samples
      ggplot_rates_vs_traits <- ggplot2::ggplot(data = data_melted_df) +

        # Plot boxplot per ranges
        ggplot2::geom_boxplot(mapping = ggplot2::aes(y = mean_rates, x = trait_value,
                                                     fill = trait_value)) +
        # Plot points for all samples
        ggplot2::geom_jitter(mapping = ggplot2::aes(y = mean_rates, x = trait_value,
                                                    fill = trait_value),
                             alpha = 0.50, size = 3, width = 0.25, shape = 21, color = "black") +

        # Adjust legend
        ggplot2::scale_fill_manual(name = "Ranges",
                                   breaks = names(colors_per_levels),
                                   values = colors_per_levels) +
        ggplot2::guides(fill = "none") +

        # Adjust axis titles
        ggplot2::ylab(y_label) +
        ggplot2::xlab("Ranges") +

        # Add title
        ggplot2::ggtitle(paste0("Mean rates vs. ranges\n",
                                "Focal time = ", focal_time)) +

        # Adjust aesthetics
        ggplot2::theme(
          plot.margin = ggplot2::margin(0.3, 0.5, 0.5, 0.5, "inches"), # trbl
          panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
          panel.background = ggplot2::element_rect(fill = NA, color = NA),
          plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                             margin = ggplot2::margin(b = 15, t = 5)),
          axis.title = ggplot2::element_text(size = 20, color = "black"),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
          axis.line = ggplot2::element_line(linewidth = 1.0),
          axis.text = ggplot2::element_text(size = 18, color = "black"),
          axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
          axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

      # Add test summary if available
      if (STRAPP_results_available)
      {
        ggplot_rates_vs_traits <- ggplot_rates_vs_traits +

          # Observed stats, Q%: Estimate,  p-value
          annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                       label = paste0(stat_name," obs = ", round(stat_estimate, digits = 3), "\n",
                                      "Q", estimate_quantile, " = ", quantile_value, "\n",
                                      "P-value = ", p_value))
      }
    }
  }

  ## Display plot if requested
  if (display_plot)
  {
    print(ggplot_rates_vs_traits)
  }

  ## Export plot if requested
  if (!is.null(PDF_file_path))
  {
    cowplot::save_plot(plot = ggplot_rates_vs_traits,
                       filename = PDF_file_path,
                       base_height = 8, base_width = 10)
  }

  ## Build output
  output <- list()

  ## Store ggplot
  output$rates_vs_trait_ggplot <- ggplot_rates_vs_traits

  ## Store melted df if requested
  if (return_mean_rates_vs_trait_data_df)
  {
    output$mean_rates_vs_trait_data_df <- as.data.frame(data_melted_df)
  }

  ## Return output
  return(invisible(output))

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

