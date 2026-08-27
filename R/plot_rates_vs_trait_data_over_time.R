
#' @title Plot mean rates vs. trait data over time-steps
#'
#' @description Plot mean rates vs. trait data of branches as extracted for all time-steps.
#'   Data are extracted from the output of a deepSTRAPP run carried out with
#'   [deepSTRAPP::run_deepSTRAPP_over_time()]) over multiple time-steps.
#'
#'   For each time-step, returns a plot showing rates vs. trait data.
#'   If the trait data are 'continuous', the plot is a scatter plot.
#'   If the trait data are 'categorical' or 'biogeographic', the plot is a boxplot.
#'
#'   If a PDF file path is provided in `PDF_file_path`, all plots will be saved directly in a unique PDF file,
#'   with one page per plot/time-step.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_over_time()]
#'    that runs the whole deepSTRAPP workflow over multiple time-steps.
#'    The list needs to include two data.frame: `$trait_data_df_over_time` and `$diversification_data_df_over_time`
#'    by setting `extract_trait_data_melted_df = TRUE` and `extract_diversification_data_melted_df = TRUE`.
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
#' @param return_mean_data_per_branches_df Logical. Whether to include in the output the data.frame of mean rates per trait values/states/ranges
#'   of all branches computed over all time-steps and used for the plots. Default is `FALSE`.
#'
#' @export
#' @importFrom cowplot save_plot
#'
#' @details The main input `deepSTRAPP_outputs` is the typical output of [deepSTRAPP::run_deepSTRAPP_over_time()].
#'   It provides information on results of a STRAPP test performed over multiple time-steps.
#'
#'   Plots are built based on both trait data and diversification data as extracted for the time-steps.
#'   Such data are recorded in the outputs of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
#'   * `extract_trait_data_melted_df` must be set to `TRUE` so that trait data are returned
#'     in a melted data.frame among the outputs under `$trait_data_df_over_time`.
#'   * `extract_diversification_data_melted_df` must be set to `TRUE` so that the diversification rates are returned
#'     among the outputs under `$diversification_data_df_over_time`.
#'   * `return_STRAPP_results` must be set to `TRUE` so that the STRAPP results objects recording the summary statistics for the STRAPP tests
#'     are returned among the outputs under `$STRAPP_results_over_time`.
#'
#'  For plotting a single `focal_time`, see [deepSTRAPP::plot_rates_vs_trait_data_for_focal_time()].
#'
#' @return The function returns a list with at least one element.
#'
#'   * `rates_vs_trait_ggplots` A list of objects of classes `gg` and `ggplot` ordered as in `$time_steps`.
#'     Each element corresponds to a ggplot for a given `focal_time`. They can be displayed on the console with `print(output$rates_vs_trait_ggplots[[i]])`.
#'     They correspond to the plots being displayed on the console one by one when the function is run, if `display_plot = TRUE`,
#'     and can be further modify for aesthetics using the ggplot2 grammar.
#'
#'   If the trait data are 'continuous', the plots are scatter plots showing how mean diversification rates varies with mean trait values across branches.
#'   If the trait data are 'categorical' or 'biogeographic', the plots are boxplots showing mean diversification rates per states/ranges per branches.
#'
#'   Each plot also displays summary statistics for the STRAPP test associated with the data displayed:
#'   * An observed statistic computed across the mean traits/ranges and rates values shown on the plot. This is not the statistic of the STRAPP test itself,
#'     which is conducted across stochastic maps X BAMM posterior samples (i.e., it is a distribution, not a unique value)..
#'   * The quantile of null statistic distribution at the significant threshold used to define test significance. The test will be considered significant
#'     (i.e., the null hypothesis is rejected) if this value is higher than zero.
#'   * The p-value of the associated STRAPP test.
#'
#'   Optional summary data.frame:
#'   * `mean_data_per_branches_df` A data.frame with four columns providing the `$mean_trait_value` (continuous) / `$states` (categorical) / `$ranges` (biogeographic)
#'     and `$mean_rates` computed along branches (`$tip_ID`) at the different `focal_time`. Rates/Traits are averaged for each stochastic maps X BAMM posterior samples
#'     used for testing, in accordance with the `uncertainty_strategy` selected when running deepSTRAPP.
#'     This is the raw data used to draw each plot for each `focal_time`. Included if `return_mean_data_per_branches_df = TRUE`.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a unique PDF file with one plot/page per `$time_steps`.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::run_deepSTRAPP_over_time()] [deepSTRAPP::plot_rates_vs_trait_data_for_focal_time()]
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  # ----- Example 1: Continuous trait ----- #
#'
#'  # Load fake trait df
#'  data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#'  # Load phylogeny with old calibration
#'  data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
#'
#'  # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#'  data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'   ## Prepare trait data
#'
#'  # Extract continuous trait data as a named vector
#'  Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#'                                      nm = Ponerinae_trait_tip_data$Taxa)
#'
#'  # Select a color scheme from lowest to highest values
#'  color_scale = c("darkgreen", "limegreen", "orange", "red")
#'
#'  # Get Ancestral Character Estimates based on a Brownian Motion model
#'  # To obtain values at internal nodes
#'  Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)
#'
#'  \donttest{ # (May take several minutes to run)
#'  # Run a Stochastic Mapping based on a Brownian Motion model
#'  # to interpolate values along branches and obtain a "contMap" object
#'  Ponerinae_contMap <- phytools::contMap(Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data,
#'                                         res = 100, # Number of time steps
#'                                         plot = FALSE)
#'  # Plot contMap = stochastic mapping of continuous trait
#'  plot_contMap(contMap = Ponerinae_contMap,
#'               color_scale = color_scale)
#'
#'  ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
#'  # nb_time_steps <- 5
#'  time_step_duration <- 5
#'  time_range <- c(0, 40)
#'
#'  ## Run deepSTRAPP on net diversification rates
#'  Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     contMap = Ponerinae_contMap,
#'     ace = Ponerinae_ACE,
#'     tip_data = Ponerinae_cont_tip_data,
#'     trait_data_type = "continuous",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     # nb_time_steps = nb_time_steps,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     uncertainty_strategy = "rates_only",
#'     return_perm_data = TRUE,
#'     # Need to be set to TRUE to save trait data
#'     extract_trait_data_melted_df = TRUE,
#'     # Need to be set to TRUE to save diversification data
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_BAMM_object = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly trait data output
#'  Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- readRDS(system.file("extdata",
#'     "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  # Explore output
#'  str(Ponerinae_deepSTRAPP_cont_old_calib_0_40, max.level = 1)
#'
#'  ## Plot for all time-steps
#'  rates_vs_trait_outputs <- plot_rates_vs_trait_data_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     color_scale = c("grey80", "purple"), # Adjust color scale
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_rates_vs_trait_0_40My.pdf",
#'     return_mean_data_per_branches_df = TRUE)
#'
#'  ## Print plot for time step 3 = 10 My
#'  print(rates_vs_trait_outputs$rates_vs_trait_ggplots[[3]])
#'
#'  ## Explore melted data.frame of rates and trait data
#'  head(rates_vs_trait_outputs$mean_data_per_branches_df)
#'
#'  # ----- Example 2: Categorical data ----- #
#'
#'  ## Load data
#'
#'  # Load trait df
#'  data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#'  # Load phylogeny
#'  data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
#'
#'  # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#'  data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Prepare trait data
#'
#'  # Extract categorical data with 3-levels
#'  Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#'                                          nm = Ponerinae_trait_tip_data$Taxa)
#'  table(Ponerinae_cat_3lvl_tip_data)
#'
#'  # Select color scheme for states
#'  colors_per_states <- c("forestgreen", "sienna", "goldenrod")
#'  names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Produce densityMaps using stochastic character mapping based on an ARD Mk model
#'  Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     phylo = Ponerinae_tree_old_calib,
#'     trait_data_type = "categorical",
#'     colors_per_levels = colors_per_states,
#'     evolutionary_models = "ARD",
#'     nb_simulations = 100,
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE,
#'     plot_map = FALSE) }
#'
#'  # Load directly trait data output
#'  data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
#'
#'  ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
#'  # nb_time_steps <- 5
#'  time_step_duration <- 5
#'  time_range <- c(0, 40)
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates across time-steps.
#'  Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
#'     ace = Ponerinae_cat_3lvl_data_old_calib$ace,
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     nb_simulations = 100,
#'     trait_data_type = "categorical",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     # nb_time_steps = nb_time_steps,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     rate_type = "net_diversification",
#'     uncertainty_strategy = "paired",
#'     seed = 1234, # Set for reproducibility
#'     alpha = 0.10, # Select a generous level of significance for the sake of the example
#'     posthoc_pairwise_tests = TRUE,
#'     return_perm_data = TRUE,
#'     # Need to be set to TRUE to save trait data
#'     extract_trait_data_melted_df = TRUE,
#'     # Need to be set to TRUE to save diversification data
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_BAMM_object = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly deepSTRAPP output
#'  Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- readRDS(system.file("extdata",
#'     "Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", package = "deepSTRAPP"))
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  # Explore output
#'  str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, max.level = 1)
#'
#'  # Adjust color scheme
#'  colors_per_states <- c("orange", "dodgerblue", "red")
#'  names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
#'
#'  ## Plot for all time-steps
#'  rates_vs_trait_outputs <- plot_rates_vs_trait_data_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#'     colors_per_levels = colors_per_states, # Adjust color scheme
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_rates_vs_trait_0_40My.pdf",
#'     return_mean_data_per_branches_df = TRUE)
#'
#'  ## Print plot for time step 3 = 10 My
#'  print(rates_vs_trait_outputs$rates_vs_trait_ggplots[[3]])
#'
#'  ## Explore melted data.frame of rates and states
#'  head(rates_vs_trait_outputs$mean_data_per_branches_df)
#' }
#'


plot_rates_vs_trait_data_over_time <- function (deepSTRAPP_outputs,
                                                rate_type = "net_diversification",
                                                select_trait_levels = "all",
                                                color_scale = NULL,
                                                colors_per_levels = NULL,
                                                display_plot = TRUE,
                                                PDF_file_path = NULL,
                                                return_mean_data_per_branches_df = FALSE)

{
  ### Check input validity
  {
    # Should all be carried within [deepSTRAPP::plot_rates_vs_trait_data_for_focal_time()]
  }

  ## Extract time-steps
  time_steps <- deepSTRAPP_outputs$time_steps

  ## Extract type of traits
  trait_data_type <- deepSTRAPP_outputs$trait_data_type

  ## Set unique color scheme with colors_per_levels for "categorical/biogeographic" traits
  if (trait_data_type %in% c("categorical", "biogeographic"))
  {
    ## Extract trait data
    trait_data_df <- deepSTRAPP_outputs$trait_data_df_over_time

    ## Get list of all states/ranges to plot

    # Extract trait levels
    states_in_trait_data_df <- unique(trait_data_df$trait_value)
    states_in_trait_data_df <- states_in_trait_data_df[order(states_in_trait_data_df)]

    ## Adjust according to select_trait_levels
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

    # Filter data to keep only the selected states/ranges
    trait_data_df <- trait_data_df[trait_data_df$trait_value %in% states_in_trait_data_df, ]

    ## Prepare colors_per_levels to use in plots
    if (is.null(colors_per_levels))
    {
      nb_groups <- length(levels(as.factor(trait_data_df$trait_value)))
      # Default: use the default ggplot palette from scales
      col_fn <- scales::hue_pal()
      colors_per_levels <- col_fn(n = nb_groups)
      names(colors_per_levels) <- levels(as.factor(trait_data_df$trait_value))
    } else {
      colors_per_levels <- colors_per_levels[states_in_trait_data_df]
    }

  }

  ## Save initial par() and reassign them on exit
  oldpar <- par(no.readonly = TRUE)
  oldpar$new <- NULL
  on.exit(par(oldpar))

  ## Loop per time-steps
  rates_vs_trait_outputs <- list()
  for (i in seq_along(time_steps))
  {
    # i <- 1

    # Extract focal time = time-step n°i
    focal_time_i <- time_steps[i]

    if (trait_data_type %in% c("categorical", "biogeographic"))
    {
      ## Adjust colors_per_levels to ensure the same color scheme is used across all time-steps
      trait_data_df_i <- trait_data_df[trait_data_df$focal_time == focal_time_i, ]
      states_in_trait_data_df_i <- levels(as.factor(trait_data_df_i$trait_value))
      colors_per_levels_i <- colors_per_levels[states_in_trait_data_df_i]

      ## Adjust select_trait_levels for the available levels
      if (!any(select_trait_levels == "all"))
      {
        select_trait_levels_i <- select_trait_levels[select_trait_levels %in% states_in_trait_data_df_i]
      } else {
        select_trait_levels_i <- "all"
      }
    } else {
      select_trait_levels_i <- NULL
      colors_per_levels_i <- NULL
    }

    ## If requested, built PDF file path for time-step n°i
    if (!is.null(PDF_file_path))
    {
      # Generate temporary PDF file name
      PDF_file_path_root <- sub(pattern = ".pdf$", replacement = "", x = PDF_file_path)
      PDF_file_path_i <- file.path(paste0(PDF_file_path_root,"_",i,".pdf"))
    } else {
      PDF_file_path_i <- NULL
    }

    ## Run plot on focal time i
    rates_vs_trait_output_i <- plot_rates_vs_trait_data_for_focal_time(
       deepSTRAPP_outputs = deepSTRAPP_outputs,
       focal_time = focal_time_i,
       rate_type = rate_type,
       select_trait_levels = select_trait_levels_i,
       color_scale = color_scale,
       colors_per_levels = colors_per_levels_i,
       display_plot = display_plot,
       PDF_file_path = PDF_file_path_i,
       return_mean_data_per_branches_df = return_mean_data_per_branches_df)

    # Store output
    rates_vs_trait_outputs[[i]] <- rates_vs_trait_output_i
  }

  ## If requested, aggregate PDF files in a unique file
  if (!is.null(PDF_file_path))
  {
    # Recreate paths to PDFs in predefined order
    all_PDF_paths <- paste0(PDF_file_path_root,"_",seq_along(time_steps),".pdf")

    # Combine PDFs in a unique PDF
    qpdf::pdf_combine(input = all_PDF_paths, output = file.path(PDF_file_path))

    # Remove temporary files
    invisible(file.remove(all_PDF_paths))
  }

  ## Extract ggplots
  rates_vs_trait_ggplots <- lapply(X = rates_vs_trait_outputs,
                                   FUN = function (x) { x$rates_vs_trait_ggplot })
  ## Build output
  output <- list()
  output$rates_vs_trait_ggplots <- rates_vs_trait_ggplots

  ## Store melted df if requested
  if (return_mean_data_per_branches_df)
  {
    # Extract melted data.frame
    mean_data_per_branches_df_list <- lapply(X = rates_vs_trait_outputs,
        FUN = function (x) { x$mean_data_per_branches_df })
    # Bind all data.frames
    mean_data_per_branches_df <- do.call(what = rbind, mean_data_per_branches_df_list)
    # Store output
    output$mean_data_per_branches_df <- as.data.frame(mean_data_per_branches_df)
  }

  ## Return output
  return(invisible(output))

}


## Possible additions

# Ensure to use the same color scheme for all time-step even if some states/ranges are missing
  # Not the case for continuous traits
# Ensure to use the same ranges for x and y even if rate/trait value ranges are changing
