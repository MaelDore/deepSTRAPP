
#' @title Plot multiple mapped phylogenies of trait/range evolution vs. diversification rates and regime shifts over time-steps
#'
#' @description Plot two mapped phylogenies with evolutionary data with branches cut off for each time-step.
#'   Main input = output of a deepSTRAPP run over time using [deepSTRAPP::run_deepSTRAPP_over_time()]).
#'
#'   * Left facet: plot the evolution of trait data/geographic ranges on the left time-calibrated phylogeny.
#'     - For continuous data: Each branch is colored according to the estimates value of the traits.
#'     - For categorical and biogeographic data: Each branch is colored according to the posterior probability
#'       of being in a given state/range. Color for each state/range are overlaid using transparency
#'       to produce a single plot for all states/ranges.
#'   * Right facet: plot the evolution of diversification rates and location of regime shits estimated from a BAMM
#'     (Bayesian Analysis of Macroevolutionary Mixtures).
#'     Each branch is colored according to the estimated rates of speciation, extinction, or net diversification
#'     stored in an object of class `bammdata`. Rates can vary along time, thus colors evolved along individual branches.
#'
#'   This function is a wrapper of [deepSTRAPP::run_deepSTRAPP_for_focal_time()] used to plot mapped phylogenies for
#'   a unique given `focal_time`.
#'
#'   Returns one plot per focal time in `$time_steps`.
#'   If a PDF file path is provided in `PDF_file_path`, the plots will be saved directly in a PDF file,
#'   with one page per focal time in `$time_steps`.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   that summarize the results of a STRAPP run over multiple time-steps in the past.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   showing the evolution of a continuous trait. From lowest values to highest values. (For continuous trait data only)
#' @param colors_per_levels Named character string. To set the colors to use to map each state/range posterior probabilities. Names = states/ranges; values = colors.
#'   If `NULL` (default), the color scale provided `densityMaps` will be used. (For categorical and biogeographic data only)
#' @param add_ACE_pies Logical. Whether to add pies of posterior probabilities of states/ranges at internal nodes on the mapped phylogeny. Default = `TRUE`.
#' @param cex_pies Numerical. To adjust the size of the ACE pies. Default = `0.5`.
#' @param rate_type A character string specifying the type of diversification rates to plot.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param keep_initial_colorbreaks Logical. Whether to keep the same color breaks as used for the most recent focal time. Typically, the current time (t = 0). Default = `FALSE`.
#' @param add_regime_shifts Logical. Whether to add the location of regime shifts on the phylogeny (Step 2). Default is `TRUE`.
#' @param configuration_type A character string specifying how to select the location of regime shifts across posterior samples.
#'   * `configuration_type = "MAP"`: Use the average locations recorded in posterior samples with the Maximum A Posteriori probability (MAP) configuration.
#'     This regime shift configuration is the most frequent configuration among the posterior samples (See [BAMMtools::getBestShiftConfiguration()]).
#'     This is the default option.
#'   * `configuration_type = "MSC"`: Use the average locations recorded in posterior samples with the Maximum Shift Credibility (MSC) configuration.
#'     This regime shift configuration has the highest product of marginal probabilities across branches (See [BAMMtools::maximumShiftCredibility()]).
#'   * `configuration_type = "index"`: Use the configuration of a unique posterior sample those index is provided in `sample_index`.
#' @param sample_index Integer. Index of the posterior samples to use to plot the location of regime shifts.
#'   Used only if `configuration_type = index`. Default = `1`.
#' @param adjust_size_to_prob Logical. Whether to scale the size of the symbols showing the location of regime shifts according to
#'   the marginal shift probability of the shift happening on each location/branch. This will only works if there is an `$MSP_tree` element
#'   summarizing the marginal shift probabilities across branches in the `BAMM_object`. Default is `TRUE`.
#' @param regimes_fill Character string. Set the color of the background of the symbols showing the location of regime shifts.
#'   Equivalent to the `bg` argument in [BAMMtools::addBAMMshifts()]. Default is `"grey"`.
#' @param regimes_size Numerical. Set the size of the symbols showing the location of regime shifts.
#'   Equivalent to the `cex` argument in [BAMMtools::addBAMMshifts()]. Default is `1`.
#' @param regimes_pch Integer. Set the shape of the symbols showing the location of regime shifts.
#'   Equivalent to the `pch` argument in [BAMMtools::addBAMMshifts()]. Default is `21`.
#' @param regimes_border_col Character string. Set the color of the border of the symbols showing the location of regime shifts.
#'   Equivalent to the `col` argument in [BAMMtools::addBAMMshifts()]. Default is `"black"`.
#' @param regimes_border_width Numerical. Set the width of the border of the symbols showing the location of regime shifts.
#'   Equivalent to the `lwd` argument in [BAMMtools::addBAMMshifts()]. Default is `1`.
#' @param ... Additional graphical arguments to pass down to [phytools::plot.contMap()], [phytools::plot.simmap()],
#'   [deepSTRAPP::plot_densityMaps_overlay()], [BAMMtools::plot.bammdata()], [BAMMtools::addBAMMshifts()], and [par()].
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#'
#' @export
#' @importFrom qpdf pdf_combine
#'
#' @details The main input `deepSTRAPP_outputs` is the typical output of [deepSTRAPP::run_deepSTRAPP_over_time()].
#'   It provides information on results of a STRAPP run performed over multiple time-steps, and can also encompass
#'   updated phylogenies with mapped trait evolution and diversification rates and regimes shifts if appropriate arguments are set.
#'
#'   * `return_updated_Maps` must be set to `TRUE` so that the updated version of mapped phylogeny (`contMap`/`densityMaps`)
#'     for the given `$time_steps` are returned among the outputs under `$updated_Maps_over_time`.
#'     The updated `contMap`/`densityMaps` consists in cutting off branches and mappings that are younger than the successive `$time_steps`.
#'     `contMap`/`densityMaps` must have been provided among the inputs to [deepSTRAPP::run_deepSTRAPP_over_time()] in the first place.
#'   * `return_updated_BAMM_objects` must be set to `TRUE` so that the list of `updated_BAMM_object` with phylogenies and mapped diversification rates
#'     cut-off at the successive `$time_steps` are returned among the outputs under `$updated_BAMM_objects_over_time`.
#'
#'   `$MAP_BAMM_object` and `$MSC_BAMM_object` elements are required in the elements of `$updated_BAMM_objects_over_time` to plot regime shift locations
#'   following the "MAP" or "MSC" `configuration_type` respectively.
#'   A `$MSP_tree` element is required to scale the size of the symbols showing the location of regime shifts according marginal shift probabilities.
#'   (If `adjust_size_to_prob = TRUE`).
#'
#'  For plotting a single time-step (i.e., `$focal_time`) at once, see [deepSTRAPP::plot_traits_vs_rates_on_phylogeny_for_focal_time()].
#'
#' @return If `display_plot = TRUE`, the function displays the successive plots with two facets in the R console:
#'  * (Left) A time-calibrated phylogeny displaying the evolution of trait/biogeographic data.
#'  * (Right) A time-calibrated phylogeny displaying diversification rates and regime shifts.
#'  A plot per time-steps is generated and displayed successively in the console.
#'
#'  If `PDF_file_path` is provided, the plots will be exported in a unique PDF file with one page per time-step.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::plot_contMap()] [phytools::plot.densityMap()] [deepSTRAPP::plot_densityMaps_overlay()] [deepSTRAPP::plot_BAMM_rates()]
#'
#' Function in deepSTRAPP needed to produce the `deepSTRAPP_outputs` as input: [deepSTRAPP::run_deepSTRAPP_over_time()]
#' Function in deepSTRAPP to a single time-step at once: [deepSTRAPP::plot_traits_vs_rates_on_phylogeny_for_focal_time()]
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  # ----- Example 1: Continuous trait ----- #
#'
#'  ## Load data
#'
#'  # Load trait df
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
#'  ## Prepare trait data
#'
#'  # Extract continuous trait data as a named vector
#'  Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#'                                      nm = Ponerinae_trait_tip_data$Taxa)
#'
#'  # Select a color scheme from lowest to highest values
#'  color_scale = c("darkgreen", "limegreen", "orange", "red")
#'
#'   # Get Ancestral Character Estimates based on a Brownian Motion model
#'  # To obtain values at internal nodes
#'  Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)
#'
#'  \donttest{ # (May take several minutes to run)
#'  # Run a Stochastic Mapping based on a Brownian Motion model
#'  # to interpolate values along branches and obtain a "contMap" object
#'  Ponerinae_contMap <- phytools::contMap(Ponerinae_tree, x = Ponerinae_cont_tip_data,
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
#'     extract_trait_data_melted_df = TRUE,
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_Maps = TRUE,
#'     return_updated_BAMM_objects = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly trait data output
#'  Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- readRDS(system.file("extdata",
#'      "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Plot updated contMap vs. updated diversification rates
#'  plot_traits_vs_rates_on_phylogeny_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     # To use the same color breaks as for t = 0 My in all BAMM plots
#'     keep_initial_colorbreaks = TRUE,
#'     color_scale = c("limegreen", "orange", "red"), # Adjust color scale on contMap
#'     legend = TRUE, labels = TRUE, # Show legend and label on BAMM plot
#'     cex = 0.7) # Adjust label size on contMap
#'     # PDF_file_path = "Updated_maps_cont_old_calib_0_40My.pdf")
#'
#'  # ----- Example 2: Biogeographic data ----- #
#'
#'  ## Load data
#'
#'  # Load phylogeny
#'  data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
#'  # Load trait df
#'  data(Ponerinae_binary_range_table, package = "deepSTRAPP")
#'
#'  # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#'  data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Prepare range data for Old World vs. New World
#'
#'  # No overlap in ranges
#'  table(Ponerinae_binary_range_table$Old_World, Ponerinae_binary_range_table$New_World)
#'
#'  Ponerinae_NO_data <- stats::setNames(object = Ponerinae_binary_range_table$Old_World,
#'                                       nm = Ponerinae_binary_range_table$Taxa)
#'  Ponerinae_NO_data <- as.character(Ponerinae_NO_data)
#'  Ponerinae_NO_data[Ponerinae_NO_data == "TRUE"] <- "O" # O = Old World
#'  Ponerinae_NO_data[Ponerinae_NO_data == "FALSE"] <- "N" # N = New World
#'  names(Ponerinae_NO_data) <- Ponerinae_binary_range_table$Taxa
#'  table(Ponerinae_NO_data)
#'
#'  colors_per_levels <- c("mediumpurple2", "peachpuff2")
#'  names(colors_per_levels) <- c("N", "O")
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run evolutionary models
#'  Ponerinae_biogeo_data <- prepare_trait_data(
#'     tip_data = Ponerinae_NO_data,
#'     trait_data_type = "biogeographic",
#'     phylo = Ponerinae_tree_old_calib,
#'     evolutionary_models = "DEC+J", # Default = "DEC" for biogeographic
#'     BioGeoBEARS_directory_path = tempdir(), # Ex: "./BioGeoBEARS_directory/"
#'     keep_BioGeoBEARS_files = FALSE,
#'     prefix_for_files = "Ponerinae_old_calib",
#'     max_range_size = 2,
#'     split_multi_area_ranges = TRUE, # Set to TRUE to use only unique areas "A" and "B"
#'     nb_simulations = 100, # Reduce to save time (Default = '1000')
#'     colors_per_levels = colors_per_levels,
#'     return_model_selection_df = TRUE,
#'     verbose = TRUE) }
#'
#'  # Load directly output
#'  data(Ponerinae_biogeo_data_old_calib, package = "deepSTRAPP")
#'
#'  ## Set for time steps of 5 My. Will generate deepSTRAPP workflows from 0 to 40 Mya.
#'  time_range <- c(0, 40)
#'  time_step_duration <- 5
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates for time-steps = 0 to 40 Mya.
#'  Ponerinae_deepSTRAPP_biogeo_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     densityMaps = Ponerinae_biogeo_data_old_calib$densityMaps,
#'     ace = Ponerinae_biogeo_data_old_calib$ace,
#'     tip_data = Ponerinae_ON_tip_data,
#'     nb_simulations = 100,
#'     trait_data_type = "biogeographic",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     uncertainty_strategy = "paired",
#'     seed = 1234, # Set seed for reproducibility
#'     alpha = 0.10, # Select a generous level of significance for the sake of the example
#'     rate_type = "net_diversification",
#'     return_perm_data = TRUE,
#'     extract_trait_data_melted_df = TRUE,
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_Maps = TRUE,
#'     return_updated_BAMM_objects = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly output
#'  Ponerinae_deepSTRAPP_biogeo_old_calib_0_40 <- readRDS(system.file("extdata",
#'      "Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rds", package = "deepSTRAPP"))
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Explore output
#'  str(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40, max.level =  1)
#'
#'  ## Plot updated contMap vs. updated diversification rates
#'  plot_traits_vs_rates_on_phylogeny_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
#'     # Adjust colors on densityMaps
#'     colors_per_levels = c("N" = "dodgerblue2", "O" = "orange"),
#'     legend = TRUE, labels = TRUE, # Show legend and label on BAMM plot
#'     cex_pies = 0.2, # Adjust size of ACE pies on densityMaps
#'     cex = 0.7) # Adjust label size on contMap
#'     # PDF_file_path = "Updated_maps_biogeo_old_calib_0_40My.pdf")
#' }
#'


plot_traits_vs_rates_on_phylogeny_over_time <- function (
    deepSTRAPP_outputs,
    color_scale = NULL,
    colors_per_levels = NULL,
    add_ACE_pies = TRUE,
    cex_pies = 0.5,
    rate_type = "net_diversification",
    keep_initial_colorbreaks = FALSE,
    add_regime_shifts = TRUE,
    configuration_type = "MAP", # MAP, MSC, or index
    sample_index = 1,
    adjust_size_to_prob = TRUE, # To adjust size of points representing regime shifts to their marginal posterior probabilities
    regimes_fill = "grey", # Replace 'bg' argument in BAMMtools::addBAMMshifts()
    regimes_size = 1, # Replace 'cex' argument in BAMMtools::addBAMMshifts()
    regimes_pch = 21, # Replace 'pch' argument in BAMMtools::addBAMMshifts()
    regimes_border_col = "black", # Replace 'col' argument in BAMMtools::addBAMMshifts()
    regimes_border_width = 1, # Replace 'lwd' argument in BAMMtools::addBAMMshifts()
    ..., # To pass down to BAMMtools::plot.bammdata(), BAMMtools::addBAMMshifts(), and par()
    display_plot = TRUE,
    PDF_file_path = NULL)
{
  ### Check input validity
  {
    ## deepSTRAPP_outputs
    # Check presence of updated trait map
    if (is.null(deepSTRAPP_outputs$updated_Maps_over_time))
    {
      stop(paste0("`deepSTRAPP_outputs` must have a `$updated_Maps_over_time` element.\n",
                  "Be sure to set `return_updated_Maps = TRUE` in [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                  "This element is needed to plot the updated `contMap`/`densityMaps` with mapped trait/range evolution."))
    }
    # Check presence of updated BAMM objects
    if (is.null(deepSTRAPP_outputs$updated_BAMM_objects_over_time))
    {
      stop(paste0("'deepSTRAPP_outputs' must have a '$updated_BAMM_objects_over_time' element.\n",
                  "Be sure to set `return_updated_BAMM_object = TRUE` in [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                  "This element is needed to plot the updated `BAMM_objects` with mapped diversification rates."))
    }

    ## display_plot OR PDF_file_path
    # Check that at least one type of output is requested
    if (!display_plot & is.null(PDF_file_path))
    {
      stop(paste0("You must request at least one option between displaying the plot (`display_plot` = TRUE), or producing a PDF (fill the `PDF_file_path` argument)."))
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

  ## Save initial par() and reassign them on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ## Extract time-steps
  time_steps <- deepSTRAPP_outputs$time_steps

  ## Loop per time-steps
  for (i in seq_along(time_steps))
  {
    # i <- 1

    # Extract focal time = time-step n°i
    focal_time_i <- time_steps[i]

    ## If requested, built PDF file path for time-step n°i
    if (!is.null(PDF_file_path))
    {
      # Generate temporary PDF file name
      PDF_file_path_root <- sub(pattern = ".pdf$", replacement = "", x = PDF_file_path)
      PDF_file_path_i <- file.path(paste0(PDF_file_path_root,"_",i,".pdf"))
    } else {
      PDF_file_path_i <- NULL
    }

    ## Plot traits vs rates for the given focal time
    plot_traits_vs_rates_on_phylogeny_for_focal_time(
      deepSTRAPP_outputs = deepSTRAPP_outputs,
      focal_time = focal_time_i,
      color_scale = color_scale,
      colors_per_levels = colors_per_levels,
      add_ACE_pies = add_ACE_pies,
      cex_pies = cex_pies,
      rate_type = rate_type,
      keep_initial_colorbreaks = keep_initial_colorbreaks,
      add_regime_shifts = add_regime_shifts,
      configuration_type = configuration_type, # MAP, MSC, or index
      sample_index = sample_index,
      adjust_size_to_prob = adjust_size_to_prob, # To adjust size of points representing regime shifts to their marginal posterior probabilities
      regimes_fill = regimes_fill, # Replace 'bg' argument in BAMMtools::addBAMMshifts()
      regimes_size = regimes_size, # Replace 'cex' argument in BAMMtools::addBAMMshifts()
      regimes_pch = regimes_pch, # Replace 'pch' argument in BAMMtools::addBAMMshifts()
      regimes_border_col = regimes_border_col, # Replace 'col' argument in BAMMtools::addBAMMshifts()
      regimes_border_width = regimes_border_width, # Replace 'lwd' argument in BAMMtools::addBAMMshifts()
      ..., # To pass down to BAMMtools::plot.bammdata(), BAMMtools::addBAMMshifts(), and par()
      display_plot = display_plot,
      PDF_file_path = PDF_file_path_i)
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
}


