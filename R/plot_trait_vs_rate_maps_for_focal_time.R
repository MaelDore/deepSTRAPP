
#' @title Plot trait/range evolution vs. diversification rates and regime shifts on phylogeny
#'
#' @description Plot two mapped phylogenies with evolutionary data with branches cut off at `focal_time`.
#'   * Left facet: plot the evolution of trait data/geographic ranges on the left time-calibrated phylogeny.
#'     - For continuous data: Each branch is colored according to the estimates value of the traits.
#'     - For categorical and biogeographic data: Each branch is colored according to the posterior probability
#'       of being in a given state/range. Color for each state/range are overlaid using transparency
#'       to produce a single plot for all states/ranges.
#'   * Right facet: plot the evolution of diversification rates and location of regime shits estimated from a BAMM
#'   (Bayesian Analysis of Macroevolutionary Mixtures).
#'   Each branch is colored according to the estimated rates of speciation, extinction, or net diversification
#'   stored in an object of class `bammdata`. Rates can vary along time, thus colors evolved along individual branches.
#'
#'   This function is a wrapper multiple plotting functions:
#'
#'   * For continuous traits: [deepSTRAPP::plot_contMap()]
#'   * For categorical and biogeographic data: [deepSTRAPP::plot_densityMaps_overlay()]
#'   * For BAMM rates and regime shifts: [deepSTRAPP::plot_BAMM_rates()]
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_for_focal_time()],
#'   that summarize the results of a STRAPP test for a specific time in the past (i.e. the `focal_time`).
#'   `deepSTRAPP_outputs` can also be extracted from the output of [deepSTRAPP::run_deepSTRAPP_over_time()] that
#'   run the whole deepSTRAPP workflow over multiple time-steps.
#' @param focal_time Numerical. (Optional) If `deepSTRAPP_outputs` comprises results over multiple time-steps
#'   (i.e., output of [deepSTRAPP::run_deepSTRAPP_over_time()], this is the time of the STRAPP test targeted for plotting.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   showing the evolution of a continuous trait. From lowest values to highest values. (For continuous trait data only)
#' @param colors_per_levels Named character string. To set the colors to use to map each state/range posterior probabilities. Names = states/ranges; values = colors.
#'   If `NULL` (default), the color scale provided `densityMaps` will be used. (For categorical and biogeographic data only)
#' @param add_ACE_pies Logical. Whether to add pies of posterior probabilities of states/ranges at internal nodes on the mapped phylogeny. Default = `TRUE`.
#' @param cex_pies Numerical. To adjust the size of the ACE pies. Default = `0.5`.
#' @param rate_type A character string specifying the type of diversification rates to plot.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
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
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with '.pdf'.
#'
#' @export
#' @importFrom graphics par
#' @importFrom grDevices pdf dev.off
#'
#' @details The main input `deepSTRAPP_outputs` is the typical output of [deepSTRAPP::run_deepSTRAPP_for_focal_time()].
#'   It provides information on results of a STRAPP test performed at a given `focal_time`, and can also encompass
#'   updated phylogenies with mapped trait evolution and diversification rates and regimes shifts if appropriate arguments are set.
#'
#'   * `return_updated_trait_data_with_Map` must be set to `TRUE` so that the trait data extracted for the given `focal_time`
#'     and the updated version of mapped phylogeny (`contMap`/`densityMaps`) are returned among the outputs under `$updated_trait_data_with_Map`.
#'     The updated `contMap`/`densityMaps` consists in cutting off branches and mappings that are younger than the `focal_time`.
#'   * `return_updated_BAMM_object` must be set to `TRUE` so that the `updated_BAMM_object` with phylogeny and mapped diversification rates
#'     cut-off at the `focal_time` are returned among the outputs under `$updated_BAMM_object`.
#'
#'   `$MAP_BAMM_object` and `$MSC_BAMM_object` elements are required in `$updated_BAMM_object` to plot regime shift locations following the
#'   "MAP" or "MSC" `configuration_type` respectively.
#'   A `$MSP_tree` element is required to scale the size of the symbols showing the location of regime shifts according marginal shift probabilities.
#'   (If `adjust_size_to_prob = TRUE`).
#'
#'   Alternatively, the main input `deepSTRAPP_outputs` can be the output of [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   providing results of STRAPP tests over multiple time-steps. In this case, you must provide a `focal_time` to select the
#'   unique time-step used for plotting.
#'   * `return_updated_trait_data_with_Map` must be set to `TRUE` so that the trait data extracted and
#'     the updated version of mapped phylogenies (`contMap`/`densityMaps`) are returned among the outputs under `$updated_trait_data_with_Map_over_time`.
#'   * `return_updated_BAMM_object` must be set to `TRUE` so that the `BAMM_objects` with phylogeny and mapped diversification rates
#'     cut-off at the specified time-steps are returned among the outputs under `$updated_BAMM_objects_over_time`.
#'
#'  For plotting all time-steps at once, see [deepSTRAPP::plot_trait_vs_rate_maps_over_time()].
#'
#' @return If `display_plot = TRUE`, the function displays a plot with two facets in the R console:
#'  * (Left) A time-calibrated phylogeny displaying the evolution of trait/biogeographic data.
#'  * (Right) A time-calibrated phylogeny displaying diversification rates and regime shifts.
#'
#'  If `PDF_file_path` is provided, the plot will be exported in a PDF file.
#'
#' @author Maël Doré
#'
#' @seealso [phytools::plot.densityMap()] [deepSTRAPP::plot_densityMaps_overlay()] [deepSTRAPP::plot_BAMM_rates()]
#'
#' Functions in deepSTRAPP needed to produce the `deepSTRAPP_outputs` as input: [deepSTRAPP::run_deepSTRAPP_for_focal_time()] [deepSTRAPP::run_deepSTRAPP_over_time()]
#' Function in deepSTRAPP to plot all time-steps at once: [deepSTRAPP::plot_trait_vs_rate_maps_over_time()]
#'
#' @examples
#' # ----- Example 1: Continuous trait ----- #
#'
#' # Load phylogeny
#' data(Ponerinae_trait_data, package = "deepSTRAPP")
#' # Load trait df
#' data(Ponerinae_tree, package = "deepSTRAPP")
#' # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#' data(Ponerinae_BAMM_object, package = "deepSTRAPP")
#'
#' ## Prepare trait data
#'
#' # Extract log(head with) data
#' Ponerinae_data_ln_HW <- setNames(object = Ponerinae_trait_data$sim_ln_HW,
#'                                  nm = Ponerinae_trait_data$Taxa)
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree, x = Ponerinae_data_ln_HW)
#'
#' # Run a Stochastic Mapping based on a Brownian Motion model
#' # to interpolate values along branches and obtain a "contMap" object
#' Ponerinae_contMap <- phytools::contMap(Ponerinae_tree, x = Ponerinae_data_ln_HW,
#'                                        res = 100, # Number of time steps
#'                                        plot = FALSE)
#' # Plot contMap = stochastic mapping of continuous trait
#' plot_contMap(Ponerinae_contMap)
#'
#' ## Set focal time to 40 Mya
#' focal_time <- 40
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#' Ponerinae_deepSTRAPP_cont_40My <- run_deepSTRAPP_for_focal_time(
#'   contMap = Ponerinae_contMap, ace = Ponerinae_ACE, tip_data = Ponerinae_data_ln_HW,
#'   trait_data_type = "continuous",
#'   BAMM_object = Ponerinae_BAMM_object,
#'   focal_time = focal_time,
#'   rate_type = "net_diversification",
#'   return_perm_data = TRUE,
#'   extract_diversification_data_melted_df = TRUE,
#'   return_updated_trait_data_with_Map = TRUE,
#'   return_updated_BAMM_object = TRUE)
#'
#' ## Explore output
#' str(Ponerinae_deepSTRAPP_cont_40My, max.level = 1)
#'
#' ## Plot updated contMap vs. updated diversification rates
#' plot_trait_vs_rate_maps_for_focal_time(
#'    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_40My,
#'    color_scale = c("limegreen", "orange", "red"), # Adjust color scale on contMap
#'    legend = TRUE, labels = TRUE, # Show legend and label on BAMM plot
#'    cex = 0.7) # Adjust label size on contMap
#'    # PDF_file_path = "Updated_maps_cont_40My.pdf") }
#'
#'
#' # ----- Example 2: Biogeographic data ----- #
#'
#' # Load phylogeny
#' data(Ponerinae_tree, package = "deepSTRAPP")
#' # Load trait df
#' data(Ponerinae_binary_range_table, package = "deepSTRAPP")
#' # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#' data(Ponerinae_BAMM_object, package = "deepSTRAPP")
#'
#' ## Prepare range data for Old World vs. New World
#'
#' # No overlap in ranges
#' table(Ponerinae_binary_range_table$Old_World, Ponerinae_binary_range_table$New_World)
#'
#' Ponerinae_ON_data <- stats::setNames(object = Ponerinae_binary_range_table$Old_World,
#'                                      nm = Ponerinae_binary_range_table$Taxa)
#' Ponerinae_ON_data <- as.character(Ponerinae_ON_data)
#' Ponerinae_ON_data[Ponerinae_ON_data == "TRUE"] <- "O" # O = Old World
#' Ponerinae_ON_data[Ponerinae_ON_data == "FALSE"] <- "N" # N = New World
#' names(Ponerinae_ON_data) <- Ponerinae_binary_range_table$Taxa
#' table(Ponerinae_ON_data)
#'
#' colors_per_levels <- c("mediumpurple2", "peachpuff2")
#' names(colors_per_levels) <- c("N", "O")
#'
#' # Load directly output of prepare_trait_data()
#' data(Ponerinae_biogeo_data, package = "deepSTRAPP")
#'
#' ## Explore output
#' str(Ponerinae_biogeo_data, 1)
#'
#' ## Set focal time to 40 Mya
#' focal_time <- 40
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#' Ponerinae_deepSTRAPP_biogeo_40My <- run_deepSTRAPP_for_focal_time(
#'     densityMaps = Ponerinae_biogeo_data$densityMaps,
#'     ace = Ponerinae_biogeo_data$ace,
#'     tip_data = Ponerinae_ON_data,
#'     trait_data_type = "biogeographic",
#'     BAMM_object = Ponerinae_BAMM_object,
#'     focal_time = focal_time,
#'     rate_type = "net_diversification",
#'     return_updated_trait_data_with_Map = TRUE,
#'     return_updated_BAMM_object = TRUE)
#'
#' ## Explore output
#' str(Ponerinae_deepSTRAPP_biogeo_10My, max.level = 1)
#'
#' ## Plot updated contMap vs. updated diversification rates
#' plot_trait_vs_rate_maps_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_40My,
#'     # Adjust colors on densityMaps
#'     colors_per_levels = c("N" = "dodgerblue2", "O" = "orange"),
#'     # Show legend and label on BAMM plot
#'     legend = TRUE, labels = TRUE,
#'     cex_pies = 0.2, # Adjust size of ACE pies on densityMaps
#'     cex = 0.7) # Adjust size of tip labels on BAMM plot
#'     # PDF_file_path = "Updated_maps_biogeo_40My.pdf") }
#'
#' # ----- Example 3: With outputs over_time ----- #
#'
#' ## Load directly outputs from run_deepSTRAPP_over_time()
#' data(Ponerinae_deepSTRAPP_biogeo_0_40, package = "deepSTRAPP")
#'
#' ## Explore output
#' str(Ponerinae_deepSTRAPP_biogeo_0_40, max.level =  1)
#'
#' ## Plot updated contMap vs. updated diversification rates for focal_time = 40My
#' plot_trait_vs_rate_maps_for_focal_time(
#'    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_0_40,
#'    focal_time = 40, # Select focal_time = 40My
#'    # Adjust colors on densityMaps
#'    colors_per_levels = c("N" = "dodgerblue2", "O" = "orange"),
#'    # Show legend and label on BAMM plot
#'    legend = TRUE, labels = TRUE,
#'    cex_pies = 0.2, # Adjust size of ACE pies on densityMaps
#'    cex = 0.7) # Adjust size of tip labels on BAMM plot
#'    # PDF_file_path = "Updated_maps_biogeo_40My.pdf")
#'


plot_trait_vs_rate_maps_for_focal_time <- function (
    deepSTRAPP_outputs,
    focal_time = NULL,
    color_scale = NULL,
    colors_per_levels = NULL,
    add_ACE_pies = TRUE,
    cex_pies = 0.5,
    rate_type = "net_diversification",
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
    if (is.null(deepSTRAPP_outputs$updated_trait_data_with_Map) & is.null(deepSTRAPP_outputs$updated_trait_data_with_Map_over_time))
    {
      stop(paste0("`deepSTRAPP_outputs` must have a `$updated_trait_data_with_Map` or `$updated_trait_data_with_Map_over_time` element.\n",
                  "Be sure to set `return_updated_trait_data_with_Map = TRUE` in [deepSTRAPP::run_deepSTRAPP_for_focal_time] or [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                  "This element is needed to plot the updated `contMap`/`densityMaps` with mapped trait/range evolution."))
    }
    # Check presence of updated BAMM_object
    if (is.null(deepSTRAPP_outputs$updated_trait_data_with_Map) & is.null(deepSTRAPP_outputs$updated_trait_data_with_Map_over_time))
    {
      stop(paste0("'deepSTRAPP_outputs' must have a '$updated_BAMM_object' or '$updated_BAMM_objects_over_time' element.\n",
                  "Be sure to set `return_updated_BAMM_object = TRUE` in [deepSTRAPP::run_deepSTRAPP_for_focal_time] or [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                  "This element is needed to plot the updated `BAMM_object` with mapped diversification rates."))
    }

    ## Identify the type of inputs
    if (is.null(deepSTRAPP_outputs$updated_trait_data_with_Map_over_time))
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
                    "For plotting all time-steps at once, please see [deepSTRAPP::plot_trait_vs_rate_maps_over_time]."))
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
      }
    }

    ## Extract updated_trait_data_with_Map & updated_BAMM_object
    if (!inputs_over_time)
    {
      # For outputs from run_deepSTRAPP_for_focal_time
      updated_trait_data_with_Map <- deepSTRAPP_outputs$updated_trait_data_with_Map
      updated_BAMM_object <- deepSTRAPP_outputs$updated_BAMM_object
    } else {
      # For outputs from run_deepSTRAPP_over_time
      focal_time_ID <- which(deepSTRAPP_outputs$time_steps == focal_time)
      updated_trait_data_with_Map <- deepSTRAPP_outputs$updated_trait_data_with_Map_over_time[[focal_time_ID]]
      updated_BAMM_object <- deepSTRAPP_outputs$updated_BAMM_objects_over_time[[focal_time_ID]]
    }

    ## Extract the type of trait
    trait_data_type <- updated_trait_data_with_Map$trait_data_type

    ## Check color_scale or colors_per_levels are provided to appropriate trait type
    if ((trait_data_type == "continuous") & !is.null(colors_per_levels))
    {
      stop(paste0("For 'continuous' traits, no `colors_per_levels`. Please provide the colors with `color_scale` if needed."))
    }
    if ((trait_data_type != "continuous") & !is.null(color_scale))
    {
      stop(paste0("For '",trait_data_type,"' data, no `color_scale`. Please provide the colors with `colors_per_levels` if needed."))
    }

    ## Check color_scale is valid
    if ((trait_data_type == "continuous") & !is.null(color_scale))
    {
      if (!all(is_color(color_scale)))
      {
        invalid_colors <- color_scale[!is_color(color_scale)]
        stop(paste0("Some color names in 'color_scale' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
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

  ## Extract focal time
  focal_time <- updated_trait_data_with_Map$focal_time

  ### Display plots
  if (display_plot)
  {
    initial_oma <- graphics::par()$oma
    graphics::par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))

    ### Plot facet A: Trait evolution

    if (trait_data_type == "continuous")
    {
      ## Case 1: Continuous traits and contMap

      plot_contMap(contMap = updated_trait_data_with_Map$contMap,
                   color_scale = color_scale,
                   ...) # May need to be filtered

    } else {

      ## Case 2: Categorical or biogeographic traits and densityMaps

      plot_densityMaps_overlay(densityMaps = updated_trait_data_with_Map$densityMaps,
                               colors_per_levels = colors_per_levels,
                               add_ACE_pies = add_ACE_pies,
                               cex_pies = cex_pies,
                               ace = NULL,
                               display_plot = TRUE,
                               PDF_file_path = NULL,
                               ...) # May need to be filtered
    }

    ### Plot facet B: BAMM rates

    plot_BAMM_rates(BAMM_object = updated_BAMM_object,
                    rate_type = rate_type,
                    method = "phylogram",
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
                    display_plot = TRUE,
                    PDF_file_path = NULL)

    ## Add a common main title across both plots
    main_title <- paste0("Focal time = ", focal_time)
    graphics::mtext(text = main_title, outer = TRUE, cex = 1.5, line = -1)

    graphics::par(mfrow = c(1, 1), oma = initial_oma)
  }

  ## Save PDF
  if (!is.null(PDF_file_path))
  {
    # Adjust width and height according to phylo
    nb_tips <- length(updated_BAMM_object$tip.label)
    width <- min(nb_tips/60*8*2, 200) # Maximum PDF size = 200 inches
    height <- width*10/8/2

    ## Open PDF
    grDevices::pdf(file = file.path(PDF_file_path),
                  width = width, height = height)

    initial_oma <- graphics::par()$oma
    top_outer_margin <- height * 0.20
    graphics::par(mfrow = c(1, 2), oma = c(0, 0, top_outer_margin, 0))

    ### Plot facet A: Trait evolution

    if (trait_data_type == "continuous")
    {
      ## Case 1: Continuous traits and contMap

      plot_contMap(contMap = updated_trait_data_with_Map$contMap,
                   color_scale = color_scale,
                   ...) # May need to be filtered

    } else {

      ## Case 2: Categorical or biogeographic traits and densityMaps

      plot_densityMaps_overlay(densityMaps = updated_trait_data_with_Map$densityMaps,
                               colors_per_levels = colors_per_levels,
                               add_ACE_pies = add_ACE_pies,
                               cex_pies = cex_pies,
                               ace = NULL,
                               display_plot = TRUE,
                               PDF_file_path = NULL,
                               ...) # May need to be filtered
    }

    ### Plot facet B: BAMM rates

    plot_BAMM_rates(BAMM_object = updated_BAMM_object,
                    rate_type = rate_type,
                    method = "phylogram",
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
                    display_plot = TRUE,
                    PDF_file_path = NULL)

    ## Add a common main title across both plots
    main_title <- paste0("Focal time = ", focal_time)
    base_cex = width/12 # Scale the size of the text relative to PDF size so it appears constant
    graphics::mtext(text = main_title, outer = TRUE, cex = base_cex*1.3, line = 0)

    graphics::par(mfrow = c(1, 1), oma = initial_oma)

    ## Close PDF
    invisible(grDevices::dev.off())
  }
}


