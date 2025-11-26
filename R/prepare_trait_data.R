## Functions to prepare trait data by mapping their evolution on the phylogeny
# One master function to select the proper pipeline according to data type
# Three sub-functions mapping trait evolution according to data type

#' @title Map trait evolution on a time-calibrated phylogeny
#'
#' @description Map trait evolution on a time-calibrated phylogeny in several steps:
#'
#'   * Step 1: Fit evolutionary models to trait data using Maximum Likelihood.
#'   * Step 2: Select the best fitting model comparing AICc.
#'   * Step 3: Infer ancestral characters estimates (ACE) at nodes.
#'   * Step 4: Run stochastic mapping simulations to generate evolutionary histories
#'     compatible with the best model and inferred ACE. (Only for categorical and biogeographic data)
#'   * Step 5: Infer ancestral states along branches.
#'     - For continuous traits: use interpolation to produce a `contMap`.
#'     - For categorical and biogeographic data: compute posterior frequencies of each state/range
#'       to produce a `densityMap` for each state/range.
#'
#' @param tip_data Named numerical or character string vector of trait values/states/ranges at tips.
#'   Names should be ordered as the tip labels in the phylogeny found in `phylo$tip.label`.
#'   For biogeographic data, ranges should follow the coding scheme of BioGeoBEARS with a unique CAPITAL letter per unique areas
#'   (ex: A, B), combined to form multi-area ranges (Ex: AB). Alternatively, you can provide tip_data as a matrix or data.frame of
#'   binary presence/absence in each area (coded as unique CAPITAL letter). In this case, columns are unique areas, rows are taxa,
#'   and values are integer (0/1) signaling absence or presence of the taxa in the area.
#' @param trait_data_type Character string. Type of trait data. Either: "continuous", "categorical" or "biogeographic".
#' @param phylo Time-calibrated phylogeny. Object of class `"phylo"` as defined in [ape].
#'   Tip labels (`phylo$tip.label`) should match names in `tip_data`.
#' @param seed Integer. Set the seed to ensure reproducibility. Default is `NULL` (a random seed is used).
#' @param evolutionary_models (Vector of) character string(s). To provide the set of evolutionary models to fit on the data.
#'   * Models available for continuous data are detailed in [geiger::fitContinuous()]. Default is `"BM"`.
#'   * Models available for categorical data are detailed in [geiger::fitDiscrete()]. Default is `"ARD"`.
#'   * Models for biogeographic data are fit with R package `BioGeoBEARS` using `BioGeoBEARS::bears_optim_run()`. Default is `"DEC"`.
#'   * See list in "Details" section.
#' @param Q_matrix Custom Q-matrix for categorical data representing transition classes between states.
#'   Transitions with similar integers are estimated with a shared rate parameter.
#'   Transitions with `0` represent rates that are fixed to zero (i.e., impossible transitions).
#'   Diagonal must be populated with `NA`. `row.names(Q_matrix)` and `col.names(Q_matrix)` are the states.
#'   Provide `"matrix"` among the model listed in 'evolutionary_models' to use the custom Q-matrix for modeling. Only for categorical data.
#' @param BioGeoBEARS_directory_path Character string. The path to the directory used to store input/output files generated
#'   for/by BioGeoBEARS during biogeographic historical inferences. Use '/' to separate directory and sub-directories. It must end with '/'.
#'   Default is `"./BioGeoBEARS_directory/"`. Only for biogeographic data.
#' @param keep_BioGeoBEARS_files Logical. Whether the `BioGeoBEARS_directory` and its content should be kept after the run. Default = `TRUE`. Only for biogeographic data.
#' @param prefix_for_files Character string. Prefix to add to all BioGeoBEARS files stored in the `BioGeoBEARS_directory_path` if `keep_BioGeoBEARS_files = TRUE`.
#'   Files will be exported such as 'prefix_*' with an underscore separating the prefix and the file name.
#'   Default is `NULL` (no prefix is added). Only for biogeographic data.
#' @param nb_cores Interger. Number of cores to use for parallel computation during BioGeoBEARS runs. Default = `1`. Only for biogeographic data.
#' @param max_range_size Integer. Maximum number of unique areas encompassed by multi-range areas. Default = `2`. Only for biogeographic data.
#' @param split_multi_area_ranges Logical. Whether to split multi-area ranges across unique areas when mapping ranges.
#'   Ex: For range EW, posterior probabilities will be split equally between Eastern Palearctic (E) and Western Palearctic (W).
#'   Default = `FALSE`. Only for biogeographic data.
#' @param ... Additional arguments to be passed down to the functions used to fit models (See `evolutionary_models`)
#'   and produce simmaps with [phytools::make.simmap()] or `BioGeoBEARS::runBSM()`.
#' @param res Integer. Define the number of time steps used to interpolate/estimate trait value/state/range in `contMap`/`densityMaps`. Default = `100`.
#' @param nb_simulations Integer. Define the number of simulations generated for stochastic mapping. Default = `1000`. Only for "categorical" and "biogeographic" data.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   showing the evolution of a continuous trait on the `contMap`. From lowest values to highest values. Only for continuous data. Default = `NULL` will use the rainbow() color palette.
#' @param colors_per_levels Named character string. To set the colors to use to map each state/range posterior probabilities. Names = states/ranges; values = colors.
#'   If `NULL` (default), the `rainbow()` color scale will be used for categorical trait; the `BioGeoBEARS::get_colors_for_states_list_0based()` will be used for biogeographic ranges.
#'   If no color is provided for multi-area ranges, they will be interpolated based on the colors provided for unique areas. Only for categorical and biogeographic data.
#' @param plot_map Logical. Whether to plot or not the phylogeny with mapped trait evolution. Default = `TRUE`.
#' @param plot_overlay Logical. If `TRUE` (default), plot a unique `densityMap` with overlapping states/ranges using transparency.
#'   If `FALSE`, plot a `densityMap` per state/range. Only for "categorical" and "biogeographic" data.
#' @param add_ACE_pies Logical. Whether to add pies of posterior probabilities of states/ranges at internal nodes on the mapped phylogeny. Default = `TRUE`.
#'   Only for categorical and biogeographic data.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#' @param return_ace Logical. Whether the named vector of ancestral characters estimates (ACE) at internal nodes should be returned in the output. Default = `TRUE`.
#' @param return_BSM Logical. (Only for Biogeographic data) Whether the summary tables of anagenetic and cladogenetic events generated during the Biogeographic Stochastic Mapping (BSM)
#'  process should be returned in the output. Default = `FALSE`.
#' @param return_simmaps Logical. Whether the evolutionary histories simulated during stochastic mapping (i.e., `simmaps`) should be returned in the output.
#'  Default = `TRUE`. Only for "categorical" and "biogeographic" data.
#' @param return_best_model_fit Logical. Whether to include the output of the best fitting model in the function output. Default = `FALSE`.
#' @param return_model_selection_df Logical. Whether to include the data.frame summarizing model comparisons used to select the best fitting model should be returned in the output. Default = `FALSE`.
#' @param verbose Logical. Should progression be displayed? A message will be printed for every steps in the process. Default is `TRUE`.
#'
#' @export
#' @importFrom stats logLik
#' @importFrom geiger fitContinuous
#' @importFrom ape all.equal.phylo
#' @importFrom phytools rescale fastAnc contMap densityMap rescaleSimmap as.Qmatrix setMap
#' @importFrom grDevices pdf dev.off rainbow colorRampPalette col2rgb rgb rgb2hsv hsv
#' @importFrom methods hasArg
#' @importFrom stats setNames
#' @importFrom BioGeoBEARS np define_tipranges_object fix_BioGeoBEARS_params_minmax check_BioGeoBEARS_run bears_optim_run get_colors_for_states_list_0based
#' @importFrom BioGeoBEARS get_inputs_for_stochastic_mapping get_LnL_from_BioGeoBEARS_results_object get_Qmat_COOmat_from_BioGeoBEARS_run_object get_sum_statetime_on_branch prt
#'
#' @details Map trait evolution on a time-calibrated phylogeny in several steps:
#'
#'  Step 1: Models are fit using Maximum Likelihood approach:
#'    * For "continuous" data models are fit with [geiger::fitContinuous()]: "BM", "OU", "EB", "rate_trend", "lambda", "kappa", "delta". Default is `"BM"`.
#'    * For "categorical" data models are fit with [geiger::fitDiscrete()]: "ER", "SYM", "ARD", "meristic", "matrix". Default is `"ARD"`.
#'    * For "biogeographic" data models are fit with R package `BioGeoBEARS`: "BAYAREALIKE", "DIVALIKE", "DEC", "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J". Default is `"DEC"`.
#'
#'  Step 2: Best model is identified among the list of `evolutionary_models` by comparing the corrected AIC (AICc)
#'    and selecting the  model with lowest AICc.
#'
#'  Step 3: For continuous traits: Ancestral characters estimates (ACE) are inferred with [phytools::fastAnc] on a tree
#'    with modified branch lengths scaled to reflect the evolutionary rates estimated from the best model using [phytools::rescale()].
#'
#'  Step 4: Stochastic Mapping.
#'
#'    For categorical and biogeographic data, stochastic mapping simulations are performed to generate evolutionary histories
#'    compatible with the best model and inferred ACE. Node states/ranges are drawn from the scaled marginal likelihoods of ACE,
#'    and states/ranges shifts along branches are simulated according to the transition matrix Q estimated from the best fitting model.
#'
#'  Step 5: Infer ancestral states along branches.
#'    * For continuous traits: ancestral trait values along branches are interpolated with [phytools::contMap()].
#'      This provides quick estimates of trait value at any point in time, but it does not provide accurate ML estimates in
#'      case of models that are time or trait-value dependent (such as "EB" or "OU") as the interpolation used to built the contMap is assuming
#'      a constant rate along each branch. However, ancestral trait values at nodes remain accurate
#'    * For categorical and biogeographic data: compute posterior frequencies of each state/range among the simulated evolutionary histories (`simmaps`)
#'      to produce a `densityMap` for each state/range that reflects the changes along branches in probability of harboring a given state/range.
#'
#'  # Note on macroevolutionary models of trait evolution
#'
#'  This function provides an easy solution to map trait evolution on a time-calibrated phylogeny
#'  and obtain the `contMap`/`densityMaps` objects needed to run the deepSTRAPP workflow ([run_deepSTRAPP_for_focal_time], [run_deepSTRAPP_over_time]).
#'  However, it does not explore the most complex options for trait evolution. You may need to explore more complex models to capture the dynamics of trait evolution.
#'  such as trait-dependent multi-rate models ([phytools::brownie.lite()], [OUwie::OUwie]), Bayesian MCMC implementations allowing a thorough exploration
#'  of location and number of regime shifts (Ex: BayesTraits, RevBayes), or RRphylo for a penalized phylogenetic ridge regression approach that allows regime shifts across all branches.
#'
#'  # Note on macroevolutionary models of biogeographic history
#'
#'  This function provides an easy solution to infer ancestral geographic ranges using the BioGeoBEARS framework.
#'  It allows to directly compare model fits across 6 models: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J.
#'  It uses a step-wise approach by using MLE estimates of previous runs as staring parameter values when increasing complexity (adding +J parameters).
#'  However, it does not explore the most complex options for historical biogeography. You may need to explore more complex models
#'  to capture the dynamics of range evolution such as time-stratification with adjacency matrices, dispersal multipliers (+W),
#'  distance-based dispersal probabilities (+X), or other features.
#'  See for instance, \href{http://phylo.wikidot.com/biogeobears}{http://phylo.wikidot.com/biogeobears}).
#'
#'  The R package `BioGeoBEARS` is needed for this function to work with biogeographic data.
#'  Please install it manually from: \href{https://github.com/nmatzke/BioGeoBEARS}{https://github.com/nmatzke/BioGeoBEARS}.
#'
#' @return The function returns a list with at least two elements.
#'
#'   * `$contMap` (For "continuous" data) Object of class `"contMap"`, typically generated with [phytools::contMap()],
#'     that contains a phylogenetic tree and associated continuous trait mapping.
#'   * `$densityMaps` (For "categorical" and "biogeographic" data) List of objects of class `"densityMap`,
#'     typically generated with [phytools::densityMap()], that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state/range along branches. The list contains one `"densityMap` per state/range found in the `tip_data`.
#'   * `$densityMaps_all_ranges` (For "biogeographic" data only, if `split_multi_area_ranges = TRUE`) Same as `$densityMaps`,
#'     but for all ranges including the multi-areas ranges (e.g., AB) while `$densityMaps` will display posterior probabilities for
#'     unique areas only (e.g., A and B), with multi-areas ranges split across the unique areas they encompass.
#'   * `$trait_data_type` Character string. Record the type of trait data. Either: "continuous", "categorical" or "biogeographic".
#'
#'   If `return_ace = TRUE`,
#'   * `$ace` For continuous traits: Named vector that record the ancestral characters estimates (ACE) at internal nodes.
#'     For categorical and biogeographic data: Matrix that record the posterior probabilities of ancestral states/ranges (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   * `$ace_all_ranges` For biogeographic data, if `split_multi_area_ranges = TRUE`: Named vector that record the ancestral characters estimates (ACE) at internal nodes,
#'     but including all ranges observed in the simmaps, including the multi-areas ranges (e.g., AB), while `$ace` will display posterior probabilities
#'     for unique areas only (e.g., A and B), with multi-areas ranges split across the unique areas they encompass.
#'
#'  If `return_BSM = TRUE`, (Only for biogeographic data)
#'   * `$BSM_output` List of two lists that contains summary information of cladogenetic (`$RES_caldo_events_tables`) and anagenetic (`$RES_ana_events_tables`) events
#'     recording across the N simulations of biogeographic histories performed during Biogeographic Stochastic Mapping (BSM).
#'     Each element of the list is a data.frame recording events occurring during one simulation.
#'
#'   If `return_simmaps = TRUE`, (Only for categorical and biogeographic data)
#'   * `$simmaps` List that contains as many objects of class `"simmap"` that `nb_simulations` were requested.
#'     Each simmap object is a phylogeny with one simulated geographic history (i.e., transitions in geographic ranges) mapped along branches.
#'
#'   If `return_best_model_fit = TRUE`,
#'   * `$best_model_fit` List that provides the output of the best fitting model.
#'
#'   If `model_selection_df = TRUE`,
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#'   For biogeographic data, the function also produce input and output files associated with BioGeoBEARS and stored in the directory
#'   specified in `BioGeoBEARS_directory_path`. The directory and its content are kept if `keep_BioGeoBEARS_files = TRUE`
#'
#' @author Maël Doré
#'
#' @seealso [geiger::fitContinuous()] [geiger::fitDiscrete()] [phytools::contMap()] [phytools::densityMap()]
#'
#' For a guided tutorial, see the associated vignettes:
#' * For continuous trait data: \code{vignette("model_continuous_trait_evolution", package = "deepSTRAPP")}
#' * For categorical trait data: \code{vignette("model_categorical_trait_evolution", package = "deepSTRAPP")}
#' * For biogeographic range data: \code{vignette("model_biogeographic_range_evolution", package = "deepSTRAPP")}
#'
#' @references For macroevolutionary models in geiger: Pennell, M. W., Eastman, J. M., Slater, G. J., Brown, J. W., Uyeda, J. C., FitzJohn, R. G., ... & Harmon, L. J. (2014).
#'  geiger v2. 0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics, 30(15), 2216-2218..
#'  \url{https://doi.org/10.1093/bioinformatics/btu181}.
#'
#'  For BioGeoBEARS: Matzke, Nicholas J. (2018). BioGeoBEARS: BioGeography with Bayesian (and likelihood) Evolutionary Analysis with R Scripts.
#'    version 1.1.1, published on GitHub on November 6, 2018. DOI: \url{http://dx.doi.org/10.5281/zenodo.1478250}. Website: \url{http://phylo.wikidot.com/biogeobears}.
#'
#' @examples
#' # ----- Example 1: Continuous data ----- #
#'
#' ## Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract body size
#' eel_data <- stats::setNames(eel.data$Max_TL_cm,
#'                             rownames(eel.data))
#'
#' ## Map trait evolution on the phylogeny
#' mapped_cont_traits <- prepare_trait_data(
#'    tip_data = eel_data,
#'    trait_data_type = "continuous",
#'    phylo = eel.tree,
#'    evolutionary_models = c("BM", "OU", "lambda", "kappa"),
#'    # Example of an additional argument ('control') that can be provided to geiger::fitContinuous()
#'    control = list(niter = 200),
#'    color_scale = c("darkgreen", "limegreen", "orange", "red"),
#'    plot_map = FALSE,
#'    return_best_model_fit = TRUE,
#'    return_model_selection_df = TRUE,
#'    verbose = TRUE)
#'
#' ## Explore output
#' plot_contMap(mapped_cont_traits$contMap) # contMap with interpolated trait values
#' mapped_cont_traits$model_selection_df # Summary of model selection
#' # Parameter estimates and optimization summary of the best model
#' # (Here, the best model is Pagel's lambda)
#' mapped_cont_traits$best_model_fit$opt
#' mapped_cont_traits$ace # Ancestral character estimates at internal nodes
#'
#' # ----- Example 2: Categorical data ----- #
#'
#' ## Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Transform feeding mode data into a 3-level factor
#' eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
#' eel_data <- as.character(eel_data)
#' eel_data[c(1, 5, 6, 7, 10, 11, 15, 16, 17, 24, 25, 28, 30, 51, 52, 53, 55, 58, 60)] <- "kiss"
#' eel_data <- stats::setNames(eel_data, rownames(eel.data))
#' table(eel_data)
#'
#' # Manually define a Q_matrix for rate classes of state transition to use in the 'matrix' model
#' # Does not allow transitions from state 1 ("bite") to state 2 ("kiss") or state 3 ("suction")
#' # Does not allow transitions from state 3 ("suction") to state 1 ("bite")
#' # Set symmetrical rates between state 2 ("kiss") and state 3 ("suction")
#' Q_matrix = rbind(c(NA, 0, 0), c(1, NA, 2), c(0, 2, NA))
#'
#' # Set colors per states
#' colors_per_states <- c("limegreen", "orange", "dodgerblue")
#' names(colors_per_states) <- c("bite", "kiss", "suction")
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run evolutionary models
#' eel_cat_3lvl_data <- prepare_trait_data(tip_data = eel_data, phylo = eel.tree,
#'     trait_data_type = "categorical",
#'     colors_per_levels = colors_per_states,
#'     evolutionary_models = c("ER", "SYM", "ARD", "meristic", "matrix"),
#'     Q_matrix = Q_matrix,
#'     nb_simulations = 1000,
#'     plot_map = TRUE,
#'     plot_overlay = TRUE,
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE) }
#'
#' # Load directly output
#' data(eel_cat_3lvl_data, package = "deepSTRAPP")
#'
#' ## Explore output
#' plot(eel_cat_3lvl_data$densityMaps[[1]]) # densityMap for state n°1 ("bite")
#' eel_cat_3lvl_data$model_selection_df # Summary of model selection
#' # Parameter estimates and optimization summary of the best model
#' # (Here, the best model is ER)
#' print(eel_cat_3lvl_data$best_model_fit)$ # Summary of the best evolutionary model
#' eel_cat_3lvl_data$ace # Posterior probabilities of each state (= ACE) at internal nodes
#'
#'
#' # ----- Example 3: Biogeographic data ----- #
#'
#' ## The R package `BioGeoBEARS` is needed for this function to work with biogeographic data.
#' Please install it manually from: \href{https://github.com/nmatzke/BioGeoBEARS}{https://github.com/nmatzke/BioGeoBEARS}.
#'
#' ## Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Transform feeding mode data into biogeographic data with ranges A, B, and AB.
#' eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
#' eel_data <- as.character(eel_data)
#' eel_data[eel_data == "bite"] <- "A"
#' eel_data[eel_data == "suction"] <- "B"
#' eel_data[c(5, 6, 7, 15, 25, 32, 33, 34, 50, 52, 57, 58, 59)] <- "AB"
#' eel_data <- stats::setNames(eel_data, rownames(eel.data))
#' table(eel_data)
#'
#' colors_per_ranges <- c("dodgerblue3", "gold")
#' names(colors_per_ranges) <- c("A", "B")
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run evolutionary models
#' eel_biogeo_data <- prepare_trait_data(
#'     tip_data = eel_data,
#'     trait_data_type = "biogeographic",
#'     phylo = eel.tree,
#'     # Default = "DEC" for biogeographic
#'     evolutionary_models = c("BAYAREALIKE", "DIVALIKE", "DEC",
#'                             "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J"),
#'     prefix_for_files = "eel",
#'     max_range_size = 2,
#'     split_multi_area_ranges = TRUE, # Set to TRUE to display the two outputs
#'     # Reduce the number of Stochastic Mapping simulations to save time (Default = '1000')
#'     nb_simulations = 100,
#'     colors_per_levels = colors_per_ranges,
#'     return_simmaps = TRUE,
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE,
#'     verbose = TRUE) }
#'
#' # Load directly output
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' ## Explore output
#' str(eel_biogeo_data, 1)
#' eel_biogeo_data$model_selection_df # Summary of model selection
#' # Parameter estimates and optimization summary of the best model
#' # (Here, the best model is DEC+J)
#' eel_biogeo_data$best_model_fit$optim_result
#'
#' # Posterior probabilities of each state (= ACE) at internal nodes
#' eel_biogeo_data$ace # Only with unique areas
#' eel_biogeo_data$ace_all_ranges # Including multi-area ranges (Here, AB)
#'
#' ## Plot densityMaps
#' # densityMap for range n°1 ("A")
#' plot(eel_biogeo_data$densityMaps[[1]])
#' # densityMaps with all unique areas overlaid
#' plot_densityMaps_overlay(eel_biogeo_data$densityMaps)
#' # densityMaps with all ranges (including multi-area ranges) overlaid
#' plot_densityMaps_overlay(eel_biogeo_data$densityMaps_all_ranges)
#'


### Master function to prepare data and select the proper test function according to data type ####

prepare_trait_data <- function (
    tip_data,
    trait_data_type,
    phylo,
    seed = NULL,
    evolutionary_models = NULL, # Default = "BM" for continuous data; "ARD" for categorical; "DEC" for biogeographic
    Q_matrix = NULL, # Custom Q-matrix for categorical data
    BioGeoBEARS_directory_path = "./BioGeoBEARS_directory/",
    keep_BioGeoBEARS_files = TRUE,
    prefix_for_files = NULL,
    nb_cores = 1,
    max_range_size = 2,
    split_multi_area_ranges = FALSE,
    ..., # To allow to pass down arguments in the functions used to fit the models
    res = 100, # Number of time steps used to interpolate trait value in the contMap
    nb_simulations = 1000, # Only for categorical and biogeographic data
    color_scale = NULL, # Only for continuous data
    colors_per_levels = NULL, # Only for categorical and biogeographic data. To set the colors to use to map each state/range posterior probabilities
    plot_map = TRUE,
    plot_overlay = TRUE, # Only for categorical and biogeographic data
    add_ACE_pies = TRUE, # Only for categorical and biogeographic data
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_BSM = FALSE, # Only for biogeographic data
    return_simmaps = FALSE, # Only for categorical and biogeographic data
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE)
{
  ### Check input validity
  {
    ## tip_data
    if (!any(c(is.matrix(tip_data), is.data.frame(tip_data))))
    {
      # tip_data must have the same names as the tip.label in the phylogeny
      if (!all((names(tip_data) %in% phylo$tip.label)))
      {
        stop(paste0("Names in 'tip_data' must match with '$tip.label' in 'phylo'."))
      }
      # Check if tip_data need to be reorder to match phylo$tip.label. If so, do it, but send a warning
      if (!all(names(tip_data) == phylo$tip.label))
      {
        tip_data <- tip_data[match(phylo$tip.label, table = names(tip_data))]
        warning(paste0("Entries in 'tip_data' were reordered to match 'phylo$tip.label."))
      }
    } else {
      # row.names(tip_data) must have the same names as the tip.label in the phylogeny
      if (!all((row.names(tip_data) %in% phylo$tip.label)))
      {
        stop(paste0("Names in 'tip_data' must match with '$tip.label' in 'phylo'."))
      }
      # Check if tip_data need to be reorder to match phylo$tip.label. If so, do it, but send a warning
      if (!all(row.names(tip_data) == phylo$tip.label))
      {
        tip_data <- tip_data[match(phylo$tip.label, table = row.names(tip_data)), ]
        warning(paste0("Entries in 'tip_data' were reordered to match 'phylo$tip.label."))
      }
    }


    ## trait_data_type
    # trait_data_type must be "continuous", categorical" or "biogeographic"
    if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
    {
      stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
    }
    # Check that trait_data_type is coherent with what is found in tip_data
    if ((trait_data_type == "continuous") & !is.numeric(tip_data))
    {
      stop("For 'trait_data_type = continuous', 'tip_data' must be a named numeric vector.")
    }
    if ((trait_data_type == "categorical") & !is.character(tip_data))
    {
      if (is.factor(tip_data))
      {
        cat("WARNING: 'tip_data' was provided as factors. It is converted to a vector of character strings.\n")

        tip_data_names <- names(tip_data)
        tip_data <- as.character(tip_data)
        names(tip_data) <- tip_data_names

      } else {
        stop("For 'trait_data_type = categorical', 'tip_data' must be a named vector of character strings providing tip states.")
      }
    }
    if ((trait_data_type == "biogeographic") & !any(c(is.character(tip_data), is.matrix(tip_data), is.data.frame(tip_data))))
    {
      stop("For 'trait_data_type = biogeographic', 'tip_data' must be a named vector of character strings providing tip ranges.\n",
           "Alternatively, you can provide a matrix/data.frame providing presence/absence binary information (0/1).\n",
           "In this case, colnames(tip_data) must be unique areas symbolized by a unique UPPERCASE letter each. row.names(tip_data) must be taxa named as in the phylogeny. Data must be numerical 0/1.\n",
           "Taxa with multi-area range should record several presences (1) in their row.")
    }

    ## phylo
    # phylo must be a "phylo" class object
    if (!("phylo" %in% class(phylo)))
    {
      stop("'phylo' must have the 'phylo' class. See ?ape::read.tree() and ?ape::read.nexus() to learn how to import phylogenies in R.")
    }
    # phylo must be rooted
    if (!(ape::is.rooted(phylo)))
    {
      stop(paste0("'phylo' must be a rooted phylogeny."))
    }
    # phylo must be fully resolved/dichotomous
    if (!(ape::is.binary(phylo)))
    {
      stop(paste0("'phylo' must be a fully resolved/dichotomous/binary phylogeny."))
    }

    ## seed
    if (!is.null(seed))
    {
      if (!is.numeric(seed))
      {
        stop(paste0("'seed' must be an interger."))
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

    ## prefix_for_files
    # If provided, prefix_for_files should be character string
    if (!is.null(prefix_for_files))
    {
      if (!is.character(prefix_for_files))
      {
        stop(paste0("'prefix_for_files' must be a character string.\n",
                    "Files will exported such as 'prefix_*' with an underscore separating the prefix and the file name."))
      }
    }

    ## Other checks are carried in dedicated sub-functions
  }

  ## Filter list of additional arguments to avoid warnings from functions
  add_args <- list(...)
  args_names_for_fitContinuous <- c("SE", "bounds", "control", "ncores")
  args_names_for_fitDiscrete <- c("transform", "bounds", "control", "ncores", "symmetric")
  args_names_for_make.simmap <- c("pi", "message", "Q", "tol", "vQ", "prior")
  args_names_to_define_BioGeoBEARS_run <- c("abbr", "description", "BioGeoBEARS_model_object", "timesfn", "distsfn",
                                             "dispersal_multipliers_fn", "area_of_areas_fn", "areas_allowed_fn",
                                             "detects_fn", "controls_fn", "states_list", "force_sparse",
                                             "use_detection_model", "tmpwd", "cluster_already_open", "use_optimx",
                                             "calc_TTL_loglike_from_condlikes_table", "calc_ancprobs", "fixnode",
                                             "fixlikes", "speedup")
  args_names_to_define_runBSM <- c("maxtries_per_branch", "save_after_every_try", "seedval", "wait_before_save", "master_nodenum_toPrint")

  ## Extract additional arguments
  args_for_fitContinuous <- add_args[names(add_args) %in% args_names_for_fitContinuous]
  args_for_fitDiscrete <- add_args[names(add_args) %in% args_names_for_fitDiscrete]
  args_for_make.simmap <- add_args[names(add_args) %in% args_names_for_make.simmap]
  args_to_define_BioGeoBEARS_run <- add_args[names(add_args) %in% args_names_to_define_BioGeoBEARS_run]
  args_to_define_runBSM <- add_args[names(add_args) %in% args_names_to_define_runBSM]

  ## Set seed
  if (!is.null(seed))
  {
    set.seed(seed = seed)
  }

  ## Compute the appropriate internal sub-function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous =   { # Case for continuous data
           # Output = contMap. Trait values are interpolated along branches.
           trait_data_output <- prepare_trait_data_for_continuous_data(
             tip_data = tip_data,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "BM" for continuous data
             # ..., # Additional arguments for geiger::fitContinuous()
             args_for_fitContinuous = args_for_fitContinuous,
             res = res,
             color_scale = color_scale,
             plot_map = plot_map,
             PDF_file_path = PDF_file_path,
             return_ace = return_ace,
             return_best_model_fit = return_best_model_fit,
             return_model_selection_df = return_model_selection_df,
             verbose = verbose
           )
         },
         categorical =  { # Case for categorical data
           # Output = densityMap (+ simmaps). Include simulations of trait states.
           trait_data_output <- prepare_trait_data_for_categorical_data(
             tip_data = tip_data,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "ARD" for categorical data
             Q_matrix = Q_matrix, # Custom Q-matrix for categorical data
             # ..., # Additional arguments for geiger::fitDiscrete() and phytools::make.simmap
             args_for_fitDiscrete = args_for_fitDiscrete, # Additional arguments for geiger::fitDiscrete()
             args_for_make.simmap = args_for_make.simmap, # Additional arguments for phytools::make.simmap()
             res = res,
             nb_simulations = nb_simulations, # Only for categorical and biogeographic data
             colors_per_levels = colors_per_levels, # Only for categorical and biogeographic data
             plot_map = plot_map,
             plot_overlay = plot_overlay, # Only for categorical and biogeographic data
             add_ACE_pies = add_ACE_pies, # Only for categorical and biogeographic data
             PDF_file_path = PDF_file_path,
             return_ace = return_ace,
             return_simmaps = return_simmaps, # Only for categorical and biogeographic data
             return_best_model_fit = return_best_model_fit,
             return_model_selection_df = return_model_selection_df,
             verbose = verbose
           )
         },
         biogeographic = { # Case for biogeographic data
           # Output = densityMap (+ simmaps). Include simulations of geographic ranges.
           trait_data_output <- prepare_trait_data_for_biogeographic_data(
             tip_data = tip_data,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "DEC" for biogeographic data
             BioGeoBEARS_directory_path = BioGeoBEARS_directory_path,
             keep_BioGeoBEARS_files = keep_BioGeoBEARS_files,
             prefix_for_files = prefix_for_files,
             nb_cores = nb_cores,
             max_range_size = max_range_size,
             split_multi_area_ranges = split_multi_area_ranges,
             # ..., # Additional arguments for BioGeoBEARS functions
             args_to_define_BioGeoBEARS_run = args_to_define_BioGeoBEARS_run, # Additional arguments for BioGeoBEARS::define_BioGeoBEARS_run()
             args_to_define_runBSM, # Additional arguments for BioGeoBEARS::runBSM()
             nb_simulations = nb_simulations, # Only for categorical and biogeographic data
             res = res,
             colors_per_levels = colors_per_levels, # Only for categorical and biogeographic data
             plot_map = plot_map,
             plot_overlay = plot_overlay, # Only for categorical and biogeographic data
             PDF_file_path = PDF_file_path,
             return_ace = return_ace,
             return_BSM = return_BSM, # Only for biogeographic data
             return_simmaps = return_simmaps, # Only for categorical and biogeographic data
             return_best_model_fit = return_best_model_fit,
             return_model_selection_df = return_model_selection_df,
             verbose = verbose
           )
         }
  )

  ## Export the output
  return(invisible(trait_data_output))
}


### Sub-function to handle continuous data ####

prepare_trait_data_for_continuous_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "BM", # Default = "BM" for continuous data
    # ..., # Additional arguments for geiger::fitContinuous()
    args_for_fitContinuous = NULL, # Additional arguments for geiger::fitContinuous()
    res = 100,
    color_scale = NULL,
    plot_map = TRUE,
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE
    )
{
  ### Check input validity
  {
    ## evolutionary_models
    if (!is.null(evolutionary_models))
    {
      if (!all(evolutionary_models %in% c("BM", "OU", "EB", "rate_trend", "lambda", "kappa", "delta")))
      {
        stop(paste0("For 'trait_data_type = continuous', 'evolutionary_models' must be selected among: 'BM', 'OU', 'EB', 'rate_trend', 'lambda', 'kappa', 'delta'.\n",
                    "See details in ?geiger::fitContinuous()."))
      }
    }

    ## res
    # res must be a positive integer.
    if ((res != abs(res)) | (res != round(res)))
    {
      stop(paste0("'res' must be a positive integer defining the number of time steps used to interpolate trait value in the contMap."))
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

  # Set default model (BM) if absent
  if (is.null(evolutionary_models))
  {
    evolutionary_models <- "BM"
    message("WARNING: No 'evolutionary_model' is provided. Use 'BM' as the default.\n")
  }

  ### Step 1 - Run all models and store their results in a list

  # ?geiger::fitContinuous

  nb_models <- length(evolutionary_models)
  if (verbose) { cat(paste0("\n", Sys.time(), " - Fit ",nb_models," evolutionary model(s): ", paste(evolutionary_models, collapse = ", "), ".\n\n")) }

  # Initiate list to store models outputs
  models_fits <- list()

  # Fit BM model
  if (any(evolutionary_models == "BM"))
  {
    BM_ID <- which(evolutionary_models == "BM")
    # BM_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "BM", ...)
    BM_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "BM"), args_for_fitContinuous))

    if (verbose) { print(BM_fit) ; cat("\n") }
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[BM_ID]] <- BM_fit
  }

  # Fit OU model
  if (any(evolutionary_models == "OU"))
  {
    OU_ID <- which(evolutionary_models == "OU")
    # OU_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "OU", ...)
    OU_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "OU"), args_for_fitContinuous))

    if (verbose) { print(OU_fit) ; cat("\n") }
    # sigsq = the drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space
      # It is not the evolutionary rate since such rate quantifying the ratio of movement in trait space per unit of time results of the combine action of the drift parameter and the strength of selection (alpha))
    # alpha = parameter describing the strength of selection pushing trait evolution toward the optimum
    # z0 = ancestral root state

    # Store output
    models_fits[[OU_ID]] <- OU_fit
  }

  # Fit EB model
  if (any(evolutionary_models == "EB"))
  {
    EB_ID <- which(evolutionary_models == "EB")
    # EB_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "EB", ...)
    EB_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "EB"), args_for_fitContinuous))

    if (verbose) { print(EB_fit) ; cat("\n") }
    # sigsq = the initial evolutionary rate at t = 0
    # a = time-dependency parameter describing the shape of the exponential decrease of sigsq though time
      # If a < 0 => the model is an EB with rates decaying in time
      # If a > 0 => the model is a late burst with rates increasing exponentially in time
    # z0 = ancestral root state

    # Store output
    models_fits[[EB_ID]] <- EB_fit
  }

  # Fit rate_trend model
  if (any(evolutionary_models == "rate_trend"))
  {
    rate_trend_ID <- which(evolutionary_models == "rate_trend")
    # rate_trend_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "rate_trend", ...)
    rate_trend_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "rate_trend"), args_for_fitContinuous))

    if (verbose) { print(rate_trend_fit) ; cat("\n") }
    # slope = coefficient of the linear trend in rates through time
      # If slope < 0 => rates decrease linearly in time
      # If slope > 0 => rates increase linearly in time
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[rate_trend_ID]] <- rate_trend_fit
  }

  # Fit Pagel's lambda model
  if (any(evolutionary_models == "lambda"))
  {
    lambda_ID <- which(evolutionary_models == "lambda")
    # lambda_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "lambda", ...)
    lambda_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "lambda"), args_for_fitContinuous))

    if (verbose) { print(lambda_fit) ; cat("\n") }
    # lambda = Internal branch length multiplier = modulates strength of the phylogenetic signal
      # Close to 0 => star-like phylogeny = no phylogenetic structure/signal
      # Close to 1 => BM = "perfect" phylogenetic signal
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[lambda_ID]] <- lambda_fit
  }

  # Fit Pagel's kappa model
  if (any(evolutionary_models == "kappa"))
  {
    kappa_ID <- which(evolutionary_models == "kappa")
    # kappa_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "kappa", ...)
    kappa_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "kappa"), args_for_fitContinuous))

    if (verbose) { print(kappa_fit) ; cat("\n") }
    # kappa = Branch length exponent = modulates effect of cladogenesis on trait evolution (punctuational/speciational) model of trait evolution by raising all branch length to power kappa
      # Close to 0 => all branches have equal length, thus only the number of speciation events affect trait evolution.
      # Close to 1 => No effect of cladogenesis = BM.
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[kappa_ID]] <- kappa_fit
  }

  # Fit Pagel's delta model
  if (any(evolutionary_models == "delta"))
  {
    delta_ID <- which(evolutionary_models == "delta")
    # delta_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "delta", ...)
    delta_fit <- do.call(what = geiger::fitContinuous, args =  c(list(phy = phylo, dat = tip_data, model = "delta"), args_for_fitContinuous))

    if (verbose) { print(delta_fit) ; cat("\n") }
    # delta = Node height exponent = time-dependent model that adjusts the relative contributions of early versus late evolution in the tree.
      # delta < 1 => relative node heights are reduced, particularly for shallower (closer to the tips) branches. Equivalent to an Early Burst model when most trait evolution happens in early times.
      # delta = 1 => No time-dependency = BM.
      # delta > 1 => relative node heights are increased, particularly for shallower (closer to the tips) branches. Equivalent to a Late Burst model when most trait evolution happens in recent times.
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[delta_ID]] <- delta_fit
  }

  # Name model fits according to the associated models
  names(models_fits) <- evolutionary_models

  ### Step 2 - Compare model fits and select best model

  # Use it to compare model fits with AICc and Akaike's weights
  models_comparison <- select_best_trait_model_from_geiger(list_model_fits = models_fits)
  best_model_name <- models_comparison$best_model_name
  best_model_fit <- models_comparison$best_model_fit

  # Display result of model comparison
  if (verbose)
  {
    cat(paste0(Sys.time(), " - Compare model fits.\n\n"))
    print(models_comparison$models_comparison_df)
  }

  ### Step 3 - Get ACE (useful in all cases to build the contMap)

  if (verbose) { cat(paste0(Sys.time(), " - Infer Ancestral Character Estimates from the best fitting model: ",best_model_name,".\n\n")) }

  ## Rescale the tree such as the relative branch length account for trait rates in the best fitting model
  # ?geiger::rescale.phylo

  # Adjust model name
  if (best_model_name == "rate_trend") { best_model_name <- "trend" }
  # Create function to rescale phylo
  rescaled_phylo_fn <- phytools::rescale(x = phylo, model = best_model_name)
  # Extract parameters from best model fit
  fun_args <- formals(fun = rescaled_phylo_fn)
  best_model_args <- best_model_fit$opt[names(fun_args)]
  # Rescale phylogeny using the estimated parameters from the best model
  rescaled_phylo <- do.call(what = rescaled_phylo_fn, args = best_model_args)

  ## Run ACE inference with BM on transformed tree
  # ?phytools::fastAnc

  ACE_output <- phytools::fastAnc(tree = rescaled_phylo,
                                  x = tip_data,
                                  vars = FALSE, # Compute the variance of ancestral state estimates
                                  CI = FALSE) # Compute the 95% CI of ancestral state estimates
  # ACE_output

  ### Step 5 - Create contMap using ACE computed with the best fitting model

  if (verbose) { cat(paste0(Sys.time(), " - Create contMap by interpolating values along branches.\n\n")) }

  contMap <- phytools::contMap(tree = phylo,
                               method = "user",
                               x = tip_data,
                               anc.states = ACE_output,
                               res = res, # Number of time steps
                               plot = FALSE)

  ## Update color scale if requested
  if (!is.null(color_scale))
  {
    # Update color palette in contMap
    contMap <- phytools::setMap(x = contMap, colors = color_scale)
    # print(contMap$cols)
  }

  ## Plot only if requested
  if (plot_map)
  {
    plot(contMap)
  }

  ## Export contMap in PDF if requested
  if (!is.null(PDF_file_path))
  {
    nb_tips <- length(phylo$tip.label)
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10

    grDevices::pdf(file = file.path(PDF_file_path),
                   width = width, height = height)

    plot_contMap(contMap)

    grDevices::dev.off()

  }

  ## Build output
  output <- list(contMap = contMap,
                 trait_data_type = "continuous")
  # Include ACE if requested
  if(return_ace) { output$ace <- ACE_output }
  # Include output of best model if requested
  if(return_best_model_fit) { output$best_model_fit <- best_model_fit }
  # Include df for model comparison if requested
  if(return_model_selection_df) { output$model_selection_df <- models_comparison$models_comparison_df }

  ## Return output
  return(invisible(output))

}


### Sub-function to handle categorical data ####

prepare_trait_data_for_categorical_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "ARD", # Default = "ARD" for categorical
    Q_matrix = NULL, # Custom Q-matrix for categorical data
    # ..., # Additional arguments for geiger::fitDiscrete() and phytools::make.simmap()
    args_for_fitDiscrete = NULL, # Additional arguments for geiger::fitDiscrete()
    args_for_make.simmap = NULL, # Additional arguments for phytools::make.simmap()
    res = 100,
    nb_simulations = 1000, # Only for categorical and biogeographic data
    colors_per_levels = NULL,
    plot_map = TRUE,
    plot_overlay = TRUE, # Only for categorical and biogeographic data
    add_ACE_pies = TRUE, # Only for categorical and biogeographic data
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_simmaps = TRUE, # Only for categorical and biogeographic data
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE)
{
  ### Check input validity
  {
    ## evolutionary_models
    if (!is.null(evolutionary_models))
    {
      if (!all(evolutionary_models %in% c("ER", "SYM", "ARD", "meristic", "matrix")))
      {
        stop(paste0("For 'trait_data_type = categorical', 'evolutionary_models' must be selected among: 'ER', 'SYM', 'ARD', 'meristic', 'matrix'.\n",
                    "See details in ?geiger::fitDiscrete()."))
      }
      ## Q_matrix
      if ("matrix" %in% evolutionary_models)
      {
        if (is.null(Q_matrix))
        {
         stop(paste0("For 'evolutionary_models = matrix', you must provide a 'Q_matrix' defining rate transition classes.\n",
              "See details in ?geiger::fitDiscrete()."))
        }
      }
    }

    ## res
    # res must be a positive integer.
    if ((res != abs(res)) | (res != round(res)))
    {
      stop(paste0("'res' must be a positive integer defining the number of time steps used to interpolate trait value in the densityMaps."))
    }

    ## nb_simulations
    # Check that nb_simulations is a positive integer.
    if ((nb_simulations != abs(nb_simulations)) | (nb_simulations != round(nb_simulations)))
    {
      stop(paste0("'nb_simulations' must be a positive integer defining the number of simulations generated for stochastic mapping."))
    }
    #  Send warnings if higher than 10000 and lower than 100
    if ((nb_simulations >= 10000))
    {
      cat(paste0("WARNING: 'nb_simulations' is set to ",nb_simulations,". High number of simulations may be time-consuming and only improve marginally the robustness of the tests.\n"))
    }
    if ((nb_simulations < 100))
    {
      cat(paste0("WARNING: 'nb_simulations' is set to ",nb_simulations,". Low number of simulations may provide biased estimates of states/ranges and affect test outputs.\n"))
    }

    # Extract list of states from tip_data
    states_list <- levels(as.factor(tip_data))

    ## colors_per_levels
    if (!is.null(colors_per_levels))
    {
      # Check that the color scale match the states
      if (!all(states_list %in% names(colors_per_levels)))
      {
        missing_states <- states_list[!(states_list %in% names(colors_per_levels))]
        stop(paste0("Not all states are found in 'colors_per_levels'.\n",
                    "Missing states: ", paste(missing_states, collapse = ", "), "."))
      }
      # Check whether all colors are valid
      if (!all(is_color(colors_per_levels)))
      {
        invalid_colors <- colors_per_levels[!is_color(colors_per_levels)]
        stop(paste0("Some color names in 'colors_per_levels' are not valid.\n",
                    "Invalid color(s): ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }

  # Set default model (ARD) if absent
  if (is.null(evolutionary_models))
  {
    evolutionary_models <- "ARD"
    message("WARNING: No 'evolutionary_model' is provided. Use 'ARD' as the default.\n")
  }

  # Get number of tips
  nb_tips <- length(phylo$tip.label)
  # Get number of nodes
  nb_nodes <- nb_tips + phylo$Nnode

  ### Step 1 - Run all models and store their results in a list
  # ?geiger::fitDiscrete

  nb_models <- length(evolutionary_models)
  if (verbose) { cat(paste0("\n", Sys.time(), " - Fit ",nb_models," evolutionary model(s): ", paste(evolutionary_models, collapse = ", "), ".\n\n")) }

  # Initiate list to store models outputs
  models_fits <- list()

  # Fit ARD model
  if (any(evolutionary_models == "ARD"))
  {
    ARD_ID <- which(evolutionary_models == "ARD")
    # ARD_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "ARD", ...)
    ARD_fit <- do.call(geiger::fitDiscrete, c(list(phy = phylo, dat = tip_data, model = "ARD"), args_for_fitDiscrete))

    if (verbose) { cat("------ ARD model ------ \n\n") ; print(ARD_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(ARD_fit$opt)
    # q12 = transition rate from state 1 to state 2
    # q21 = transition rate from state 2 to state 1
    # q13 = transition rate from state 1 to state 3
    # All transition rates can differ

    # Store output
    models_fits[[ARD_ID]] <- ARD_fit
  }

  # Fit ER model
  if (any(evolutionary_models == "ER"))
  {
    ER_ID <- which(evolutionary_models == "ER")
    # ER_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "ER", ...)
    ER_fit <- do.call(geiger::fitDiscrete, c(list(phy = phylo, dat = tip_data, model = "ER"), args_for_fitDiscrete))

    if (verbose) { cat("------ ER model ------ \n\n") ; print(ER_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(ER_fit$opt)
    # q12 = q21 = q13 = q31 = q23 = q32
    # Transition rates between states are all equals

    # Store output
    models_fits[[ER_ID]] <- ER_fit
  }

  # Fit SYM model
  if (any(evolutionary_models == "SYM"))
  {
    SYM_ID <- which(evolutionary_models == "SYM")
    # SYM_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "SYM", ...)
    SYM_fit <- do.call(geiger::fitDiscrete, c(list(phy = phylo, dat = tip_data, model = "SYM"), args_for_fitDiscrete))

    if (verbose) { cat("------ SYM model ------ \n\n") ; print(SYM_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(SYM_fit$opt)
    # q12 = q21
    # Transition rates from state 1 to state 2 are equal in both direction
    # But transition between different states can differ: q12 ≠ q13

    # Store output
    models_fits[[SYM_ID]] <- SYM_fit
  }

  # Fit 'meristic' model = step-wise transitions => 1 <-> 2 <-> 3
  if (any(evolutionary_models == "meristic"))
  {
    meristic_ID <- which(evolutionary_models == "meristic")
    # meristic_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "meristic", ...)
    meristic_fit <- do.call(geiger::fitDiscrete, c(list(phy = phylo, dat = tip_data, model = "meristic"), args_for_fitDiscrete))
    # Can define if the Q-matrix is symmetrical or not with 'symmetric = TRUE/FALSE'. Default is TRUE.

    if (verbose) { cat("------ Meristic model ------ \n\n") ; print(meristic_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(meristic_fit$opt)
    # q13 = q31 = 0
    # No transition possible between non-sequential states

    # Store output
    models_fits[[meristic_ID]] <- meristic_fit
  }

  # Fit 'matrix' model = User-defined Q-matrix defining state transition classes
  if (any(evolutionary_models == "matrix"))
  {
    matrix_ID <- which(evolutionary_models == "matrix")
    # matrix_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = Q_matrix, ...)
    matrix_fit <- do.call(geiger::fitDiscrete, c(list(phy = phylo, dat = tip_data, model = Q_matrix), args_for_fitDiscrete))
    # Can define if the Q-matrix is symmetrical or not with symmetric = TRUE/FALSE. Default is TRUE.

    if (verbose) { cat("------ User-defined matrix model ------ \n\n") ; print(matrix_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(matrix_fit$opt)

    # Store output
    models_fits[[matrix_ID]] <- matrix_fit
  }

  # Name model fits according to the associated models
  names(models_fits) <- evolutionary_models

  ### Step 2 - Compare model fits and select best model

  # Compare model fits with AICc and Akaike's weights
  models_comparison <- select_best_trait_model_from_geiger(list_model_fits = models_fits)

  # Extract best model
  best_model_name <- models_comparison$best_model_name
  best_model_fit <- models_comparison$best_model_fit

  # Display result of model comparison
  if (verbose)
  {
    cat(paste0(Sys.time(), " - Compare model fits.\n\n"))
    print(models_comparison$models_comparison_df)
  }

  ### Step 4 - Run stochastic mapping = obtain simmaps

  # Display steps
  if (verbose)
  {
    cat(paste0(Sys.time(), " - Run simulations for stochastic mapping.\n\n"))
  }

  # Extract Q-matrix from best model
  best_Q_matrix <- phytools::as.Qmatrix(best_model_fit)

  # all_simmaps <- phytools::make.simmap(tree = phylo, x = tip_data, model = best_model_name,
  #                                  nsim = nb_simulations, # Run at least 100 simulations because you want to observe the mean trend
  #                                  Q = best_Q_matrix,
  #                                  ...) # Can provide additional arguments for root states

  all_simmaps <- do.call(what = phytools::make.simmap,
                         args = c(list(tree = phylo, x = tip_data, model = best_model_name,
                                       nsim = nb_simulations, # Run at least 100 simulations because you want to observe the mean trend
                                       Q = best_Q_matrix),
                                  args_for_make.simmap)) # Can provide additional arguments for root states

  ### Step 3 - Get ACE at nodes based on posterior sampling of the simmaps (if return_ace)

  # Display steps
  if (verbose)
  {
    cat(paste0(Sys.time(), " - Extract ACE as posterior sampling from stochastic mapping.\n"))
  }

  # Use phytools summary function
  simmaps_summary_obj <- phytools::describe.simmap(tree = all_simmaps, plot = FALSE)
  ace_matrix <- simmaps_summary_obj$ace

  ### Step 5 - Create densityMaps from simmaps

  if (verbose) { cat(paste0(Sys.time(), " - Create densityMaps by summarizing simulations of evolutionary history (simmaps).\n\n")) }

  # Works only for binary traits.
  # Need to binarize every state as 0/1 with 0 = Not the focal state; 1 = Focal state.

  ## If not provided, define a color per state
  if (is.null(colors_per_levels))
  {
    colors_per_levels <- grDevices::rainbow(n = length(states_list), s = 0.8) # With HSV color space
    # colors_per_levels <- colorspace::rainbow_hcl(n = length(states_list), c = 90, l = 65) # With HLC color space
    names(colors_per_levels) <- states_list
  } else {
    # If provided, ensure it is properly ordered
    colors_per_levels <- colors_per_levels[match(x = states_list, table = names(colors_per_levels))]
  }

  # Initiate list of densityMaps
  densityMaps_all_states <- list()

  ## Loop per state
  for (i in seq_along(states_list))
  {
    # i <- 1

    # Extract state
    state_i <- states_list[i]
    # Define other states by contrast
    other_states <- states_list[states_list != state_i]

    ## Binarize states in simmaps
    simmaps_binary_i <- all_simmaps
    simmaps_binary_i <- lapply(X = simmaps_binary_i, FUN = phytools::mergeMappedStates, old.states = other_states, new.state = "0")
    simmaps_binary_i <- lapply(X = simmaps_binary_i, FUN = phytools::mergeMappedStates, old.states = state_i, new.state = "1")

    class(simmaps_binary_i) <- c("list", "multiSimmap", "multiPhylo")

    # Check that remaining states are all binary
    # unique(unlist(lapply(X = simmaps_binary_i, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
    # simmaps_binary_i[[1]]$maps

    ## Estimate the posterior probabilities of states along all branches (from the set of simulated maps)

    densityMap_state_i <- densityMap_custom(trees = simmaps_binary_i,
                                            tol = 1e-5, verbose = verbose,
                                            col_scale = NULL,
                                            plot = FALSE)

    ## Update color gradient

    # Set color gradient from grey to focal color
    focal_color <- colors_per_levels[i]
    col_fn <- grDevices::colorRampPalette(colors = c("grey90", focal_color))
    col_scale <- col_fn(n = 1001)

    # Update color gradient
    densityMap_state_i <- phytools::setMap(densityMap_state_i, c("grey90", focal_color))

    ## Update state names
    densityMap_state_i$states <- c(paste0("Not ", state_i), state_i)

    # plot(densityMap_state_i)

    ## Store in final object with all density maps
    densityMaps_all_states[[i]] <- densityMap_state_i

    ## Print progress for each state
    if (verbose)
    {
      cat(paste0(Sys.time(), " - Posterior probabilities computed for State = ",state_i," - n\u00B0", i, "/", length(states_list),"\n"))
    }
  }
  names(densityMaps_all_states) <- paste0("Density_map_", states_list)

  ## Plot densityMaps (if plot_map)
  if (plot_map)
  {
    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    if (!plot_overlay)
    {
      ## Print progress for each state
      if (verbose)
      {
        cat(paste0("\n", Sys.time(), " - Plot one densityMap for each state.\n"))
      }

      ## Plot one densityMap per state
      for (i in seq_along(densityMaps_all_states))
      {
        plot(densityMaps_all_states[[i]])
      }
    } else {

      ## Print progress for each state
      if (verbose)
      {
        cat(paste0("\n", Sys.time(), " - Plot a unique densityMap with for all states overlaid.\n"))
      }

      ## Plot the overlay of densityMaps with alpha
      plot_densityMaps_overlay(densityMaps = densityMaps_all_states,
                               add_ACE_pies = TRUE,
                               ace = ace_matrix[1:(nrow(ace_matrix)-nb_tips), ])
    }

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Export densityMap(s) in PDF if requested
  if (!is.null(PDF_file_path))
  {
    ## Print progress for each state
    if (verbose)
    {
      cat(paste0("\n", Sys.time(), " - Export PDF.\n"))
    }

    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    # Adjust width/height according to the nb of tips
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10

    if (!plot_overlay)
    {
      ## Plot one densityMap per state
      grDevices::pdf(file = file.path(PDF_file_path),
                     width = width, height = height)

      for (i in seq_along(densityMaps_all_states))
      {
        plot(densityMaps_all_states[[i]])
      }

      grDevices::dev.off()

    } else {

      ## Plot the overlay of densityMaps with alpha
      grDevices::pdf(file = file.path(PDF_file_path),
                     width = width, height = height)

      plot_densityMaps_overlay(densityMaps = densityMaps_all_states,
                               add_ACE_pies = TRUE,
                               ace = ace_matrix[1:(nrow(ace_matrix)-nb_tips), ])

      grDevices::dev.off()
    }

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Build output
  output <- list(densityMaps = densityMaps_all_states,
                 trait_data_type = "categorical")
  # Include simmaps if requested
  if(return_simmaps) { output$simmaps <- all_simmaps } # Filter to include only internal nodes (to be consistent with continuous trait)
  # Include ACE if requested
  if(return_ace) { output$ace <- ace_matrix[1:(nrow(ace_matrix)-nb_tips), ] } # Filter to include only internal nodes (to be consistent with continuous trait)
  # Include output of best model if requested
  if(return_best_model_fit) { output$best_model_fit <- best_model_fit }
  # Include df for model comparison if requested
  if(return_model_selection_df) { output$model_selection_df <- models_comparison$models_comparison_df }

  ## Return output
  return(invisible(output))
}


### Sub-function to handle biogeographic data ####

prepare_trait_data_for_biogeographic_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "DEC", # Default = "DEC" for biogeographic
    BioGeoBEARS_directory_path = "./BioGeoBEARS_directory/",
    keep_BioGeoBEARS_files = TRUE,
    prefix_for_files = NULL,
    nb_cores = 1,
    max_range_size = 2,
    split_multi_area_ranges = FALSE,
    # ..., # Additional arguments for BioGeoBEARS::define_BioGeoBEARS_run() and BioGeoBEARS::runBSM()
    args_to_define_BioGeoBEARS_run = args_to_define_BioGeoBEARS_run, # Additional arguments for BioGeoBEARS::define_BioGeoBEARS_run()
    args_to_define_runBSM, # Additional arguments for BioGeoBEARS::runBSM()
    res = 100,
    colors_per_levels = NULL,
    nb_simulations = 1000, # Only for categorical and biogeographic data
    plot_map = TRUE,
    plot_overlay = TRUE, # Only for categorical and biogeographic data
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_BSM = FALSE, # Only for biogeographic data
    return_simmaps = FALSE, # Only for categorical and biogeographic data
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE)
{
  ## Control for BioGeoBEARS install
  if (!requireNamespace("BioGeoBEARS", quietly = TRUE))
  {
    stop("Package 'BioGeoBEARS' is needed for this function to work.
       Please install it manually from: https://github.com/nmatzke/BioGeoBEARS")
  }

  ## Convert tip_data into vector of character strings if needed
  if (any(c(is.matrix(tip_data), is.data.frame(tip_data))))
  {
    # Check initial format
    PA_test <- apply(X = tip_data, MARGIN = 2, FUN = function (x) { all(x %in% c(0, 1)) })
    if (!all(PA_test))
    {
      stop("You provided 'tip_data' for 'biogeographic' data as a matrix/data.frame.\n",
           "Please ensure presences/absences are recorded as numerical 0/1 values.\n")
    }
    unique_letters_test <- all(nchar(colnames(tip_data)) == 1)
    if (!all(unique_letters_test))
    {
      stop("You provided 'tip_data' for 'biogeographic' data as a matrix/data.frame.\n",
           "Please ensure colnames(tip_data) are unique areas symbolized by a unique UPPERCASE letter each.\n")
    }
    uppercase_test <- colnames(tip_data) == toupper(colnames(tip_data))
    if (!all(uppercase_test))
    {
      stop("You provided 'tip_data' for 'biogeographic' data as a matrix/data.frame.\n",
           "Please ensure colnames(tip_data) are unique areas symbolized by a unique UPPERCASE letter each.\n")
    }
    # Convert to vector of character strings
    tip_data_df <- as.data.frame(tip_data)
    tip_data <- setNames(object = rep(NA, times = nrow(tip_data_df)), nm = row.names(tip_data_df))
    for (i in 1:nrow(tip_data_df))
    {
      # i <- 1

      # Extract unique areas
      unique_areas_i <- names(tip_data_df)[tip_data_df[i, ] == 1]
      unique_areas_i <- unique_areas_i[order(unique_areas_i)]

      # Collapse into range
      range_i <- paste(unique_areas_i, collapse = "")
      # Store range
      tip_data[i] <- range_i
    }
  }

  ## Extract and order ranges
  all_ranges <- unique(tip_data)
  all_ranges <- all_ranges[order(all_ranges)]
  unique_areas <- all_ranges[nchar(all_ranges) == 1]
  multi_area_ranges <- setdiff(all_ranges, unique_areas)
  all_ranges <- c(unique_areas, multi_area_ranges)

  ### Check input validity
  {
    ## tip_data
    # Check that tip_data follow the coding scheme of BioGeoBEARS
    if (!identical(toupper(tip_data), tip_data))
    {
      stop(paste0("For 'trait_data_type = biogeographic', 'tip_data' must be only in CAPITALS following BioGeoBEARS scheme.\n",
                  "See details in documentation from R package `BioGeoBEARS`."))
    }
    # Check that all multi-area ranges are covered by unique areas
    unique_areas_in_multi_area_ranges <- unique(unlist(strsplit(x = multi_area_ranges, split = "")))
    if (!all(unique_areas_in_multi_area_ranges %in% unique_areas))
    {
      missing_unique_areas <- unique_areas_in_multi_area_ranges[!(unique_areas_in_multi_area_ranges %in% unique_areas)]
      detect_problematic_multi_area_ranges <- unlist(lapply(X = strsplit(x = multi_area_ranges, split = ""), FUN = function (x) { any(missing_unique_areas %in% x) }))
      problematic_multi_area_ranges <- multi_area_ranges[detect_problematic_multi_area_ranges]
      stop(paste0("For 'trait_data_type = biogeographic', 'tip_data' must following BioGeoBEARS scheme.\n",
                  "Currently, range(s) ",paste(problematic_multi_area_ranges, collapse = ", ")," is/are not covered by the unique areas found in tip_data.\n",
                  "Unique areas: ",paste(unique_areas, collapse = ", "),"."))
    }

    ## evolutionary_models
    if (!is.null(evolutionary_models))
    {
      if (!all(evolutionary_models %in% c("BAYAREALIKE", "DIVALIKE", "DEC", "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J")))
      {
        stop(paste0("For 'trait_data_type = biogeographic', 'evolutionary_models' must be selected among: 'BAYAREALIKE', 'DIVALIKE', 'DEC', 'BAYAREALIKE+J', 'DIVALIKE+J', 'DEC+J'.\n",
                    "See details in documentation from R package `BioGeoBEARS`."))
      }
    }

    ## BioGeoBEARS_directory_path
    # BioGeoBEARS_directory_path must be a directory, so it must end with '/'
    if (!stringr::str_detect(BioGeoBEARS_directory_path, pattern = "/$"))
    {
      stop(paste0("'BioGeoBEARS_directory_path' must end with '/'"))
    }
    # Create the directory if it does not exist yet
    if (!dir.exists(BioGeoBEARS_directory_path))
    {
      dir.create(path = BioGeoBEARS_directory_path, recursive = T)
      message(paste0("WARNING: '",BioGeoBEARS_directory_path,"' does not exist. The directory was created to be able to run BioGeoBEARS."))
    }

    ## max_range_size
    if ((max_range_size != abs(max_range_size)) | (max_range_size != round(max_range_size)))
    {
      stop(paste0("For 'trait_data_type = biogeographic', 'max_range_size' represents the maximum number of unique areas encompassed by multi-area ranges.\n",
                  "This must be a positive integer."))
    }
    if (max_range_size > length(unique_areas))
    {
      stop(paste0("For 'trait_data_type = biogeographic', 'max_range_size' represents the maximum number of unique areas encompassed by multi-area ranges.\n",
                  "It cannot be larger than the number of unique areas found in 'tip_data'.\n",
                  "Currently, 'max_range_size' = ", max_range_size,". Number of unique areas = ",length(unique_areas),"."))
    }

    ## split_multi_area_ranges
    # Check if ranges need to be split
    if (split_multi_area_ranges & max_range_size == 1)
    {
      stop(paste0("You requested to split multi-area ranges, but no r 'trait_data_type = biogeographic', 'evolutionary_models' must be selected among: 'BAYAREALIKE', 'DIVALIKE', 'DEC', 'BAYAREALIKE+J', 'DIVALIKE+J', 'DEC+J'.\n",
                  "See details in documentation from R package `BioGeoBEARS`."))
    }

    ## res
    # res must be a positive integer.
    if ((res != abs(res)) | (res != round(res)))
    {
      stop(paste0("'res' must be a positive integer defining the number of time steps used to interpolate trait value in the densityMaps."))
    }

    ## nb_simulations
    # Check that nb_simulations is a positive integer.
    if ((nb_simulations != abs(nb_simulations)) | (nb_simulations != round(nb_simulations)))
    {
      stop(paste0("'nb_simulations' must be a positive integer defining the number of simulations generated for stochastic mapping."))
    }
    #  Send warnings if higher than 10000 and lower than 100
    if ((nb_simulations >= 10000))
    {
      warning(paste0("'nb_simulations' is set to ",nb_simulations,". High number of simulations may be time-conusming and only improve marginally the robustness of the tests."))
    }
    if ((nb_simulations < 100))
    {
      warning(paste0("'nb_simulations' is set to ",nb_simulations,". Low number of simulations may provide bias estimates of states/ranges and affect test outputs."))
    }

    ## colors_per_levels
    if (!is.null(colors_per_levels))
    {
      # Check that the color scale match at least the unique areas
      if (!all(unique_areas %in% names(colors_per_levels)))
      {
        missing_unique_areas <- unique_areas[!(unique_areas %in% names(colors_per_levels))]
        stop(paste0("Not all unique areas are found in 'colors_per_levels'.\n",
                    "You must provide at least a color per unique areas. (Multi-area ranges can be interpolated)\n",
                    "Missing unique area(s): ", paste(missing_unique_areas, collapse = ", "), ".\n\n",
                    "If set to 'NULL', all colors will be automatically generated with 'BioGeoBEARS::get_colors_for_states_list_0based()'."))
      }
      # Check whether all colors are valid
      if (!all(is_color(colors_per_levels)))
      {
        invalid_colors <- colors_per_levels[!is_color(colors_per_levels)]
        stop(paste0("Some color names in 'colors_per_levels' are not valid.\n",
                    "Invalid color(s): ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }

  # Set default model (DEC) if absent
  if (is.null(evolutionary_models))
  {
    message("WARNING: No 'evolutionary_model' is provided. Use 'DEC' as the default.\n")
    evolutionary_models <- "DEC"
  }

  # Get number of tips
  nb_tips <- length(phylo$tip.label)
  # Get number of nodes
  nb_nodes <- nb_tips + phylo$Nnode

  ### Step 1 - Run all models and store their results in a list

  ### Prepare data for BioGeoBEARS

  ## Prepare phylogeny for BioGeoBEARS

  # Build path to phylogeny
  if (is.null(prefix_for_files))
  {
    path_to_phylo <- BioGeoBEARS::np(paste0(BioGeoBEARS_directory_path, "phylo.tree"))
  } else {
    path_to_phylo <- BioGeoBEARS::np(paste0(BioGeoBEARS_directory_path, prefix_for_files,"_phylo.tree"))
  }
  # Export phylogeny in Newick format
  write.tree(phylo, file = path_to_phylo)

  ## Prepare tip ranges for BioGeoBEARS

  # Convert to df
  ranges_df <- as.data.frame(tip_data)
  # Get list of all unique areas
  unique_areas_in_ranges_list <- strsplit(x = tip_data, split = "")
  # Loop per unique areas
  for (i in seq_along(unique_areas))
  {
    # i <- 1

    # Extract unique area
    unique_area_i <- unique_areas[i]
    # Detect presence in ranges
    binary_match_i <- unlist(lapply(X = unique_areas_in_ranges_list, FUN = function (x) { unique_area_i %in% x } ))

    # Add to ranges_df
    ranges_df <- cbind(ranges_df, binary_match_i)
  }
  # Extract binary df of presence/absence
  binary_df <- ranges_df[, -1]
  # Convert character strings into numerical factors
  binary_df_num <- as.data.frame(apply(X = binary_df, MARGIN = 2, FUN = as.numeric))
  row.names(binary_df_num) <- names(tip_data)
  names(binary_df_num) <- unique_areas

  # Produce tipranges object from numeric df
  Taxa_bioregions_tipranges_obj <- BioGeoBEARS::define_tipranges_object(tmpdf = binary_df_num)

  # Set path to tip ranges object
  if (is.null(prefix_for_files))
  {
    path_to_tip_ranges <- BioGeoBEARS::np(paste0(BioGeoBEARS_directory_path,"tip_ranges.data"))
  } else {
    path_to_tip_ranges <- BioGeoBEARS::np(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_tip_ranges.data"))
  }

  # Export tip ranges in Lagrange/PHYLIP format
  BioGeoBEARS::save_tipranges_to_LagrangePHYLIP(
     tipranges_object = Taxa_bioregions_tipranges_obj,
     lgdata_fn = path_to_tip_ranges,
     areanames = colnames(Taxa_bioregions_tipranges_obj@df))

  ### Prepare models

  if (verbose) { cat(paste0("# ----------- Step 0: Set BioGeoBEARS models ----------- #\n\n")) }

  ## Set the default DEC run

  ## Properties of DEC model

  # Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
  # No jump dispersal (j = 0)
  # Here, relative weight of y,v,s are fixed to 1.
  # Model as only 2 free parameters = d (range extension) and e (range contraction)

  # DEC_run <- BioGeoBEARS::define_BioGeoBEARS_run(
  #    num_cores_to_use = nb_cores,
  #    max_range_size = max_range_size, # To set the maximum number of bioregion encompassed by a lineage range at any time
  #    trfn = path_to_phylo, # To provide path to the input tree file
  #    geogfn = path_to_tip_ranges, # To provide path to the LagrangePHYLIP file with binary ranges
  #    return_condlikes_table = TRUE, # To ask to obtain all marginal likelihoods computed by the model and used to display ancestral states
  #    print_optim = verbose,
  #    ...) # Additional arguments

  DEC_run <- do.call(what = BioGeoBEARS::define_BioGeoBEARS_run,
                     args = c(list(num_cores_to_use = nb_cores,
                                   max_range_size = max_range_size, # To set the maximum number of bioregion encompassed by a lineage range at any time
                                   trfn = path_to_phylo, # To provide path to the input tree file
                                   geogfn = path_to_tip_ranges, # To provide path to the LagrangePHYLIP file with binary ranges
                                   return_condlikes_table = TRUE, # To ask to obtain all marginal likelihoods computed by the model and used to display ancestral states
                                   print_optim = verbose),
                              args_to_define_BioGeoBEARS_run)) # Additional arguments

  # Check that starting parameter values are inside the min/max
  DEC_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_run)
  # Check validity of set-up before run
  invisible(BioGeoBEARS::check_BioGeoBEARS_run(DEC_run))

  ## Inspect the parameter table

  # DEC_run$BioGeoBEARS_model_object@params_table

  # d = dispersal rate = anagenetic range extension. Ex: A -> AB
  # e = extinction rate = anagenetic range contraction. Ex: AB -> A
  # a = range-switching rate = anagenetic jump dispersal. Ex: A -> B

  # x = used in DEC+X to account for effect of geographic distances between bioregions
  # n = used to account for effect of environmental distances between bioregions
  # w = used in DEC+W to account for effect of manual dispersal multipliers
  # u = used to account for bioregion-specific risk of extinction based on areas of bioregions

  # j = Jump-dispersal relative weight = cladogenetic founder-event. Ex: A -> (A),(B)
  # y = Non-transitional speciation relative weight = cladogenetic inheritence. Ex: A -> (A),(A)
  # s = Subset speciation relative weight = cladogenetic sympatric speciation. Ex: AB -> (AB),(A)
  # v = Vicariance relative weight = cladogenetic vicariance. Ex: AB -> (A),(B)

  ## Set the DIVALIKE run
  if (any(c("DIVALIKE", "DIVALIKE+J") %in% evolutionary_models))
  {
    # DIVALIKE is a likelihood interpretation of parsimony DIVA, and it is "like DIVA"
    # Similar to, but not identical to, parsimony DIVA.

    # Allows narrow sympatric speciation (y), and narrow AND WIDE vicariance (v)
    # No subset sympatric speciation (s = 0)
    # No jump dispersal (j = 0)
    # Here, relative weight of y and v are fixed to 1.
    # Model has 2 free parameters = d (range extension), e (range contraction)

    # Use the DEC model as template
    DIVALIKE_run <- DEC_run

    # Remove subset sympatric speciation
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

    # Adjust the fixed values of ysv, ys, y, and v to account for the absence of s
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j" # Instead of 3-j in DEC
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2" # Instead of ysv*2/3 in DEC
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2" # Instead of ysv*1/3 in DEC
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2" # Instead of ysv*1/3 in DEC

    # Allow widespread vicariance; Narrow vicariance and Wide-range vicariance are set equiprobable
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    DIVALIKE_run$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

    # Check that starting parameter values are inside the min/max
    DIVALIKE_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DIVALIKE_run) # Check that starting parameter values are inside the min/max
    # Check validity of set-up before run
    invisible(BioGeoBEARS::check_BioGeoBEARS_run(DIVALIKE_run))
  }

  ## Set the BAYAREALIKE run
  if (any(c("BAYAREALIKE", "BAYAREALIKE+J") %in% evolutionary_models))
  {
    # BAYAREALIKE is a likelihood interpretation of the Bayesian implementation of BAYAREA model, and it is "like BAYAREA"
    # Similar to, but not identical to, Bayesian BAYAREA.

    # Nothing is happening during cladogenesis. Only inheritance of previous range (y = 1)
    # Allows widespread sympatric speciation: AB => (AB),(AB)
    # No narrow or wide vicariance (v = 0)
    # No subset sympatric speciation (s = 0)
    # No jump dispersal (j = 0)
    # Only changes are happening along branches through range extension and contraction
    # Model has 2 free parameters = d (range extension), e (range contraction)

    # Use the DEC model as template
    BAYAREALIKE_run <- DEC_run

    # Remove subset sympatric speciation
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

    # Remove vicariance
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

    # Adjust the fixed values of ysv, ys, y  to account for the absence of s and v
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j" # Instead of 3-j in DEC
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1" # Instead of ysv*2/3 in DEC
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/1" # Instead of ysv*1/3 in DEC

    # Force the exact copying of ancestral range during cladogenesis.
    # As a consequence: allows widespread sympatric speciation
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BAYAREALIKE_run$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

    # Check that starting parameter values are inside the min/max
    BAYAREALIKE_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BAYAREALIKE_run) # Check that starting parameter values are inside the min/max
    # Check validity of set-up before run
    invisible(BioGeoBEARS::check_BioGeoBEARS_run(BAYAREALIKE_run))
  }

  ## Set the DEC+J run
  if ("DEC+J" %in% evolutionary_models)
  {
    # Allows (A) -> (A),(B) # Cladogenetic transition as founder-event speciation = jump-speciation

    # Allows narrow sympatric speciation (y), narrow vicariance (v), and subset sympatric speciation (s)
    # Allows jump dispersal (j)
    # Here, relative weight of y,v,s are fixed to 1.
    # Model as only 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

    # Use the DEC model as template
    DEC_J_run <- DEC_run

    # Update status of jump speciation parameter to be estimated
    DEC_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
    # Set initial value of J for optimization to an arbitrary low non-null value
    j_start <- 0.0001
    DEC_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
    DEC_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization

    # Check validity of set-up before run
    DEC_J_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_J_run) # Check that starting parameter values are inside the min/max
    invisible(BioGeoBEARS::check_BioGeoBEARS_run(DEC_J_run))

  }

  ## Set the DIVALIKE+J run
  if ("DIVALIKE+J" %in% evolutionary_models)
  {
    # Allows (A) -> (A),(B) # Cladogenetic transition as founder-event speciation = jump-speciation

    # Allows narrow and WIDE vicariance (v)
    # Allows jump dispersal (j)
    # Here, relative weight of y,v are fixed to 1.
    # Model as only 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

    # Use the DIVALIKE model as template
    DIVALIKE_J_run <- DIVALIKE_run

    # Update status of jump speciation parameter to be estimated
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
    # Set initial value of J for optimization to an arbitrary low non-null value
    j_start <- 0.0001
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999 # Under DIVALIKE+J, the max of "j" should be 2, because y+v+s = 2.

    # Check validity of set-up before run
    DIVALIKE_J_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DIVALIKE_J_run) # Check that starting parameter values are inside the min/max
    invisible(BioGeoBEARS::check_BioGeoBEARS_run(DIVALIKE_J_run))
  }

  ## Set the BAYAREALIKE+J run
  if ("BAYAREALIKE+J" %in% evolutionary_models)
  {
    # Allows (A) -> (A),(B) # Cladogenetic transition as founder-event speciation = jump-speciation

    # Allows widespread sympatric speciation: (AB) => (AB, AB)
    # No narrow or wide vicariance (v = 0)
    # No subset sympatric speciation (s = 0)
    # Allows jump dispersal (j)
    # Model as only 3 free parameters = d (range extension), e (range contraction), and j (jump dispersal)

    # Use the BAYAREALIKE model as template
    BAYAREALIKE_J_run <- BAYAREALIKE_run

    # Update status of jump speciation parameter to be estimated
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
    # Set initial value of J for optimization to an arbitrary low non-null value
    j_start <- 0.0001
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start # MLE will evolved after optimization
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999 # Under BAYAREALIKE+J, the max of "j" should be 1, because y+v+s = 1.

    # Check validity of set-up before run
    BAYAREALIKE_J_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BAYAREALIKE_J_run) # Check that starting parameter values are inside the min/max
    invisible(BioGeoBEARS::check_BioGeoBEARS_run(BAYAREALIKE_J_run))
  }

  ### Run BioGeoBEARS models

  if (verbose) { cat(paste0("# ----------- Step 1: Run BioGeoBEARS models ----------- #\n\n")) }

  # Initiate list to store models outputs
  models_fits <- list()

  ## Run DEC model
  if (any(c("DEC", "DEC+J") %in% evolutionary_models))
  {
    if (verbose) { cat(paste0(Sys.time(), "- Run DEC model\n\n")) }

    # Open log file
    if (is.null(prefix_for_files))
    {
      log_file <- file(paste0(BioGeoBEARS_directory_path,"DEC_run_log.txt"), open = "wt")
    } else {
      log_file <- file(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DEC_run_log.txt"), open = "wt")
    }

    # Run model while saving logs
    sink(file = log_file, type = "output", split = TRUE)
    cat(paste0(Sys.time(), " - Start DEC model run\n"))
    Start_time <- Sys.time()
    DEC_fit <- BioGeoBEARS::bears_optim_run(DEC_run)
    End_time <- Sys.time()
    cat(paste0(Sys.time(), " - End of DEC model run\n"))
    cat(paste0("Total run time of:\n"))
    print(End_time - Start_time)
    sink()
    closeAllConnections()

    if ("DEC" %in% evolutionary_models)
    {
      # Store output
      DEC_ID <- which(evolutionary_models == "DEC")
      models_fits[[DEC_ID]] <- DEC_fit

      # Save model fit
      if (is.null(prefix_for_files))
      {
        saveRDS(object = DEC_fit, file = paste0(BioGeoBEARS_directory_path,"DEC_fit.rds"))
      } else {
        saveRDS(object = DEC_fit, file = paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DEC_fit.rds"))
      }

    } else {
      # Remove log if DEC model not requested in outputs
      if (is.null(prefix_for_files))
      {
        file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DEC_run_log.txt")))
      } else {
        file.remove(file.path(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DEC_run_log.txt")))
      }
    }
  }

  ## Run DEC+J model
  if ("DEC+J" %in% evolutionary_models)
  {
    if (verbose) { cat(paste0(Sys.time(), "- Run DEC+J model\n\n")) }

    ## Set starting values for optimization based on MLE in the DEC model
    d_start <- DEC_fit$outputs@params_table["d","est"]
    e_start <- DEC_fit$outputs@params_table["e","est"]

    DEC_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
    DEC_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
    DEC_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
    DEC_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

    # Open log file
    if (is.null(prefix_for_files))
    {
      log_file <- file(paste0(BioGeoBEARS_directory_path,"DEC_J_run_log.txt"), open = "wt")
    } else {
      log_file <- file(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DEC_J_run_log.txt"), open = "wt")
    }

    # Run model while saving logs
    sink(file = log_file, type = "output", split = TRUE)
    cat(paste0(Sys.time(), " - Start DEC+J model run\n"))
    Start_time <- Sys.time()
    DEC_J_fit <- BioGeoBEARS::bears_optim_run(DEC_J_run)
    End_time <- Sys.time()
    cat(paste0(Sys.time(), " - End of DEC+J model run\n"))
    cat(paste0("Total run time of:\n"))
    print(End_time - Start_time)
    sink()
    closeAllConnections()

    # Save model fit
    if (is.null(prefix_for_files))
    {
      saveRDS(object = DEC_J_fit, file = paste0(BioGeoBEARS_directory_path,"DEC_J_fit.rds"))
    } else {
      saveRDS(object = DEC_J_fit, file = paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DEC_J_fit.rds"))
    }

    # Store output
    DEC_J_ID <- which(evolutionary_models == "DEC+J")
    models_fits[[DEC_J_ID]] <- DEC_J_fit
  }

  ## Run DIVALIKE model
  if (any(c("DIVALIKE", "DIVALIKE+J") %in% evolutionary_models))
  {
    if (verbose) { cat(paste0(Sys.time(), "- Run DIVALIKE model\n\n")) }

    # Open log file
    if (is.null(prefix_for_files))
    {
      log_file <- file(paste0(BioGeoBEARS_directory_path,"DIVALIKE_run_log.txt"), open = "wt")
    } else {
      log_file <- file(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DIVALIKE_run_log.txt"), open = "wt")
    }

    # Run model while saving logs
    sink(file = log_file, type = "output", split = TRUE)
    cat(paste0(Sys.time(), " - Start DIVALIKE model run\n"))
    Start_time <- Sys.time()
    DIVALIKE_fit <- BioGeoBEARS::bears_optim_run(DIVALIKE_run)
    End_time <- Sys.time()
    cat(paste0(Sys.time(), " - End of DIVALIKE model run\n"))
    cat(paste0("Total run time of:\n"))
    print(End_time - Start_time)
    sink()
    closeAllConnections()

    if ("DIVALIKE" %in% evolutionary_models)
    {
      # Store output
      DIVALIKE_ID <- which(evolutionary_models == "DIVALIKE")
      models_fits[[DIVALIKE_ID]] <- DIVALIKE_fit

      # Save model fit
      if (is.null(prefix_for_files))
      {
        saveRDS(object = DIVALIKE_fit, file = paste0(BioGeoBEARS_directory_path,"DIVALIKE_fit.rds"))
      } else {
        saveRDS(object = DIVALIKE_fit, file = paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DIVALIKE_fit.rds"))
      }

    } else {
      # Remove log if DIVALIKE model not requested in outputs
      if (is.null(prefix_for_files))
      {
        file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DIVALIKE_run_log.txt")))
      } else {
        file.remove(file.path(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DIVALIKE_run_log.txt")))
      }
    }
  }

  ## Run DIVALIKE+J model
  if ("DIVALIKE+J" %in% evolutionary_models)
  {
    if (verbose) { cat(paste0(Sys.time(), "- Run DIVALIKE+J model\n\n")) }

    ## Set starting values for optimization based on MLE in the DIVALIKE model
    d_start <- DIVALIKE_fit$outputs@params_table["d","est"]
    e_start <- DIVALIKE_fit$outputs@params_table["e","est"]

    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
    DIVALIKE_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

    # Open log file
    if (is.null(prefix_for_files))
    {
      log_file <- file(paste0(BioGeoBEARS_directory_path,"DIVALIKE_J_run_log.txt"), open = "wt")
    } else {
      log_file <- file(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DIVALIKE_J_run_log.txt"), open = "wt")
    }

    # Run model while saving logs
    sink(file = log_file, type = "output", split = TRUE)
    cat(paste0(Sys.time(), " - Start DIVALIKE+J model run\n"))
    Start_time <- Sys.time()
    DIVALIKE_J_fit <- BioGeoBEARS::bears_optim_run(DIVALIKE_J_run)
    End_time <- Sys.time()
    cat(paste0(Sys.time(), " - End of DIVALIKE+J model run\n"))
    cat(paste0("Total run time of:\n"))
    print(End_time - Start_time)
    sink()
    closeAllConnections()

    # Save model fit
    if (is.null(prefix_for_files))
    {
      saveRDS(object = DIVALIKE_J_fit, file = paste0(BioGeoBEARS_directory_path,"DIVALIKE_J_fit.rds"))
    } else {
      saveRDS(object = DIVALIKE_J_fit, file = paste0(BioGeoBEARS_directory_path,prefix_for_files,"_DIVALIKE_J_fit.rds"))
    }

    # Store output
    DIVALIKE_J_ID <- which(evolutionary_models == "DIVALIKE+J")
    models_fits[[DIVALIKE_J_ID]] <- DIVALIKE_J_fit
  }

  ## Run BAYAREALIKE model
  if (any(c("BAYAREALIKE", "BAYAREALIKE+J") %in% evolutionary_models))
  {
    if (verbose) { cat(paste0(Sys.time(), "- Run BAYAREALIKE model\n\n")) }

    # Open log file
    if (is.null(prefix_for_files))
    {
      log_file <- file(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_run_log.txt"), open = "wt")
    } else {
      log_file <- file(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_BAYAREALIKE_run_log.txt"), open = "wt")
    }

    # Run model while saving logs
    sink(file = log_file, type = "output", split = TRUE)
    cat(paste0(Sys.time(), " - Start BAYAREALIKE model run\n"))
    Start_time <- Sys.time()
    BAYAREALIKE_fit <- BioGeoBEARS::bears_optim_run(BAYAREALIKE_run)
    End_time <- Sys.time()
    cat(paste0(Sys.time(), " - End of BAYAREALIKE model run\n"))
    cat(paste0("Total run time of:\n"))
    print(End_time - Start_time)
    sink()
    closeAllConnections()

    if ("BAYAREALIKE" %in% evolutionary_models)
    {
      # Store output
      BAYAREALIKE_ID <- which(evolutionary_models == "BAYAREALIKE")
      models_fits[[BAYAREALIKE_ID]] <- BAYAREALIKE_fit

      # Save model fit
      if (is.null(prefix_for_files))
      {
        saveRDS(object = BAYAREALIKE_fit, file = paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_fit.rds"))
      } else {
        saveRDS(object = BAYAREALIKE_fit, file = paste0(BioGeoBEARS_directory_path,prefix_for_files,"_BAYAREALIKE_fit.rds"))
      }

    } else {
      # Remove log if DIVALIKE model not requested in outputs
      if (is.null(prefix_for_files))
      {
        file.remove(file.path(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_run_log.txt")))
      } else {
        file.remove(file.path(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_BAYAREALIKE_run_log.txt")))
      }
    }
  }

  ## Run BAYAREALIKE+J model
  if ("BAYAREALIKE+J" %in% evolutionary_models)
  {
    if (verbose) { cat(paste0("# ----------- Run BAYAREALIKE+J model ----------- #\n\n")) }

    ## Set starting values for optimization based on MLE in the BAYAREA model
    d_start <- BAYAREALIKE_fit$outputs@params_table["d","est"]
    e_start <- BAYAREALIKE_fit$outputs@params_table["e","est"]

    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start # MLE will evolved after optimization
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
    BAYAREALIKE_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start # MLE will evolved after optimization

    # Open log file
    if (is.null(prefix_for_files))
    {
      log_file <- file(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_J_run_log.txt"), open = "wt")
    } else {
      log_file <- file(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_BAYAREALIKE_J_run_log.txt"), open = "wt")
    }

    # Run model while saving logs
    sink(file = log_file, type = "output", split = TRUE)
    cat(paste0(Sys.time(), " - Start BAYAREALIKE+J model run\n"))
    Start_time <- Sys.time()
    BAYAREALIKE_J_fit <- BioGeoBEARS::bears_optim_run(BAYAREALIKE_J_run)
    End_time <- Sys.time()
    cat(paste0(Sys.time(), " - End of BAYAREALIKE+J model run\n"))
    cat(paste0("Total run time of:\n"))
    print(End_time - Start_time)
    sink()
    closeAllConnections()

    # Save model fit
    if (is.null(prefix_for_files))
    {
      saveRDS(object = BAYAREALIKE_J_fit, file = paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_J_fit.rds"))
    } else {
      saveRDS(object = BAYAREALIKE_J_fit, file = paste0(BioGeoBEARS_directory_path,prefix_for_files,"_BAYAREALIKE_J_fit.rds"))
    }

    # Store output
    BAYAREALIKE_J_ID <- which(evolutionary_models == "BAYAREALIKE+J")
    models_fits[[BAYAREALIKE_J_ID]] <- BAYAREALIKE_J_fit
  }

  # Name model fits according to the associated models
  names(models_fits) <- evolutionary_models

  ### Step 2 - Compare model fits and select best model

  if (verbose) { cat(paste0("# ----------- Step 2: Compare model fits ----------- #\n\n")) }

  # Compare model fits with AICc and Akaike's weights
  models_comparison <- select_best_model_from_BioGeoBEARS(list_model_fits = models_fits)

  # Extract best model
  best_model_name <- models_comparison$best_model_name
  best_model_fit <- models_comparison$best_model_fit

  # Display result of model comparison
  if (verbose)
  {
    print(models_comparison$models_comparison_df)
  }

  if (verbose) { cat(paste0("# ----------- Step 3: Extract Ancestral Ranges ----------- #\n\n")) }

  if (return_ace)
  {
    # Retrieve list of ranges in the model
    all_ranges_in_models <- generate_list_ranges(areas_list = unique_areas,
                                                 max_range_size = max_range_size,
                                                 include_null_range = TRUE)

    # Extract ACE as scaled marginal likelihood of each state per tips/nodes
    ace_matrix <- best_model_fit$ML_marginal_prob_each_state_at_branch_top_AT_node
    row.names(ace_matrix) <- c(phylo$tip.label, (nb_tips+1):nb_nodes)
    colnames(ace_matrix) <- all_ranges_in_models

    # print(ace_matrix)
  }

  if (verbose) { cat(paste0("# ----------- Step 4: Run simulations for stochastic mapping ----------- #\n\n")) }

  ## Extract inputs needed for BSM from model fit object
  BSM_inputs <- BioGeoBEARS::get_inputs_for_stochastic_mapping(res = best_model_fit)

  ## Run BSM

  # BSM_output <- BioGeoBEARS::runBSM(
  #    res = best_model_fit,  # Model fit object
  #    stochastic_mapping_inputs_list = BSM_inputs,
  #    maxnum_maps_to_try = 100, # Maximum number of stochastic maps to try to simulate before giving up
  #    nummaps_goal = nb_simulations, # Number of accepted stochastic maps to target
  #    savedir = BioGeoBEARS_directory_path,
  #    ...)

  BSM_output <- do.call(what = BioGeoBEARS::runBSM, args = c(list(res = best_model_fit,  # Model fit object
                                                                  stochastic_mapping_inputs_list = BSM_inputs,
                                                                  maxnum_maps_to_try = 100, # Maximum number of stochastic maps to try to simulate before giving up
                                                                  nummaps_goal = nb_simulations, # Number of accepted stochastic maps to target
                                                                  savedir = BioGeoBEARS_directory_path),
                                                             args_to_define_runBSM))

  ## Rename BSM outputs
  if (!is.null(prefix_for_files))
  {
    file.rename(from = file.path(paste0(BioGeoBEARS_directory_path,"RES_ana_events_tables.Rdata")),
                to = file.path(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_RES_ana_events_tables.Rdata")))
    file.rename(from = file.path(paste0(BioGeoBEARS_directory_path,"RES_ana_events_tables_PARTIAL.Rdata")),
                to = file.path(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_RES_ana_events_tables_PARTIAL.Rdata")))
    file.rename(from = file.path(paste0(BioGeoBEARS_directory_path,"RES_clado_events_tables.Rdata")),
                to = file.path(paste0(BioGeoBEARS_directory_path,prefix_for_files,"_RES_clado_events_tables.Rdata")))
  }

  ## Convert BioGeoBEARS BSM output to simmaps
  all_simmaps <- BSMs_to_phytools_simmaps(model_fit = best_model_fit, phylo = phylo, BSM_output = BSM_output)

  # Print simmap n°1
  # plot(all_simmaps[1])

  if (verbose) { cat(paste0("# ----------- Step 5: Create densityMaps per ranges ----------- #\n\n")) }

  ## Define and/or interpolate the colors per ranges

  # Extract list of effective ranges used in the model
  returned_mats <- BioGeoBEARS::get_Qmat_COOmat_from_BioGeoBEARS_run_object(
    BioGeoBEARS_run_object = best_model_fit$inputs,
    include_null_range = best_model_fit$inputs$include_null_range)
  effective_ranges_list <- returned_mats$ranges_list

  # If no colors are provided, define automatically
  if (is.null(colors_per_levels))
  {
    # Extract list of 0based ranges
    effective_0based_states <- returned_mats$states_list

    # Extract max number of areas combined in a range/state
    effective_max_range_size <- max(unlist(lapply(X = effective_0based_states, FUN = length)))

    # Is there a null range?
    null_range_check <- any(unlist(lapply(X = effective_0based_states, FUN = function (x) any(is.na(x)))))

    # Produce list of colors for each possible state
    colors_list_for_ranges <- BioGeoBEARS::get_colors_for_states_list_0based(
        areanames = unique_areas,
        states_list_0based = effective_0based_states,
        max_range_size = effective_max_range_size,
        plot_null_range = null_range_check)
    names(colors_list_for_ranges) <- effective_ranges_list
  } else {

    ## If at least unique areas are provided, interpolate colors for missing multi-area ranges

    # Identify missing multi-area ranges
    missing_ranges <- effective_ranges_list[!(effective_ranges_list %in% names(colors_per_levels))]
    missing_ranges <- setdiff(missing_ranges, "_")

    # Identify unique areas encompassed by the range
    unique_areas_in_missing_ranges <- strsplit(x = missing_ranges, split = "")

    ## Loop per missing multi-area ranges
    colors_for_missing_ranges <- c()
    for (i in seq_along(missing_ranges))
    {
      # i <- 1

      # Extract unique areas
      unique_areas_in_missing_range_i <- unique_areas_in_missing_ranges[[i]]
      # Get colors of all unique areas
      colors_in_missing_range_i <- colors_per_levels[unique_areas_in_missing_range_i]
      # Get RGB values of all unique areas
      RGB_in_missing_range_i <- col2rgb(colors_in_missing_range_i)
      # Get HSV values of all unique areas
      HSV_in_missing_range_i <- grDevices::rgb2hsv(RGB_in_missing_range_i)
      # Interpolate as the mean of all unique areas
      HSV_range_i <- apply(X = HSV_in_missing_range_i, MARGIN = 1, FUN = mean)
      # Hue value is an angle, interpolation is not a simple mean
      HSV_range_i[1] <- mean_angle(a = HSV_in_missing_range_i[1, ]*360, unit = "degrees")/360
      col_range_i <- grDevices::hsv(h = HSV_range_i[1], s = HSV_range_i[2], v = HSV_range_i[3])
      # Store interpolated color
      colors_for_missing_ranges[i] <- col_range_i
    }
    names(colors_for_missing_ranges) <- missing_ranges

    # Initiate vector of colors for all ranges found in the model
    colors_list_for_ranges <- setNames(object = effective_ranges_list, nm = effective_ranges_list)
    # Replace already provided colors
    colors_list_for_ranges[match(x = names(colors_per_levels), table = names(colors_list_for_ranges))] <- colors_per_levels
    # Replace interpolated colors for missing range
    colors_list_for_ranges[match(x = names(colors_for_missing_ranges), table = names(colors_list_for_ranges))] <- colors_for_missing_ranges
    # Provide color for null range
    colors_list_for_ranges[names(colors_list_for_ranges) == "_"] <- "#000000"
  }

  ## Identify ranges with a non-zero residence time

  total_residence_times_df <- data.frame()
  # Loop per simmaps to compute residence time
  for (i in seq_along(all_simmaps))
  {
    # i <- 1

    # Extract simmap
    simmap_i <- all_simmaps[[i]]
    # Compute residence time per ranges
    total_residence_times_i <- apply(X = simmap_i$mapped.edge, MARGIN = 2, FUN = sum)
    # Store in total_residence_times_d
    total_residence_times_df <- rbind(total_residence_times_df, total_residence_times_i)
  }
  names(total_residence_times_df) <- names(total_residence_times_i)

  # Sum residence times across simmaps
  total_residence_times_per_ranges <- apply(total_residence_times_df, MARGIN = 2, FUN = sum)

  # Select ranges with non-null residence time
  selected_ranges <- names(total_residence_times_per_ranges)[total_residence_times_per_ranges > 0]

  ## Create densityMaps from simmaps

  # Initiate list of densityMaps
  densityMaps_all_ranges <- list()

  ## Loop per selected ranges
  for (i in seq_along(selected_ranges))
  {
    # i <- 1

    # Extract range
    range_i <- selected_ranges[i]
    # Define other ranges by contrast
    other_ranges <- selected_ranges[selected_ranges != range_i]

    ## Binarize ranges in simmaps
    simmaps_binary_i <- all_simmaps
    simmaps_binary_i <- lapply(X = simmaps_binary_i, FUN = phytools::mergeMappedStates, old.states = other_ranges, new.state = "0")
    simmaps_binary_i <- lapply(X = simmaps_binary_i, FUN = phytools::mergeMappedStates, old.states = range_i, new.state = "1")

    class(simmaps_binary_i) <- c("list", "multiSimmap", "multiPhylo")

    # Check that remaining ranges are all binary
    # unique(unlist(lapply(X = simmaps_binary_i, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
    # simmaps_binary_i[[1]]$maps

    ## Estimate the posterior probabilities of ranges along all branches (from the set of simulated maps)

    densityMap_range_i <- densityMap_custom(trees = simmaps_binary_i,
                                            tol = 1e-5, verbose = verbose,
                                            col_scale = NULL,
                                            plot = FALSE)

    ## Update color gradient

    # Set color gradient from grey to focal color
    focal_color <- colors_list_for_ranges[selected_ranges][i]
    col_fn <- grDevices::colorRampPalette(colors = c("grey90", focal_color))
    col_scale <- col_fn(n = 1001)

    # Update color gradient
    densityMap_range_i <- phytools::setMap(densityMap_range_i, c("grey90", focal_color))

    ## Update range names
    densityMap_range_i$states <- c(paste0("Not ", range_i), range_i)

    # plot(densityMap_range_i)

    ## Store in final object with all density maps
    densityMaps_all_ranges[[i]] <- densityMap_range_i

    ## Print progress for each range
    if (verbose)
    {
      cat(paste0(Sys.time(), " - Posterior probabilities computed for Range = ",range_i," - n\u00B0", i, "/", length(selected_ranges),"\n"))
    }
  }
  names(densityMaps_all_ranges) <- paste0("Density_map_", selected_ranges)

  ## Aggregate multi-are ranges in densityMaps and ace_matrix if requested
  if (split_multi_area_ranges)
  {
    # Initiate new densityMaps
    densityMaps_unique_areas <- densityMaps_all_ranges

    # Detect unique areas vs. multi-area ranges
    selected_areas <- selected_ranges[nchar(selected_ranges) == 1]
    ranges_to_split <- setdiff(selected_ranges, selected_areas)
    areas_list_for_ranges_to_split <- strsplit(x = ranges_to_split, split = "")

    ## Loop per ranges to split
    for (i in seq_along(ranges_to_split))
    {
      # i <- i

      # Extract multi-area range to split
      range_to_split_i <- ranges_to_split[i]
      # Extract associated areas to blend-in
      areas_list_for_range_to_split_i <- areas_list_for_ranges_to_split[[i]]

      # Update densityMaps
      densityMaps_unique_areas <- split_state_densityMaps(densityMaps = densityMaps_unique_areas,
                                                          state_to_split = range_to_split_i,
                                                          states_to_blend_in = areas_list_for_range_to_split_i)
    }

    # Compare before and after splitting
    # plot_densityMaps_overlay(densityMaps_all_ranges)
    # plot_densityMaps_overlay(densityMaps_unique_areas)

    ## Updata ace_matrix

    # Initiate new ace_matrix
    ace_matrix_unique_areas <- ace_matrix

    ## Loop per ranges to split
    for (i in seq_along(ranges_to_split))
    {
      # i <- 1

      # Extract multi-area range to split
      range_to_split_i <- ranges_to_split[i]
      # Extract associated areas to blend-in
      areas_list_for_range_to_split_i <- areas_list_for_ranges_to_split[[i]]

      # Identify nb of areas to blend-in
      nb_areas_to_blend_in <- length(areas_list_for_range_to_split_i)

      # Extract PP for range_to_split
      PP_to_split_i <- ace_matrix_unique_areas[, range_to_split_i]
      # Split by nb of unique areas to blend-in
      PP_split_i <- PP_to_split_i/nb_areas_to_blend_in

      # Add PP to unique areas to blend-in
      ace_matrix_unique_areas[, areas_list_for_range_to_split_i] <- ace_matrix_unique_areas[, areas_list_for_range_to_split_i] + PP_split_i
    }
    # Avoid trespassing PP > 1 due to rounding
    ace_matrix_unique_areas[ace_matrix_unique_areas[] > 1] <- 1

    # Remove split ranges
    ace_matrix_unique_areas <- ace_matrix_unique_areas[, !(colnames(ace_matrix_unique_areas) %in% ranges_to_split)]

    # View(ace_matrix_unique_areas)
  }

  ## Plot densityMaps (if plot_map)
  if (plot_map)
  {
    # Allows plotting outside of figure margins
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    if (!plot_overlay)
    {
      ## Print progress for each range
      if (verbose)
      {
        cat(paste0("\n", Sys.time(), " - Plot one densityMap for each range.\n"))
      }

      ## Plot one densityMap per range
      for (i in seq_along(densityMaps_all_ranges))
      {
        plot(densityMaps_all_ranges[[i]])
      }
    } else {

      ## Print progress for each range
      if (verbose)
      {
        cat(paste0("\n", Sys.time(), " - Plot a unique densityMap with for all ranges overlaid.\n"))
      }

      ## Filter ace_matrix to display only slected ranges for internal nodes
      ace_matrix_for_plot <- ace_matrix[(nb_tips+1):nrow(ace_matrix), selected_ranges]

      ## Plot the overlay of densityMaps with alpha
      plot_densityMaps_overlay(densityMaps = densityMaps_all_ranges,
                               add_ACE_pies = TRUE,
                               ace = ace_matrix_for_plot)
    }

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Export densityMap(s) in PDF if requested
  if (!is.null(PDF_file_path))
  {
    ## Print progress for each state
    if (verbose)
    {
      cat(paste0("\n", Sys.time(), " - Export PDF.\n"))
    }

    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    # Adjust width/height according to the nb of tips
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10

    if (!plot_overlay)
    {
      ## Plot one densityMap per state
      grDevices::pdf(file = file.path(PDF_file_path),
                     width = width, height = height)

      for (i in seq_along(densityMaps_all_ranges))
      {
        plot(densityMaps_all_ranges[[i]])
      }

      grDevices::dev.off()

    } else {

      ## Filter ace_matrix to display only selected ranges for internal nodes
      ace_matrix_for_plot <- ace_matrix[(nb_tips+1):nrow(ace_matrix), selected_ranges]

      ## Plot the overlay of densityMaps with alpha
      grDevices::pdf(file = file.path(PDF_file_path),
                     width = width, height = height)

      plot_densityMaps_overlay(densityMaps = densityMaps_all_ranges,
                               add_ACE_pies = TRUE,
                               ace = ace_matrix_for_plot)

      grDevices::dev.off()
    }

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Clean BioGeoBEARS output files if needed
  if (!keep_BioGeoBEARS_files)
  {
    if (is.null(prefix_for_files))
    {
      files_path <- BioGeoBEARS_directory_path
    } else {
      files_path <- paste0(BioGeoBEARS_directory_path, prefix_for_files, "_")
    }

    # Remove input files
    file.remove(file.path(paste0(BioGeoBEARS_directory_path,"phylo.tree")))
    file.remove(file.path(paste0(BioGeoBEARS_directory_path,"tip_ranges.data")))

    # Remove BSM files
    file.remove(file.path(paste0(BioGeoBEARS_directory_path,"RES_ana_events_tables.Rdata")))
    file.remove(file.path(paste0(BioGeoBEARS_directory_path,"RES_ana_events_tables_PARTIAL.Rdata")))
    file.remove(file.path(paste0(BioGeoBEARS_directory_path,"RES_clado_events_tables.Rdata")))

    # Remove models outputs
    if ("DEC" %in% evolutionary_models)
    {
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DEC_run_log.txt")))
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DEC_fit.rds")))
    }
    if ("DEC+J" %in% evolutionary_models)
    {
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DEC_J_run_log.txt")))
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DEC_J_fit.rds")))
    }
    if ("DIVALIKE" %in% evolutionary_models)
    {
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DIVALIKE_run_log.txt")))
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DIVALIKE_fit.rds")))
    }
    if ("DIVALIKE+J" %in% evolutionary_models)
    {
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DIVALIKE_J_run_log.txt")))
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"DIVALIKE_J_fit.rds")))
    }
    if ("BAYAREALIKE" %in% evolutionary_models)
    {
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_run_log.txt")))
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_fit.rds")))
    }
    if ("BAYAREALIKE+J" %in% evolutionary_models)
    {
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_J_run_log.txt")))
      file.remove(file.path(paste0(BioGeoBEARS_directory_path,"BAYAREALIKE_J_fit.rds")))
    }

    # If empty, remove the BioGeoBEARS_directory
    if (length(list.files(path = file.path(BioGeoBEARS_directory_path))) == 0)
    {
      unlink(x = file.path(BioGeoBEARS_directory_path), force = TRUE)
    }
  }

  ## Build output
  if (!split_multi_area_ranges)
  {
    # With all ranges
    output <- list(densityMaps = densityMaps_all_ranges,
                   trait_data_type = "biogeographic")
  } else {
    # With muti-area ranges split into unique areas + all ranges
    output <- list(densityMaps = densityMaps_unique_areas,
                   densityMaps_all_ranges = densityMaps_all_ranges,
                   trait_data_type = "biogeographic")
  }

  # Include ACE if requested
  if(return_ace)
  {
    if (!split_multi_area_ranges)
    {
      output$ace <- ace_matrix[(nb_tips+1):nrow(ace_matrix), selected_ranges] # Filter to include only internal nodes (to be consistent with continuous trait)
    } else {
      output$ace <- ace_matrix_unique_areas[(nb_tips+1):nrow(ace_matrix), selected_areas] # Filter to include only internal nodes (to be consistent with continuous trait)
      output$ace_all_ranges <- ace_matrix[(nb_tips+1):nrow(ace_matrix), selected_ranges] # Filter to include only internal nodes (to be consistent with continuous trait)
    }
  }

  # Include BSM_output if requested
  if(return_BSM) { output$BSM_output <- BSM_output }

  # Include simmaps if requested
  if(return_simmaps) { output$simmaps <- all_simmaps }

  # Include output of best model if requested
  if(return_best_model_fit) { output$best_model_fit <- best_model_fit }

  # Include df for model comparison if requested
  if(return_model_selection_df) { output$model_selection_df <- models_comparison$models_comparison_df }

  ## Return output
  return(invisible(output))
}


### Helper function to compare trait model fits from geiger::fitContinuous or geiger::fitDiscrete with AICc and Akaike's weights ####
# Input = list of models fit with geiger::fitContinuous or geiger::fitDiscrete

#' @title Compare trait evolutionary model fits with AICc and Akaike's weights
#'
#' @description Compare evolutionary models fit with [geiger::fitContinuous()] or [geiger::fitDiscrete()]
#'   using AICc and Akaike's weights. Generate a data.frame summarizing information.
#'   Identify the best model and extract its results.
#'
#' @param list_model_fits Named list with the results of a model fit with [geiger::fitContinuous()] or [geiger::fitDiscrete()] in each element.
#'
#' @importFrom phytools aic.w
#' @export
#'
#' @return The function returns a list with three elements.
#' * `$model_comparison_df` Data.frame summarizing information to compare model fits. It includes the model name (`$model`),
#'   the log-likelihood (`$logLik`), the number of free-parameters (`$k`), the AIC (`$AICc`), the corrected AIC (`$AICc`),
#'   the delta to the best/lowest AICc (`$delta_AICc`), the Akaike weights (`$Akaike_weights`), and their rank based on AICc (`$rank`).
#' * `$best_model_name` Character string. Name of the best model.
#' * `$best_model_fit` List containing the output of [geiger::fitContinuous()] or [geiger::fitDiscrete()] for the model with the best fit.
#'
#' @author Maël Doré
#'
#' @seealso [geiger::fitContinuous()] [geiger::fitDiscrete()]
#'
#' @examples
#'
#' # ----- Example 1: Continuous data ----- #
#'
#' # Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract body size
#' eel_data <- stats::setNames(eel.data$Max_TL_cm,
#'                             rownames(eel.data))
#'
#' # Fit BM model
#' BM_fit <- geiger::fitContinuous(phy = eel.tree, dat = eel_data, model = "BM")
#' # Fit EB model
#' EB_fit <- geiger::fitContinuous(phy = eel.tree, dat = eel_data, model = "EB")
#' # Fit OU model
#' OU_fit <- geiger::fitContinuous(phy = eel.tree, dat = eel_data, model = "OU")
#'
#' # Store models
#' list_model_fits <- list(BM = BM_fit, EB = EB_fit, OU = OU_fit)
#'
#' # Compare models
#' model_comparison_output <- select_best_trait_model_from_geiger(list_model_fits = list_model_fits)
#'
#' # Explore output
#' str(model_comparison_output, max.level = 2)
#'
#' # Print comparison
#' print(model_comparison_output$models_comparison_df)
#'
#' # Print best model fit
#' print(model_comparison_output$best_model_fit)
#'
#' # ----- Example 2: Categorical data ----- #
#'
#' # Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract feeding mode data
#' eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
#' table(eel_data)
#'
#' # Fit ER model
#' ER_fit <- geiger::fitDiscrete(phy = eel.tree, dat = eel_data, model = "ER")
#' # Fit ARD model
#' ARD_fit <- geiger::fitDiscrete(phy = eel.tree, dat = eel_data, model = "ARD")
#'
#' # Store models
#' list_model_fits <- list(ER = ER_fit, ARD = ARD_fit)
#'
#' # Compare models
#' model_comparison_output <- select_best_trait_model_from_geiger(list_model_fits = list_model_fits)
#'
#' # Explore output
#' str(model_comparison_output, max.level = 2)
#'
#' # Print comparison
#' print(model_comparison_output$models_comparison_df)
#'
#' # Print best model fit
#' print(model_comparison_output$best_model_fit)
#'

select_best_trait_model_from_geiger <- function (list_model_fits)
{
  # Homemade function to extract the number of free-parameters from geiger model outputs
  extract_k <- function (x)
  {
    k <- x$opt$k
    return(k)
  }

  # Homemade function to extract AICc from geiger model outputs
  extract_AIC <- function (x)
  {
    AIC <- x$opt$aic
    return(AIC)
  }

  # Homemade function to extract AICc from geiger model outputs
  extract_AICc <- function (x)
  {
    AICc <- x$opt$aicc
    return(AICc)
  }

  # Extract ln-likelihood, number of parameters, and AICc from models' list
  models_comparison_df <- data.frame(model = names(list_model_fits),
                                     logL = sapply(X = list_model_fits, FUN = stats::logLik),
                                     k = sapply(X = list_model_fits, FUN = extract_k),
                                     AIC = sapply(X = list_model_fits, FUN = extract_AIC),
                                     AICc = sapply(X = list_model_fits, FUN = extract_AICc))

  # Compute Delta AICc
  best_AICc <- min(models_comparison_df$AICc)
  models_comparison_df$delta_AICc <- models_comparison_df$AICc - best_AICc

  # Compute Akaike's weights based on AICc
  models_comparison_df$Akaike_weights <- round(phytools::aic.w(models_comparison_df$AICc), 3) * 100

  # Compute ranks based on AICc (The lowest = the best)
  models_comparison_df$rank <- rank(models_comparison_df$AICc, na.last = "keep", ties.method = "first")

  # Extract name of best model
  best_model_ID <- which(models_comparison_df$rank == 1)
  best_model_name <- names(list_model_fits)[best_model_ID]

  # Extract output of best model
  best_model_fit <- list_model_fits[[best_model_ID]]

  # Export output without printing best model fit
  return(invisible(list(models_comparison_df = models_comparison_df,
                        best_model_name = best_model_name,
                        best_model_fit = best_model_fit)))
}


### Helper functions to extract info from BioGeoBEARS model outputs ####
# Input = list of models fit with BioGeoBEARS::bear_optim_run()

#' @title Compare model fits with AICc and Akaike's weights
#'
#' @description Compare models fit with `BioGeoBEARS::bear_optim_run()` using AICc and Akaike's weights.
#'   Generate a data.frame summarizing information. Identify the best model and extract its results.
#'
#' @param list_model_fits Named list with the results of a model fit with `BioGeoBEARS::bear_optim_run()` in each element.
#'
#' @importFrom phytools aic.w
#' @export
#'
#' @return The function returns a list with three elements.
#' * `$model_comparison_df` Data.frame summarizing information to compare model fits. It includes the model name (`$model`),
#'   the log-likelihood (`$logLik`), the number of free-parameters (`$k`), the AIC (`$AIC`), the corrected AIC (`$AICc`),
#'   the delta to the best/lowest AICc (`$delta_AICc`), the Akaike weights (`$Akaike_weights`), and their rank based on AICc (`$rank`).
#' * `$best_model_name` Character string. Name of the best model.
#' * `$best_model_fit` List containing the output of `BioGeoBEARS::bear_optim_run()` for the model with the best fit.
#'
#' @author Maël Doré
#'
#' @seealso `BioGeoBEARS::bear_optim_run()` `BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object()` `BioGeoBEARS::AICstats_2models()`
#'
#' @examples
#'
#' # Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Transform feeding mode data into biogeographic data with ranges A, B, and AB.
#' eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
#' eel_data <- as.character(eel_data)
#' eel_data[eel_data == "bite"] <- "A"
#' eel_data[eel_data == "suction"] <- "B"
#' eel_data[c(5, 6, 7, 15, 25, 32, 33, 34, 50, 52, 57, 58, 59)] <- "AB"
#' eel_data <- stats::setNames(eel_data, rownames(eel.data))
#' table(eel_data)
#'
#' colors_per_levels <- c("dodgerblue3", "gold")
#' names(colors_per_levels) <- c("A", "B")
#'
#' ### Uncomment to run ###
#'
#' # ## The R package `BioGeoBEARS` is needed for this function to work with biogeographic data.
#' # Please install it manually from: \href{https://github.com/nmatzke/BioGeoBEARS}{https://github.com/nmatzke/BioGeoBEARS}.
#' #
#' # ## Prepare phylo
#' #
#' # # Set path to BioGeoBEARS directory
#' # BioGeoBEARS_directory_path = "./BioGeoBEARS_directory/"
#' #
#' # # Export phylo
#' # path_to_phylo <- BioGeoBEARS::np(paste0(BioGeoBEARS_directory_path, "phylo.tree"))
#' # write.tree(phy = eel.tree, file = path_to_phylo)
#' #
#' # ## Prepare range data
#' #
#' # tip_data <- eel_data
#' # # Convert tip_data to df
#' # ranges_df <- as.data.frame(tip_data)
#' # # Get list of all unique areas
#' # unique_areas_in_ranges_list <- strsplit(x = tip_data, split = "")
#' # unique_areas <- unique(unlist(unique_areas_in_ranges_list))
#' # unique_areas <- unique_areas[order(unique_areas)]
#' #
#' # # Loop per unique areas
#' # for (i in seq_along(unique_areas))
#' # {
#' #   # i <- 1
#' #
#' #   # Extract unique area
#' #   unique_area_i <- unique_areas[i]
#' #   # Detect presence in ranges
#' #   binary_match_i <- unlist(lapply(X = unique_areas_in_ranges_list,
#' #      FUN = function (x) { unique_area_i %in% x } ))
#' #
#' #   # Add to ranges_df
#' #   ranges_df <- cbind(ranges_df, binary_match_i)
#' # }
#' # # Extract binary df of presence/absence
#' # binary_df <- ranges_df[, -1]
#' # Convert character strings into numerical factors
#' # binary_df_num <- as.data.frame(apply(X = binary_df, MARGIN = 2, FUN = as.numeric))
#' # row.names(binary_df_num) <- names(tip_data)
#' # names(binary_df_num) <- unique_areas
#' #
#' # # Produce tipranges object from numeric df
#' # Taxa_bioregions_tipranges_obj <- BioGeoBEARS::define_tipranges_object(tmpdf = binary_df_num)
#' #
#' # # Set path to tip ranges object
#' # path_to_tip_ranges <- BioGeoBEARS::np(paste0(BioGeoBEARS_directory_path,"tip_ranges.data"))
#' #
#' # # Export tip ranges in Lagrange/PHYLIP format
#' # BioGeoBEARS::save_tipranges_to_LagrangePHYLIP(
#' #   tipranges_object = Taxa_bioregions_tipranges_obj,
#' #   lgdata_fn = path_to_tip_ranges,
#' #   areanames = colnames(Taxa_bioregions_tipranges_obj@df))
#' #
#' # ## Prepare models
#' #
#' # # Prepare DEC model run
#' # DEC_run <- BioGeoBEARS::define_BioGeoBEARS_run(
#' #    num_cores_to_use = 1,
#' #    max_range_size = 2, # To set the maximum number of bioregion encompassed
#' #       # by a lineage range at any time
#' #    trfn = path_to_phylo, # To provide path to the input tree file
#' #    geogfn = path_to_tip_ranges, # To provide path to the LagrangePHYLIP file with binary ranges
#' #    return_condlikes_table = TRUE, # To ask to obtain all marginal likelihoods
#' #       # computed by the model and used to display ancestral states
#' #    print_optim = TRUE)
#' #
#' # # Check that starting parameter values are inside the min/max
#' # DEC_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_run)
#' # # Check validity of set-up before run
#' # BioGeoBEARS::check_BioGeoBEARS_run(DEC_run)
#' #
#' # # Use the DEC model as template for DEC+J model
#' # DEC_J_run <- DEC_run
#' #
#' # # Update status of jump speciation parameter to be estimated
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
#' # # Set initial value of J for optimization to an arbitrary low non-null value
#' # j_start <- 0.0001
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["j","init"] <- j_start
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["j","est"] <- j_start
#' #
#' # # Check validity of set-up before run
#' # DEC_J_run <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = DEC_J_run)
#' # invisible(BioGeoBEARS::check_BioGeoBEARS_run(DEC_J_run))
#' #
#' # ## Run models
#' #
#' # # Run DEC model
#' # DEC_fit <- BioGeoBEARS::bears_optim_run(DEC_run)
#' #
#' # # Set starting values for optimization of DEC+J model based on MLE in the DEC model
#' # d_start <- DEC_fit$outputs@params_table["d","est"]
#' # e_start <- DEC_fit$outputs@params_table["e","est"]
#' #
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["d","init"] <- d_start
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["d","est"] <- d_start
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["e","init"] <- e_start
#' # DEC_J_run$BioGeoBEARS_model_object@params_table["e","est"] <- e_start
#' #
#' # # Run DEC+J model
#' # DEC_J_fit <- BioGeoBEARS::bears_optim_run(DEC_J_run)
#' #
#' # ## Store model outputs
#' # list_model_fits <- list(DEC = DEC_fit, DEC_J = DEC_J_fit)
#' #
#' # ## Compare models
#' # model_comparison_output <- select_best_model_from_BioGeoBEARS(list_model_fits = list_model_fits)
#' #
#' # # Explore output
#' # str(model_comparison_output, max.level = 1)
#' # # Print comparison
#' # print(model_comparison_output$models_comparison_df)
#' # # Print best model fit
#' # print(model_comparison_output$best_model_fit$outputs)
#' #

select_best_model_from_BioGeoBEARS <- function (list_model_fits)
{
  ## Control for BioGeoBEARS install
  if (!requireNamespace("BioGeoBEARS", quietly = TRUE))
  {
    stop("Package 'BioGeoBEARS' is needed for this function to work.
       Please install it manually from: https://github.com/nmatzke/BioGeoBEARS")
  }

  # Homemade function to extract the number of observations from BioGeoBEARS model outputs
  extract_nobs_from_BioGeoBEARS_output <- function (x)
  {
    nobs <- (length(x$computed_likelihoods_at_each_node) + 1)/2
    names(nobs) <- "nobs"
    return(nobs)
  }

  # Homemade function to extract the number of free-parameters from BioGeoBEARS model outputs
  extract_k_from_BioGeoBEARS_output <- function (x)
  {
    k <- sum(x$inputs$BioGeoBEARS_model_object@params_table$type == "free")
    names(k) <- "k"
    return(k)
  }

  # Home-made function to compute AIC from logLk and number of parameters (k)
  compute_AIC <- function (logLk, k)
  {
    AIC <- - 2 * (logLk - k)
  }

  # Home-made function to compute AICc from AIC, number of observations, and number of parameters (k)
  compute_AICc <- function (AIC, nobs, k)
  {
    AICc <- AIC + (2 * k) * (k - 1) / (nobs - k - 1)
  }

  # Extract nb of observation from the first model
  nb_obs <- extract_nobs_from_BioGeoBEARS_output(list_model_fits[[1]])

  # Extract ln-likelihood, number of parameters, and AICc from models' list
  models_comparison_df <- data.frame(model = names(list_model_fits),
                                     logLk = sapply(X = list_model_fits, FUN = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object),
                                     k = sapply(X = list_model_fits, FUN = extract_k_from_BioGeoBEARS_output))

  # Compute AIC from logLk and k
  models_comparison_df$AIC <- unlist(purrr::map2(.x = models_comparison_df$logLk, .y = models_comparison_df$k, .f = compute_AIC))

  # Compute AICc from AIC, k and nobs
  models_comparison_df$AICc <- unlist(purrr::map2(.x = models_comparison_df$AIC, nobs = nb_obs, .y = models_comparison_df$k, .f = compute_AICc))

  # Compute Delta AICc
  best_AICc <- min(models_comparison_df$AICc)
  models_comparison_df$delta_AICc <- models_comparison_df$AICc - best_AICc

  # Compute Akaike's weights based on AICc
  models_comparison_df$Akaike_weights <- round(phytools::aic.w(models_comparison_df$AICc), 3) * 100

  # Compute ranks based on AICc (The lowest = the best)
  models_comparison_df$rank <- rank(models_comparison_df$AICc, na.last = "keep", ties.method = "first")

  # Extract name of best model
  best_model_ID <- which(models_comparison_df$rank == 1)
  best_model_name <- names(list_model_fits)[best_model_ID]

  # Extract output of best model
  best_model_fit <- list_model_fits[[best_model_ID]]

  # Export output without printing best model fit
  return(invisible(list(models_comparison_df = models_comparison_df,
                        best_model_name = best_model_name,
                        best_model_fit = best_model_fit)))
}


### Helper function to generate list of ranges from unique areas ####

#'
#' @importFrom cladoRcpp rcpp_areas_list_to_states_list
#'
#' @noRd
#'

generate_list_ranges <- function (areas_list, max_range_size, include_null_range = TRUE)
{
  states_list_0based = cladoRcpp::rcpp_areas_list_to_states_list(areas = areas_list, maxareas = max_range_size, include_null_range = include_null_range)
  ranges_list <- NULL
  for (i in 1:length(states_list_0based))
  {
    if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
    {
      tmprange = "_"
    } else {
      tmprange = paste(areas_list[states_list_0based[[i]]+1], collapse="")
    }
    ranges_list = c(ranges_list, tmprange)
  }
  return(ranges_list)
}

### Helper function to produce densityMap from simmaps ####
# Original function written by Liam Revell, 2012
# Input = simmaps

#' @title Plot posterior density of stochastic mapping on a tree
#'
#' @description Visualize posterior probability density from stochastic mapping using a color gradient on the tree.
#'   Original function written by Liam Revell, 2012 in the [phytools] package: [phytools::densityMap()].
#'
#' @inheritParams phytools::densityMap
#' @param tol Positive numerical. To set the tolerance used to match node ages and time steps (i.e., onsider them equal). Default = 1e-5.
#' @param verbose Logical. To display or progress every 100 edges. Default = `TRUE`.
#' @param col_scale Character string vector. To set the color scale manually. Need to provide 1001 colors for the scale.
#'   If `NULL` (the default), the `rainbow()` color scale will be used.
#'
#' @return The function plots a tree with mapped trait probability densities and returns an object of class `densityMap` invisibly.
#'   A `densityMap` is a list with three elements.
#'     * `$tree` List of at least 8 elements. Includes the phylogeny, the trait evolution model data from the simmaps, and the newly mapped trait posterior densities.
#'       * `$maps` List of N elements, one per edge. Each list comprises a named numerical vector that represent changes in posterior probability density of the focal state along segments of equal time.
#'         Named are posterior probabilities scaled from 0 to 1000. Values are length of the segments. Segments are ordered from root to tips.
#'       * `$mapped.edge` Matrix of edge per posterior probability summarizing the overall length of each edge attributed to a specific posterior probabiliy value.
#'       * `$Q` Numerical square matrix summarizing instantaneous transition rates between states as estimated from the evolutionary model.
#'       Rows = initial states. Cols = final states.
#'       * `$logL` Numerical. Log-likelihood of the data as optimized when estimated model parameters.
#'     * `$col` Named character string vector. Color scale used to map posterior probabilities. Names are the posterior probabilities scaled from 0 to 1000. Values are the colors.
#'     * `$states` Character string. The name of the states.
#'
#' @details Wrapped function of [phytools::densityMap()].
#'   Additions to the initial function:
#'   * Can modify manually the tolerance to handle issue with mismatch between node ages and time steps used.
#'   * Can print progress across egdes
#'   * Can provide a manual color scale to replace the default rainbow scale. The color scale must have 1001 colors.
#'
#' @author Maël Doré. Initial function by Liam Revell, 2012 in the [phytools] package.
#'
#' @seealso [phytools::densityMap()]
#'
#' @noRd
#'

densityMap_custom <- function (trees, res = 100, fsize = NULL, ftype = NULL, lwd = 3,
                               tol = 1e-5, verbose = T, col_scale = NULL,
                               check = FALSE, legend = NULL, outline = FALSE, type = "phylogram",
                               direction = "rightwards", plot = TRUE, ...)
{
  # Set graphical parameters
  if (methods::hasArg(mar))
    mar <- list(...)$mar
  else mar <- rep(0.3, 4)
  if (methods::hasArg(offset))
    offset <- list(...)$offset
  else offset <- NULL
  if (methods::hasArg(states))
    states <- list(...)$states
  else states <- NULL
  if (methods::hasArg(hold))
    hold <- list(...)$hold
  else hold <- TRUE
  if (length(lwd) == 1)
    lwd <- rep(lwd, 2)
  else if (length(lwd) > 2)
    lwd <- lwd[1:2]

  # Adjust tolerance
  tol <- tol

  # Check validity of the class of the input object
  if (!inherits(trees, "multiPhylo") && inherits(trees, "phylo"))
  {
    stop("trees not \"multiPhylo\" object; just use plotSimmap.")
  }
  if (!inherits(trees, "multiPhylo"))
  {
    stop("trees should be an object of class \"multiPhylo\".")
  }

  # Extract root age
  h <- sapply(unclass(trees), function(x) max(phytools::nodeHeights(x)))

  # Define time steps
  steps <- 0:res/res * max(h)

  # Rescale trees to ensure they all have the same root age
  trees <- phytools::rescaleSimmap(trees, totalDepth = max(h))

  # Check that phylogeny topology and branch length are equal
  if (check)
  {
    X <- matrix(FALSE, length(trees), length(trees))
    for (i in 1:length(trees)) X[i, ] <- sapply(trees, ape::all.equal.phylo,
                                                current = trees[[i]])
    if (!all(X))
      stop("some of the trees don't match in topology or relative branch lengths")
  }

  # Extract first tree as reference
  tree <- trees[[1]]

  # Remove class
  trees <- unclass(trees)

  # Extract all states from the first tree (dangerous if some states are not present in this tree but in other!)
  if (is.null(states))
  {
    ss <- sort(unique(c(phytools::getStates(tree, "nodes"), phytools::getStates(tree, "tips"))))
  }  else {
    ss <- states
  }

  # If states are not binary, rename the first two states as "0" and "1" (dangerous as if there are more states, will lead to errors)
  if (!all(ss == c("0", "1")))
  {
    c1 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    c2 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    trees <- lapply(trees, phytools::mergeMappedStates, ss[1], c1)
    trees <- lapply(trees, phytools::mergeMappedStates, ss[2], c2)
    trees <- lapply(trees, phytools::mergeMappedStates, c1, "0")
    trees <- lapply(trees, phytools::mergeMappedStates, c2, "1")
  }

  # Extract all node ages per edge
  H <- phytools::nodeHeights(tree)
  # message("sorry - this might take a while; please be patient")

  # Reinitiate the map of the reference tree with NULL data
  tree$maps <- vector(mode = "list", length = nrow(tree$edge))

  # Loop per edge/item in tree$maps
  for (i in 1:nrow(tree$edge))
  {
    # i <- 5
    # i <- 983

    # YY = Matrix of ages to use as time step along edge i
    # Include start (0) and end (edge length)
    # One raw = one interval
    # Columns = Start and End relative age
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]),
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) - H[i, 1]

    # Initiate vector of final time step values for edge i
    ZZ <- rep(0, nrow(YY))

    # Loop per trees/simmaps
    for (j in 1:length(trees))
    {
      # j <- 1

      # XX = a matrix of the states detected for edge i on simmap j
      # One row = one state
      # Colmuns = start and end relative age
      XX <- matrix(data = 0,
                   nrow = length(trees[[j]]$maps[[i]]),
                   ncol = 2,
                   dimnames = list(names(trees[[j]]$maps[[i]]),
                                   c("start", "end")))
      # Fill the first raw with information on start and end relative age of the first state
      XX[1, 2] <- trees[[j]]$maps[[i]][1]

      # Case with multiple states: fill information for other states
      if (length(trees[[j]]$maps[[i]]) > 1)
      {
        for (k in 2:length(trees[[j]]$maps[[i]]))
        {
          XX[k, 1] <- XX[k - 1, 2]
          XX[k, 2] <- XX[k, 1] + trees[[j]]$maps[[i]][k]
        }
      }

      # Loop per time interval wanted for the density mapping
      for (k in 1:nrow(YY))
      {
        # k <- 1

        # Detect which state start before the k time step
        lower <- which(XX[, 1] <= YY[k, 1])
        lower <- lower[length(lower)] # Take the last one as the last state recorded before the beginning of the time step

        # Detect which state end after the k time step
        upper <- which(XX[, 2] >= (YY[k, 2] - tol))[1]


        AA <- 0
        names(lower) <- names(upper) <- NULL
        if (!all(XX == 0))
        {
          # Case for internal edge (end time > 0)

          # Loop per states on the time interval
          for (l in lower:upper)
          {
            # Compute weighted mean of the time interval
            AA <- AA + (min(XX[l, 2], YY[k, 2]) - max(XX[l, 1], YY[k, 1]))/(YY[k, 2] - YY[k, 1]) * as.numeric(rownames(XX)[l])
          }
        } else {
          # Case for tips (or null branches) (start and end time = 0)
          AA <- as.numeric(rownames(XX)[1]) # Use the tip state
        }
        # Increment the final value of the edge i
        ZZ[k] <- ZZ[k] + AA/length(trees)
      }
    }
    # Record time steps length
    tree$maps[[i]] <- YY[, 2] - YY[, 1]

    # Record time steps continuous value as names of the maps of edge i
    names(tree$maps[[i]]) <- round(ZZ * 1000) # Convert proportion in a scale from 0 to 1000

    # Print progress every 100 edges
    if ((i %% 100 == 0) & verbose)
    {
      cat(paste0(Sys.time(), " - Posterior probability computed for edge n\u00B0", i, "/", nrow(tree$edge),"\n"))
    }
  }

  # Create color scale
  if (is.null(col_scale))
  {
    # cols <- grDevices::rainbow(n =1001, start = 0.7, end = 0) # Initial HSV color scale from phytools::contMap
    cols <- grDevices::rainbow(n = 1001, s = 0.8, start = 0.7, end = 0) # HSV color space with less saturation
    # cols <- colorspace::rainbow_hcl(n = 1001, c = 90, l = 65, start = 240) # With HCL color space
    names(cols) <- 0:1000
  } else {
    cols <- col_scale
  }

  # Recreate map using continuous values
  tree$mapped.edge <- makeMappedEdge(tree$edge, tree$maps)
  tree$mapped.edge <- tree$mapped.edge[, order(as.numeric(colnames(tree$mapped.edge)))]

  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  attr(tree, "map.order") <- "right-to-left"
  x <- list(tree = tree, cols = cols, states = ss)
  class(x) <- "densityMap"

  # Plot
  if (plot)
    plot(x, fsize = fsize, ftype = ftype, lwd = lwd,
         legend = legend, outline = outline, type = type,
         mar = mar, direction = direction, offset = offset,
         hold = hold)

  # Return final simmap
  invisible(x)
}


### Helper functions to convert BioGeoBEARS BSM outputs into phytools.simmaps ####

#' @title Convert Biogeographic Stochastic Map (BSM) to phytools SIMMAP stochastic map (SM) format
#'
#' @description These functions converts a Biogeographic Stochastic Map (BSM) output from BioGeoBEARS into
#'  a `simmap` object from R package [phytools] (See [phytools::make.simmap()]).
#'
#'  They require a model fit with `BioGeoBEARS::bears_optim_run()` and the output of a Biogeographic Stochastic Mapping
#'  performed with `BioGeoBEARS::runBSM()` to produce `simmap` objects as phylogenies with the associated
#'  mapping of range evolution along branches across simulations.
#'
#'   * [BSM_to_phytools_simmap()]: Produce one `simmap` for the required simulation (index of the simulation provided with `sim_index`).
#'   * [BSMs_to_phytools_simmaps()]: Produce all `simmap` objects for all simulations stored in a unique `multiSimmap` object.
#'
#'  Initial functions in R package BioGeoBEARS by Nicholas J. Matzke:
#'  * `BioGeoBEARS::BSM_to_phytools_SM()`
#'  * `BioGeoBEARS::BSMs_to_phytools_SMs()`
#'
#' # Notes on BioGeoBEARS
#'
#'  The R package `BioGeoBEARS` is needed for this function to work with biogeographic data.
#'  Please install it manually from: \href{https://github.com/nmatzke/BioGeoBEARS}{https://github.com/nmatzke/BioGeoBEARS}.
#'
#' # Notes on using the resulting simmap object in phytools (adapted from Nicholas J. Matzke)
#'
#'  The phytools functions, like [phytools::countSimmap()], will only count the anagenetic events
#'  (range transitions occurring along branches) as it was written assuming purely anagenetic models.
#'
#'  It remains possible to extract cladogenetic events (range transitions occurring at speciation)
#'  by comparing the last-state-below-a-node with the descendant-pairs-above-a-node.
#'  However, it is recommended to use the built-in functions from BioGeoBEARS to summarize
#'  the biogeographic history based on the tables of cladogenetic and anagenetic events obtained
#'  from `BioGeoBEARS::runBSM()`. `simmap` objects should primarily be considered as a tool for visualization.
#'
#'  Associated functions in R package `BioGeoBEARS`:
#'
#'  * `BioGeoBEARS::simulate_source_areas_ana_clado()`: To select randomly a unique area source for transition from a multi-area state to a single area.
#'  * `BioGeoBEARS::get_dmat_times_from_res()`: To generate matrices of range expansion from source area to destination area.
#'  * `BioGeoBEARS::count_ana_clado_events()`: To count the number and type of events from BSM tables.
#'  * `BioGeoBEARS::hist_event_counts()`: To plot histograms of event counts across BSM tables.
#'
#'  Please note carefully that area-to-area dispersal events are not identical with the state transitions.
#'  For example, a state can be a geographic range with multipleareas, but under the logic of DEC-type models,
#'  a range-expansion event like  ABC->ABCD actually means that a dispersal happened from some specific area (A, B, or C)
#'  to the new area. BSMs track this area-to-area sourcing in its cladogenetic and anagenetic event tables,
#'  at least if `BioGeoBEARS::simulate_source_areas_ana_clado()` has been run on the output of `BioGeoBEARS::runBSM()`.
#'
#' @param model_fit A BioGeoBEARS results object, produced by ML inference via `BioGeoBEARS::bears_optim_run()`.
#' @param phylo Time-calibrated phylogeny used in the BioGeoBEARS analyses to produce the historical biogeographic inference
#'  and run the Biogeographic Stochastic Mapping. Object of class `"phylo"` as defined in [ape].
#' @param BSM_output A list with two objects, a cladogenetic events table and an anagenetic events table, as the result of
#'   Biogeographic Stochastic Mapping conducted with `BioGeoBEARS::runBSM()`.
#' @param sim_index Integer. Index of the biogeographic simulation targeted to produce the `simmap` with [BSM_to_phytools_simmap()].
#'
#' @export
#' @importFrom BioGeoBEARS get_Qmat_COOmat_from_BioGeoBEARS_run_object get_sum_statetime_on_branch prt
#' @importFrom plyr ldply
#'
#' @details These functions are slight adaptations of original functions from the R Package `BioGeoBEARS` by N. Matzke.
#'  Initial functions: `BioGeoBEARS::BSM_to_phytools_SM()` `BioGeoBEARS::BSMs_to_phytools_SMs()`
#'  Changes:
#'  * Solves issue with differences in ranges allowed across time-strata.
#'  * Requires directly the output of `BioGeoBEARS::runBSM()` instead of separated cladogenetic and anagenetic event tables.
#'  * Update the documentation.
#'
#' @return The [BSM_to_phytools_simmap()] function returns a list with two elements:
#'   * `$simmap` A unique `simmap` for a given biogeographic simulation as an object of classes `c("simmap", "phylo")`.
#'     This is a modified [ape] tree with additional elements to report range mapping, model parameters and likelihood.
#'     - `$maps` A list of named numerical vectors. Provides the mapping of ranges along each remaining edge.
#'       Names are the ranges. Values are residence times in each state across segments
#'     - `$mapped.edge` A numerical matrix. Provides the evolutionary time spent across ranges (columns) along the edges (rows).
#'       row.names() are the node ID at the rootward and tipward ends of each edge.
#'     - `$Q` Numerical matrix. The transition rates across ranges calculated from the ML parameter estimates of the model.
#'     - `$logL` Numeric. The log-likelihood of the data under the ML model.
#'   * `$residence_times` Data.frame with two rows. Summarizes the residence time spent in each range along all branches,
#'    in (raw) evolutionary time (i.e., branch lengths), and in percentage (perc).
#'
#' The [BSMs_to_phytools_simmaps()] function loop around the [BSM_to_phytools_simmap()] function to aggregate all `simmaps`
#'   from all biogeographic simulations in a unique list of classes `c("multiSimmap", "multiPhylo")`.
#'   * Each element in the `$simmap` of a biogeographic simulation obtained with [BSM_to_phytools_simmap()].
#'   * `$residence_times` summary data.frames are not preserved.
#'
#' @seealso [phytools::countSimmap()] [phytools::make.simmap()]
#' `BioGeoBEARS::simulate_source_areas_ana_clado()` `BioGeoBEARS::get_dmat_times_from_res()`
#' `BioGeoBEARS::count_ana_clado_events()` `BioGeoBEARS::hist_event_counts()`
#'
#' @author Nicholas J. Matzke. Contact: \email{matzke@@berkeley.edu}
#' @author Changes by Maël Doré (see Details)
#'
#' @references For BioGeoBEARS: Matzke, Nicholas J. (2018). BioGeoBEARS: BioGeography with Bayesian (and likelihood) Evolutionary Analysis with R Scripts.
#'    version 1.1.1, published on GitHub on November 6, 2018. DOI: \url{http://dx.doi.org/10.5281/zenodo.1478250}. Website: \url{http://phylo.wikidot.com/biogeobears}.
#'
#' @examples
#'
#' # Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#'
#' # Load directly output of prepare_trait_data() run on biogeographic data
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' ## Convert BSM output into a unique simmap, including residence times
#' simmap_1 <- BSM_to_phytools_simmap(model_fit = eel_biogeo_data$best_model_fit,
#'                                    phylo = eel.tree,
#'                                    BSM_output = eel_biogeo_data$BSM_output,
#'                                    sim_index = 1)
#' # Explore output
#' str(simmap_1, max.level = 1)
#' # Print residence times in each range
#' simmap_1$residence_times
#' # Plot simmap
#' plot(simmap_1$simmap)
#'
#'
#' \dontrun{  (May take several minutes to run)
#' ## Convert BSM output into all simmaps in a multiSimmap/multiPhylo object
#' all_simmaps <- BSMs_to_phytools_simmaps(model_fit = eel_biogeo_data$best_model_fit,
#'                                         phylo = eel.tree,
#'                                         BSM_output = eel_biogeo_data$BSM_output)
#' # Explore output
#' str(all_simmaps, max.level = 1)
#' # Plot simmap n°1
#' plot(all_simmaps[[1]]) }
#'


BSM_to_phytools_simmap <- function(model_fit, phylo, BSM_output, sim_index)
{
  if (!requireNamespace("BioGeoBEARS", quietly = TRUE))
  {
    stop("Package 'BioGeoBEARS' is needed for this function to work.
       Please install it manually from: https://github.com/nmatzke/BioGeoBEARS")
  }

  # Extract the tables of cladogenetic and anagenetic events
  clado_events_table <- BSM_output$RES_clado_events_tables[[sim_index]]
  ana_events_table <- BSM_output$RES_ana_events_tables[[sim_index]]

  run_parent_brs_TF = TRUE
  tr <- phylo

  # Error check
  if ((inherits(ana_events_table, 'data.frame')) && (dim(ana_events_table)[1] > 0))
  {
    run_parent_brs_TF = TRUE
  } else {
    ana_events_table = NA
    run_parent_brs_TF = FALSE
  }

  if (is.null(ana_events_table) == TRUE)
  {
    ana_events_table = NA
    run_parent_brs_TF = FALSE
  }

  # Is it time-stratified?
  stratTF = (length(model_fit$inputs$timeperiods) > 0)

  returned_mats = BioGeoBEARS::get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = model_fit$inputs, include_null_range = model_fit$inputs$include_null_range)
  returned_mats
  areanames = returned_mats$areanames
  ranges_list = returned_mats$ranges_list

  # Fix for issue with using only reduced range list
  reduced_ranges_list = returned_mats$ranges_list
  max_range_size <- max(nchar(reduced_ranges_list))
  ranges_list <- generate_list_ranges(areas_list = areanames, max_range_size = max_range_size, include_null_range = model_fit$inputs$include_null_range)

  # Get list of edges:
  trtable = BioGeoBEARS::prt(tr, printflag=FALSE)
  trtable

  # View(trtable)
  # View(clado_events_table)

  # Convert pruningwise edge numbers to cladewise edgenums
  # The clado_events_table (non-stratified) has the edge numbers in
  # pruningwise order
  if (stratTF == FALSE)
  {
    pruningwise_edgenums = clado_events_table$parent_br
    cladewise_edgenums = trtable$parent_br
    translation_pruning_to_clade_edgenums = as.data.frame(cbind(pruningwise_edgenums, cladewise_edgenums), stringsAsFactors=FALSE)
    translation_pruning_to_clade_edgenums

    # Error trap for when there are no anagenetic events
    if (run_parent_brs_TF == TRUE)
    {
      ana_events_edgenums_indexes_in_clado_events_table = match(x = ana_events_table$parent_br, table = clado_events_table$parent_br)
      ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[ana_events_edgenums_indexes_in_clado_events_table]
    }
    clado_events_table$parent_br = trtable$parent_br
    # } else {
    #   # If same issue with time-stratified clado_events_table, also have duplicates due to branches split across time-strata!
    #   pruningwise_edgenums = unique(clado_events_table$parent_br)
    #   cladewise_edgenums = trtable$parent_br
    #   translation_pruning_to_clade_edgenums = as.data.frame(cbind(pruningwise_edgenums, cladewise_edgenums), stringsAsFactors=FALSE)
    #   translation_pruning_to_clade_edgenums
    #   # Error trap for when there are no anagenetic events
    #   if (run_parent_brs_TF == TRUE)
    #   {
    #     # ana_events_edgenums_indexes_in_clado_events_table = match(x = ana_events_table$parent_br, table = clado_events_table$parent_br)
    #     # ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[ana_events_edgenums_indexes_in_clado_events_table]
    #     ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[match(x = ana_events_table$parent_br, table = translation_pruning_to_clade_edgenums$pruningwise_edgenums)]
    #   }
    #   clado_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[match(x = clado_events_table$parent_br, table = translation_pruning_to_clade_edgenums$pruningwise_edgenums)]
  }


  # Get the edgenums, exclude the "NA" for the root branch
  # Order edgenums from smallest to largest
  nonroot_TF = trtable$node.type != "root"
  edgenums = trtable$parent_br[nonroot_TF]
  edgenums_order = order(edgenums)
  edgenums = edgenums[edgenums_order]

  # Get the "ancestor node, descendant node"
  ancnodenums = trtable$ancestor[nonroot_TF][edgenums_order]
  decnodenums = trtable$node[nonroot_TF][edgenums_order]

  # rownames for mapped.edge
  rownames_for_mapped_edge = paste0(ancnodenums, ",", decnodenums)
  rownames_for_mapped_edge

  # instantiate "maps" for phytools (a list, with array of state residence times
  maps = list()

  if (stratTF == TRUE)
  {
    time_tops = sort(unique(clado_events_table$time_top))
    time_bots = sort(unique(clado_events_table$time_bot))
  }


  # Loop through the edges, record any anagenetic events on the branches
  # fill in "maps"
  for (i in 1:length(edgenums))
  {
    # i <- 410
    edgenum = edgenums[i]

    # Trap for if ana_events_table is NA (common, if there are no anagenetic events in
    # the tree at all)
    if ( (length(ana_events_table) == 1) && (is.na(ana_events_table)) )
    {
      edgefound_TF = FALSE
    } else {
      edgefound_TF = ana_events_table$parent_br == edgenum
    }

    # If no anagenetic events are found, the whole branchlength is in the
    # starting state, as specified in "clado_events_table"
    if (sum(edgefound_TF) == 0)
    {
      # The states should be the same at the branch bottom and top:
      clado_row_TF = clado_events_table$parent_br == edgenum
      clado_row_TF[is.na(clado_row_TF)] = FALSE

      # 			if (stratTF == TRUE)
      # 				{
      # 				match_time_tops_TF = as.numeric(tmptable_rows$abs_event_time[nnr]) >= as.numeric(trtable$time_top)
      # 				match_time_bots_TF = as.numeric(tmptable_rows$abs_event_time[nnr]) < as.numeric(trtable$time_bot)
      # 				}


      ## NJM 2019-03-12_ fix: doubles can be found in time-strat, FIX
      clado_events_table[clado_row_TF,]


      # Error check
      if (sum(clado_row_TF) < 1)
      {
        txt = paste0("STOP ERROR #1 in BSM_to_phytools_SM(): 0 rows in clado_events_table match edge number/branch number (parent_br==", edgenum, ").")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }

      if (sum(clado_row_TF) == 1)
      {
        bottom_state_num_1based = clado_events_table$sampled_states_AT_brbots[clado_row_TF]
        top_state_num_1based = clado_events_table$sampled_states_AT_nodes[clado_row_TF]

        # Error handle for the root
        if (clado_events_table$node.type[clado_row_TF] == "root")
        {
          bottom_state_num_1based <- top_state_num_1based
        }
      }

      if (sum(clado_row_TF) > 1)
      {
        bottom_state_num_1based = unique(clado_events_table$sampled_states_AT_brbots[clado_row_TF])
        top_state_num_1based = unique(clado_events_table$sampled_states_AT_nodes[clado_row_TF])

        # Error check: there should be only 1 unique state corresponding to this
        # node (because we are in the section where no anagenetic histories were
        # found on the branch).
        if (length(top_state_num_1based) != 1)
        {
          txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'top_state_num_1based' corresponding to the node specified. Printing the matching rows of 'clado_events_table'."
          cat("\n\n")
          cat(txt)
          cat("\n\n")

          print("clado_events_table[clado_row_TF,]")
          print(clado_events_table[clado_row_TF,])

          print("top_state_num_1based")
          print(top_state_num_1based)

          stop(txt)
        } # END if (length(top_state_num_1based) != 1)

        if (length(bottom_state_num_1based) != 1)
        {
          txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'bottom_state_num_1based' corresponding to the node specified. Printing the matching rows of 'clado_events_table'."
          cat("\n\n")
          cat(txt)
          cat("\n\n")

          print("clado_events_table[clado_row_TF,]")
          print(clado_events_table[clado_row_TF,])

          print("bottom_state_num_1based")
          print(bottom_state_num_1based)

          stop(txt)
        } # END if (length(bottom_state_num_1based) != 1)
      } # END if (sum(clado_row_TF) > 1)

      # Error check
      if (bottom_state_num_1based != top_state_num_1based)
      {
        txt = paste0("STOP ERROR #2 in BSM_to_phytools_SM(): the top_state_num_1based (", top_state_num_1based, "), and bottom_state_num_1based (", bottom_state_num_1based, ") have to match at edge number/branch number (parent_br==", edgenum, "), because no anagenetic events were recorded on this branch.")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }

      # No events detected, so put in the states_array just one state
      # But, since there are NO events on this branch, don't add to it if
      # there is already something there (e.g. don't add the same branch for
      # multiple branch segments)
      names_of_states_array = c(ranges_list[bottom_state_num_1based])
      times_in_each_state_array = unique(c(clado_events_table$edge.length[clado_row_TF]))

      if (length(times_in_each_state_array) != 1)
      {
        txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'times_in_each_state_array' corresponding to the node specified. There should only be one time, because no anagenetic events were detected on this branch. Printing the matching rows of 'clado_events_table'."
        cat("\n\n")
        cat(txt)
        cat("\n\n")

        print("clado_events_table[clado_row_TF,]")
        print(clado_events_table[clado_row_TF,])

        print("times_in_each_state_array")
        print(times_in_each_state_array)

        stop(txt)
      } # END if (length(times_in_each_state_array) != 1)


      names(times_in_each_state_array) = names_of_states_array
      times_in_each_state_array
    } # END if (sum(edgefound_TF) == 0)

    # If some anagenetic events are found, the whole branchlength is in the
    # starting state, as specified in "clado_events_table"
    if (sum(edgefound_TF) > 0)
    {
      # The states will be listed in the ana_events_table
      rows_matching_edgenum_TF = ana_events_table$parent_br == edgenum
      tmp_ana_events_table = ana_events_table[rows_matching_edgenum_TF,]

      # Make sure the tmp_ana_events_table is sorted by REVERSE event_time (along branch)
      # (Do REVERSE, because older events have larger event ages)
      tmp_ana_events_table = tmp_ana_events_table[rev(order(tmp_ana_events_table$abs_event_time)),]
      tmp_ana_events_table

      numevents = sum(rows_matching_edgenum_TF)

      # 1 event, 2 states on branch
      if (numevents == 1)
      {
        first_state_name = c(tmp_ana_events_table$current_rangetxt[1])
        abs_time_bp_at_branch_top = tmp_ana_events_table$time_bp[1]
        abs_time_bp_at_branch_bot = abs_time_bp_at_branch_top + tmp_ana_events_table$edge.length[1]
        first_state_timelength = c(abs_time_bp_at_branch_bot - tmp_ana_events_table$abs_event_time[1])

        # The rest of the branch is the 2nd state
        further_state_name = c(tmp_ana_events_table$new_rangetxt[1])
        further_state_time = c(tmp_ana_events_table$edge.length[1] - first_state_timelength)

        times_in_each_state_array = c(first_state_timelength, further_state_time)
        names_of_states_array = c(first_state_name, further_state_name)
      } # END if (numevents == 1)

      # 2+ events, 3+ states on branch
      if (numevents >= 2)
      {
        first_state_name = c(tmp_ana_events_table$current_rangetxt[1])
        #first_state_time = c(tmp_ana_events_table$event_time[1])
        abs_time_bp_at_branch_top = tmp_ana_events_table$time_bp[1]
        abs_time_bp_at_branch_bot = abs_time_bp_at_branch_top + tmp_ana_events_table$edge.length[1]
        first_state_timelength = c(abs_time_bp_at_branch_bot - tmp_ana_events_table$abs_event_time[1])

        nonfirst_rows = 2:numevents
        nonlast_rows = 1:(numevents-1)

        further_state_names = c(tmp_ana_events_table$new_rangetxt)
        further_state_times = tmp_ana_events_table$abs_event_time[nonlast_rows] - tmp_ana_events_table$abs_event_time[nonfirst_rows]

        # How much of the branch is left
        last_time = c(tmp_ana_events_table$edge.length[numevents] - sum(c(first_state_timelength,further_state_times)))

        times_in_each_state_array = c(first_state_timelength, further_state_times, last_time)
        names_of_states_array = c(first_state_name, further_state_names)
      } # END if (numevents >= 2)

      names(times_in_each_state_array) = names_of_states_array
    } # END if (sum(edgefound_TF) > 0)

    # Store
    maps[[i]] = times_in_each_state_array
  } # END for (i in 1:length(edgenums))

  maps


  # Check that the sums of state residence times add up to the branch lengths
  simmap_times <- unlist(lapply(X=maps, FUN=sum))
  tree_times <- tr$edge.length
  all(simmap_times == tree_times)
  cbind(simmap_times, tree_times)

  # Check if issue with order of edges
  simmap_times <- simmap_times[order(simmap_times)]
  tree_times <- tree_times[order(tree_times)]
  cbind(simmap_times, tree_times)

  # Compute residence times per states
  edge_states_matrix <- plyr::ldply(.data = maps, .fun = base::rbind)
  state_times <- apply(X = edge_states_matrix, MARGIN = 2, FUN = sum, na.rm = T)
  state_times_reordered <- state_times[ranges_list]
  names(state_times_reordered) <- ranges_list
  state_times_reordered[is.na(state_times_reordered)] <- 0

  # Add total
  state_times_df <- as.data.frame(t(state_times_reordered))
  state_times_df$total <- sum(state_times_reordered)

  # Add percentage of times
  state_times_perc <- state_times_df / state_times_df$total * 100
  state_times_df <- rbind(state_times_df, state_times_perc)
  row.names(state_times_df) <- c("raw", "perc")

  # # Remove empty states
  # state_times_df <- state_times_df[, state_times_perc != 0]

  # Make the mapped.edge output
  mapped.edge = matrix(data=0.0, nrow=length(edgenums), ncol=length(ranges_list))
  row.names(mapped.edge) = rownames_for_mapped_edge
  colnames(mapped.edge) = ranges_list
  mapped.edge

  # For each branch,
  # 1. Get the list of observed states
  i = 3
  observed_states = sort(unique(names(maps[[i]])))
  observed_states


  # sapply to get the sum of each
  sapply(X = observed_states, FUN = BioGeoBEARS::get_sum_statetime_on_branch, branch_history_map = maps[[i]])

  # Fill in the table for each branch
  for (i in 1:nrow(mapped.edge))
  {
    observed_states = sort(unique(names(maps[[i]])))
    total_residence_times = sapply(X = observed_states, FUN = BioGeoBEARS::get_sum_statetime_on_branch, branch_history_map = maps[[i]])
    names_observed_states = names(total_residence_times)
    mapped.edge[i,names_observed_states] = unname(total_residence_times)
  }
  mapped.edge


  # Get the transition matrix and logL
  Q = returned_mats$Qmat
  # row.names(Q) = ranges_list
  # colnames(Q) = ranges_list
  row.names(Q) = reduced_ranges_list
  colnames(Q) = reduced_ranges_list
  logL = model_fit$total_loglikelihood

  tr_wSimmap = tr
  tr_wSimmap$maps = maps
  tr_wSimmap$mapped.edge = mapped.edge
  tr_wSimmap$Q = Q
  tr_wSimmap$logL = logL
  class(tr_wSimmap) = c("simmap", "phylo")
  tr_wSimmap

  tr_wSimmap$maps
  tr_wSimmap$mapped.edge

  return(list(simmap = tr_wSimmap, residence_times = state_times_df))
}

### Same but for list of multiple BSM maps

#'
#' @rdname BSM_to_phytools_simmap
#'

BSMs_to_phytools_simmaps <- function(model_fit, phylo, BSM_output)
{
  # Initiate final output
  simmaps_list = list()

  for (i in 1:length(BSM_output$RES_clado_events_tables))
  {
    simmaps_list[[i]] = BSM_to_phytools_simmap(model_fit = model_fit, phylo = phylo, BSM_output = BSM_output, sim_index = i)$simmap
  }

  class(simmaps_list) = c("multiSimmap", "multiPhylo")
  return(simmaps_list)
}


### Helper function to split states in densityMaps ####

split_state_densityMaps <- function (densityMaps, state_to_split, states_to_blend_in)
{
  # Initiate output
  densityMaps_updated <- densityMaps

  # Extract all states/ranges
  all_states <- unname(unlist(lapply(X = densityMaps_updated, FUN = function (x) { x$states[2] })))
  # Identify state to split
  state_to_split_ID <- which(all_states == state_to_split)
  # Identify state to blend in
  states_to_blend_in_ID <- which(all_states %in% states_to_blend_in)

  # Get number of states to blend in
  nb_states_to_blend_in <- length(states_to_blend_in)

  # Get mappings from state to split
  state_to_split_maps <- densityMaps_updated[[state_to_split_ID]]$tree$maps

  # Split posterior probabilities by number of states to blend in
  state_to_split_maps_split <- lapply(X = state_to_split_maps, FUN = function (x) { names(x) <- as.character(round(as.numeric(names(x))/nb_states_to_blend_in, 0)) ; return(x) })

  # Add posterior probabilities to states to blend in
  for (i in seq_along(states_to_blend_in))
  {
    # i <- 1

    # Extract state to blend in
    state_to_blend_in_i <- states_to_blend_in[i]
    state_to_blend_in_ID_i <- states_to_blend_in_ID[i]

    # Get mappings from state to split
    state_to_blend_in_maps_i <- densityMaps_updated[[state_to_blend_in_ID_i]]$tree$maps

    # Add posterior probabilities to $maps
    state_to_blend_in_maps_i <- purrr::map2(.x = state_to_blend_in_maps_i, .y = state_to_split_maps_split,
        .f = function (x, y) { names(x) <- as.character(as.numeric(names(x)) + as.numeric(names(y))) ; return(x) })
    # Avoid trespassing 1000 due to rounding issues
    state_to_blend_in_maps_i <- lapply(X = state_to_blend_in_maps_i, FUN = function (x) { names(x)[as.numeric(names(x)) > 1000] <- as.character(1000) ; return(x) })
    # Update $maps
    densityMaps_updated[[state_to_blend_in_ID_i]]$tree$maps <- state_to_blend_in_maps_i

    # Update posterior probabilities in $mapped.edge
    updated_mapped.edge <- makeMappedEdge(edge = densityMaps_updated[[state_to_blend_in_ID_i]]$tree$edge,
        maps = densityMaps_updated[[state_to_blend_in_ID_i]]$tree$maps)
    updated_mapped.edge <- updated_mapped.edge[, colnames(updated_mapped.edge)[order(as.numeric(colnames(updated_mapped.edge)))]]
    densityMaps_updated[[state_to_blend_in_ID_i]]$tree$mapped.edge <- updated_mapped.edge
  }

  # Remove split states
  densityMaps_updated <- densityMaps_updated[-state_to_split_ID]

  # Plot initial and updated densityMaps
  # plot_densityMaps_overlay(densityMaps)
  # plot_densityMaps_overlay(densityMaps_updated)

  # # Check that sum of PP still equals roughly 1000
  # test1 <- lapply(X = densityMaps[[1]]$tree$maps, FUN = function (x) { as.numeric(names(x))})
  # test2 <- lapply(X = densityMaps[[2]]$tree$maps, FUN = function (x) { as.numeric(names(x))})
  # test3 <- purrr::map2(.x = test1, .y = test2,
  #                      .f = function (x, y) { x + y })
  # test3

  ## Export updated densityMaps
  return(densityMaps_updated)
}

### Utility function from phytools used to produce $mapped.edge in densityMap ####
# Original function written by Liam Revell, 2012

makeMappedEdge <- function (edge, maps)
{
  st <- sort(unique(unlist(sapply(maps, function(x) names(x)))))
  mapped.edge <- matrix(0, nrow(edge), length(st))
  rownames(mapped.edge) <- apply(edge, 1, function(x) paste(x, collapse = ","))
  colnames(mapped.edge) <- st

  for (i in 1:length(maps))
  {
    for (j in 1:length(maps[[i]]))
    {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }

  return(mapped.edge)
}

### Utility function to check that an object is a valid color according to grDevices ####

is_color <- function(x, null_ok = FALSE)
{
  if (is.null(x) && null_ok)
  {
    return(TRUE)
  }
  vapply(x, function(i)
  {
    tryCatch({
      is.matrix(grDevices::col2rgb(i))
    }, error = function(e) FALSE)
  }, TRUE)
}

### Utility function to compute means of a series of angles

mean_angle <- function(a, unit = c("degrees", "radians"))
{
  a <- a[!is.na(a)]

  deg2rad <- function(x) { x * pi/180}
  rad2deg <- function(x) { x * 180/pi }

  deg2vec <- function(x, unit = c("degrees", "radians"))
  {
    if(unit == "degrees")
    {
      a <- c(sin(deg2rad(x)), cos(deg2rad(x)))
    } else if (unit == "radians") {
      a <- c(sin(x), cos(x))
    }
    return(a)
  }

  vec2deg <- function(x)
  {
    res <- rad2deg(atan2(x[1], x[2]))
    if (res < 0) { res <- 360 + res }
    return(res)
  }

  mean_vec <- function(x)
  {
    y <- lapply(x, deg2vec, unit = unit)
    Reduce(`+`, y)/length(y)
  }
  return( vec2deg(mean_vec(a)) )
}


# ### Utility functions for tree transformation ####
#
# Used to copy the S3 method rescale.phylo() from geiger
#
# # Wrapper function for tree transformation from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# rescale.phylo <- function(phy,
#                           model = c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"),
#                           ...)
# {
#
#   model = match.arg(model, c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda",
#                              "kappa", "delta", "white", "depth"))
#
#   if (!("phylo" %in% class(phy))) stop("supply 'phy' as a 'phylo' object")
#
#   FUN = switch(model,
#                BM = .bm.phylo(phy),
#                OU = .ou.phylo(phy),
#                EB = .eb.phylo(phy),
#                nrate = .nrate.phylo(phy),
#                lrate = .lrate.phylo(phy),
#                trend = .trend.phylo(phy),
#                lambda = .lambda.phylo(phy),
#                kappa = .kappa.phylo(phy),
#                delta = .delta.phylo(phy),
#                white = .white.phylo(phy),
#                depth = .depth.phylo(phy)
#   )
#   class(FUN) = c("transformer", "function")
#   dots = list(...)
#   if(!missing(...))
#   {
#     if(!all(names(dots) %in% argn(FUN)))
#       stop(paste("The following parameters are expected:\n\t", paste(argn(FUN),
#                                                                      collapse = "\n\t", sep = ""), sep = ""))
#     return(FUN(...))
#   } else {
#     return(FUN)
#   }
# }
#
# # Tree transformation for BM from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .bm.phylo = function(phy)
# {
#   el = phy$edge.length
#   z = function(sigsq)
#   {
#     phy$edge.length = el * sigsq
#     phy
#   }
#   attr(z, "argn") = "sigsq"
#   return(z)
# }
#
# # Tree transformation for OU from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .ou.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t1 = Tmax - ht$end[phy$edge[mm,1]]
#   ht$t2 = ht$start - ht$end + ht$t1
#   z = function(alpha, sigsq=1)
#   {
#     if(alpha<0) stop("'alpha' must be positive valued")
#     bl = (1/(2*alpha)) * exp(-2*alpha * (Tmax-ht$t2)) * (1 - exp(-2 * alpha * ht$t2)) - (1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t1)) * (1 - exp(-2 * alpha * ht$t1))
#     phy$edge.length = bl[phy$edge[,2]]
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("alpha", "sigsq")
#   return(z)
# }
#
# # Tree transformation for EB from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .eb.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t1 = Tmax - ht$end[phy$edge[mm,1]]
#   ht$t2 = ht$start - ht$end + ht$t1
#
#   z = function(a, sigsq = 1)
#   {
#     if(a == 0) return(phy)
#     bl = (exp(a*ht$t2) - exp(a*ht$t1)) / (a)
#     phy$edge.length = bl[phy$edge[,2]]
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("a", "sigsq")
#   return(z)
# }
#
# # Tree transformation for multi-regime BM with rates split by time from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .nrate.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t = Tmax - ht$end
#   ht$e = ht$start - ht$end
#   ht$a = ht$t - ht$e
#   ht$rS = ht$a/Tmax
#   ht$rE = ht$t/Tmax
#
#   dd = phy$edge[,2]
#
#   relscale.brlen = function(start, end, len, dat)
#   {
#     ss = start <= dat[,"time"]
#     strt = min(which(ss))
#
#     ee = dat[,"time"] < end
#     etrt = max(which(ee)) + 1
#
#     bl = numeric()
#     fragment = numeric()
#     marker = start
#     for(i in strt:etrt)
#     {
#       fragment = c(fragment, (nm <- (min(c(end, dat[i, "time"])))) - marker)
#       bl = c(bl,dat[i, "rate"])
#       marker = nm
#     }
#     fragment = fragment/(sum(fragment))
#     sclbrlen = numeric()
#     for (i in 1:length(bl)) sclbrlen = c(sclbrlen, len * fragment[i] * bl[i])
#     sc = structure(as.numeric(sclbrlen), names = strt:etrt)
#     return(sc)
#   }
#
#
#   z = function(time, rate, sigsq=1)
#   {
#     if(any(time > 1) | any(time < 0)) stop("supply 'time' as a vector of relative time:\n\tvalues should be in the range 0 (root) to 1 (present)")
#     if(any(rate < 0)) stop("'rate' must consist of positive values")
#     if(length(time) != length(rate)) stop("'time' and 'rate' must be of equal length")
#     phy$edge.length = phy$edge.length * sigsq
#     ordx = order(time)
#     time = time[ordx]
#     rate = rate[ordx]
#     dat = cbind(time=c(0,time, 1), rate = (c(1, 1, rate)))
#     rs = sapply(dd, function(x) as.numeric(sum(relscale.brlen(ht$rS[x], ht$rE[x], ht$e[x], dat))))
#     phy$edge.length = rs
#     phy
#   }
#   attr(z, "argn") = c("time", "rate", "sigsq")
#   return(z)
# }
#
# # Tree transformation for multi-regime BM with rates split by clades from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .lrate.phylo = function(phy)
# {
#   N = Ntip(phy)
#   n = Nnode(phy)
#   cache = .cache.tree(phy)
#   cache$phy = phy
#   vv = rep(1, N+n-1)
#
#   z = function(node, rate, sigsq = 1)
#   {
#     shifts = c(sort(node[node > N]), node[node <= N])
#     mm = match(shifts, node)
#     rates = rate[mm]
#
#     phy$edge.length = phy$edge.length * sigsq
#     for(i in 1:length(shifts)) vv = .assigndescendants(vv, shifts[i], rates[i], exclude = shifts, cache = cache)
#     phy$edge.length = vv * phy$edge.length
#     phy
#   }
#
#   attr(z, "argn") = c("node", "rate", "sigsq")
#   return(z)
# }
#
# # Tree transformation for "rate_trend" model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .trend.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$head = Tmax - ht$end[phy$edge[mm,1]] # age
#   ht$tail = ht$head + (ht$start - ht$end)
#
#   z = function(slope, sigsq = 1)
#   {
#     # begin (age): head
#     # end: tail
#     ht$br = 1 + ht$head * slope
#     ht$er = 1 + ht$tail * slope
#     scl = sapply(1:nrow(ht),
#                  function(idx)
#                  {
#                    if(idx == N+1) return(NA)
#                    if(ht$br[idx] > 0 & ht$er[idx] > 0)
#                    {
#                      return((ht$br[idx]+ht$er[idx])/2)
#                    } else if (ht$br[idx] < 0 & ht$er[idx] < 0) {
#                      return(0)
#                    } else {
#                      si = -1/slope
#                      return(ht$br[idx] * (si-ht$head[idx]) / (2 * (ht$tail[idx] - ht$head[idx])))
#                    }
#                  })
#
#     phy$edge.length = phy$edge.length * scl[phy$edge[,2]]
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("slope", "sigsq")
#   return(z)
# }
#
# # Tree transformation for Pagel's lambda model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .lambda.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:N, phy$edge[,2])
#   ht$e = ht$start - ht$end
#   paths = .paths.phylo(phy)
#
#   z = function(lambda, sigsq = 1)
#   {
#     if(lambda < 0) stop("'lambda' must be positive valued")
#
#     bl = phy$edge.length * lambda
#     bl[mm] = bl[mm] + (paths - (paths * lambda))
#     phy$edge.length = bl
#     if(any(phy$edge.length < 0))
#     {
#       warning("negative branch lengths encountered:\n\tlambda may be too large")
#     }
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("lambda", "sigsq")
#   return(z)
# }
#
# # Tree transformation for Pagel's kappa model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .kappa.phylo = function(phy)
# {
#   z = function(kappa, sigsq = 1)
#   {
#     if(kappa < 0) stop("'kappa' must be positive valued")
#
#     phy$edge.length = phy$edge.length^kappa
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("kappa", "sigsq")
#   return(z)
# }
#
# # Tree transformation for Pagel's delta model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .delta.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t = Tmax - ht$end
#   ht$e = ht$start - ht$end
#   ht$a = ht$t - ht$e
#   if(sum(ht$a < -0.1) > 0) stop("Calculation error; contact developers (LJ Harmon and JM Eastman).")
#   ht$a[ht$a < 0] <- 0
#
#   z = function(delta, sigsq = 1, rescale=TRUE)
#   {
#     if(delta < 0) stop("'delta' must be positive valued")
#     bl = (ht$a + ht$e)^delta - ht$a^delta
#     phy$edge.length = bl[phy$edge[,2]]
#
#     if(rescale)
#     {
#       scl = Tmax^delta
#       phy$edge.length = (phy$edge.length/scl) * Tmax
#     }
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("delta", "sigsq")
#   return(z)
# }
#
# # Tree transformation for a non-phylogenetic "white" model from R package geiger
# # that assumes there is no covariance structure among species
# # Internal branches are set to a length of zero
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .white.phylo = function(phy)
# {
#   N = Ntip(phy)
#   phy$edge.length[] = 0
#   phy$edge.length[phy$edge[,2] <= N] = 1
#   el = phy$edge.length
#   z = function(sigsq = 1)
#   {
#     phy$edge.length = el * sigsq
#     phy
#   }
#   attr(z, "argn") = "sigsq"
#   return(z)
# }
#
# # Tree transformation "depth" model from R package geiger
# # Adjusts the total depth of the tree
# # which in return increase/decrease rates of evolution
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .depth.phylo = function(phy)
# {
#   orig = max(heights.phylo(phy))
#   z = function(depth)
#   {
#     phy$edge.length <- (phy$edge.length/orig) * depth
#     if(!is.null(phy$root.edge)) phy$root.edge = (phy$root.edge/orig) * depth
#     phy
#   }
#   attr(z, "argn") = "depth"
#   z
# }
#
# # Internal function from R package geiger
# # Compute path length from root to tip
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .paths.phylo = function(phy, ...)
# {
#
#   ## much from ape:::vcv.phylo()
#   phy <- reorder(phy, "postorder")
#
#   FUN = function(vcv = FALSE)
#   {
#     n <- length(phy$tip.label)
#     pp <- .cache.descendants(phy)$tips
#     e1 <- phy$edge[, 1]
#     e2 <- phy$edge[, 2]
#     EL <- phy$edge.length
#     xx <- numeric(n + phy$Nnode)
#     if(vcv) vmat = matrix(0, n, n)
#     for (i in length(e1):1)
#     {
#       var.cur.node <- xx[e1[i]]
#       xx[e2[i]] <- var.cur.node + EL[i]
#       if(vcv)
#       {
#         j <- i - 1L
#         while (e1[j] == e1[i] && j > 0)
#         {
#           left = pp[[e2[j]]]
#           right = pp[[e2[i]]]
#           vmat[left, right] <- vmat[right, left] <- var.cur.node
#           j <- j - 1L
#         }
#       }
#     }
#     if(vcv)
#     {
#       diags <- 1 + 0:(n - 1) * (n + 1)
#       vmat[diags] <- xx[1:n]
#       colnames(vmat) <- rownames(vmat) <- phy$tip.label
#       return(vmat)
#     } else {
#       return(xx[1:n])
#     }
#   }
#
#   FUN(...)
# }
#
# ## Internal function from R package geiger
# # Associated C scripts in /src/utilities.cpp
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
#
# # Declare the use of compiled C code in the package
# # Add an import in NAMESPACE
# # Need to add things listed in Rcpp::Rcpp.package.skeleton() to get the necessary changes for the package to include Rcpp
# # Need to add CallEntries in init.c such as R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
# #' @useDynLib deepSTRAPP
# #' @importFrom Rcpp evalCpp
#
# .cache.descendants = function(phy)
# {
#   # Fetches all tips subtended by each internal node
#
#   N = as.integer(Ntip(phy))
#   n = as.integer(Nnode(phy))
#
#   phy = reorder(phy, "postorder")
#
#   zz = list(N=N,
#             MAXNODE=N+n,
#             ANC=as.integer(phy$edge[,1]),
#             DES=as.integer(phy$edge[,2])
#   )
#
#   res = .Call("cache_descendants", phy = zz, PACKAGE = "geiger")
#   return(res)
# }
#
#
# # Internal function from R package geiger
# # Compute heights of nodes in phylo
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# heights.phylo = function(x)
# {
#   phy = x
#   phy <- reorder(phy, "postorder")
#   n <- length(phy$tip.label)
#   n.node <- phy$Nnode
#   xx <- numeric(n + n.node) # ending times
#   for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
#   root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
#   labs = c(phy$tip.label, phy$node.label)
#   depth = max(xx)
#   tt = depth - xx # time to 'present day' of branch starts
#   idx = 1:length(tt)
#   #dd = phy$edge.length[idx]
#   mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
#   dd = c(phy$edge.length, root)[mm] # reordered bls
#   ss = tt + dd
#   res = cbind(ss, tt)
#   rownames(res) = idx
#   colnames(res) = c("start", "end")
#   res = data.frame(res)
#   res
# }

### Dummy function to force the use of Rcpp ####
#' @useDynLib deepSTRAPP
#' @importFrom Rcpp evalCpp
dummy_Rcpp <- function ()
{
  Rcpp::evalCpp(code = "__cplusplus" )
}
