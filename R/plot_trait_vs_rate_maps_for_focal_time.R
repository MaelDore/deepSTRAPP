
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
#'   * For continuous traits: [phytools::plot.contMap()]
#'   * For categorical and biogeographic data: [deepSTRAPP::plot_densityMaps_overlay()]
#'   * For BAMM rates and regme shifts: [deepSTRAPP::plot_BAMM_rates()]
#'
#' #' @param STRAPP_results List of elements generated with [deepSTRAPP::compute_STRAPP_test_for_focal_time()],
#'   that summarize the results of a STRAPP test for a specific time in the past (i.e. the `focal_time`).
#'   `STRAPP_results` can also be extracted from the output of [deepSTRAPP::run_deepSTRAPP_over_time()] that
#'   run the whole deepSTRAPP workflow and store results inside `$STRAPP_results`.
#'
#' @param rate_type A character string specifying the type of diversification rates to plot.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param method A character string indicating the method for plotting the phylogenetic tree.
#'   * `method = "phylogram"` (default) plots the phylogenetic tree using rectangular coordinates.
#'   * `method = "polar"` plots the phylogenetic tree using polar coordinates.
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
#'
#'
#' @param ... Additional graphical arguments to pass down to [phytools::plot.contMap()], [deepSTRAPP::plot_densityMaps_overlay()],
#'   [BAMMtools::plot.bammdata()], [BAMMtools::addBAMMshifts()], and [par()].
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with '.pdf'.
#'
#' @export
#' @importFrom grDevices pdf dev.off gray
#' @importFrom BAMMtools plot.bammdata addBAMMshifts
#'
#' @details The main input `BAMM_object` is the typical output of [deepSTRAPP::prepare_diversification_data()].
#'   It provides information on rates and regimes shifts across the posterior samples of a BAMM.
#'
#'   `$MAP_BAMM_object` and `$MSC_BAMM_object` elements are required to plot regime shift locations following the
#'   "MAP" or "MSC" `configuration_type` respectively.
#'   A `$MSP_tree` element is required to scale the size of the symbols showing the location of regime shifts according marginal shift probabilities.
#'   (If `adjust_size_to_prob = TRUE`).
#'
#'   The default option to display regime shift is to use the average locations from the posterior samples with the Maximum A Posteriori probability (MAP) configuration.
#'   However, sometimes, multiple configurations have similarly high frequency in the posterior samples (See [BAMMtools::credibleShiftSet()].
#'   An alternative is to use the average locations from posterior samples with the Maximum Shift Credibility (MSC) configuration instead.
#'   This regime shift configuration has the highest product of marginal probabilities across branches where a shift is estimated.
#'   It may differ from the MAP configuration. (See [BAMMtools::maximumShiftCredibility()]).
#'
#' @return The function returns (invisibly) a list with three three elements similarly to [BAMMtools::plot.bammdata()].
#'  * `$coords`: A matrix of plot coordinates. Rows correspond to branches. Columns 1-2 are starting (x,y) coordinates of each branch and columns 3-4 are ending (x,y) coordinates of each branch. If method = "polar" a fifth column gives the angle(in radians) of each branch.
#'  * `$colorbreaks`: A vector of percentiles used to group macroevolutionary rates into color bins.
#'  * `$colordens`: A matrix of the kernel density estimates (column 2) of evolutionary rates (column 1) and the color (column 3) corresponding to each rate value.
#'
#' @author Maël Doré
#'
#' @seealso Initial functions in BAMMtools: [BAMMtools::plot.bammdata()] [BAMMtools::addBAMMshifts()]
#'
#' Associated functions in deepSTRAPP: [deepSTRAPP::prepare_diversification_data()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()] [deepSTRAPP::run_deepSTRAPP_for_focal_time] [deepSTRAPP::run_deepSTRAPP_over_time()]
#'
#' @examples
#' # Load BAMM output
#' data(whale_BAMM_object, package = "deepSTRAPP")
#'
#' ## Plot overall mean rates with MAP configuration for regime shifts
#' # (rates are averaged only all posterior samples)
#' plot_BAMM_rates(whale_BAMM_object, add_regime_shifts = TRUE,
#'                 configuration_type = "MAP", bg = "black",
#'                 regimes_size = 3)
#' ## Plot overall mean rates with MSC configuration for regime shifts
#' # (rates are averaged only all posterior samples)
#' plot_BAMM_rates(whale_BAMM_object, add_regime_shifts = TRUE,
#'                 configuration_type = "MSC", bg = "black",
#'                 regimes_size = 3)
#'
#' ## Plot mean MAP rates with regime shifts
#' # (rates are averaged only across MAP samples)
#' plot_BAMM_rates(whale_BAMM_object$MAP_BAMM_object, add_regime_shifts = TRUE,
#'                 configuration_type = "index",
#'                 # Set to index to use the regime shift data from the '$MAP_BAMM_object'
#'                 regimes_size = 3, bg = "black")
#' ## Plot mean MSC rates (rates averaged only across MSC samples) with regime shifts
#' # (rates averaged only across MSC samples)
#' plot_BAMM_rates(whale_BAMM_object$MSC_BAMM_object, add_regime_shifts = TRUE,
#'                 configuration_type = "index",
#'                 # Set to index to use the regime shift data from the '$MSC_BAMM_object'
#'                 regimes_size = 3, bg = "black")
#'
