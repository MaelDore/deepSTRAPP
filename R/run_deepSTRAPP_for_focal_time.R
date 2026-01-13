
#' @title Run deepSTRAPP to test for a relationship between diversification rates and trait data at a given focal time
#'
#' @description Wrapper function to run deepSTRAPP workflow for a given point in the past (i.e. the `focal_time`).
#'   It starts from traits mapped on a phylogeny (trait data) and BAMM output (diversification data)
#'   and carries out the appropriate statistical method to test for a relationship between diversification rates and trait data.
#'   Tests are based on block-permutations: rates data are randomized across tips following blocks
#'   defined by the diversification regimes identified on each tip (typically from a BAMM).
#'
#'   Such tests are called STructured RAte Permutations on Phylogenies (STRAPP) as described in
#'   Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification.
#'   Systematic biology, 65(2), 181-193. \doi{10.1093/sysbio/syv066}.
#'
#'   See the original [BAMMtools::traitDependentBAMM()] function used to
#'   carry out STRAPP test on extant time-calibrated phylogenies.
#'
#'   Tests can be carried out on speciation, extinction and net diversification rates.
#'
#' @param contMap For continuous trait data. Object of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()] or [phytools::contMap()],
#'   that contains a phylogenetic tree and associated continuous trait mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param densityMaps For categorical trait or biogeographic data. List of objects of class `"densityMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state/range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state/range. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param ace (Optional) Ancestral Character Estimates (ACE) at the internal nodes.
#'   Obtained with [deepSTRAPP::prepare_trait_data()] as output in the `$ace` slot.
#'   * For continuous trait data: Named numerical vector typically generated with [phytools::fastAnc()], [phytools::anc.ML()], or [ape::ace()].
#'     Names are nodes_ID of the internal nodes. Values are ACE of the trait.
#'   * For categorical trait or biogeographic data: Matrix that record the posterior probabilities of ancestral states/ranges.
#'     Rows are internal nodes_ID. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   Needed in all cases to provide accurate estimates of trait values.
#' @param tip_data (Optional) Named vector of tip values of the trait.
#'   * For continuous trait data: Named numerical vector of trait values.
#'   * For categorical trait or biogeographic data: Character string vector of states/ranges
#'   Names are nodes_ID of the internal nodes. Needed to provide accurate tip values.
#'   * For biogeographic data, ranges should follow the coding scheme of BioGeoBEARS with a unique CAPITAL letter per unique areas
#'   (ex: A, B), combined to form multi-area ranges (Ex: AB). Alternatively, you can provide tip_data as a matrix or data.frame of
#'   binary presence/absence in each area (coded as unique CAPITAL letter). In this case, columns are unique areas, rows are taxa,
#'   and values are integer (0/1) signaling absence or presence of the taxa in the area.
#' @param trait_data_type Character string. Specify the type of trait data. Must be one of "continuous", "categorical", "biogeographic".
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::prepare_diversification_data()],
#'   that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#'   The phylogenetic tree must the same as the one associated with the `contMap`/`densityMaps`, `ace` and `tip_data`.
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   at which data must be extracted and the phylogeny and mappings must be cut.
#'   It must be smaller than the root age of the phylogeny.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated phylogeny. Default is `TRUE`.
#' @param rate_type A character string specifying the type of diversification rates to use. Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param seed Integer. Set the seed to ensure reproducibility. Default is `NULL` (a random seed is used).
#' @param nb_permutations Integer. To select the number of random permutations to perform during the tests.
#'   If NULL (default), all posterior samples will be used once.
#' @param replace_samples Logical. To specify whether to allow 'replacement' (i.e., multiple use) of a posterior sample
#'   when drawing samples used to carry out the STRAPP test. Default is `FALSE`.
#' @param alpha Numerical. Significance level to use to compute the `estimate` corresponding to the values of the test statistic used to assess significance of the test.
#'   This does NOT affect p-values. Default is `0.05`.
#' @param two_tailed Logical. To define the type of tests. If `TRUE` (default), tests for correlations/differences in rates will be carried out with a null hypothesis
#'   that rates are not correlated with trait values (continuous data) or equals between trait states (categorical and biogeographic data).
#'   If `FALSE`, one-tailed tests are carried out.
#'   * For continuous data, it involves defining a `one_tailed_hypothesis` testing for
#'   either a "positive" or "negative" correlation under the alternative hypothesis.
#'   * For binary data (two states), it involves defining a `one_tailed_hypothesis` indicating which states have higher rates under the alternative hypothesis.
#'   * For multinominal data (more than two states), it defines the type of post hoc pairwise tests to carry out between pairs of states.
#'     If `posthoc_pairwise_tests = TRUE`, all two-tailed (if `two_tailed = TRUE`) or one-tailed (if `two_tailed = FALSE`) tests are automatically carried out.
#' @param one_tailed_hypothesis A character string specifying the alternative hypothesis in the one-tailed test.
#'   For continuous data, it is either "negative" or "positive" correlation.
#'   For binary data, it lists the trait states with states ordered in increasing rates under the alternative hypothesis, separated by a greater-than such as c('A > B').
#' @param posthoc_pairwise_tests Logical. Only for multinominal data (with more than two states). If `TRUE`, all possible post hoc pairwise (Dunn) tests
#'   will be computed across all pairs of states. This is a way to detect which pairs of states have significant differences in rates
#'   if the overall test (Kruskal-Wallis) is significant. Default is `FALSE`.
#' @param p.adjust_method A character string. Only for multinominal data (with more than two states). It specifies the type of correction to apply to the p-values
#'  in the post hoc pairwise tests to account for multiple comparisons. See [stats::p.adjust()] for the available methods. Default is `none`.
#' @param return_perm_data Logical. Whether to return the stats data computed from the posterior samples for observed and permuted data in the output.
#'  This is needed to plot the histogram of the null distribution used to assess significance of the test with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'  Default is `FALSE`.
#' @param nthreads Integer. Number of threads to use for paralleled computing of the STRAPP tests across the permutations.
#'  The R package `parallel` must be loaded for `nthreads > 1`. Default is `1`.
#' @param print_hypothesis Logical. Whether to print information on what test is carried out, detailing the null and alternative hypotheses,
#'  and what significant level is used to rejected or not the null hypothesis. Default is `TRUE`.
#' @param extract_diversification_data_melted_df Logical. Specify whether diversification data (regimes ID and tip rates) must be extracted from the `updated_BAMM_object`
#'   and returned in a melted data.frame. Default is `FALSE`.
#' @param return_updated_trait_data_with_Map Logical. Specify whether the `trait_data` extracted
#'   for the given `focal_time` and the updated version of mapped phylogeny (`contMap`/`densityMaps`) provided as input
#'   should be returned among the outputs. The updated `contMap`/`densityMaps` consists in cutting off branches and mapping
#'   that are younger than the `focal_time`. Default is `FALSE`.
#' @param return_updated_BAMM_object Logical. Specify whether the `updated_BAMM_object` with phylogeny and
#'   mapped diversification rates cut-off at the `focal_time` should be returned among the outputs.
#' @param verbose Logical. Should progression be displayed? A message will be printed at each stepof the deepSTRAPP workflow,
#'   and for every batch of 100 BAMM posterior samples whose rates are regimes are updated, and optionally extracted in a melted data.frame
#'   (if `extract_diversification_data_melted_df = TRUE`). Default is `TRUE`.
#'
#' @export
#'
#' @details The function encapsulates several functions carrying out each step of the deepSTRAPP workflow:
#'
#'   ## Extract trait data
#'
#'   [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] extracts the most likely trait values
#'   found along branches at the `focal_time`.
#'   Optionally, the function can update the mapped phylogeny (`contMap`/`densityMaps`) such as
#'   branches overlapping the `focal_time` are shorten to the `focal_time`, and
#'   the trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#'   ## Extract diversification data
#'
#'   [deepSTRAPP::update_rates_and_regimes_for_focal_time()] updates the `BAMM_object` to obtain
#'   the diversification rates/regimes found along branches the `focal_time`.
#'   Optionally, the function can update the `BAMM_object` to display a mapped phylogeny
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`
#'
#'   ## Extract diversification data in a melted df
#'
#'   If requested (`extract_diversification_data_melted_df = TRUE`), [deepSTRAPP::extract_diversification_data_melted_df_for_focal_time()]
#'   will be used to extract regimes ID and tip rates from the `updated_BAMM_object` and provide a melted data.frame summarizing the diversification data
#'   as found on the phylogeny for the `focal_time`.
#'
#'   ## Compute STRAPP test
#'
#'   [deepSTRAPP::compute_STRAPP_test_for_focal_time()] carries out the appropriate statistical method to test for
#'   a relationship between diversification rates and trait data for a given point in the past (i.e. the `focal_time`).
#'   It can handle three types of statistical tests depending on the type of trait data provided:
#'   * Continuous trait data: Test for correlations with the Spearman's rank correlation test (See [stats::cor.test]).
#'   * Binary trait data (two states): Test for differences in rates between states with the Mann-Whitney-Wilcoxon rank-sum test (See [stats::wilcox.test]).
#'   * Multinominal trait data (More than two states): Test for differences in rates across all states with the Kruskal-Wallis H test (See [stats::kruskal.test]).
#'     If `posthoc_pairwise_tests = TRUE`, Dunn's post hoc pairwise rank-sum tests between pairs of states will be carried out too (See [dunn.test::dunn.test]).
#'
#' @return The function returns a list with at least two elements.
#'
#'   * `$STRAPP_results` List with at least eight elements summarizing the results of the STRAPP tests.
#'     See [deepSTRAPP::compute_STRAPP_test_for_focal_time()] for a detailed description of the output.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the data were extracted and the STRAPP test carried out.
#'
#'   Optional formatted output:
#'   * `$diversification_data_df` A data.frame with six columns summarizing the diversification data as found on the phylogeny for the `focal_time`.
#'     See [extract_diversification_data_melted_df_for_focal_time()] for a detailed description of the output.
#'
#'   Optional data updated for the `focal_time`:
#'   * `$updated_trait_data_with_Map` A list with four elements that contains trait data found at the `focal_time` and an updated `contMap` or `densityMaps`
#'     that can be used as input of [deepSTRAPP::plot_contMap()] or [deepSTRAPP::plot_densityMaps_overlay()] to display a phylogeny mapped with trait values/states/ranges with branches cut at the `focal_time`.
#'     See [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] for a detailed description of the output.
#'   * `$updated_BAMM_object` An updated `BAMM_object` of class `"bammdata"` that contains rates and regimes ID found at the `focal_time`.
#'     Can be used as input of [deepSTRAPP::plot_BAMM_rates()] to display a phylogeny mapped with diversification rates with branches cut at the `focal_time`.
#'     See [deepSTRAPP::update_rates_and_regimes_for_focal_time()] for a detailed description of the output.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#' [deepSTRAPP::extract_diversification_data_melted_df_for_focal_time()] [deepSTRAPP::compute_STRAPP_test_for_focal_time()]
#'
#' For a guided tutorial on complete deepSTRAPP workflow, see the associated vignettes:
#' * For continuous trait data: \code{vignette("deepSTRAPP_continuous_data", package = "deepSTRAPP")}
#' * For categorical trait data: \code{vignette("deepSTRAPP_categorical_3lvl_data", package = "deepSTRAPP")}
#' * For biogeographic range data: \code{vignette("deepSTRAPP_biogeographic_data", package = "deepSTRAPP")}
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  # ----- Example 1: Continuous trait ----- #
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
#'  \donttest{ # (May take several minutes to run)
#'  # Get Ancestral Character Estimates based on a Brownian Motion model
#'  # To obtain values at internal nodes
#'  Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)
#'
#'  # Run a Stochastic Mapping based on a Brownian Motion model
#'  # to interpolate values along branches and obtain a "contMap" object
#'  Ponerinae_contMap <- phytools::contMap(Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data,
#'                                         res = 100, # Number of time steps
#'                                         plot = FALSE)
#'  # Plot contMap = stochastic mapping of continuous trait
#'  plot_contMap(contMap = Ponerinae_contMap,
#'              color_scale = color_scale)
#'
#'  ## Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#'  deepSTRAPP_output <- run_deepSTRAPP_for_focal_time(
#'    contMap = Ponerinae_contMap,
#'    ace = Ponerinae_ACE,
#'    tip_data = Ponerinae_cont_tip_data,
#'    trait_data_type = "continuous",
#'    BAMM_object = Ponerinae_BAMM_object_old_calib,
#'    focal_time = focal_time,
#'    rate_type = "net_diversification",
#'    return_perm_data = TRUE,
#'    extract_diversification_data_melted_df = TRUE,
#'    return_updated_trait_data_with_Map = TRUE,
#'    return_updated_BAMM_object = TRUE)
#'
#'  ## Explore output
#'  str(deepSTRAPP_output, max.level = 1)
#'
#'  # Access deepSTRAPP results
#'  str(deepSTRAPP_output$STRAPP_results)
#'
#'  # Access trait data
#'  head(deepSTRAPP_output$updated_trait_data_with_Map$trait_data)
#'
#'  # Access the diversification data in a melted data.frame
#'  head(deepSTRAPP_output$diversification_data_df)
#'
#'  # Plot rates vs. trait values across branches
#'  plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_output)
#'
#'  # Plot updated contMap
#'  plot_contMap(deepSTRAPP_output$updated_trait_data_with_Map$contMap)
#'  ape::nodelabels(text =
#'    deepSTRAPP_output$updated_trait_data_with_Map$contMap$tree$initial_nodes_ID)
#'
#'  # Plot diversification rates on updated phylogeny
#'  plot_BAMM_rates(deepSTRAPP_output$updated_BAMM_object, labels = TRUE)
#'
#'  # Plot histogram of test stats
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_output) }
#'
#'  # ----- Example 2: Categorical trait ----- #
#'
#'  ## Load data
#'
#'  # Load phylogeny
#'  data(Ponerinae_tree, package = "deepSTRAPP")
#'  # Load trait df
#'  data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
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
#'  ## Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#'  Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     phylo = Ponerinae_tree_old_calib,
#'     trait_data_type = "categorical",
#'     colors_per_levels = colors_per_states,
#'     evolutionary_models = "ARD", # Use default ARD model
#'     nb_simulations = 100, # Reduce number of simulations to save time
#'     seed = 1234, Seet seed for reproducibility
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE,
#'     plot_map = FALSE) }
#'
#'  # Load directly output
#'  data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
#'
#'  ## Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#'  deepSTRAPP_output <- run_deepSTRAPP_for_focal_time(
#'      densityMaps = Ponerinae_cat_data_old_calib$densityMaps,
#'      ace = Ponerinae_cat_data_old_calib$ace,
#'      tip_data = Ponerinae_cat_3lvl_tip_data,
#'      trait_data_type = "categorical",
#'      BAMM_object = Ponerinae_BAMM_object_old_calib,
#'      focal_time = focal_time,
#'      rate_type = "net_diversification",
#'      posthoc_pairwise_tests = TRUE,
#'      return_perm_data = TRUE,
#'      extract_diversification_data_melted_df = TRUE,
#'      return_updated_trait_data_with_Map = TRUE,
#'      return_updated_BAMM_object = TRUE)
#'
#'  ## Explore output
#'  str(deepSTRAPP_output, max.level = 1)
#'
#'  # Access deepSTRAPP results
#'  str(deepSTRAPP_output$STRAPP_results, max.level = 2)
#'  # Result for overall Kruskal-Wallis test
#'  deepSTRAPP_output$STRAPP_results[1:3]
#'  # Results for posthoc pairwise Dunn's tests
#'  deepSTRAPP_output$STRAPP_results$posthoc_pairwise_tests$summary_df
#'
#'  # Access trait data
#'  head(deepSTRAPP_output$updated_trait_data_with_Map$trait_data)
#'
#'  # Access the diversification data in a melted data.frame
#'  head(deepSTRAPP_output$diversification_data_df)
#'
#'  # Plot rates vs. states across branches
#'  plot_rates_vs_trait_data_for_focal_time(
#'      deepSTRAPP_outputs = deepSTRAPP_output,
#'      colors_per_levels = colors_per_states)
#'
#'  # Plot updated densityMaps cut at focal time
#'  plot_densityMaps_overlay(deepSTRAPP_output$updated_trait_data_with_Map$densityMaps)
#'
#'  # Plot diversification rates on updated phylogeny
#'  plot_BAMM_rates(BAMM_object = deepSTRAPP_output$updated_BAMM_object, legend = TRUE, labels = FALSE,
#'     colorbreaks = deepSTRAPP_output$updated_BAMM_object$initial_colorbreaks$net_diversification)
#'
#'  # Plot histogram of Kruskal-Wallis overall test stats
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_output)
#'
#'  # Plot histograms of posthoc pairwise Dunn's test stats
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_output,
#'     plot_posthoc_tests = TRUE) }
#'
#'  # ----- Example 3: Biogeographic ranges ----- #
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
#'  colors_per_ranges <- c("mediumpurple2", "peachpuff2")
#'  names(colors_per_ranges) <- c("N", "O")
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run evolutionary models
#'  Ponerinae_biogeo_data <- prepare_trait_data(
#'      tip_data = Ponerinae_NO_data,
#'      trait_data_type = "biogeographic",
#'      phylo = Ponerinae_tree_old_calib,
#'      evolutionary_models = "DEC+J", # Default = "DEC" for biogeographic
#'      BioGeoBEARS_directory_path = tempdir(), # Ex: "./BioGeoBEARS_directory/"
#'      keep_BioGeoBEARS_files = FALSE,
#'      prefix_for_files = "Ponerinae_old_calib",
#'      max_range_size = 2,
#'      split_multi_area_ranges = TRUE, # Set to TRUE to display the two outputs
#'      nb_simulations = 100, # Reduce to save time (Default = '1000')
#'      colors_per_levels = colors_per_ranges,
#'      return_model_selection_df = TRUE,
#'      verbose = TRUE) }
#'
#'  # Load directly output
#'  data(Ponerinae_biogeo_data_old_calib, package = "deepSTRAPP")
#'
#'  ## Explore output
#'  str(Ponerinae_biogeo_data_old_calib, 1)
#'
#'  ## Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#'  deepSTRAPP_output <- run_deepSTRAPP_for_focal_time(
#'     densityMaps = Ponerinae_biogeo_data_old_calib$densityMaps,
#'     ace = Ponerinae_biogeo_data_old_calib$ace,
#'     tip_data = Ponerinae_NO_data,
#'     trait_data_type = "biogeographic",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     focal_time = focal_time,
#'     rate_type = "net_diversification",
#'     return_perm_data = TRUE,
#'     extract_diversification_data_melted_df = TRUE,
#'     return_updated_trait_data_with_Map = TRUE,
#'     return_updated_BAMM_object = TRUE)
#'
#'  ## Explore output
#'  str(deepSTRAPP_output, max.level = 1)
#'
#'  # Access deepSTRAPP results
#'  str(deepSTRAPP_output$STRAPP_results, max.level = 2)
#'  # Result for Mann-Whitney-Wilcoxon test
#'  deepSTRAPP_output$STRAPP_results[1:3]
#'
#'  # Access trait data
#'  head(deepSTRAPP_output$updated_trait_data_with_Map$trait_data)
#'
#'  # Access the diversification data in a melted data.frame
#'  head(deepSTRAPP_output$diversification_data_df)
#'
#'  # Plot rates vs. ranges across branches
#'  plot_rates_vs_trait_data_for_focal_time(
#'      deepSTRAPP_outputs = deepSTRAPP_output,
#'      colors_per_levels = colors_per_ranges)
#'
#'  # Plot updated densityMaps cut at focal time
#'  plot_densityMaps_overlay(deepSTRAPP_output$updated_trait_data_with_Map$densityMaps)
#'
#'  # Plot diversification rates on updated phylogeny
#'  plot_BAMM_rates(BAMM_object = deepSTRAPP_output$updated_BAMM_object, legend = TRUE, labels = FALSE,
#'     colorbreaks = deepSTRAPP_output$updated_BAMM_object$initial_colorbreaks$net_diversification)
#'
#'  # Plot histogram of Mann-Whitney-Wilcoxon test
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     STRAPP_results = deepSTRAPP_output$STRAPP_results) }
#' }
#'


run_deepSTRAPP_for_focal_time <- function (contMap = NULL,
                                           densityMaps = NULL,
                                           ace = NULL,
                                           tip_data = NULL,
                                           trait_data_type,
                                           BAMM_object,
                                           focal_time,
                                           keep_tip_labels = TRUE,
                                           rate_type = "net_diversification",
                                           seed = NULL,
                                           nb_permutations = NULL,
                                           replace_samples = FALSE,
                                           alpha = 0.05,
                                           two_tailed = TRUE,
                                           one_tailed_hypothesis = NULL,
                                           posthoc_pairwise_tests = FALSE,
                                           p.adjust_method = "none",
                                           return_perm_data = FALSE,
                                           nthreads = 1,
                                           print_hypothesis = TRUE,
                                           extract_diversification_data_melted_df = FALSE,
                                           return_updated_trait_data_with_Map = FALSE,
                                           return_updated_BAMM_object = FALSE,
                                           verbose = TRUE)
{
  ### Check input validity
  {
    ## Checks should all already be included in the wrapped functions

    ## Convert tip_data into vector of character strings if needed
    if (any(c(is.matrix(tip_data), is.data.frame(tip_data))) & trait_data_type == "biogeographic")
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
  }

  # ------ Extract trait data ------ #

  if (verbose)
  {
    cat(paste0(Sys.time(), " - Extract trait data for focal-time = ", focal_time, "\n\n"))
  }

  trait_data_list <- extract_most_likely_trait_values_for_focal_time(
    contMap = contMap,
    densityMaps = densityMaps,
    ace = ace,
    tip_data = tip_data,
    trait_data_type = trait_data_type,
    focal_time = focal_time,
    update_map = return_updated_trait_data_with_Map,
    keep_tip_labels = keep_tip_labels)

  # ------ Extract diversification data ------ #

  if (verbose)
  {
    cat(paste0(Sys.time(), " - Extract diversification data for focal-time = ", focal_time, "\n"))
  }

  updated_BAMM_object <- update_rates_and_regimes_for_focal_time(
    BAMM_object = BAMM_object,
    focal_time = focal_time,
    update_rates = TRUE,
    update_regimes = TRUE,
    update_tree = return_updated_BAMM_object,
    update_plot = return_updated_BAMM_object,
    update_all_elements = return_updated_BAMM_object,
    keep_tip_labels = keep_tip_labels,
    verbose = verbose)

  # ------ Convert diversification data in a melted df ------ #
  if (extract_diversification_data_melted_df)
  {
    if (verbose)
    {
      cat(paste0("\n", Sys.time(), " - Convert diversification data in a melted data.frame\n"))
    }

    diversification_data_df <- extract_diversification_data_melted_df_for_focal_time(
      BAMM_object = updated_BAMM_object,
      verbose = verbose)
  }

  # ------ Compute STRAPP test ------ #

  if (verbose)
  {
    cat(paste0("\n", Sys.time(), " - Compute STRAPP test for focal-time = ", focal_time, "\n\n"))
  }

  STRAPP_results <- compute_STRAPP_test_for_focal_time(
    BAMM_object = updated_BAMM_object,
    trait_data_list = trait_data_list,
    rate_type = rate_type,
    seed = seed,
    nb_permutations = nb_permutations,
    replace_samples = replace_samples,
    alpha = alpha,
    two_tailed = two_tailed,
    one_tailed_hypothesis = one_tailed_hypothesis,
    posthoc_pairwise_tests = posthoc_pairwise_tests,
    p.adjust_method = p.adjust_method,
    return_perm_data = return_perm_data,
    nthreads = nthreads,
    print_hypothesis = print_hypothesis)

  # ------ Export outputs ------ #

  deepSTRAPP_output <- list(
    STRAPP_results = STRAPP_results,
    focal_time = focal_time)

  # Add the melted df of diversification data if requested
  if (extract_diversification_data_melted_df)
  {
    deepSTRAPP_output$diversification_data_df <- diversification_data_df
  }

  # Add the updated trait data with updated contMap/densityMaps if requested
  if (return_updated_trait_data_with_Map)
  {
    deepSTRAPP_output$updated_trait_data_with_Map <- trait_data_list
  }

  # Add the updated BAMM object if requested
  if (return_updated_BAMM_object)
  {
    deepSTRAPP_output$updated_BAMM_object <- updated_BAMM_object
  }

  ## Return final object
  return(deepSTRAPP_output)

}

