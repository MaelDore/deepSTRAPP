
#' @title Run deepSTRAPP to test for a relationship between diversification rates and trait data over multiple time steps
#'
#' @description Wrapper function to run deepSTRAPP workflows over multiple time steps in the past.
#'   It starts from traits mapped on a phylogeny (trait data) and BAMM output (diversification data)
#'   and carries out the appropriate statistical method to test for a relationship between diversification rates and trait data.
#'   The workflow is repeated over multiple points in time (i.e. the `time_steps`) and results are summarized in a data.frame.
#'   The function can also provide summaries of trait values and diversification rates
#'   extracted along branches over the different `time_steps`.
#'
#'   Statistical tests are based on block-permutations: rates data are randomized across tips following blocks
#'   defined by the diversification regimes identified on each tip (typically from a BAMM).
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
#' @param contMaps For continuous trait data. List of objects of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains multiple continuous stochastic maps that represent possible trait evolutionary histories
#'   conditioned to the observed trait values and model fit.
#' @param densityMaps For categorical trait or biogeographic data. List of objects of class `"densityMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state/range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state/range. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param simmaps For categorical trait or biogeographic data. List of objects of class `"simmap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()] or [phytools::make.simmap()],
#'   that represent discrete character/geographic evolutionary history
#'   (i.e., transitions in character states/geographic ranges) mapped along branches.
#'   This is needed to be able to track which simulated history provided which trait data in downstream analyses
#'   employing the "paired" or "full" strategies to account for uncertainty in ancestral trait estimates.
#' @param nb_simulations For categorical trait or biogeographic data. Integer. The number of stochastic maps used to simulate
#'   trait evolution. This is needed for the "paired" and "full" strategies to account for trait estimate uncertainty,
#'   if only densityMaps summarizing posterior state/range density are provided, but not the simmaps representing all evolutionary histories.
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
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated phylogeny. Default is `TRUE`.
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::prepare_diversification_data()],
#'   that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#'   The phylogenetic tree must the same as the one associated with the `contMap`, `ace` and `tip_data`.
#' @param rate_type A character string specifying the type of diversification rates to use. Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param time_steps Numerical vector. Time steps at which the STRAPP tests should be carried out. If `NULL` (the default),
#'  `time_steps` will be generated from a combination of two arguments among `time_range`, `nb_time_steps`, and/or `time_step_duration`.
#' @param time_range Vector of two numerical values. Time boundaries within with the `time_steps` must be defined if not provided.
#'   If `NULL` (the default), and `time_range` is needed to generate the `time_steps`, the depth of the tree is used by default: `c(0, root_age)`.
#'   However, no time step will be generated for the 'root_age'.
#' @param nb_time_steps,time_step_duration Numerical. Number of time steps and duration of each time step used to generate `time_steps` if not provided.
#'   You must provide at least one of those two arguments to be able to generate `time_steps`.
#' @param uncertainty_strategy Character string. To select the strategy used to account for uncertainty in estimates.
#'   * `"rates_only"`: Only accounts for diversification-rate uncertainty across BAMM posterior samples. Uses ML estimates for continuous traits
#'                     and the most frequent state/range observed across stochastic maps for categorical and biogeographic data.
#'   * `"paired"`: Default option. Accounts for both diversification-rate and ancestral trait/range reconstruction uncertainty by pairing BAMM posterior samples with stochastic maps.
#'                 When the number of BAMM samples and stochastic maps differ, random pairing with replacement from the smaller set is used
#'                 so that all posterior samples and stochastic maps contribute to the analysis.
#'   * `"full"`: Exhaustive option that accounts for trait/range- and rate- uncertainty by crossing all BAMM posterior samples with all stochastic maps.
#'               Accounts for both diversification-rate and ancestral reconstruction uncertainty by evaluating every combination of BAMM posterior sample and stochastic map.
#'               WARNING: This exhaustive approach can substantially increase computation time and memory requirements and is therefore recommended only for moderate-sized analyses.
#' @param trait_maps_vs_BAMM_samples_list (Optional) List of two elements manually providing the names to associate stochastic maps (`$trait_map_ID`)
#'   with BAMM samples (`$BAMM_posterior_sample_ID`). This is typically used to ensure the same stochastic maps and BAMM samples are used to test across multiple time-steps.
#'   Values are the names of the objects such as "Map_X" and "BAMM_X". Default = `NULL`.
#'     * For uncertainty_strategy == 'rates_only', the `$trait_map_ID` must be "Map_ML" as only the ML estimates of trait values/states/ranges are used.
#'     * For uncertainty_strategy == 'paired', each pair of stochastic map and BAMM sample will be used once.
#'     * For uncertainty_strategy == 'full', all stochastic maps will be matched with all BAMM samples.
#'   Those may partly differ from the actual maps and BAMM samples used for the tests as recorded in `$STRAPP_results_over_time[[X]]$perm_data_df`
#'   because invalid maps with not enough states/ranges are discarded at each time-step.
#' @param seed Integer. Set the seed to ensure reproducibility. Default is `NULL` (a random seed is used).
#' @param nb_permutations Integer. To select the number of random permutations to perform during the tests.
#'   If NULL (default), all posterior samples will be used once.
#' @param alpha Numerical. Significance level to use to compute the `estimate` corresponding to the values of the test statistic used to assess significance of the test.
#'   This does NOT affect p-values. Default is `0.05`.
#' @param two_tailed Logical. To define the type of tests. If `TRUE` (default), tests for correlations/differences in rates will be carried out with a null hypothesis
#'   that rates are not correlated with trait values (continuous data) or equals between trait states (categorical and biogeographic data).
#'   If `FALSE`, one-tailed tests are carried out.
#'   * For continuous data, it involves defining a `one_tailed_hypothesis` testing for either a "positive" or "negative" correlation under the alternative hypothesis.
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
#'  This is needed to plot the histograms of the null distribution used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'  (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`). Default is `FALSE`.
#' @param nthreads Integer. Number of threads to use for paralleled computing of the STRAPP tests across the permutations.
#'  The R package `parallel` must be loaded for `nthreads > 1`. Default is `1`.
#' @param print_hypothesis Logical. Whether to print information on what test is carried out, detailing the null and alternative hypotheses,
#'  what significant level is used to rejected or not the null hypothesis, and how uncertainty in trait estimates is handled. Default is `TRUE`.
#' @param extract_trait_data_melted_df Logical. Specify whether trait data must be extracted from the mapped phylogenies at each time step
#'  and returned in a melted data.frame. Default is `FALSE`.
#' @param extract_diversification_data_melted_df Logical. Specify whether diversification data (regimes ID and tip rates) must be extracted from the `updated_BAMM_object`
#'   at each time step and returned in a melted data.frame. Default is `FALSE`.
#' @param return_STRAPP_results Logical. Specify whether the `STRAPP_results` objects summarizing the results of the STRAPP tests carried out at each time step
#'   should be returned among the outputs in addition to the `$pvalues_summary_df` already providing test stat estimates and p-values obtained across all `time_steps`.
#' @param return_updated_Maps Logical. Specify whether the updated versions of the mapped phylogenies (`contMap(s)`/`densityMaps`/`simmaps`),
#'   should be returned among the outputs. The updated `contMap(s)`/`densityMaps`/`simmaps` consist in cutting off branches and mapping
#'   that are younger than the `time_steps` used to perform the tests. Default is `FALSE`.
#' @param return_updated_BAMM_objects Logical. Specify whether the `updated_BAMM_objects` with phylogeny and
#'   mapped diversification rates cut-off at the `time_steps` used to perform the tests should be returned among the outputs.
#' @param verbose Logical. Should progression per `time_steps` be displayed? Default is `TRUE`.
#' @param verbose_extended Should progression per `time_steps` AND within each deepSTRAPP workflow de displayed?
#'   In addition to printing progress along `time_steps`, a message will be printed at each step of the deepSTRAPP workflow,
#'   and for every batch of 100 BAMM posterior samples whose rates are regimes are updated. If `extract_diversification_data_melted_df = TRUE`,
#'   a message for will also be printed when rates are extracted. Default is `FALSE`.
#'
#' @export
#' @importFrom phytools nodeHeights
#'
#' @details The function is a wrapper of [deepSTRAPP::run_deepSTRAPP_for_focal_time()] that runs the
#'   deepSTRAPP workflow over multiple `time_steps`.
#'
#'   The deepSTRAPP workflow is described step by step in the [deepSTRAPP::run_deepSTRAPP_for_focal_time()] documentation.
#'
#'   Its main output is the `$pvalues_summary_df`: a data.frame providing test stat estimates and p-values obtained across all `time_steps`,
#'   that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
#'   If using multinominal data (with more than two states) and `posthoc_pairwise_tests = TRUE`, the output will also contain
#'   a data.frame providing test stat estimates and p-values for post hoc pairwise tests in `$pvalues_summary_df_for_posthoc_pairwise_tests`.
#'
#'   The function offers options to generate summary data.frames of the data extracted across `time_steps`:
#'   * If `extract_trait_data_melted_df = TRUE`, a data.frame of trait values found along branches at each time step
#'     is provided in `$trait_data_df_over_time`.
#'   * If `extract_diversification_data_melted_df = TRUE`, a data.frame of diversification data (regimes ID and tip rates)
#'     found along branches at each time step is provided in `$diversification_data_df_over_time`.
#'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
#'     the evolution diversification rates across trait values over time.
#'
#'   The function also allows to keep records of the intermediate objects generated during the STRAPP workflow:
#'   * If `return_STRAPP_results = TRUE`, a list of STRAPP test outputs is provided in `$STRAPP_results_over_time`.
#'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
#'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
#'   * If `return_updated_Maps = TRUE`, a list of objects containing the updated `contMap(s)`/`densityMaps`/`simmaps` provided as inputs
#'     is provided in `$updated_Maps_over_time`.
#'     The updated mapped phylogenies can be visualized with [deepSTRAPP::plot_contMap()] or [deepSTRAPP::plot_densityMaps_overlay()],
#'     to display a phylogeny mapped with trait values/states/ranges with branches cut at each `focal_time`.
#'   * If `return_updated_BAMM_objects = TRUE`, a list of updated `BAMM_objects` of class `"bammdata"` that contains rates and regimes ID
#'     found at each `focal_time` is provided in `$updated_BAMM_objects_over_time`.
#'     Updated `BAMM_objects` can be plotted with [deepSTRAPP::plot_BAMM_rates()] to display
#'     a phylogeny mapped with diversification rates with branches cut at each `focal_time`.
#'
#' @return The function returns a list with at least five elements.
#'
#'   * `$pvalues_summary_df` Data.frame with three columns providing test stat `$estimate` and `$p_value` obtained for each time step (i.e., `$focal_time`),
#'     that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
#'   * `$time_steps` Numerical vector. Time steps at which the STRAPP tests were carried out in the same order as the objects returned in the output lists.
#'   * `$trait_data_type` Character string. Specify the type of trait data. Possible values are: "continuous", "categorical", "biogeographic".
#'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. One of 'continuous', 'binary', or 'multinominal'.
#'   * `$rate_type` Character string. The type of diversification rates used in the tests: 'speciation', 'extinction' or 'net_diversification'.
#'   * `$uncertainty_strategy` Character string. The strategy used to account for uncertainty in estimates. One of 'rates_only', 'paired', or 'full'.
#'   * `$trait_maps_vs_BAMM_samples_list` List of two elements recording the stochastic maps (`$trait_map_ID`) and BAMM samples (`$BAMM_posterior_sample_ID`) chosen for testing across time-steps.
#'     Those may partly differ from the actual maps and BAMM samples used for the tests as recorded in `$STRAPP_results_over_time[[X]]$perm_data_df` because invalid maps with not enough states/ranges are discarded.
#'
#'   Optional summary df for multinominal data, if `posthoc_pairwise_tests = TRUE`:
#'   * `$pvalues_summary_df_for_posthoc_pairwise_tests` Data.frame with four or five columns providing test stat `$estimate`, `$p_value`, and `$p_value_adjusted`
#'     (if `p.adjust_method` used is not "none") for each `$pair` of states involved in post hoc Dunn's tests obtained for each time step (i.e., `$focal_time`).
#'      This data.frame can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the post hoc test results across time.
#'
#'   Optional melted data.frames:
#'   * `$trait_data_df_over_time` Data.frame with five columns providing `$trait_value` associated with each `$tip_ID` found along each time step (i.e., `$focal_time`)
#'     across all stochastic maps (`$Map_ID`). Set `extract_trait_data_melted_df = TRUE` to include it in the output.
#'   * `$diversification_data_df_over_time` Data.frame with six columns providing diversification regimes (`$regime_ID`) and `$rates` sorted by `$rate_type` along tips (`$tip_ID`)
#'     found across all posterior samples (`$BAMM_sample_ID`) over each time step (i.e., `$focal_time`).
#'     Set `extract_diversification_data_melted_df = TRUE` to include it in the output.
#'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
#'     the evolution diversification rates across trait values over time.
#'
#'   Optional objects generated for each time step (i.e., `focal_time`) and ordered as in `$time_steps`:
#'   * `$STRAPP_results_over_time` List of objects summarizing the results of the STRAPP tests
#'     See [compute_STRAPP_test_for_focal_time()] for a detailed description of the elements in each object.
#'     Set `return_STRAPP_results = TRUE` to include it in the output.
#'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
#'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
#'   * `$updated_Maps_over_time` List of objects containing updated `contMap(s)`/`densityMaps`/`simmaps`.
#'     The updated mapped phylogenies can be plotted with [deepSTRAPP::plot_contMap()] or [deepSTRAPP::plot_densityMaps_overlay()],
#'     to display a phylogeny mapped with trait values with branches cut at each `focal_time`.
#'   * `$updated_BAMM_objects_over_time` List of objects containing rates and regimes ID mapped on phylogeny.
#'     Updated `BAMM_objects` can be plotted with [deepSTRAPP::plot_BAMM_rates()] to display a phylogeny mapped with
#'     diversification rates with branches cut at each `focal_time`.
#'
#' @author Maël Doré
#'
#' @seealso To run the deepSTRAPP workflow for a single `focal_time`: [deepSTRAPP::run_deepSTRAPP_for_focal_time()]
#' [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
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
#' # ----- Example 1: Continuous trait ----- #
#'  ## Load data
#'
#'  # Load trait df
#' data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
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
#'  # Map trait evolution across multiple simulations (i.e., continuous stochastic maps)
#'  Ponerinae_cont_data_old_calib <- prepare_trait_data(
#'      tip_data = Ponerinae_cont_tip_data,
#'      trait_data_type = "continuous",
#'      phylo = Ponerinae_tree_old_calib,
#'      seed = 1234,
#'      evolutionary_models = "BM",
#'      plot_map = FALSE,
#'      verbose = TRUE)
#'
#'  ## Load directly trait data output
#'  data(Ponerinae_cont_data_old_calib, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  # Plot contMap = ML estimates of continuous trait evolution
#'  plot_contMap(contMap = Ponerinae_cont_data_old_calib$contMap,
#'               color_scale = color_scale)
#'
#'  ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
#'  # nb_time_steps <- 5
#'  time_step_duration <- 5
#'  time_range <- c(0, 40)
#'
#'  ## Run deepSTRAPP on net diversification rates
#'  Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     contMap = Ponerinae_cont_data_old_calib$contMap, # Provide ML estimates as a contMap
#'     ace = Ponerinae_cont_data_old_calib$ace,
#'     tip_data = Ponerinae_cont_tip_data,
#'     trait_data_type = "continuous",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     # nb_time_steps = nb_time_steps,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     uncertainty_strategy = "rates_only",
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
#'  ## Load directly deepSTRAPP output
#'  data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#' # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Explore output
#' str(Ponerinae_deepSTRAPP_cont_old_calib_0_40, max.level = 1)
#'  deepSTRAPP_outputs <- Ponerinae_deepSTRAPP_cont_old_calib_0_40
#'
#'  # Display test summary
#'  # Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
#'  # showing the evolution of the test results across time.
#'  deepSTRAPP_outputs$pvalues_summary_df
#'
#'  # Access trait data in a melted data.frame
#'  head(deepSTRAPP_outputs$trait_data_df_over_time)
#'  # Trait data includes only the ML estimates
#'  table(deepSTRAPP_outputs$trait_data_df_over_time$Map_ID)
#'
#'  # Access the diversification data in a melted data.frame
#'  head(deepSTRAPP_outputs$diversification_data_df_over_time)
#'  # Diversification data includes 1000 BAMM posteriors
#'  table(deepSTRAPP_outputs$diversification_data_df$BAMM_sample_ID)
#'
#'  # Access deepSTRAPP results ordered by time-steps
#'  str(deepSTRAPP_outputs$STRAPP_results, max.level = 2)
#'
#'  ## Visualize updated phylogenies
#'
#'  \donttest{ # (May take time to plot)
#'  # Plot updated contMap for time step n°2
#'  contMap_step2 <- deepSTRAPP_outputs$updated_Maps_over_time[[2]]$contMap
#'  plot_contMap(contMap_step2, ftype = "off", color_scale = color_scale)
#'
#'  # Plot diversification rates on updated phylogeny for time step n°2
#'  BAMM_object_step2 <- deepSTRAPP_outputs$updated_BAMM_objects_over_time[[2]]
#'  plot_BAMM_rates(BAMM_object = BAMM_object_step2,
#'    legend = TRUE, labels = FALSE,
#'    colorbreaks = BAMM_object_step2$initial_colorbreaks$net_diversification) }
#'
#'  ## Visualize test results
#'
#'  \donttest{ # (May take time to plot)
#'  # Plot p-values of Spearman tests across all time-steps
#'  plot_STRAPP_pvalues_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     alpha = 0.05)
#'
#'  # Plot evolution of mean rates through time
#'  plot_rates_through_time(deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40)
#'
#'  # Plot rates vs. trait values across branches for time step = 10 My
#'  plot_rates_vs_trait_data_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     focal_time = 10)
#'
#' # Plot histogram of Spearman test stats for time step = 10 My
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     focal_time = 10)
#'
#'  # Plot histograms of Spearman test results (One plot per time-step)
#'  plot_histograms_STRAPP_tests_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     display_plots = TRUE)  }
#'
#'  # ----- Example 2: Categorical trait ----- #
#'
#' ## Load data
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
#'     evolutionary_models = "ER", # Use default ER model
#'     nb_simulations = 100, # Reduce number of simulations to save time
#'     seed = 1234, # Set seed for reproducibility
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE,
#'     plot_map = FALSE) }
#'
#'   # Load directly trait data output
#'  data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
#'
#' ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
#'  # nb_time_steps <- 5
#'  time_step_duration <- 5
#'  time_range <- c(0, 40)
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates across time-steps.
#'  Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
#'     nb_simulations = 100,
#'     ace = Ponerinae_cat_3lvl_data_old_calib$ace,
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     trait_data_type = "categorical",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     # nb_time_steps = nb_time_steps,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     uncertainty_strategy = "paired",
#'     rate_type = "net_diversification",
#'     seed = 1234,
#'     alpha = 0.10, # Select a generous level of significance for the sake of the example
#'     posthoc_pairwise_tests = TRUE,
#'     return_perm_data = TRUE,
#'     extract_trait_data_melted_df = TRUE,
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_Maps = TRUE,
#'     return_updated_BAMM_objects = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly deepSTRAPP output
#'  data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, package = "deepSTRAPP")
#'  deepSTRAPP_outputs <- Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Explore output
#'  str(deepSTRAPP_outputs, max.level = 1)
#'
#'  # Display test summaries
#'  # Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
#'  # showing the evolution of the test results across time.
#'  deepSTRAPP_outputs$pvalues_summary_df
#'  # Results for posthoc pairwise Dunn's tests over time-steps
#'  deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests
#'
#'  # Access trait data in a melted data.frame
#'  # Trait data were extracted from densityMaps recording only frequency of states.
#'  # Thus, trait states have been distributed accordingly across 'Dummy_maps'
#'  # to reproduce the recorded frequencies.
#'  head(deepSTRAPP_outputs$trait_data_df_over_time)
#'  table(deepSTRAPP_outputs$trait_data_df_over_time$Map_ID)
#'
#'  # Access the diversification data in a melted data.frame
#'  head(deepSTRAPP_outputs$diversification_data_df_over_time)
#'  # Diversification data includes 1000 BAMM posteriors
#'  table(deepSTRAPP_outputs$diversification_data_df_over_time$BAMM_sample_ID)
#'
#'  # Access details of deepSTRAPP results across time-steps
#'  str(deepSTRAPP_outputs$STRAPP_results, max.level = 2)
#'
#'  ## Visualize updated phylogenies
#'
#'  \donttest{ # (May take time to plot)
#'  # Plot updated densityMaps for time step n°2 = 10My
#'  densityMaps_10My <- deepSTRAPP_outputs$updated_Maps_over_time[[2]]
#'  plot_densityMaps_overlay(densityMaps_10My$densityMaps, colors_per_levels = colors_per_states)
#'
#'  # Plot diversification rates on updated phylogeny for time step n°2
#'  BAMM_object_10My <- deepSTRAPP_outputs$updated_BAMM_objects_over_time[[2]]
#'  plot_BAMM_rates(BAMM_object = BAMM_object_10My,
#'     legend = TRUE, labels = FALSE,
#'    colorbreaks = BAMM_object_10My$initial_colorbreaks$net_diversification) }
#'
#'  ## Visualize test results
#'
#'  \donttest{ # (May take time to plot)
#'  # Plot p-values of overall Kruskal-Wallis test across all time-steps
#'  plot_STRAPP_pvalues_over_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'    alpha = 0.10)
#'
#'  # Plot p-values of post hoc pairwise Dunn's tests between pairs of tests across all time-steps
#'  plot_STRAPP_pvalues_over_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     alpha = 0.10,
#'     plot_posthoc_tests = TRUE)
#'
#'  # Plot evolution of mean rates through time
#'  plot_rates_through_time(deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     colors_per_levels = colors_per_states)
#'
#'  # Plot rates vs. trait values across branches for time step n°2 = 10 My
#'  plot_rates_vs_trait_data_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     focal_time = 10,
#'     colors_per_levels = colors_per_states)
#'
#'  # Plot histogram of overall Kruskal-Wallis test for time step n°2 = 10 My
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     focal_time = 10)
#'
#'  # Plot histograms of overall Kruskal-Wallis test results across all time-steps
#'  # (One plot per time-step)
#'  plot_histograms_STRAPP_tests_over_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     plot_posthoc_tests = FALSE)
#'
#'  # Plot histograms of post hoc pairwise Dunn's test results across all time-steps
#'  # (One multifaceted plot per time-step)
#'  plot_histograms_STRAPP_tests_over_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     plot_posthoc_tests = TRUE) }
#'
#' # ----- Example 3: Biogeographic ranges ----- #
#'
#' ## Load data
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
#'  Ponerinae_NO_tip_data <- stats::setNames(object = Ponerinae_binary_range_table$Old_World,
#'                                       nm = Ponerinae_binary_range_table$Taxa)
#'  Ponerinae_NO_tip_data <- as.character(Ponerinae_NO_tip_data)
#'  Ponerinae_NO_tip_data[Ponerinae_NO_tip_data == "TRUE"] <- "O" # O = Old World
#'  Ponerinae_NO_tip_data[Ponerinae_NO_tip_data == "FALSE"] <- "N" # N = New World
#'  names(Ponerinae_NO_tip_data) <- Ponerinae_binary_range_table$Taxa
#'  table(Ponerinae_NO_tip_data)
#'
#'  colors_per_ranges <- c("mediumpurple2", "peachpuff2")
#'  names(colors_per_ranges) <- c("N", "O")
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run evolutionary models
#'  Ponerinae_biogeo_data <- prepare_trait_data(
#'    tip_data = Ponerinae_NO_tip_data,
#'     trait_data_type = "biogeographic",
#'     phylo = Ponerinae_tree_old_calib,
#'     evolutionary_models = "DEC+J", # Default = "DEC" for biogeographic
#'     BioGeoBEARS_directory_path = tempdir(), # Ex: "./BioGeoBEARS_directory/"
#'     keep_BioGeoBEARS_files = FALSE,
#'     prefix_for_files = "Ponerinae_old_calib",
#'     max_range_size = 2,
#'     split_multi_area_ranges = TRUE, # Set to TRUE to use only unique areas "A" and "B"
#'     nb_simulations = 100, # Reduce to save time (Default = '1000')
#'     colors_per_levels = colors_per_ranges,
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
#'    densityMaps = Ponerinae_biogeo_data_old_calib$densityMaps,
#'     nb_simulations = 100,
#'     ace = Ponerinae_biogeo_data_old_calib$ace,
#'     tip_data = Ponerinae_NO_tip_data,
#'     trait_data_type = "biogeographic",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     seed = 1234, # Set seed for reproducibility
#'     alpha = 0.10, # Select a generous level of significance for the sake of the example
#'     uncertainty_strategy = "paired",
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
#'  data(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40, package = "deepSTRAPP")
#'  deepSTRAPP_outputs <- Ponerinae_deepSTRAPP_biogeo_old_calib_0_40
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Explore output
#'  str(deepSTRAPP_outputs, max.level = 1)
#'
#'  # Display test summaries
#'  # Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
#'  # showing the evolution of the test results across time.
#'  deepSTRAPP_outputs$pvalues_summary_df
#'
#'  # Access bioregeographic range data in a melted data.frame
#'  head(deepSTRAPP_outputs$trait_data_df_over_time)
#'  # Ranges were extracted from densityMaps recording only frequency of states.
#'  # Thus, ranges have been distributed accordingly across 'Dummy_maps'
#'  # to reproduce the recorded frequencies.
#' table(deepSTRAPP_outputs$trait_data_df_over_time$Map_ID)
#'
#'  # Access the diversification data in a melted data.frame
#'  head(deepSTRAPP_outputs$diversification_data_df_over_time)
#'  # Diversification data includes 1000 BAMM posteriors
#'  table(deepSTRAPP_outputs$diversification_data_df_over_time$BAMM_sample_ID)
#'
#'  # Access details of deepSTRAPP results
#'  str(deepSTRAPP_outputs$STRAPP_results, max.level = 2)
#'
#'  ## Visualize updated phylogenies
#'
#'  \donttest{ # (May take time to plot)
#'  # Plot updated densityMaps for time step n°2 = 10My
#'  densityMaps_10My <- deepSTRAPP_outputs$updated_Maps_over_time[[2]]
#'  plot_densityMaps_overlay(densityMaps_10My$densityMaps)
#'
#'  # Plot diversification rates on updated phylogeny for time step n°2
#'  BAMM_object_10My <- deepSTRAPP_outputs$updated_BAMM_objects_over_time[[2]]
#'  plot_BAMM_rates(BAMM_object = BAMM_object_10My,
#'    legend = TRUE, labels = FALSE,
#'    colorbreaks = BAMM_object_10My$initial_colorbreaks$net_diversification) }
#'
#'  ## Visualize test results
#'
#'  \donttest{ # (May take time to plot)
#'  # Plot p-values of Mann-Whitney-Wilcoxon tests across all time-steps
#'  plot_STRAPP_pvalues_over_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     alpha = 0.10) # Set a generous significance threshold for the example
#'
#' # Plot evolution of mean rates through time
#'  plot_rates_through_time(deepSTRAPP_outputs = deepSTRAPP_outputs,
#'    colors_per_levels = colors_per_ranges)
#'
#'  # Plot rates vs. trait values across branches for time step n°2 = 10 My
#'  plot_rates_vs_trait_data_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     focal_time = 10,
#'     colors_per_levels = colors_per_ranges)
#'
#'  # Plot histogram of Mann-Whitney-Wilcoxon test stats for time step n°2 = 10My
#'  plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'    focal_time = 10)
#'
#' # Plot histograms of Mann-Whitney-Wilcoxon test stats for all time-steps (One plot per time-step)
#'  plot_histograms_STRAPP_tests_over_time(
#'     deepSTRAPP_outputs = deepSTRAPP_outputs,
#'     display_plots = TRUE,
#'     plot_posthoc_tests = FALSE) }
#' }
#'


run_deepSTRAPP_over_time <- function (contMap = NULL,
                                      contMaps = NULL,
                                      densityMaps = NULL,
                                      simmaps = NULL,
                                      nb_simulations = NULL,
                                      ace = NULL,
                                      tip_data = NULL,
                                      trait_data_type,
                                      keep_tip_labels = TRUE,
                                      BAMM_object,
                                      rate_type = "net_diversification",
                                      time_steps = NULL,
                                      time_range = NULL,
                                      nb_time_steps = NULL,
                                      time_step_duration = NULL,
                                      uncertainty_strategy = "paired",
                                      trait_maps_vs_BAMM_samples_list = NULL,
                                      seed = NULL,
                                      nb_permutations = NULL,
                                      alpha = 0.05,
                                      two_tailed = TRUE,
                                      one_tailed_hypothesis = NULL,
                                      posthoc_pairwise_tests = FALSE,
                                      p.adjust_method = "none",
                                      return_perm_data = FALSE,
                                      nthreads = 1,
                                      print_hypothesis = TRUE,
                                      extract_trait_data_melted_df = FALSE,
                                      extract_diversification_data_melted_df = FALSE,
                                      return_STRAPP_results = FALSE,
                                      return_updated_Maps = FALSE,
                                      return_updated_BAMM_objects = FALSE,
                                      verbose = TRUE,
                                      verbose_extended = FALSE)
{
  ### Check input validity
  {
    ## contMap(s) OR densityMaps/simmaps
    if (((!is.null(contMap) | !is.null(contMaps)) & (!is.null(densityMaps) | !is.null(simmaps))))
    {
      stop(paste0("You must provide a 'contMap'/'contMaps' (for continuous traits) OR a 'densityMaps'/'simmaps' (for categorical and biogeographic traits) but not both, according to the type of trait data.\n",
                  "'ace' and tip_data' can also be provided in complement to use accurate ML estimates of trait values at nodes and tips."))
    }

    ## trait_data_type
    # trait_data_type must be "continuous", categorical" or "biogeographic"
    if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
    {
      stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
    }

    ## Check that what is provided contMap(s) OR densityMaps/simmaps match the trait_data_type
    if ((!is.null(contMap) & trait_data_type != "continuous"))
    {
      stop(paste0("You provided a 'contMap' but selected '",trait_data_type,"' as 'trait_data_type'. contMaps are used to map continuous traits.\n",
                  "If you wish to extract trait values for a continuous trait, provide 'trait_data_type = continuous'.\n",
                  "If you wish to extract trait states/ranges for ",trait_data_type," data, provide 'densityMaps' or 'simmaps' as input instead of a 'contMap'"))
    }
    if ((!is.null(contMaps) & trait_data_type != "continuous"))
    {
      stop(paste0("You provided 'contMaps' but selected '",trait_data_type,"' as 'trait_data_type'. contMaps are used to map continuous traits.\n",
                  "If you wish to extract trait values for a continuous trait, provide 'trait_data_type = continuous'.\n",
                  "If you wish to extract trait states/ranges for ",trait_data_type," data, provide 'densityMaps' or 'simmaps' as input instead of 'contMaps'"))
    }
    if (((is.null(contMap) & is.null(contMaps)) & trait_data_type == "continuous"))
    {
      stop(paste0("You selected 'trait_data_type = continuous' but did not provide 'contMap' or 'contMaps' as input.\n",
                  "contMap(s) are needed to extract continuous trait data.\n",
                  "See ?deepSTRAPP::prepare_trait_data(), ?phytools::contMap(), ?contsimmap::make.contsimmap(), and ?deepSTRAPP::convert_contsimmap_to_contMaps_list() to learn how to generate those objects."))
    }

    if ((!is.null(densityMaps) & !(trait_data_type %in% c("categorical", "biogeographic"))))
    {
      stop(paste0("You provided 'densityMaps' but selected '",trait_data_type,"' as 'trait_data_type'. densityMaps are used to map categorical or biogeographic data.\n",
                  "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                  "If you wish to extract trait values for a continuous trait, provide 'contMaps' as input instead of 'densityMaps'.\n"))
    }
    if ((!is.null(simmaps) & !(trait_data_type %in% c("categorical", "biogeographic"))))
    {
      stop(paste0("You provided 'simmaps' but selected '",trait_data_type,"' as 'trait_data_type'. simmaps are used to map categorical or biogeographic data.\n",
                  "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                  "If you wish to extract trait values for a continuous trait, provide 'contMaps' as input instead of 'simmaps'.\n"))
    }
    if ((is.null(densityMaps) & is.null(simmaps)) & (trait_data_type %in% c("categorical", "biogeographic")))
    {
      stop(paste0("You selected 'trait_data_type = ",trait_data_type,"' but did not provide 'densityMaps' or 'simmaps' as input.\n",
                  "'simmaps' are needed to extract state/range data from stochastic maps.\n",
                  "Alternatively, 'densityMaps' and 'nb_simulations' can be provided to reconstruct states/ranges from posterior probabilities across dummy stochastic maps."))
    }

    ## Extract root_age from a valid input mapped phylogeny
    if (!is.null(contMap))
    {
      root_age <- max(phytools::nodeHeights(contMap$tree)[,2])
    }
    if (!is.null(contMaps))
    {
      root_age <- max(phytools::nodeHeights(contMaps[[1]]$tree)[,2])
    }
    if (!is.null(densityMaps))
    {
      root_age <- max(phytools::nodeHeights(densityMaps[[1]]$tree)[,2])
    }
    if (!is.null(simmaps))
    {
      root_age <- max(phytools::nodeHeights(simmaps[[1]])[,2])
    }

    ## If provided, check that time_range is two positive numerical values, in increasing order, not above root_age.
    if (!is.null(time_range))
    {
      if (!all(time_range >= 0))
      {
        stop("'time_range' must be two positive numerical values.")
      } else {
        if (time_range[1] > time_range[2])
        {
          stop("'time_range' values must be provided in increasing order.")
        } else {
          if (time_range[2] > root_age)
          {
            stop("'time_range' older boundary must be equal or lower than the age of the root of the tree.\n",
                 "Max time_range = ",max(time_range),". Root age = ",root_age, ".")
          }
          if (time_range[2] == root_age)
          {
            cat(paste0("WARNING: 'time_range' older boundary is equal to the age of the root of the tree.\n",
                       "\t Although, no time step is generated for the root age.\n"))
          }
        }
      }
    }

    ## If provided, check that nb_time_steps is a positive integer (warning over 100 against time consuming option).
    if (!is.null(nb_time_steps))
    {
      if (!(round(nb_time_steps) == nb_time_steps) | !(nb_time_steps > 0))
      {
        stop("'nb_time_steps' must be a positive integer.")
      } else {
        if (nb_time_steps > 100)
        {
          cat("WARNING: Setting 'nb_time_steps' higher than 100 may take a long-time and produce very large outputs.\n")
        }
      }
    }

    ## If provided, check that time_step_duration is a positive numerical, smaller than root_age.
    if (!is.null(time_step_duration))
    {
      if (!(time_step_duration > 0))
      {
        stop("'time_step_duration' must be a positive numerical.\n")
      } else {
        if (time_step_duration >= root_age)
        {
          stop("'time_step_duration' must be smaller than the root of the tree.\n",
               "'time_step_duration' = ",time_step_duration,". Root age = ",root_age,".")
        }
      }
    }

    ## If time_step_duration AND time_range provided, check that time_step_duration is smaller or equal to the range of time_range.
    if (!is.null(time_step_duration) && !is.null(time_range))
    {
      if (time_step_duration > diff(time_range))
      {
        stop("'time_step_duration' must be smaller or equal to the range of 'time_range'.\n",
             "'time_step_duration' = ",time_step_duration,". 'time_range' = ",time_range,".")
      }
    }

    ## Check validity of time steps arguments while generating the time steps vector
    if (is.null(time_steps))
    {
      cat(paste0("WARNING: 'time_steps' were not provided. Instead, they are tentatively generated from a combination of 'time_range', 'nb_time_steps', and/or 'time_step_duration'.\n"))

      if (!is.null(time_range) && !is.null(nb_time_steps) && !is.null(time_step_duration))
      {
        stop("In order to generate the 'time_steps', please provide only two of the three arguments: 'time_range', 'nb_time_steps' and/or 'time_step_duration'.")
      }

      if (is.null(nb_time_steps) && is.null(time_step_duration))
      {
        stop("In order to generate the 'time_steps', you must provide at least one of those two arguments: 'nb_time_steps' or 'time_step_duration'.")
      }

      # Case with 'nb_time_steps' AND 'time_step_duration'
      if (!is.null(nb_time_steps) && !is.null(time_step_duration))
      {
        # Produce time_steps from nb_time_steps & time_step_duration
        time_steps <- seq(from = 0, by = time_step_duration, length.out = nb_time_steps)

        # Check that generated time steps are within the range of the tree age.
        if (max(time_steps) > root_age)
        {
          stop(paste0("Some generated 'time_steps' are older than the age of the root of the tree.\n",
                      "Max time step = ",max(time_steps),". Root age = ",root_age, ".\n",
                      "You must provide 'nb_time_steps' and 'time_step_duration' arguments that produce time steps that are compatible with the root age."))
        }
        # Case if the oldest generated time step equals the root age => ignore it.
        if (max(time_steps) == root_age)
        {
          cat(paste0("WARNING: The oldest time step generated from 'nb_time_steps' and 'time_step_duration' equals the root age.\n",
                     "\t This time step will be ignored.\n"))
          time_steps <- time_steps[-length(time_steps)]
        }

      } else {
        # Cases with 'nb_time_steps' OR 'time_step_duration'

        # Check time_range
        # If absent, use c(0, root_age) and send warning
        if (is.null(time_range))
        {
          cat(paste0("WARNING: 'time_range' was not provided. The depth of the tree is used by default: c(0, root_age).\n",
                     "\t Although, no time step is generated for the root age.\n"))
          time_range <- c(0, root_age)
        }

        # Case with 'nb_time_steps' and no 'time_step_duration'
        if (!is.null(nb_time_steps))
        {
          # Produce time_steps from time_range & nb_time_steps
          time_steps <- seq(from = time_range[1], to = time_range[2], length.out = nb_time_steps)
          # Case if the oldest generated time step equals the root age => ignore it.
          if (max(time_steps) == root_age)
          {
            cat(paste0("WARNING: The oldest time step generated from 'time_range' and 'nb_time_steps' equals the root age.\n",
                       "\t This time step will be ignored.\n"))
            time_steps <- time_steps[-length(time_steps)]
          }
        } else {
          # Case with no 'nb_time_steps' and 'time_step_duration'
          # Produce time_steps from time_range & time_step_duration
          time_steps <- seq(from = time_range[1], to = time_range[2], by = time_step_duration)
          # Case if the oldest generated time step equals the root age => ignore it.
          if (max(time_steps) == root_age)
          {
            cat(paste0("WARNING: The oldest time step generated from 'time_range' and 'time_step_duration' equals the root age.\n",
                       "\t This time step will be ignored.\n"))
            time_steps <- time_steps[-length(time_steps)]
          }
        }
      }

      # Print time steps generated
      if (verbose | verbose_extended)
      {
        cat(paste0("Time steps generated from a combination of 'time_range', 'nb_time_steps', and/or 'time_step_duration': ",paste(time_steps, collapse = ", "),".\n\n"))
      }
    }

    ## Check if there is more than one time step
    if (length(time_steps) < 2)
    {
      stop(paste0("Only one time step was provided/generated.\n",
                  "Please use run_deepSTRAPP_for_focal_time() to run the deepSTRAPP workflow for a unique time step/focal time."))
    }

    ## If return_perm_data = TRUE, but return_STRAPP_results = FALSE, the distribution of the test statistics will not be returned in the output
    if (return_perm_data && !return_STRAPP_results)
    {
      stop("'return_perm_data' = TRUE but 'return_STRAPP_results' = FALSE.\n",
           "The distribution of the test statistics will not be returned in the output if 'return_STRAPP_results' = FALSE.\n",
           "Set 'return_STRAPP_results' = TRUE if you want to retrieve the test statistics in order to plot histograms of the STRAPP test results.")
    }
    # Only useful if return_STRAPP_results = TRUE

    ## Other validity checks should all already be included in the wrapped functions
  }

  ### Run deepSTRAPP workflow per time steps

  deepSTRAPP_outputs_over_time <- list()
  for (i in seq_along(time_steps))
  {
    # i <- 1

    # Extract current focal_time
    focal_time_i <- time_steps[i]

    # Print progress per time steps
    if (verbose | verbose_extended)
    {
      cat(paste0(Sys.time(), " - deepSTRAPP running for focal_time = ",focal_time_i," - n\u00B0 ", i, "/", length(time_steps),"\n"))
    }

    ## Run deepSTRAPP workflow for the current focal_time
    deepSTRAPP_output_i <- run_deepSTRAPP_for_focal_time(
      contMap = contMap,
      contMaps = contMaps,
      densityMaps = densityMaps,
      simmaps = simmaps,
      nb_simulations = nb_simulations,
      ace = ace,
      tip_data = tip_data,
      trait_data_type = trait_data_type,
      keep_tip_labels = keep_tip_labels,
      BAMM_object = BAMM_object,
      rate_type = rate_type,
      focal_time = focal_time_i,
      uncertainty_strategy = uncertainty_strategy,
      trait_maps_vs_BAMM_samples_list = trait_maps_vs_BAMM_samples_list,
      seed = seed,
      nb_permutations = nb_permutations,
      alpha = alpha,
      two_tailed = two_tailed,
      one_tailed_hypothesis = one_tailed_hypothesis,
      posthoc_pairwise_tests = posthoc_pairwise_tests,
      p.adjust_method = p.adjust_method,
      return_perm_data = return_perm_data,
      nthreads = nthreads,
      print_hypothesis = print_hypothesis,
      extract_trait_data_melted_df = extract_trait_data_melted_df,
      extract_diversification_data_melted_df = extract_diversification_data_melted_df,
      return_updated_Maps = return_updated_Maps,
      return_updated_BAMM_object = return_updated_BAMM_objects,
      verbose = verbose_extended)

    ## Store output for current focal time
    deepSTRAPP_outputs_over_time[[i]] <- deepSTRAPP_output_i

    ## Special case for the first step
    if (i == 1)
    {
      # Retrieve the maps/BAMM samples used for testing and update the input such as the same are used for each time-step
      trait_maps_vs_BAMM_samples_list <- deepSTRAPP_output_i$STRAPP_results$trait_maps_vs_BAMM_samples_list
      # Remove nb_permutations as pairing has been fixed already for all subsequent steps
      nb_permutations <- NULL
    }
  }

  ## Initiate list for final output
  final_ouput <- list()

  ## Extract time-steps with STRAPP results
  trait_data_type_for_stats_list <- unlist(lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$STRAPP_results$trait_data_type_for_stats }))
  time_steps_with_STRAPP_results_ID <- which(trait_data_type_for_stats_list != "none")
  time_steps_with_STRAPP_results <- time_steps[time_steps_with_STRAPP_results_ID]

  ## Extract p_values and estimates in summary_df
  p_values <- unlist(lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$STRAPP_results$p_value }))
  estimates <- unlist(lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$STRAPP_results$estimate }))

  # Store values for all steps, including NA for time-steps without STRAPP results
  estimates_all_steps <- p_values_all_steps <- rep(x = NA, times = length(time_steps))
  estimates_all_steps[time_steps_with_STRAPP_results_ID] <- estimates
  p_values_all_steps[time_steps_with_STRAPP_results_ID] <- p_values

  pvalues_summary_df <- data.frame(focal_time = time_steps,
                                   estimate = estimates_all_steps,
                                   p_value = p_values_all_steps)

  # Store pvalues_summary_df in final output
  final_ouput$pvalues_summary_df <- pvalues_summary_df

  ## Store time steps in final output
  final_ouput$time_steps <- time_steps

  ## Store type of trait data tested
  final_ouput$trait_data_type <- trait_data_type

  ## Store the type of trait data used to select the appropriate statistical method
  trait_data_type_for_stats_list <- unlist(lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$STRAPP_results$trait_data_type_for_stats } ))
  trait_data_type_for_stats_list <- setdiff(trait_data_type_for_stats_list, "none")
  final_ouput$trait_data_type_for_stats <- unique(trait_data_type_for_stats_list)
  if (length(final_ouput$trait_data_type_for_stats) > 1) # Deal with cases with a mix of binary and multinominal tests due to the absence of states at some time-steps
  {
    final_ouput$trait_data_type_for_stats <- "multinominal"
  }

  ## Store type of diversification rates tested
  final_ouput$rate_type <- rate_type

  ## Store type of uncertainty_strategy used
  final_ouput$uncertainty_strategy <- uncertainty_strategy

  ## Record trait maps and BAMM samples chsoen for testing
  final_ouput$trait_maps_vs_BAMM_samples_list <- trait_maps_vs_BAMM_samples_list

  ## Extract p_values and estimates in summary_df for post hoc tests
  if (posthoc_pairwise_tests)
  {
    # Extract summary_df from $STRAPP_results
    pvalues_summary_df_for_posthoc_pairwise_tests <- lapply(
      X = deepSTRAPP_outputs_over_time,
      FUN = function (x) { x$STRAPP_results$posthoc_pairwise_tests$summary_df } )

    # Detect empty slots (because less than two levels where present at focal time)
    time_steps_with_posthoc_tests_ID <- which(!unlist(lapply(X = pvalues_summary_df_for_posthoc_pairwise_tests, FUN = is.null)))

    # Add focal_time and rename columns to match with $pvalues_summary_df
    for (i in seq_along(time_steps_with_posthoc_tests_ID))
    {
      pvalues_summary_df_i <- pvalues_summary_df_for_posthoc_pairwise_tests[[i]]
      pvalues_summary_df_i$focal_time <- deepSTRAPP_outputs_over_time[[i]]$STRAPP_results$focal_time
      pvalues_summary_df_i <- pvalues_summary_df_i[, c("focal_time", "pairs", "estimates", "p_values", "p_values_adjusted")]
      names(pvalues_summary_df_i) <- c("focal_time", "pair", "estimate", "p_value", "p_value_adjusted")
      pvalues_summary_df_for_posthoc_pairwise_tests[[i]] <- pvalues_summary_df_i
    }
    # Bind across all time steps
    pvalues_summary_df_for_posthoc_pairwise_tests <- do.call(what = rbind, args = pvalues_summary_df_for_posthoc_pairwise_tests)

    # Store pvalues_summary_df in final output
    final_ouput$pvalues_summary_df_for_posthoc_pairwise_tests <- pvalues_summary_df_for_posthoc_pairwise_tests
  }

  ## Extract melted data.frame with trait data if requested
  # Can be used to plot rates through time ~ trait values
  if (extract_trait_data_melted_df)
  {
    trait_data_over_time <- lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$trait_data_df } )
    # Bind across all time steps
    trait_data_df_over_time <- do.call(what = rbind, args = trait_data_over_time)

    # Store melted data.frame with trait data in final output
    final_ouput$trait_data_df_over_time <- trait_data_df_over_time
  }

  ## Extract melted data.frame with rates/regimes data if requested
  # Can be used to plot rates through time ~ trait values
  if (extract_diversification_data_melted_df)
  {
    diversification_data_over_time <- lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$diversification_data_df } )
    # Bind across all time steps
    diversification_data_df_over_time <- do.call(what = rbind, args = diversification_data_over_time)

    # Store melted data.frame with rates/regimes data in final output
    final_ouput$diversification_data_df_over_time <- diversification_data_df_over_time
  }

  ## Extract STRAPP_results if requested
  # Can be used to plot histograms of STRAPP test statistics for each focal_time
  if (return_STRAPP_results)
  {
    STRAPP_results_over_time <- lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$STRAPP_results } )

    # Store STRAPP_results in final output
    final_ouput$STRAPP_results_over_time <- STRAPP_results_over_time
  }

  ## Extract updated Maps if requested
  # Can be used to plot an updated mapped phylogenies for each focal_time
  if (return_updated_Maps)
  {
    updated_Maps_over_time <- lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$updated_Maps } )

    # Store updated contMap(s)/densityMaps/simmaps in final output
    final_ouput$updated_Maps_over_time <- updated_Maps_over_time
  }

  ## Extract updated BAMM_object with associated rates and regimes if requested
  # Can be used to plot diversification rates on updated an phylogeny for each focal_time
  if (return_updated_BAMM_objects)
  {
    updated_BAMM_objects_over_time <- lapply(X = deepSTRAPP_outputs_over_time, FUN = function (x) { x$updated_BAMM_object } )

    # Store updated BAMM objects in final output
    final_ouput$updated_BAMM_objects_over_time <- updated_BAMM_objects_over_time
  }

  ## Export final output
  return(final_ouput)

}

## Adjust all examples/tutorials in run_deepSTRAPP to include the new uncertainty_strategy argument and set it as non-default if needed !!!!! #############

##### Update all tutorials accordingly to advertise the multiple options!
# Clarify what inputs are available for what strategy
# Possible to provide anything (contMap, contMaps; densityMaps, simmaps) for "rates_only"
# Requires contMaps, densityMaps (but no trackdown of stochastic map ID), simmaps for "paired' and "full" strategies
# Advertise changes as modularity in the manuscript and response to reviewers

## Make a tutorial dedicated to loading external results into deepSTRAPP
 # contMap / contMaps / densityMaps / simmaps
 # BAMM object with build_BAMM_object
 # Add reference to it into the vignettes/tutorials for each trait evolution / BAMM

## Once tutorials are ready, push the update on GitHub and ask Bonnie for review/testing

# Then work the second example to include


