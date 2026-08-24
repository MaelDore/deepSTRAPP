## Functions to compute STRAPP tests
# One master function to prepare data and select the proper test function according to data type
# There sub-functions carrying out tests according to data type

#' @title Compute STRAPP to test for a relationship between diversification rates and trait data
#'
#' @description Carries out the appropriate statistical method to test for a relationship between
#'   diversification rates and trait data for a given point in the past (i.e. the `focal_time`).
#'   Tests are based on block-permutations: rates data are randomized across tips following blocks
#'   defined by the diversification regimes identified on each tip (typically from a BAMM).
#'
#'   Such tests are called STructured RAte Permutations on Phylogenies (STRAPP) as described in
#'   Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification.
#'   Systematic biology, 65(2), 181-193. \doi{10.1093/sysbio/syv066}.
#'
#'   The function is an extension of the original [BAMMtools::traitDependentBAMM()] function used to
#'   carry out STRAPP test on extant time-calibrated phylogenies.
#'
#'   Tests can be carried out on speciation, extinction and net diversification rates.
#'
#'   `deepSTRAPP::compute_STRAPP_test_for_focal_time()` can handle three types of statistical tests depending on the type of trait data provided:
#'
#'   ## Continuous trait data
#'
#'   Tests for correlations between trait and rates carried out with `deepSTRAPP::compute_STRAPP_test_for_continuous_data()`.
#'   The associated test is the Spearman's rank correlation test (See [stats::cor.test]).
#'
#'   ## Binary trait data
#'
#'   For categorical and biogeographic trait data that have only two states (ex: 'Nearctic' vs. 'Neotropics').
#'   Tests for differences in rates between states are carried out with `deepSTRAPP::compute_STRAPP_test_for_binary_data()`.
#'   The associated test is the Mann-Whitney-Wilcoxon rank-sum test (See [stats::wilcox.test]).
#'
#'   ## Multinominal trait data
#'
#'   For categorical and biogeographic trait data with more than two states (ex: 'No leg' vs. 'Two legs' vs. 'Four legs').
#'   Tests for differences in rates between states are carried out with `deepSTRAPP::compute_STRAPP_test_for_multinominal_data()`.
#'   The associated test for all states is the Kruskal-Wallis H test (See [stats::kruskal.test]).
#'   If `posthoc_pairwise_tests = TRUE`, post hoc pairwise tests between pairs of states will be carried out too.
#'   The associated test for post hoc pairwise tests is the Dunn's post hoc pairwise rank-sum test (See [dunn.test::dunn.test]).
#'
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'   that contains a phylogenetic tree and associated diversification rates
#'   across selected posterior samples updated to a specific time in the past (i.e. the `focal_time`).
#' @param trait_data_list List obtained from [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'   or [deepSTRAPP::extract_all_trait_values_for_focal_time()] that contains at least
#'   a `$trait_data` element, a `$focal_time` element, and a `$trait_data_type`.
#'   `$trait_data` is a named vector with the trait data found on the phylogeny at `focal_time`.
#'   `$focal_time` informs on the time in the past at which the trait and rates data will be tested.
#'   `$trait_data_type` informs on the type of trait data: 'continuous', 'categorical', or 'biogeographic'.
#' @param rate_type A character string specifying the type of diversification rates to use. Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param uncertainty_strategy Character string. To select the strategy used to account for uncertainty in estimates.
#'   * `"rates_only"`: Only accounts for diversification-rate uncertainty across BAMM posterior samples. Uses ML estimates for continuous traits
#'                     and the most frequent state/range observed across stochastic maps for categorical and biogeographic data.
#'   * `"paired"`: Default option. Accounts for both diversification-rate and ancestral trait/range reconstruction uncertainty by pairing BAMM posterior samples with stochastic maps.
#'                 When the number of BAMM samples and stochastic maps differ, random pairing with replacement from the smaller set is used
#'                 so that all posterior samples and stochastic maps contribute to the analysis.
#'   * `"full"`: Exhaustive option that accounts for trait/range- and rate- uncertainty by crossing all BAMM posterior samples with all stochastic maps.
#'               Accounts for both diversification-rate and ancestral reconstruction uncertainty by evaluating every combination of BAMM posterior sample and stochastic map.
#'               WARNING: This exhaustive approach can substantially increase computation time and memory requirements and is therefore recommended only for moderate-sized analyses.
#' @param trait_maps_vs_BAMM_samples_list (Optional) Data.frame of two variables manually providing the names to associate stochastic maps (`$trait_map_ID`)
#'   with BAMM samples (`$BAMM_posterior_sample_ID`). This is typically used to ensure the same stochastic maps and BAMM samples are used to test across multiple time-steps.
#'   Values are the names of the objects such as "Map_X" and "BAMM_X". Default = `NULL`.
#'     * For uncertainty_strategy == 'rates_only', the `$trait_map_ID` must be "Map_ML" as only the ML estimates of trait values/states/ranges are used.
#'     * For uncertainty_strategy == 'paired', each pair of stochastic map and BAMM sample will be used once.
#'     * For uncertainty_strategy == 'full', all stochastic maps will be matched with all BAMM samples.
#'   Those may partly differ from the actual maps and BAMM samples used for the test, as recorded in `$perm_data_df`, because invalid maps with not enough states/ranges are discarded.
#' @param seed Integer. Set the seed to ensure reproducibility. Default is `NULL` (a random seed is used).
#' @param nb_permutations Integer. To select the number of random permutations to perform during the tests.
#'   If NULL (default), all BAMM posterior samples will be used once.
#' @param alpha Numerical. Significance level to use to compute the `estimate` corresponding to the values of the test statistic used to assess significance of the test. This does NOT affect p-values. Default is `0.05`.
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
#' @param posthoc_pairwise_tests Logical. Only for multinominal data (with more than two states). If `TRUE`, all possible post hoc pairwise (Dunn) tests will be computed across all pairs of states.
#'   This is a way to detect which pairs of states have significant differences in rates if the overall test (Kruskal-Wallis) is significant. Default is `FALSE`.
#' @param p.adjust_method A character string. Only for multinominal data (with more than two states). It specifies the type of correction to apply to the p-values
#'  in the post hoc pairwise tests to account for multiple comparisons. See [stats::p.adjust()] for the available methods. Default is `none`.
#' @param return_perm_data Logical. Whether to return the stats data computed from the posterior samples for observed and permuted data in the output.
#'  This is needed to plot the histogram of the null distribution used to assess significance of the test with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'  Default is `FALSE`.
#' @param nthreads Integer. Number of threads to use for paralleled computing of the tests across the permutations. The R package `parallel` must be loaded for `nthreads > 1`. Default is `1`.
#' @param print_hypothesis Logical. Whether to print information on what test is carried out, detailing the null and alternative hypotheses,
#'  what significant level is used to rejected or not the null hypothesis, and how uncertainty in trait estimates is handled. Default is `TRUE`.
#'
#' @export
#' @importFrom stats wilcox.test cor.test kruskal.test p.adjust median qchisq qnorm quantile sd
#' @importFrom dunn.test dunn.test
#' @importFrom utils capture.output
#'
#' @details These set of functions carries out the STructured RAte Permutations on Phylogenies (STRAPP) test as defined in
#'   Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification.
#'   Systematic biology, 65(2), 181-193.
#'
#'   It is an extension of the original [BAMMtools::traitDependentBAMM()] function used to
#'   carry out STRAPP test on extant time-calibrated phylogenies, but allowing here to test for
#'   differences/correlations at any point in the past (i.e. the `focal_time`).
#'
#'   It takes an object of class `"bammdata"` (`BAMM_object`) that was updated such as
#'   its diversification rates (`$tipLambda` and `$tipMu`) and regimes (`$tipStates`) are reflecting
#'   values observed at at a specific time in the past (i.e. the `$focal_time`).
#'   Similarly, it takes a list (`trait_data_list`) that provides `$trait_data` as observed on branches
#'   at the same `focal_time` than the diversification rates and regimes.
#'
#'   A STRAPP test is carried out by drawing a random set of posterior samples from the `BAMM_object`, then randomly permuting rates
#'   across blocks of tips defined by the macroevolutionary regimes. Test statistics are then computed across the initial observed data
#'   and the permuted data for each sample.
#'   In a two-tailed test, the p-value is the proportion of posterior samples in which the test stats is as extreme in the permuted than in the observed data.
#'   In a one-tailed test, the p-value is the proportion of posterior samples in which the test stats is higher in the permuted than in the observed data.
#'
#'   ----------  Major changes compared to [BAMMtools::traitDependentBAMM()]  ----------
#'
#'   * Allow to account for uncertainty in trait estimates by pairing/mapping multiple trait data extracted from stochastic maps
#'     across the BAMM samples. See the `uncertainty_strategy` argument for details.
#'   * Add post hoc pairwise tests (Dunn test) for multinominal data. Use `posthoc_pairwise_tests = TRUE`.
#'   * Provide outputs tailored for histogram plots [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()]
#'     and p-value time-series plots [deepSTRAPP::plot_STRAPP_pvalues_over_time()].
#'   * Add prints detailing what test is carried out, what are the null and alternative hypotheses,
#'     and what significant level is used to rejected or not the null hypothesis. (Enabled with `print_hypothesis = TRUE`).
#'   * Split the function in multiple sub-functions according to the type of data (`$trait_data_type`).
#'   * Prevent using Pearson's correlation tests and applying log-transformation for continuous data.
#'     The rationale is that there is no reason to assume that tip rates are distributed normally or log-normally.
#'     Thus, a Spearman's rank correlation test is favored.
#'
#' @return The function returns a list with at least eleven elements.
#'
#'   Summary elements for the main test:
#'   * `$estimate` Named numerical. Value of the test statistic used to assess significance of the test
#'   according to the significance level provided (`alpha`). The test is significant if `$estimate` is higher than zero.
#'   * `$stats_median` Numerical. Median value of the distribution of test statistics.
#'   * `$nb_test_stats` Integer. Number of test stats in the distribution.
#'   * `$p-value` Numerical. P-value of the test. The test is considered significant if `$p-value` is lower than `alpha`.
#'   * `$method` Character string. The statistical method used to carry out the test.
#'   * `$rate_type` Character string. The type of diversification rates tested. One of 'speciation', 'extinction' or 'net_diversification'.
#'   * `$trait_data_type` Character string. The type of trait data as found in 'trait_data_list$trait_data_type'. One of 'continuous', 'categorical', or 'biogeographic'.
#'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. One of 'continuous', 'binary', or 'multinominal'.
#'   * `$uncertainty_strategy` Character string. The strategy used to account for uncertainty in estimates.
#'   * `$trait_maps_vs_BAMM_samples_list` List of two elements recording the stochastic maps (`$trait_map_ID`) and BAMM samples (`$BAMM_posterior_sample_ID`) chosen for testing.
#'      Those may partly differ from the actual maps and BAMM samples used for the test as recorded in `$perm_data_df` because invalid maps with not enough states/ranges are discarded.
#'   * `$focal_time` The time in the past at which the trait and rates data were tested.
#'
#'   If using continuous or binary data:
#'   * `$two-tailed` Logical. Record the type of test used: two-tailed if `TRUE`, one-tailed if `FALSE`.
#'   If `one_tailed_hypothesis` is provided (only for continuous and binary trait data):
#'   * `$one_tailed_hypothesis` Character string. Record of the alternative hypothesis used for the one-tailed tests.
#'
#'   If `posthoc_pairwise_tests = TRUE` (only for multinomial trait data):
#'   * `$posthoc_pairwise_tests` List of at least 3 sub-elements:
#'     + `$summary_df` Data.frame of six variables providing the summary results of post hoc pairwise tests
#'     + `$method` Character string. The statistical method used to carry out the test. Here, "Dunn".
#'     + `$two-tailed` Logical. Record the type of post hoc pairwise tests used: two-tailed if `TRUE`, one-tailed if `FALSE`.
#'
#'   If `return_perm_data = TRUE`, the stats data computed from the posterior samples for observed and permuted data are provided.
#'   This is needed to plot the histogram of the null distribution used to assess significance of the test with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'   * `$perm_data_df` A data.frame with four variables summarizing the data generated during the STRAPP test:
#'     + `$trait_map_ID` Integer. ID of the stochastic map from which trait data were extracted and used for the STRAPP test.
#'     + `$BAMM_posterior_sample_ID` Integer. ID of the posterior samples randomly drawn and used for the STRAPP test.
#'     + `$*_obs` Numerical. Test stats computed from the observed data in the posterior samples. Name depends on the test used.
#'     + `$*_perm` Numerical. Test stats computed from the permuted data in the posterior samples. Name depends on the test used.
#'     + `$delta_*` OR `$abs_delta_*` Numerical. Test stats computed for the STRAPP test comparing observed stats and permuted stats.
#'       Name depends on the test used and the type of tests (two-tailed compare absolute values; one-tailed compare raw values).
#'   Combined with `posthoc_pairwise_tests = TRUE`, the stats data are also provided for the post hoc pairwise tests:
#'   * `$posthoc_pairwise_tests$perm_data_array` A 3D array containing stats data for all post hoc pairwise tests in a similar format that `$perm_data_df`.
#'
#'   If no STRAPP test was performed in the case of categorical/biogeographic data with a single state/range at `focal_time`,
#'   only the `$trait_data_type`, `$trait_data_type_for_stats` = "none", and `$focal_time` are returned.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#'
#'   Original function in BAMMtools: [BAMMtools::traitDependentBAMM()]
#'
#'   Statistical tests: [stats::cor.test()] [stats::wilcox.test()] [stats::kruskal.test()] [dunn.test::dunn.test()]
#'
#'   For a guided tutorial, see this vignette: \code{vignette("explore_STRAPP_test_types", package = "deepSTRAPP")}
#'
#' @references For STRAPP: Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification.
#'   Systematic biology, 65(2), 181-193. \doi{10.1093/sysbio/syv066}.
#'
#'   For STRAPP in deep times: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'   Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity, Nature Communications.
#'   \doi{10.1038/s41467-025-63709-3}
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  # ------ Prepare data ------ #
#'
#'  ## Load the BAMM_object summarizing 1000 posterior samples of BAMM with diversification rates
#'  # for ponerine ants extracted for 10My ago.
#'  data(Ponerinae_BAMM_object_10My, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#' # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  # Plot the associated phylogeny with mapped rates
#'  plot_BAMM_rates(Ponerinae_BAMM_object_10My)
#'
#'  ## Load the object containing head width trait data for ponerine ants extracted for 10My ago.
#'  data(Ponerinae_trait_cont_tip_data_10My, package = "deepSTRAPP")
#'
#'  # Plot the associated contMap (continuous trait stochastic map)
#'  plot_contMap(Ponerinae_trait_cont_tip_data_10My$contMap)
#'
#'  # Check that objects are ordered in the same fashion
#'  identical(names(Ponerinae_BAMM_object_10My$tipStates[[1]]),
#'            names(Ponerinae_trait_cont_tip_data_10My$trait_data))
#'
#'  # Save continuous data
#'  trait_data_continuous <- Ponerinae_trait_cont_tip_data_10My
#'
#'  ## Transform trait data into binary and multinominal data
#'
#'  # Binarize data into two states
#'  trait_data_binary <- trait_data_continuous
#'  trait_data_binary$trait_data[trait_data_continuous$trait_data < 0.5] <- "state_A"
#'  trait_data_binary$trait_data[trait_data_continuous$trait_data >= 0.5] <- "state_B"
#'  trait_data_binary$trait_data_type <- "categorical"
#'
#'  table(trait_data_binary$trait_data)
#'
#'  # Categorize data into three states
#'  trait_data_multinominal <- trait_data_continuous
#'  trait_data_multinominal$trait_data[trait_data_continuous$trait_data < 0.6] <- "state_B"
#'  trait_data_multinominal$trait_data[trait_data_continuous$trait_data < 0.4] <- "state_A"
#'  trait_data_multinominal$trait_data[trait_data_continuous$trait_data >= 0.6] <- "state_C"
#'  trait_data_multinominal$trait_data_type <- "categorical"
#'
#'  table(trait_data_multinominal$trait_data)
#'
#'  ## Duplicate trait data as if extracted from multiple stochastic maps
#'
#'  trait_data_continuous_multimaps <- Ponerinae_trait_cont_tip_data_10My
#'  trait_data_initial <- trait_data_continuous_multimaps$trait_data
#'
#'  trait_data_multimaps <- list()
#'  for (i in 1:10)
#'  {
#'    # Add a bit of randomness across all stochastic maps data
#'    trait_data_multimaps[[i]] <- trait_data_initial +
#'        rnorm(n = length(trait_data_initial), mean = 0, sd = sd(trait_data_initial) / 10)
#'  }
#'  names(trait_data_multimaps) <- paste0("Map_", 1:10)
#'  trait_data_continuous_multimaps$trait_data <- trait_data_multimaps
#'
#'
#'  \donttest{ # (May take several minutes to run)
#'  # ------ Compute STRAPP test for continuous data ------ #
#'
#'  plot(x = trait_data_continuous$trait_data, y = Ponerinae_BAMM_object_10My$tipLambda[[1]])
#'
#'  # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
#'  # between "net_diversification" rates and trait data
#'  STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'     BAMM_object = Ponerinae_BAMM_object_10My,
#'     trait_data_list = trait_data_continuous,
#'     uncertainty_strategy = "rates_only",
#'     two_tailed = FALSE,
#'     one_tailed_hypothesis = "negative",
#'     return_perm_data = TRUE)
#'  str(STRAPP_results, max.level = 2)
#'  # Data from the posterior samples is available in STRAPP_results$perm_data_df
#'  head(STRAPP_results$perm_data_df)
#'
#'  # ------ Compute STRAPP test for binary data ------ #
#'
#'  # Compute STRAPP test under the alternative hypothesis that "state_A" is associated
#'  # with higher "net_diversification" that "state_B"
#'  STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'     BAMM_object = Ponerinae_BAMM_object_10My,
#'     trait_data_list = trait_data_binary,
#'     uncertainty_strategy = "rates_only",
#'     two_tailed = FALSE,
#'     one_tailed_hypothesis = c("state_A > state_B"))
#'  str(STRAPP_results, max.level = 1)
#'
#'  # Compute STRAPP test under the alternative hypothesis that "state_B" is associated
#'  # with higher "net_diversification" that "state_A"
#'  STRAPP_results <- compute_STRAPP_test_for_focal_time(BAMM_object = Ponerinae_BAMM_object_10My,
#'     trait_data_list = trait_data_binary,
#'     uncertainty_strategy = "rates_only",
#'     two_tailed = FALSE,
#'     one_tailed_hypothesis = c("state_B > state_A"))
#'  str(STRAPP_results, max.level = 1)
#'
#'  # ------ Compute STRAPP test for multinominal data ------ #
#'
#'  # Compute STRAPP test between all three states, and compute post hoc tests
#'  # for differences in rates between all possible pairs of states
#'  # with a p-value adjusted for multiple comparison using Bonferroni's correction
#'  STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'     BAMM_object = Ponerinae_BAMM_object_10My,
#'     trait_data_list = trait_data_multinominal,
#'     uncertainty_strategy = "rates_only",
#'     posthoc_pairwise_tests = TRUE,
#'     two_tailed = TRUE,
#'     p.adjust_method = "bonferroni",
#'     return_perm_data = TRUE)
#'  str(STRAPP_results, max.level = 2)
#'  # All post hoc pairwise test summaries are available in $summary_df
#'  STRAPP_results$posthoc_pairwise_tests$summary_df
#'
#'  # ------ Compute STRAPP test with the 'paired' strategy ------ #
#'
#'  # Account for uncertainty in trait estimates
#'  # by pairing stochastic maps with BAMM posterior samples
#'
#'  # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
#'  # between "net_diversification" rates and trait data
#'  STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_continuous_multimaps,
#'    uncertainty_strategy = "paired",
#'    two_tailed = FALSE,
#'    one_tailed_hypothesis = "negative",
#'    return_perm_data = TRUE)
#'  str(STRAPP_results, max.level = 2)
#'  # Data from the paired stochastic maps and BAMM posterior samples
#'  # is available in STRAPP_results$perm_data_df
#'  head(STRAPP_results$perm_data_df)
#'  # Tests were performed across 1000 iterations since the stochastic maps
#'  # have been paired with the 1000 BAMM samples
#'  nrow(STRAPP_results$perm_data_df)
#'  # Each of the 10 maps was paired ca. 100 times with a unique BAMM sample
#'  table(STRAPP_results$perm_data_df$trait_map_ID)
#'
#'  # ------ Compute STRAPP test with the 'full' strategy ------ #
#'
#'  # Account for uncertainty in trait estimates
#'  # by crossing stochastic maps with BAMM posterior samples
#'
#'  # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
#'  # between "net_diversification" rates and trait data
#'  STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_continuous_multimaps,
#'    uncertainty_strategy = "full",
#'    two_tailed = FALSE,
#'    one_tailed_hypothesis = "negative",
#'    return_perm_data = TRUE)
#'  str(STRAPP_results, max.level = 2)
#'  # Data from the combined stochastic maps X BAMM posterior samples
#'  # is available in STRAPP_results$perm_data_df
#'  head(STRAPP_results$perm_data_df)
#'  # Tests were performed across 10000 iterations since the 10 stochastic maps
#'  # have been combined with the 1000 BAMM samples
#'  nrow(STRAPP_results$perm_data_df)
#'  # Each of the 10 maps was paired with all 1000 BAMM samples
#'  table(STRAPP_results$perm_data_df$trait_map_ID)
#'  }
#' }
#'


### Master function to prepare data and select the proper test function according to data type ####

compute_STRAPP_test_for_focal_time <- function (BAMM_object, trait_data_list,
                                                rate_type = "net_diversification",
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
                                                print_hypothesis = TRUE)
{

  ### Check input validity
  {
    ## Check if trait data is recorded for ML estimates (as a named vector) or across multiple stochastic maps (as a list of named vectors)
    if (is.list(trait_data_list$trait_data))
    {
      trait_data_is_ML_estimates <- FALSE
    } else {
      trait_data_is_ML_estimates <- TRUE
    }

    ## BAMM_object
    # BAMM_object must be a 'bammdata' object
    if (!("bammdata" %in% class(BAMM_object)))
    {
      stop("'BAMM_object' must have the 'bammdata' class. See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects.")
    }
    # Number of posterior sample data must be equal between $tipStates, $tipLambda and $tipMu
    posterior_samples_length <- c(length(BAMM_object$tipStates), length(BAMM_object$tipLambda), length(BAMM_object$tipMu))
    if (length(unique(posterior_samples_length)) != 1)
    {
      stop("Number of posterior samples in 'BAMM_object' must be equal between $tipStates, $tipLambda and $tipMu.\n",
           "Please check the structure of your 'BAMM_object' with str(BAMM_object, 1).\n",
           "See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects.")
    }
    # Number of branches in each posterior sample must be equal within $tipStates, $tipLambda and $tipMu
    tipStates_data_length <- unlist(lapply(X = BAMM_object$tipStates, FUN = length))
    if (length(unique(tipStates_data_length)) != 1)
    {
      stop("Number of branches in each posterior sample of 'BAMM_object$tipStates' must be equal.\n",
           "Please check the structure of your 'BAMM_object' with str(BAMM_object$tipStates, 1).\n",
           "See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects.")
    }
    tipLambda_data_length <- unlist(lapply(X = BAMM_object$tipLambda, FUN = length))
    if (length(unique(tipLambda_data_length)) != 1)
    {
      stop("Number of branches in each posterior sample of 'BAMM_object$tipLambda' must be equal.\n",
           "Please check the structure of your 'BAMM_object' with str(BAMM_object$tipLambda, 1)\n",
           "See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects.")
    }
    tipMu_data_length <- unlist(lapply(X = BAMM_object$tipMu, FUN = length))
    if (length(unique(tipMu_data_length)) != 1)
    {
      stop("Number of branches in each posterior sample of 'BAMM_object$tipMu' must be equal.\n",
           "Please check the structure of your 'BAMM_object' with str(BAMM_object$tipMu, 1)\n",
           "See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects.")
    }
    # Number of branches in each posterior sample must be equal between $tipStates, $tipLambda and $tipMu
    posterior_samples_data_length <- c(unique(tipStates_data_length), unique(tipLambda_data_length), unique(tipMu_data_length))
    if (length(unique(posterior_samples_data_length)) != 1)
    {
      stop(paste0("Number of branches in posterior samples of 'BAMM_object$tipMu', 'BAMM_object$tipLambda', and 'BAMM_object$tipMu' must be equal.\n",
                  "There respective number of branches is: ",paste(posterior_samples_data_length, collapse = ", "),".\n",
                  "Please check the structure of your 'BAMM_object' with str(BAMM_object, 2)\n",
                  "See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects."))
    }

    ## trait_data_list
    # trait_data_list must be a list with $trait_data and $trait_data_type
    if (is.null(trait_data_list$trait_data) | is.null(trait_data_list$trait_data_type))
    {
      stop("'trait_data_list' must be a list with $trait_data and $trait_data_type elements.")
    }

    # trait_data_list$trait_data_type can only be "continuous", categorical" or "biogeographic"
    if (!(trait_data_list$trait_data_type %in% c("continuous", "categorical", "biogeographic")))
    {
      stop("'trait_data_list$trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
    }

    ## Check trait_data_list$trait_data for ML_estimates
    if (trait_data_is_ML_estimates)
    {
      # trait_data_list$trait_data must be a named vector (can be numerical or character string)
      if (is.null(names(trait_data_list$trait_data)))
      {
        stop(paste0("'trait_data_list$trait_data' must be a named vector with names matching those found in BAMM_object$tipStates, BAMM_object$tipLambda, and BAMM_object$tipMu.\n",
                    "Names are either tip.label or tipward_node_ID of the branches cut at 'trait_data_list$focal_time' with deepSTRAPP::extract_most_likely_trait_values_for_focal_time."))
      }

      # Check compatibility with $trait_data_type
      if (trait_data_list$trait_data_type == "continuous" & !is.numeric(trait_data_list$trait_data))
      {
        stop("'trait_data_list$trait_data' must be numeric if 'trait_data_list$trait_data_type' is 'continuous'.")
      }
      if (trait_data_list$trait_data_type == "categorical" & !is.character(trait_data_list$trait_data))
      {
        stop("'trait_data_list$trait_data' must be a vector of character strings if 'trait_data_list$trait_data_type' is 'categorical'.")
      }
      if (trait_data_list$trait_data_type == "biogeographic" & !is.character(trait_data_list$trait_data))
      {
        stop("'trait_data_list$trait_data' must be a vector of character strings if 'trait_data_list$trait_data_type' is 'biogeographic'.")
      }

      # Length of $trait_data should match length of $tipStates, $tipLambda and $tipMu (for each posterior sample)
      if (length(trait_data_list$trait_data) != unique(posterior_samples_data_length))
      {
        stop("Number of branches in 'trait_data_list$trait_data' must be equal to number of branches in posterior samples in 'BAMM_object$tipStates', 'BAMM_object$tipLambda' and 'BAMM_object$tipMu'.\n",
        "Please check the structure of your 'BAMM_object' with str(BAMM_object, 2)")
      }
      # Names of $trait_data should match names in $tipStates, $tipLambda and $tipMu (for each posterior sample)
      if (!all(names(trait_data_list$trait_data) %in% names(BAMM_object$tipStates[[1]])))
      {
        stop("Names of 'trait_data_list$trait_data' should match names in 'BAMM_object$tipStates', 'BAMM_object$tipLambda' and 'BAMM_object$tipMu' (for each posterior sample).")
      }
      # Names of $trait_data should be ordered similarly as in in $tipStates, $tipLambda and $tipMu (for each posterior sample)
      if (!all(names(trait_data_list$trait_data) == names(BAMM_object$tipStates[[1]])))
      {
        warning(paste0("Branch data in 'trait_data_list$trait_data' should be ordered similarly as in 'BAMM_object$tipStates', 'BAMM_object$tipLambda' and 'BAMM_object$tipMu' (for each posterior sample).\n",
                       "This was not the case, likely because the initial phylogenies used to model trait evolution and diversification dynamics were not ordered in the same fashion.\n",
                       "See attr(x = phylo, which = 'order') to detect the order used in each phylogeny.\n",
                       "'trait_data_list$trait_data' have been reordered to match data in 'BAMM_object'."))
        trait_data_list$trait_data <- trait_data_list$trait_data[match(x = names(BAMM_object$tipStates[[1]]), table = names(trait_data_list$trait_data))]
      }
    }
    if (!trait_data_is_ML_estimates)
    {
      # trait_data_list$trait_data must be a named vector (can be numerical or character string)
      if (is.null(names(trait_data_list$trait_data[[1]])))
      {
        stop(paste0("'trait_data_list$trait_data' must be (a list of) named vector(s) with names matching those found in BAMM_object$tipStates, BAMM_object$tipLambda, and BAMM_object$tipMu.\n",
                    "Names are either tip.label or tipward_node_ID of the branches cut at 'trait_data_list$focal_time' with deepSTRAPP::extract_most_likely_trait_values_for_focal_time."))
      }

      # Check compatibility with $trait_data_type
      if (trait_data_list$trait_data_type == "continuous" & !all(unlist(lapply(X = trait_data_list$trait_data, FUN = is.numeric))))
      {
        stop("All elements in 'trait_data_list$trait_data' must be numeric if 'trait_data_list$trait_data_type' is 'continuous'.")
      }
      if (trait_data_list$trait_data_type == "categorical" & !all(unlist(lapply(X = trait_data_list$trait_data, FUN = is.character))))
      {
        stop("All elements in 'trait_data_list$trait_data' must be a vector of character strings if 'trait_data_list$trait_data_type' is 'categorical'.")
      }
      if (trait_data_list$trait_data_type == "biogeographic" & !all(unlist(lapply(X = trait_data_list$trait_data, FUN = is.character))))
      {
        stop("All elements in 'trait_data_list$trait_data' must be a vector of character strings if 'trait_data_list$trait_data_type' is 'biogeographic'.")
      }

      # Length of $trait_data should match length of $tipStates, $tipLambda and $tipMu (for each posterior sample)
      if (unique(unlist(lapply(X = trait_data_list$trait_data, FUN = length))) != unique(posterior_samples_data_length))
      {
        stop("Number of branches in elements of 'trait_data_list$trait_data' must be equal to number of branches in posterior samples in 'BAMM_object$tipStates', 'BAMM_object$tipLambda' and 'BAMM_object$tipMu'.\n",
        "Please check the structure of your 'BAMM_object' with str(BAMM_object, 2)")
      }
      # Names of $trait_data should match names in $tipStates, $tipLambda and $tipMu (for each posterior sample)
      if (!all(names(trait_data_list$trait_data[[1]]) %in% names(BAMM_object$tipStates[[1]])))
      {
        stop("Names in elements of 'trait_data_list$trait_data' should match names in 'BAMM_object$tipStates', 'BAMM_object$tipLambda' and 'BAMM_object$tipMu' (for each posterior sample).")
      }
      # Names of $trait_data should be ordered similarly as in in $tipStates, $tipLambda and $tipMu (for each posterior sample)
      if (!all(names(trait_data_list$trait_data[[1]]) == names(BAMM_object$tipStates[[1]])))
      {
        warning(paste0("Branch data in elements of 'trait_data_list$trait_data' should be ordered similarly as in 'BAMM_object$tipStates', 'BAMM_object$tipLambda' and 'BAMM_object$tipMu' (for each posterior sample).\n",
                       "This was not the case, likely because the initial phylogenies used to model trait evolution and diversification dynamics were not ordered in the same fashion.\n",
                       "See attr(x = phylo, which = 'order') to detect the order used in each phylogeny.\n",
                       "'trait_data_list$trait_data' have been reordered to match data in 'BAMM_object'."))
        for (i in seq_along(trait_data_list$trait_data))
        {
          trait_data_list$trait_data[[i]] <- trait_data_list$trait_data[[i]][match(x = names(BAMM_object$tipStates[[1]]), table = names(trait_data_list$trait_data[[i]]))]
        }
      }
    }

    ## focal_time
    # $focal_time in $trait_data_list and in BAMM_object must be equal.
    if (BAMM_object$focal_time != trait_data_list$focal_time)
    {
      stop(paste0("$focal_time should be the same in 'BAMM_object$focal_time' and 'trait_data_list$focal_time'.\n",
                  "You provided a 'BAMM_object' with '$focal_time' = ",BAMM_object$focal_time,".\n",
                  "You provided a 'trait_data_list' with '$focal_time' = ",trait_data_list$focal_time,"."))
    }

    ## rate_type must be either "speciation", "extinction" or "net_diversification"
    if (!(rate_type %in% c("speciation", "extinction", "net_diversification")))
    {
      stop("'rate_type' can only be 'speciation', 'extinction', or 'net_diversification'.")
    }

    ## uncertainty_strategy
    # uncertainty_strategy must be either "rates_only", "paired", or "full"
    if (!uncertainty_strategy %in% c("rates_only", "paired", "full"))
    {
      stop(paste0("'uncertainty_strategy' must be either 'rates_only', 'paired', or 'full'."))
    }
    # Ensure the trait_data available is compatible with the requested 'uncertainty_strategy'
    if (!trait_data_is_ML_estimates & (uncertainty_strategy == "rates_only"))
    {
      stop(paste0("You requested to compute a STRAPP test with the 'rates_only' strategy to account for uncertainty in estimates.\n",
                  "Yet you provided trait data recorded across multiple stochastic maps in 'trait_data_list$trait_data'.\n",
                  "For the 'rates_only' strategy, only ML estimates of trait values/states/ranges must be provided.\n",
                  "See ?deepSTRAPP::extract_most_likely_trait_values_for_focal_time() to learn how to obtain such data."))
    }
    if (trait_data_is_ML_estimates & (uncertainty_strategy %in% c("paired", "full")))
    {
      stop(paste0("You requested to compute a STRAPP test with the '",uncertainty_strategy,"' strategy to account for uncertainty in estimates.\n",
                  "Yet you only provided ML estimates of trait data (or trait data as extracted from a single stochastic map) in 'trait_data_list$trait_data'.\n",
                  "For the '",uncertainty_strategy,"' strategy, trait values/states/ranges as recorded across multiple stochastic maps must be provided.\n",
                  "See ?deepSTRAPP::extract_all_trait_values_for_focal_time() to learn how to obtain such data."))
    }

    ## trait_maps_vs_BAMM_samples_list
    if (!is.null(trait_maps_vs_BAMM_samples_list))
    {
      if (!all(names(trait_maps_vs_BAMM_samples_list) == c("trait_map_ID", "BAMM_posterior_sample_ID")))
      {
        stop("'trait_maps_vs_BAMM_samples_list' must be a list of two elements named 'trait_map_ID' and 'BAMM_posterior_sample_ID'.")
      }
      if (!all(grepl(pattern = "^BAMM_", x = trait_maps_vs_BAMM_samples_list$BAMM_posterior_sample_ID)))
      {
        stop("'trait_maps_vs_BAMM_samples_list$BAMM_posterior_sample_ID' must be a character vector with all values starting with 'BAMM_' to assign BAMM samples for testing.")
      }
      if (uncertainty_strategy == "rates_only")
      {
        if (!all(trait_maps_vs_BAMM_samples_list$trait_map_ID == "Map_ML"))
        {
          stop("For uncertainty_strategy = 'rates_only', 'trait_maps_vs_BAMM_samples_list$trait_map_ID' must be a unique character string = 'Map_ML' to reflect the use of the ML trait estimates for all tests.")
        }
      } else {
        if ((!all(grepl(pattern = "^Map_", x = trait_maps_vs_BAMM_samples_list$trait_map_ID))) & (!all(grepl(pattern = "^Dummy_map_", x = trait_maps_vs_BAMM_samples_list$trait_map_ID))))
        {
          stop(paste0("For uncertainty_strategy = ",uncertainty_strategy,", 'trait_maps_vs_BAMM_samples_list$trait_map_ID' must be a character vector with all values starting with 'Map_' or 'Dummy_map_' to assign trait stochastic maps for testing."))
        }
      }
      if ((uncertainty_strategy == "paired") & (length(trait_maps_vs_BAMM_samples_list$trait_map_ID) != length(trait_maps_vs_BAMM_samples_list$BAMM_posterior_sample_ID)))
      {
        stop(paste0("For uncertainty_strategy = 'paired', 'trait_maps_vs_BAMM_samples_list$trait_map_ID' and 'trait_maps_vs_BAMM_samples_list$BAMM_posterior_sample_ID' must be the same length, so they can be paired."))
      }
      if (!is.null(nb_permutations))
      {
        warning(paste0("'nb_permutations' is ignored when the pairing between stochastic maps and BAMM samples is provided manually with 'trait_maps_vs_BAMM_samples_list'."))
      }
    }

    ## seed
    if (!is.null(seed))
    {
      if (!is.numeric(seed))
      {
        stop(paste0("'seed' must be an integer."))
      }
    }

    ## nb_permutations
    if (is.null(nb_permutations) & is.null(trait_maps_vs_BAMM_samples_list))
    {
      # If NULL, set to the number of BAMM posterior samples
      nb_permutations <- length(BAMM_object$tipStates)
    }

    ## alpha
    # alpha must be set between 0 and 1.
    if ((alpha < 0) | (alpha > 1))
    {
      stop(paste0("'alpha' reflects the quantile used to extract 'estimate' values and assess significance of the test. It must be between 0 and 1.\n",
                  "Current value of 'alpha' is ",alpha,"."))
    }

    ## p.adjust_method. Check that it is one of the available option. See [stats::p.adjust()] for the available methods.
    if (!(p.adjust_method %in% c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")))
    {
      stop(paste0("'p.adjust_method' specifies the type of correction to apply to the p-values. See ?stats::p.adjust for the available methods.\n"))
    }

    ## nthreads
    if (nthreads > 1) {
      if (!"package:parallel" %in% search())
      {
        stop("Please load package 'parallel' for using the multi-thread option\n")
      }
    }
  }

  ## Set seed
  if (!is.null(seed))
  {
    set.seed(seed = seed)
  }

  ## Extract BAMM rates and regimes data, and assign BAMM_ID
  BAMM_data <- list(tipStates = BAMM_object$tipStates, tipLambda = BAMM_object$tipLambda, tipMu = BAMM_object$tipMu)
  for (i in seq_along(BAMM_data))
  {
    BAMM_data_i <- BAMM_data[[i]]
    names(BAMM_data_i) <- paste0("BAMM_", 1:length(BAMM_data_i))
    BAMM_data[[i]] <- BAMM_data_i
  }

  ## Extract trait data
  trait_data <- trait_data_list$trait_data

  ## Extract type of trait data
  trait_data_type <- trait_data_list$trait_data_type

  ## Filter data to keep only the designated Map x BAMM samples ('$trait_maps_vs_BAMM_samples_list')
  if (!is.null(trait_maps_vs_BAMM_samples_list))
  {
    # Filter BAMM samples
    BAMM_data <- lapply(X = BAMM_data, FUN = function (x) {x[trait_maps_vs_BAMM_samples_list$BAMM_posterior_sample_ID]} )
    BAMM_data$BAMM_posterior_sample_ID <- trait_maps_vs_BAMM_samples_list$BAMM_posterior_sample_ID

    # Filter trait maps
    if (!trait_data_is_ML_estimates)
    {
      trait_data <- trait_data[trait_maps_vs_BAMM_samples_list$trait_map_ID]
    }
  }

  ## Check test validity for case with only ML estimates trait data = a unique (stochastic) map
  if (trait_data_is_ML_estimates)
  {
    run_test <- TRUE

    ## If trait_data_type is "categorical" or "biogeographic", reclassify according to the number of states
    if (trait_data_type %in% c("categorical", "biogeographic"))
    {
      nb_levels <- nlevels(as.factor(trait_data))

      if (nb_levels == 1) # Case with a single state
      {
        # Send warning
        warning(paste0("There is only a single state/range found at focal time = ",BAMM_object$focal_time,": '",levels(as.factor(trait_data)),"'.\n",
                       "No STRAPP test for difference in rates can be performed with a single level.\n",
                       "'STRAPP_results' were provided without test results and 'trait_data_type_for_stats' was tagged as 'none'.\n"))

        ## Create empty STRAPP_results object with no test
        STRAPP_results <- list(trait_data_type = trait_data_type,
                               trait_data_type_for_stats = "none")
        run_test <- FALSE
      }

      if (nb_levels == 2) # Case with two states
      {
        trait_data_type_for_stats <- "binary"
      } else { # Case with more than two states
        trait_data_type_for_stats <- "multinominal"
      }
    } else {
      trait_data_type_for_stats <- "continuous"
    }

    ## posthoc_pairwise_tests = TRUE. Only makes sense if more than 2 states/ranges in categorical/biogeographic data
    if (trait_data_type_for_stats == "continuous" & posthoc_pairwise_tests == TRUE)
    {
      stop(paste0("'posthoc_pairwise_tests = TRUE' does not make sense for a continuous trait.\n",
                  "Please set 'posthoc_pairwise_tests = FALSE' or provide categorical/biogeographic data with more than two states/ranges."))
    }
    if (trait_data_type_for_stats == "binary" & posthoc_pairwise_tests == TRUE)
    {
      warning(paste0("There are only two states/ranges found at focal time = ",BAMM_object$focal_time,": ",paste(levels(as.factor(trait_data)), collapse = ", "),".\n",
                     "'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                     "If you want to test specific hypotheses with continuous or categorical binary data, use 'two_tailed = FALSE' and provide the 'one_tailed_hypothesis'.\n"))
    }
  }

  ## Check test validity for case with multiple trait data = multiple stochastic maps
  # Run test if at least one trait data have valid data, but send a warning to inform on the number of stochastic maps discarded

  if (!trait_data_is_ML_estimates)
  {
    run_test <- TRUE

    ## If trait_data_type is "categorical" or "biogeographic", reclassify according to the number of states
    if (trait_data_type %in% c("categorical", "biogeographic"))
    {
      nb_levels_list <- unlist(lapply(X = trait_data, FUN = function (x) { nlevels(as.factor(x)) } ))

      if (all(nb_levels_list == 1)) # Case with a single state across all maps
      {
        # Send warning
        warning(paste0("There is only a single state/range found at focal time = ",BAMM_object$focal_time," across all stochastic maps.\n",
                       "No STRAPP test for difference in rates can be performed with a single level.\n",
                       "'STRAPP_results' were provided without test results and 'trait_data_type_for_stats' was tagged as 'none'.\n"))

        ## Create empty STRAPP_results object with no test
        STRAPP_results <- list(trait_data_type = trait_data_type,
                               trait_data_type_for_stats = "none")
        run_test <- FALSE
      }

      ## Remove trait data from stochastic maps with a single level
      trait_data_to_remove_ID <- which(nb_levels_list == 1)
      if (length(trait_data_to_remove_ID) > 0)
      {
        # Send warning
        nb_maps_removed <- length(trait_data_to_remove_ID)
        nb_maps_total <- length(trait_data)
        nb_maps_kept <- nb_maps_total - nb_maps_removed
        perc_maps_kept <- round(nb_maps_kept/nb_maps_total*100, 1)
        warning(paste0("There are ",length(trait_data_to_remove_ID)," cases of stochastic maps with only a single state/range found at focal time = ",BAMM_object$focal_time,".\n",
                       "No STRAPP test for difference in rates can be performed with a single level.\n",
                       "Trait data from those maps have been removed, but STRAPP results must be interpreted with caution\n",
                       "as results are drawn from ",nb_maps_kept," out of ",nb_maps_total," = ",perc_maps_kept," % of the initial trait data.\n"))

        # Remove problematic trait data
        trait_data <- trait_data[-trait_data_to_remove_ID]
        # If uncertainty_strategy = 'paired', and pairing have been provided, we also need to remove the associated BAMM samples
        if ((uncertainty_strategy == "paired") & !is.null(trait_maps_vs_BAMM_samples_list))
        {
          BAMM_data_names_to_keep <- names(BAMM_data$tipStates)[setdiff(1:length(BAMM_data$tipStates), trait_data_to_remove_ID)]
          BAMM_data <- lapply(X = BAMM_data, FUN = function (x) {x[BAMM_data_names_to_keep]} )
          BAMM_data$BAMM_posterior_sample_ID <- BAMM_data_names_to_keep
        }
      }

      if (max(nb_levels_list) == 2) # Case with maximum two states
      {
        trait_data_type_for_stats <- "binary"
      } else { # Case with more than two states in at least one stochastic map
        trait_data_type_for_stats <- "multinominal"
      }
    } else {
      trait_data_type_for_stats <- "continuous"
    }

    ## posthoc_pairwise_tests = TRUE. Only makes sense if more than 2 states/ranges in categorical/biogeographic data
    if (trait_data_type_for_stats == "continuous" & posthoc_pairwise_tests == TRUE)
    {
      stop(paste0("'posthoc_pairwise_tests = TRUE' does not make sense for a continuous trait.\n",
                  "Please set 'posthoc_pairwise_tests = FALSE' or provide categorical/biogeographic data with more than two states/ranges."))
    }
    if (trait_data_type_for_stats == "binary" & posthoc_pairwise_tests == TRUE)
    {
      warning(paste0("There are only two states/ranges found at focal time = ",BAMM_object$focal_time," across all stochastic maps.\n",
                     "'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                     "If you want to test specific hypotheses with continuous or categorical binary data, use 'two_tailed = FALSE' and provide the 'one_tailed_hypothesis'.\n"))
    }
  }

  ### Assign stochastic maps to BAMM samples if association in not provided manually
  if (is.null(trait_maps_vs_BAMM_samples_list))
  {
    ## Filter BAMM_data to keep only the requested number

    # Case with just enough BAMM samples
    if (length(BAMM_data$tipStates) == nb_permutations)
    {
      BAMM_data$BAMM_posterior_sample_ID <- sample(x = names(BAMM_data$tipStates), size = nb_permutations, replace = FALSE)
    } else {
      # Case with not enough BAMM samples
      if (length(BAMM_data$tipStates) < nb_permutations)
      {
        # Replicate BAMM samples
        BAMM_posterior_sample_ID <- sample(x = names(BAMM_data$tipStates), size = nb_permutations, replace = TRUE)
        warning(paste0("There are only ",length(BAMM_data$tipStates)," BAMM posterior samples.\n",
                       "They were resampled with replacement to match the requested number of permutations for the test: ",nb_permutations,".\n"))

      }
      # Case with too many BAMM samples
      if (length(BAMM_data$tipStates) > nb_permutations)
      {
        # Replicate trait data
        BAMM_posterior_sample_ID <- sample(x = names(BAMM_data$tipStates), size = nb_permutations, replace = FALSE)
        warning(paste0("There are too many BAMM posterior samples: ",length(BAMM_data$tipStates),".\n",
                       "They were subsampled to match the requested number of permutations for the test: ",nb_permutations,".\n"))

      }
      BAMM_data <- lapply(X = BAMM_data, FUN = function (x) { x[BAMM_posterior_sample_ID] } )
      BAMM_data$BAMM_posterior_sample_ID <- BAMM_posterior_sample_ID
    }

    ## Filter trait data to keep only the requested number

    # For 'rates_only' strategy, there is only one trait data, so no need to resample
    # For 'full' strategy, all trait data are matched with all BAMM samples, so no need to resample

    Map_sample_ID <- names(trait_data)

    # For 'paired' strategy, apply sub-sampling/replication step to match nb_permutations only on the valid trait data / stochastic maps !
    if (uncertainty_strategy == "paired" & (length(trait_data) != nb_permutations))
    {
      # Case with not enough trait data to pair with BAMM samples
      if (length(trait_data) < nb_permutations)
      {
        # Replicate trait data
        Map_sample_ID <- sample(x = names(trait_data), size = nb_permutations, replace = TRUE)
        warning(paste0("There are only ",length(trait_data)," trait datasets valid for testing found at focal time = ",BAMM_object$focal_time," across all stochastic maps.\n",
                       "They were resampled with replacement to match the requested number of permutations/BAMM samples = ",nb_permutations,", under the 'paired' strategy for uncertainty estimates.\n"))

      }
      # Case with too many trait data to pair with BAMM samples
      if (length(trait_data) > nb_permutations)
      {
        # Replicate trait data
        Map_sample_ID <- sample(x = names(trait_data), size = nb_permutations, replace = FALSE)
        warning(paste0("There are ",length(trait_data)," trait datasets valid for testing found at focal time = ",BAMM_object$focal_time," across all stochastic maps.\n",
                       "Since you requested ",nb_permutations," permutations/BAMM samples, the stochastic maps have been subsampled to match BAMM samples, under the 'paired' strategy for uncertainty estimates.\n"))

      }
      trait_data <- trait_data[Map_sample_ID]
    }

    ## Build trait_maps_vs_BAMM_samples_list to record the association

    if (uncertainty_strategy == "rates_only")
    {
      trait_maps_vs_BAMM_samples_list <- list(trait_map_ID = "Map_ML", BAMM_posterior_sample_ID = BAMM_data$BAMM_posterior_sample_ID)
    }
    if (uncertainty_strategy == "paired")
    {
      trait_maps_vs_BAMM_samples_list <- list(trait_map_ID = Map_sample_ID, BAMM_posterior_sample_ID = BAMM_data$BAMM_posterior_sample_ID)
    }
    if (uncertainty_strategy == "full")
    {
      trait_maps_vs_BAMM_samples_list <- list(trait_map_ID = names(trait_data), BAMM_posterior_sample_ID = BAMM_data$BAMM_posterior_sample_ID)
    }
    # if (uncertainty_strategy == "full")
    # {
    #   trait_maps_vs_BAMM_samples_list <- expand.grid(trait_map_ID = names(trait_data), BAMM_posterior_sample_ID = BAMM_data$BAMM_posterior_sample_ID)
    # }
  }

  if (run_test)
  {
    ## Compute the appropriate internal function depending on the type of trait data

    switch(EXPR = trait_data_type_for_stats,
           continuous =   { # Case for continuous data
             # Stat test = Spearman's rank Rho test
             STRAPP_results <- compute_STRAPP_test_for_continuous_data(
               BAMM_data = BAMM_data,
               trait_data = trait_data,
               trait_data_type = trait_data_type,
               rate_type = rate_type,
               uncertainty_strategy = uncertainty_strategy,
               nb_permutations = nb_permutations,
               alpha = alpha,
               two_tailed = two_tailed,
               one_tailed_hypothesis = one_tailed_hypothesis,
               return_perm_data = return_perm_data,
               nthreads = nthreads,
               print_hypothesis = print_hypothesis)
           },
           binary =       { # Case for binary data (Special case of categorical/biogeographic data with only two states)
             # Stat test = Mann-Whitney U test
             STRAPP_results <- compute_STRAPP_test_for_binary_data(
               BAMM_data = BAMM_data,
               trait_data = trait_data,
               trait_data_type = trait_data_type,
               rate_type = rate_type,
               uncertainty_strategy = uncertainty_strategy,
               nb_permutations = nb_permutations,
               alpha = alpha,
               two_tailed = two_tailed,
               one_tailed_hypothesis = one_tailed_hypothesis,
               return_perm_data = return_perm_data,
               nthreads = nthreads,
               print_hypothesis = print_hypothesis)
           },
           multinominal = { # Case for multinominal data (Case of categorical/biogeographic data with more than two states)
             # Stat test = Kruskal-Wallis H test
             # Can define the post hoc pairwise tests to compute
             STRAPP_results <- compute_STRAPP_test_for_multinominal_data(
               BAMM_data = BAMM_data,
               trait_data = trait_data,
               trait_data_type = trait_data_type,
               rate_type = rate_type,
               uncertainty_strategy = uncertainty_strategy,
               nb_permutations = nb_permutations,
               alpha = alpha,
               posthoc_pairwise_tests = posthoc_pairwise_tests,
               two_tailed = two_tailed,
               p.adjust_method = p.adjust_method,
               return_perm_data = return_perm_data,
               nthreads = nthreads,
               print_hypothesis = print_hypothesis)
           }
    )
  }

  ## Include selected stochastic maps X BAMM samples for testing in the output
  # Those may differ from the actual maps X BAMM samples used for the test as recorded in '$perm_data_df' because invalid maps
  # with not enough states will be discarded.
  STRAPP_results$trait_maps_vs_BAMM_samples_list <- trait_maps_vs_BAMM_samples_list

  ## Include focal_time in the output
  STRAPP_results$focal_time <- BAMM_object$focal_time

  ## Export the STRAPP test output
  return(STRAPP_results)
}


### Sub-function to handle continuous data ####

compute_STRAPP_test_for_continuous_data <- function (
    BAMM_data, trait_data,
    trait_data_type = "continuous",
    rate_type = "net_diversification",
    uncertainty_strategy = "paired",
    nb_permutations = NULL,
    alpha = 0.05,
    two_tailed = TRUE,
    one_tailed_hypothesis = NULL,
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{
  ### Check input validity
  {
    ## one_tailed_hypothesis
    # For continuous data, it is either "negative" or "positive" correlation.
    if (!is.null(one_tailed_hypothesis))
    {
      if (!(one_tailed_hypothesis %in% c("negative", "positive")))
      {
        stop(paste0("The 'one_tailed_hypothesis' must be either 'negative' or 'positive' for continuous trait data."))
      }
    }

    ## two_tailed & one_tailed_hypothesis
    if (!two_tailed & is.null(one_tailed_hypothesis))
    {
      stop(paste0("You selected a one-tailed test ('two_tailed' = FALSE), but 'one_tailed_hypothesis' is not specified.\n",
                  "You must specify the alternative hypothesis for a 'negative' or 'positive' correlation ",
                  "between trait values and diversification rates using the 'one_tailed_hypothesis' argument."))
    }

    if (two_tailed & !is.null(one_tailed_hypothesis))
    {
      stop(paste0("You selected a two-tailed test ('two_tailed' = TRUE), but also specified a 'one_tailed_hypothesis': '",one_tailed_hypothesis,"'.\n",
                  "If you want to test that hypothesis, please select a one-tailed test ('two_tailed' = FALSE).\n",
                  "If you want to compute a two-tailed test, remove the 'one_tailed_hypothesis' or replace it with 'NULL'."))
    }
  }

  ## Extract rates data
  if (rate_type == "speciation")
  {
    rates_data <- BAMM_data$tipLambda
  }
  else if (rate_type == "extinction")
  {
    rates_data <- BAMM_data$tipMu
  }
  else if (rate_type == "net_diversification")
  {
    rates_data <- lapply(X = 1:length(BAMM_data$tipLambda), FUN = function (i) { BAMM_data$tipLambda[[i]] - BAMM_data$tipMu[[i]] })
  }

  ## Extract regime data
  regimes_data <- BAMM_data$tipStates

  ## Set number of permutations
  if (is.null(nb_permutations))
  {
    # If NULL, set to the number of posterior samples
    nb_permutations <- length(rates_data)
  }

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:nb_permutations)
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[l]],
                                                           regimes = regimes_data[[l]], stringsAsFactors = FALSE)
  }

  # Permute tip rates on tips using blocks defined by regime membership
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run permutations in parallel
    posterior_samples_permuted_rates_data <- parallel::parLapply(
      cl = cl, X = posterior_samples_random_rates_data, fun = block_permute_rates_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    posterior_samples_permuted_rates_data <- lapply(posterior_samples_random_rates_data, block_permute_rates_data)
  }
  names(posterior_samples_permuted_rates_data) <- BAMM_data$BAMM_posterior_sample_ID

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })
  names(posterior_samples_obs_rates_data) <- BAMM_data$BAMM_posterior_sample_ID

  ## Print what is tested
  if (print_hypothesis)
  {
    if (two_tailed) # For two-tailed test
    {
      cat(paste0("Selected two-tailed Spearman's rank correlation test:\n\n",
                 "Null hypothesis: no correlation between trait data and diversification rates.\n\n",
                 "Alternative hypothesis: negative or positive correlation between trait data diversification rates.\n\n",
                 "'Estimate' stats is the ",alpha*100,"% quantile of differences in absolute rho-stats between observed and permuted data.\n",
                 "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))

    } else { # For one-tailed test

      if (one_tailed_hypothesis == "positive")
      {
        cat(paste0("Selected one-tailed positive Spearman's rank correlation test:\n\n",
                   "Null hypothesis: negative or no correlation between trait data and diversification rates.\n\n",
                   "Alternative hypothesis: positive correlation between trait data diversification rates.\n",
                   "Low trait values associated with low diversification rates, and conversely.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of rho differences between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      } else {
        cat(paste0("Selected one-tailed negative Spearman's rank correlation test:\n\n",
                   "Null hypothesis: positive or no correlation between trait data and diversification rates.\n\n",
                   "Alternative hypothesis: negative correlation between trait data diversification rates.\n",
                   "Low trait values associated with high diversification rates, and conversely.\n\n",
                   "'Estimate' stats is the ",(1-alpha)*100,"% quantile of rho differences between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is lower than zero / p-value lower than ",alpha,".\n\n"))
      }
    }
  }

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

  ### Compute test for 'rates_only' strategy: run tests for each BAMM posterior with a unique trait dataset
  if (uncertainty_strategy == "rates_only")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed with a unique trait dataset of ML estimates matched across rates from ",nb_permutations," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",nb_permutations,".\n\n"))
    }

    ## Compute correlation test on each permutation. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Compute test on observed data
      rho_obs <- parallel::parLapply(cl = cl,
                                     X = posterior_samples_obs_rates_data,
                                     fun = spearman_test,
                                     trait_data = trait_data)
      # Compute test on permuted data
      rho_perm <- parallel::parLapply(cl = cl,
                                      X = posterior_samples_permuted_rates_data,
                                      fun = spearman_test,
                                      trait_data = trait_data)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      rho_obs <- lapply(X = posterior_samples_obs_rates_data,
                        FUN = spearman_test,
                        trait_data = trait_data)
      rho_perm <- lapply(X = posterior_samples_permuted_rates_data,
                         FUN = spearman_test,
                         trait_data = trait_data)
    }
    # Unlist outputs
    rho_obs <- unlist(rho_obs)
    rho_perm <- unlist(rho_perm)
    # Assign names to track origin of data
    names(rho_obs) <- names(rho_perm) <- paste0("Map_ML_",names(posterior_samples_obs_rates_data))
  }

  ### Compute test for 'paired' strategy: run tests for each (BAMM posterior, trait dataset) pair
  if (uncertainty_strategy == "paired")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed across trait data extracted from ",nb_permutations," stochastic maps paired with rates from ",nb_permutations," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",nb_permutations,".\n\n"))
    }

    ## Compute correlation test on each pair of trait map X BAMM posterior. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Add the custom spearman_test function to the environment of each worker
      parallel::clusterExport(cl, varlist = c("spearman_test"), envir = environment())
      # Compute test on observed data
      rho_obs <- parallel::parLapply(
         cl = cl,
         X = seq_along(posterior_samples_obs_rates_data), # Loop along pair indices
         # Function is pairing rates and traits
         fun = function(i, rates, traits)
           {
             spearman_test(trait_data = traits[[i]], rates = rates[[i]])
           },
         rates = posterior_samples_obs_rates_data,
         traits = trait_data)
      # Compute test on permuted data
      rho_perm <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_permuted_rates_data), # Loop along pair indices
        # Function is pairing rates and traits
        fun = function(i, rates, traits)
        {
          spearman_test(trait_data = traits[[i]], rates = rates[[i]])
        },
        rates = posterior_samples_permuted_rates_data,
        traits = trait_data)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      rho_obs <- purrr::map2(.x = posterior_samples_obs_rates_data,
                             .y = trait_data,
                             .f = spearman_test)
      rho_perm <- purrr::map2(.x = posterior_samples_permuted_rates_data,
                              .y = trait_data,
                              .f = spearman_test)
    }
    # Unlist outputs
    rho_obs <- unlist(rho_obs)
    rho_perm <- unlist(rho_perm)
    # Assign names to track origin of data
    names(rho_obs) <- names(rho_perm) <- paste0(names(trait_data),"_",names(posterior_samples_obs_rates_data))
  }

  ### Compute test for 'full' strategy: run tests for each combination of BAMM posterior × trait dataset
  if (uncertainty_strategy == "full")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed across trait data extracted from ",length(trait_data)," stochastic maps crossed with rates from ",length(rates_data)," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",length(trait_data)*length(rates_data),".\n\n"))
    }

    ## Compute correlation test on each nested combination. For observed data and permuted data.
    # Nested with BAMM samples > trait_data
    # This prevent to have to distribute the BAMM_data object across multiple workers to increase efficency
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Add the custom spearman_test function to the environment of each worker
      parallel::clusterExport(cl, varlist = c("spearman_test"), envir = environment())
      # Compute test on observed data
      rho_obs <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_obs_rates_data), # Loop along BAMM rates
        # Function is mapping rates X traits
        function(i, rates, traits)
        {
          lapply(X = traits,
                 # Function is matching traits with a given BAMM rates
                 FUN = function (tr)
                 {
                   spearman_test(trait_data = tr, rates = rates[[i]])
                 })
        },
        rates = posterior_samples_obs_rates_data,
        traits = trait_data)
      # Compute test on permuted data
      rho_perm <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_permuted_rates_data), # Loop along BAMM rates
        # Function is mapping rates X traits
        function(i, rates, traits)
        {
          lapply(X = traits,
                 # Function is matching traits with a given BAMM rates
                 FUN = function (tr)
                 {
                   spearman_test(trait_data = tr, rates = rates[[i]])
                 })
        },
        rates = posterior_samples_permuted_rates_data,
        traits = trait_data)
      # Close cluster
      parallel::stopCluster(cl)

      # Add names of BAMM samples
      names(rho_obs) <- names(rho_perm) <- names(posterior_samples_obs_rates_data)
      # Unlist nested results
      rho_obs <- unlist(rho_obs)
      rho_perm <- unlist(rho_perm)

    } else { # In series
      rho_obs <- purrr::map(
        .x = posterior_samples_obs_rates_data,
        .f = function (rate_sample)
        {
          purrr::map(
            .x = trait_data,
            .f = function (trait_sample)
            {
              spearman_test(trait_data = trait_sample, rates = rate_sample)
            }
          )
        })
      rho_obs <- unlist(rho_obs)
      rho_perm <- purrr::map(
        .x = posterior_samples_permuted_rates_data,
        .f = function (rate_sample)
        {
          purrr::map(
            .x = trait_data,
            .f = function (trait_sample)
            {
              spearman_test(trait_data = trait_sample, rates = rate_sample)
            }
          )
        })
     rho_perm <- unlist(rho_perm)
    }
    # Assign names to track origin of data
    names(rho_obs) <- sub(x = names(rho_obs), pattern = "\\.rho$", replacement = "")
    names_df <- do.call(rbind, strsplit(names(rho_obs), split = "\\."))
    names(rho_obs) <- paste0(names_df[,2],"_",names_df[,1])

    names(rho_perm) <- sub(x = names(rho_perm), pattern = "\\.rho$", replacement = "")
    names_df <- do.call(rbind, strsplit(names(rho_perm), split = "\\."))
    names(rho_perm) <- paste0(names_df[,2],"_",names_df[,1])
  }

  ## Compute p-value for two-tailed test
  if (two_tailed)
  {
    # Ho: correlation is equal in observed data than in permuted data
    # Ha: correlation is lower or higher in observed data than in permuted data
    # P-value = frequency of cases where observed stats is less extreme
    # (closer from null hypothesis) than the permuted stats
    p_value <- sum(abs(rho_obs) <= abs(rho_perm)) / length(rho_perm)
  } else {

    ## Compute p-value for one-tailed tests

    if (one_tailed_hypothesis == "positive")
    { # Test for positive correlation
      # Ho: correlation is lower in observed data than in permuted data
      # Ha: correlation is higher in observed data than in permuted data
      # P-value = frequency of cases where observed stats is lower than the permuted stats
      p_value <- sum(rho_obs <= rho_perm) / length(rho_perm)
    } else { # Test for negative correlation
      # Ho: correlation is higher in observed data than in permuted data
      # Ha: correlation is lower in observed data than in permuted data
      # P-value = frequency of cases where observed stats is higher than the permuted stats
      p_value <- sum(rho_obs >= rho_perm) / length(rho_perm)
    }
  }

  ## Save test stats
  if (two_tailed)
  {
    # If two-tailed test, need to compare the abs_delta_rho with alpha % quantile to see if higher than zero.
    STRAPP_results <- list(
      estimate = stats::quantile(abs(rho_obs) - abs(rho_perm), p = alpha))
  } else {

    if (one_tailed_hypothesis == "positive")
    {
      # If one-tailed test for positive correlation, need to compare the delta_rho with alpha % quantile to see if higher than zero.
      STRAPP_results <- list(
        estimate = stats::quantile(as.numeric(rho_obs) - as.numeric(rho_perm), p = alpha))
    } else {
      # If one-tailed test for negative correlation, need to compare the delta_rho with (1-alpha) % quantile to see if lower than zero.
      STRAPP_results <- list(
        estimate = stats::quantile(as.numeric(rho_obs) - as.numeric(rho_perm), p = 1 - alpha))
    }
  }

  ## Save test summary results
  if (two_tailed)
  {
    # For two-tailed test, distribution based on difference in absolute correlations
    STRAPP_results$stats_median <- stats::median(abs(rho_obs) - abs(rho_perm))
  } else {
    # For one-tailed test, distribution based on difference in correlations
    STRAPP_results$stats_median <- stats::median(as.numeric(rho_obs) - as.numeric(rho_perm))
  }
  STRAPP_results$nb_test_stats <- length(rho_obs) # Nb of test stats in the distribution
  STRAPP_results$p_value <- p_value # P-value of the test
  STRAPP_results$method <- "Spearman" # Stats method
  STRAPP_results$two_tailed <- two_tailed # Type of test: two-tailed or not
  STRAPP_results$one_tailed_hypothesis <- one_tailed_hypothesis # Type of hypothesis if one-tailed test
  STRAPP_results$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification
  STRAPP_results$trait_data_type <- trait_data_type # Type of trait data: continuous, categorical, or biogeographic
  STRAPP_results$trait_data_type_for_stats <- "continuous" # Type of trait data used to select statistical method: continuous, binary, or multinominal
  STRAPP_results$uncertainty_strategy <- uncertainty_strategy # Type of strategy employed to account for uncertainty in estimates

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    # Extract Map_ID and BAMM_ID from names
    trait_map_ID <- sub(x = names(rho_obs), pattern = "_BAMM_.*", replacement = "")
    BAMM_posterior_sample_ID <- paste0("BAMM_",sub(x = names(rho_obs), pattern = ".*_BAMM_", replacement = ""))

    # Build perm_data_df
    perm_data_df <- data.frame(trait_map_ID = trait_map_ID, # Stochastic map ID
                               BAMM_posterior_sample_ID = BAMM_posterior_sample_ID, # As filtered in the master function
                               rho_obs = as.numeric(rho_obs),
                               rho_perm = as.numeric(rho_perm))

    if (two_tailed)
    { # For two-tailed test, distribution based on difference in absolute correlations
      perm_data_df$abs_delta_rho <- abs(rho_obs) - abs(rho_perm)
    } else { # For one-tailed test, distribution based on difference in correlations
      perm_data_df$delta_rho <- as.numeric(rho_obs) - as.numeric(rho_perm)
    }

    # Reorder per Trait_map_ID
    # perm_data_df <- perm_data_df[order(perm_data_df$trait_map_ID), ]
    # Store in output
    STRAPP_results$perm_data_df <- perm_data_df
  }

  ## Export output
  return(STRAPP_results)
}


### Sub-function to handle binary data ####

compute_STRAPP_test_for_binary_data <- function (
    BAMM_data, trait_data,
    trait_data_type,
    rate_type = "net_diversification",
    uncertainty_strategy = "paired",
    nb_permutations = NULL,
    alpha = 0.05,
    two_tailed = TRUE,
    one_tailed_hypothesis = NULL,
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{

  ### Check input validity
  {
    ## two_tailed & one_tailed_hypothesis
    if (!two_tailed & is.null(one_tailed_hypothesis))
    {
      stop(paste0("You selected a one-tailed test ('two_tailed' = FALSE), but 'one_tailed_hypothesis' is not specified.\n",
                  "You must specify the hypothesis by providing a character string vector with states ordered in increasing rates under the alternative hypothesis, separated by a greater-than such as c('A > B').\n"))
    }

    if (two_tailed & !is.null(one_tailed_hypothesis))
    {
      stop(paste0("You selected a two-tailed test ('two_tailed' = TRUE), but also specified a 'one_tailed_hypothesis': '",one_tailed_hypothesis,"'.\n",
                  "If you want to test that hypothesis, please select a one-tailed test ('two_tailed' = FALSE).\n",
                  "If you want to compute a two-tailed test, remove the 'one_tailed_hypothesis' or replace it with 'NULL'."))
    }
  }

  ## Extract rates data
  if (rate_type == "speciation")
  {
    rates_data <- BAMM_data$tipLambda
  }
  else if (rate_type == "extinction")
  {
    rates_data <- BAMM_data$tipMu
  }
  else if (rate_type == "net_diversification")
  {
    rates_data <- lapply(X = 1:length(BAMM_data$tipLambda), FUN = function (i) { BAMM_data$tipLambda[[i]] - BAMM_data$tipMu[[i]] })
  }

  ## Extract regime data
  regimes_data <- BAMM_data$tipStates

  ## Set number of permutations
  if (is.null(nb_permutations))
  {
    # If NULL, set to the number of posterior samples
    nb_permutations <- length(rates_data)
  }

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:nb_permutations)
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[l]],
                                                           regimes = regimes_data[[l]], stringsAsFactors = FALSE)
  }

  # Permute tip rates on tips using blocks defined by regime membership
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run permutations in parallel
    posterior_samples_permuted_rates_data <- parallel::parLapply(
      cl = cl, X = posterior_samples_random_rates_data, fun = block_permute_rates_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    posterior_samples_permuted_rates_data <- lapply(posterior_samples_random_rates_data, block_permute_rates_data)
  }
  names(posterior_samples_permuted_rates_data) <- BAMM_data$BAMM_posterior_sample_ID

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })
  names(posterior_samples_obs_rates_data) <- BAMM_data$BAMM_posterior_sample_ID

  ## Display hypothesis

  if (uncertainty_strategy == "rates_only") # Case with only one trait dataset
  {
    obs_trait_states <- unique(trait_data)[order(unique(trait_data))]
  } else { # Case with multiple trait datasets
    obs_trait_states <- unique(unlist(trait_data))
    obs_trait_states <- obs_trait_states[order(obs_trait_states)]
  }

  trait_states <- NA
  if (two_tailed)
  { # Case for two-tailed test

    if (print_hypothesis)
    {
      cat(paste0("Selected two-tailed Mann-Whitney-Wilcoxon rank-sum test:\n\n",
                 "Null hypothesis: taxa with state/range '",
                 obs_trait_states[1], "' have equal ",rate_type," rates than those with state/range '",
                 obs_trait_states[2], "'.\n",
                 "Alternative hypothesis: taxa with state/range '",
                 obs_trait_states[1], "' have higher or lower ",rate_type," rates than those with state/range '",
                 obs_trait_states[2], "'.\n\n",
                 "'Estimate' stats is the ",alpha*100,"% quantile of differences in absolute U-stats between observed and permuted data.\n",
                 "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
    }

  } else { # Case for one-tailed test

    # one_tailed_hypothesis <- c("state_A > state_B")

    # Parse one_tailed_hypothesis
    one_tailed_hypothesis_parsed <- gsub(pattern = " ", replacement = "", x = one_tailed_hypothesis)
    trait_states <- as.character(unlist(strsplit(x = one_tailed_hypothesis_parsed, split = ">")))

    if (!identical(trait_states, obs_trait_states) & !identical(rev(trait_states), obs_trait_states))
    {
      stop(paste0("States specified in the 'one_tailed_hypothesis' do not match with observed states in 'trait_data'.\n",
                  "States in 'one_tailed_hypothesis' = ", paste(trait_states, collapse = ", "), ". States in 'trait_data' = ", paste(obs_trait_states, collapse = ", "),".\n",
                  "Please use this format to provide the 'one_tailed_hypothesis': Two states separated by greater-than sign, with the state that is expected to have higher ",rate_type," rates in first.\n",
                  "Example: 'A > B'.\n"))
    } else {
      if (print_hypothesis)
      {
        cat(paste0("Selected one-tailed Mann-Whitney-Wilcoxon rank-sum test:\n\n",
                   "Null hypothesis: taxa with state/range '",
                   trait_states[1], "' have lower or equal ",rate_type," rates than those with state/range '",
                   trait_states[2], "'.\n",
                   "Alternative hypothesis: taxa with state/range '",
                   trait_states[1], "' have higher ",rate_type," rates than those with state/range '",
                   trait_states[2],"'.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of differences in absolute U-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      }
    }
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

  ### Compute test for 'rates_only' strategy: run tests for each BAMM posterior with a unique trait dataset
  if (uncertainty_strategy == "rates_only")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed with a unique trait dataset of ML estimates matched across rates from ",nb_permutations," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",nb_permutations,".\n\n"))
    }

    ## Compute MWW test on each permutation. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Compute test on observed data
      U_obs <- parallel::parLapply(cl = cl,
                                   X = posterior_samples_obs_rates_data,
                                   fun = mann_whitney_wilcoxon_test,
                                   trait_data = trait_data,
                                   two_tailed = two_tailed,
                                   trait_states = trait_states)
      # Compute test on permuted data
      U_perm <- parallel::parLapply(cl = cl,
                                    X = posterior_samples_permuted_rates_data,
                                    fun = mann_whitney_wilcoxon_test,
                                    trait_data = trait_data,
                                    two_tailed = two_tailed,
                                    trait_states = trait_states)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      U_obs <- lapply(X = posterior_samples_obs_rates_data,
                      FUN = mann_whitney_wilcoxon_test,
                      trait_data = trait_data,
                      two_tailed = two_tailed,
                      trait_states = trait_states)
      U_perm <- lapply(X = posterior_samples_permuted_rates_data,
                       FUN = mann_whitney_wilcoxon_test,
                       trait_data = trait_data,
                       two_tailed = two_tailed,
                       trait_states = trait_states)
    }

    ## Unlist outputs
    U_obs <- unlist(U_obs)
    U_perm <- unlist(U_perm)
    # Assign names to track origin of data
    names(U_obs) <- names(U_perm) <- paste0("Map_ML_",names(posterior_samples_obs_rates_data))

    ## Center stats around location shift of the null hypothesis (mu)
    # Null hypothesis is that ranks of the values of the two groups are random
    # Compute location shift (mu) from state frequencies as average of the products of frequencies
    trait_data_counts <- table(trait_data)
    trait_data_counts <- trait_data_counts[!is.na(names(trait_data_counts))] # Remove NA
    stat_mu <- prod(trait_data_counts)/2
    # Center U-stats to get an estimate of how greater/lower (far away) than the null hypothesis (mu) are the calculated U-stats
    U_obs <- U_obs - stat_mu
    U_perm <- U_perm - stat_mu
  }

  ### Compute test for 'paired' strategy: run tests for each (BAMM posterior, trait dataset) pair
  if (uncertainty_strategy == "paired")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed across trait data extracted from ",nb_permutations," stochastic maps paired with rates from ",nb_permutations," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",nb_permutations,".\n\n"))
    }

    ## Compute Mann-Whitney-Wilcoxon test on each pair of trait map X BAMM posterior. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Add the custom mann_whitney_wilcoxon_test function to the environment of each worker
      parallel::clusterExport(cl, varlist = c("mann_whitney_wilcoxon_test", "two_tailed", "trait_states"), envir = environment())
      # Compute test on observed data
      U_obs <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_obs_rates_data), # Loop along pair indices
        # Function is pairing rates and traits
        fun = function(i, rates, traits)
        {
          mann_whitney_wilcoxon_test(trait_data = traits[[i]], rates = rates[[i]],
                                     two_tailed = two_tailed,
                                     trait_states = trait_states)
        },
        rates = posterior_samples_obs_rates_data,
        traits = trait_data)
      # Compute test on permuted data
      U_perm <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_permuted_rates_data), # Loop along pair indices
        # Function is pairing rates and traits
        fun = function(i, rates, traits)
        {
          mann_whitney_wilcoxon_test(trait_data = traits[[i]], rates = rates[[i]],
                                     two_tailed = two_tailed,
                                     trait_states = trait_states)
        },
        rates = posterior_samples_permuted_rates_data,
        traits = trait_data)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      U_obs <- purrr::map2(.x = posterior_samples_obs_rates_data,
                             .y = trait_data,
                             .f = mann_whitney_wilcoxon_test,
                             two_tailed = two_tailed,
                             trait_states = trait_states)
      U_perm <- purrr::map2(.x = posterior_samples_permuted_rates_data,
                              .y = trait_data,
                              .f = mann_whitney_wilcoxon_test,
                              two_tailed = two_tailed,
                              trait_states = trait_states)
    }
    # Unlist outputs
    U_obs <- unlist(U_obs)
    U_perm <- unlist(U_perm)
    # Assign names to track origin of data
    names(U_obs) <- names(U_perm) <- paste0(names(trait_data),"_",names(posterior_samples_obs_rates_data))

    ## Center stats around location shift of the null hypothesis (mu)
    # Null hypothesis is that ranks of the values of the two groups are random
    # Compute location shift (mu) from state frequencies as average of the products of frequencies

    # Loop per trait_data as their distribution of states differ
    for (i in seq_along(trait_data))
    {
      trait_data_counts_i <- table(trait_data[[i]])
      trait_data_counts_i <- trait_data_counts_i[!is.na(names(trait_data_counts_i))] # Remove NA
      stat_mu_i <- prod(trait_data_counts_i)/2
      # Center U-stats to get an estimate of how greater/lower (far away) than the null hypothesis (mu) are the calculated U-stats
      U_obs[[i]] <- U_obs[[i]] - stat_mu_i
      U_perm[[i]] <- U_perm[[i]] - stat_mu_i
    }
  }

  ### Compute test for 'full' strategy: run tests for each combination of BAMM posterior × trait dataset
  if (uncertainty_strategy == "full")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed across trait data extracted from ",length(trait_data)," stochastic maps crossed with rates from ",length(rates_data)," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",length(trait_data)*length(rates_data),".\n\n"))
    }

    ## Compute Mann-Whitney-Wilcoxon test on each nested combination. For observed data and permuted data.
    # Nested with BAMM samples > trait_data
    # This prevent to have to distribute the BAMM_data object across multiple workers to increase efficency
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Add the custom mann_whitney_wilcoxon_test function to the environment of each worker
      parallel::clusterExport(cl, varlist = c("mann_whitney_wilcoxon_test", "two_tailed", "trait_states"), envir = environment())
      # Compute test on observed data
      U_obs <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_obs_rates_data), # Loop along BAMM rates
        # Function is mapping rates X traits
        function(i, rates, traits)
        {
          lapply(X = traits,
                 # Function is matching traits with a given BAMM rates
                 FUN = function (tr)
                 {
                   mann_whitney_wilcoxon_test(trait_data = tr, rates = rates[[i]],
                                              two_tailed = two_tailed,
                                              trait_states = trait_states)
                 })
        },
        rates = posterior_samples_obs_rates_data,
        traits = trait_data)
      # Compute test on permuted data
      U_perm <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_permuted_rates_data), # Loop along BAMM rates
        # Function is mapping rates X traits
        function(i, rates, traits)
        {
          lapply(X = traits,
                 # Function is matching traits with a given BAMM rates
                 FUN = function (tr)
                 {
                   mann_whitney_wilcoxon_test(trait_data = tr, rates = rates[[i]],
                                              two_tailed = two_tailed,
                                              trait_states = trait_states)
                 })
        },
        rates = posterior_samples_permuted_rates_data,
        traits = trait_data)
      # Close cluster
      parallel::stopCluster(cl)

      # Add names of BAMM samples
      names(U_obs) <- names(U_perm) <- names(posterior_samples_obs_rates_data)
      # Unlist nested results
      U_obs <- unlist(U_obs)
      U_perm <- unlist(U_perm)

    } else { # In series
      U_obs <- purrr::map(
        .x = posterior_samples_obs_rates_data,
        .f = function (rate_sample)
        {
          purrr::map(
            .x = trait_data,
            .f = function (trait_sample)
            {
              mann_whitney_wilcoxon_test(trait_data = trait_sample, rates = rate_sample,
                                         two_tailed = two_tailed, trait_states = trait_states)
            }
          )
        })
      U_obs <- unlist(U_obs)
      U_perm <- purrr::map(
        .x = posterior_samples_permuted_rates_data,
        .f = function (rate_sample)
        {
          purrr::map(
            .x = trait_data,
            .f = function (trait_sample)
            {
              mann_whitney_wilcoxon_test(trait_data = trait_sample, rates = rate_sample,
                                         two_tailed = two_tailed, trait_states = trait_states)
            }
          )
        })
      U_perm <- unlist(U_perm)
    }

    # Assign names to track origin of data
    names(U_obs) <- sub(x = names(U_obs), pattern = "\\.W$", replacement = "")
    names_df <- do.call(rbind, strsplit(names(U_obs), split = "\\."))
    names(U_obs) <- paste0(names_df[,2],"_",names_df[,1])

    names(U_perm) <- sub(x = names(U_perm), pattern = "\\.W$", replacement = "")
    names_df <- do.call(rbind, strsplit(names(U_perm), split = "\\."))
    names(U_perm) <- paste0(names_df[,2],"_",names_df[,1])

    ## Center stats around location shift of the null hypothesis (mu)
    # Null hypothesis is that ranks of the values of the two groups are random
    # Compute location shift (mu) from state frequencies as average of the products of frequencies

    # Loop per trait_data as their distribution of states differ
    for (i in seq_along(trait_data))
    {
      trait_data_counts_i <- table(trait_data[[i]])
      trait_data_counts_i <- trait_data_counts_i[!is.na(names(trait_data_counts_i))] # Remove NA
      stat_mu_i <- prod(trait_data_counts_i)/2
      # Center U-stats to get an estimate of how greater/lower (far away) than the null hypothesis (mu) are the calculated U-stats
      U_obs[[i]] <- U_obs[[i]] - stat_mu_i
      U_perm[[i]] <- U_perm[[i]] - stat_mu_i
    }
  }

  ## Compute p-value for two-tailed test
  if (two_tailed)
  {
    # Ho: rate differences in ranks between states are equal in observed data and permuted data
    # Ha: rate differences in ranks between states are lower or higher in observed data than in permuted data
    # P-value = frequency of cases where observed stats is less extreme
    # (closer from null hypothesis) than the permuted stats
    p_value <- sum(abs(U_obs) <= abs(U_perm)) / length(U_perm)

  } else {

    ## Compute p-value for one-tailed tests

    # Ho: rate differences in ranks between states are lower or equal in observed data than in permuted data
    # Ha: rate differences in ranks between states are higher in observed data than in permuted data
    # P-value = frequency of cases where observed stats is lower than the permuted stats
    p_value <- sum(U_obs <= U_perm) / length(U_perm)
  }

  ## Save test stats
  if (two_tailed)
  {
    # If two-tailed test, need to compare the abs_delta_U with alpha % quantile to see if higher than zero.
    STRAPP_results <- list(
      estimate = stats::quantile(abs(U_obs) - abs(U_perm), p = alpha))
  } else {

    # If one-tailed test, need to compare the delta_U with alpha % quantile to see if higher than zero.
    STRAPP_results <- list(
      estimate = stats::quantile(as.numeric(U_obs) - as.numeric(U_perm), p = alpha))
  }

  ## Save test summary results
  if (two_tailed)
  {
    # For two-tailed test, distribution based on difference in absolute U-stats
    STRAPP_results$stats_median <- stats::median(abs(U_obs) - abs(U_perm))
  } else {
    # For one-tailed test, distribution based on difference in U-stats
    STRAPP_results$stats_median <- stats::median(as.numeric(U_obs) - as.numeric(U_perm))
  }
  STRAPP_results$nb_test_stats <- length(U_obs) # Nb of test stats in the distribution
  STRAPP_results$p_value <- p_value # P-value of the test
  STRAPP_results$method <- "Mann-Whitney-Wilcoxon" # Stats method
  STRAPP_results$two_tailed <- two_tailed # Type of test: two-tailed or not
  STRAPP_results$one_tailed_hypothesis <- one_tailed_hypothesis # Type of hypothesis if one-tailed test
  STRAPP_results$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification
  STRAPP_results$trait_data_type <- trait_data_type # Type of trait data: continuous, categorical, or biogeographic
  STRAPP_results$trait_data_type_for_stats <- "binary" # Type of trait data used to select statistical method: continuous, binary, or multinominal
  STRAPP_results$uncertainty_strategy <- uncertainty_strategy # Type of strategy employed to account for uncertainty in estimates

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    # Extract Map_ID and BAMM_ID from names
    trait_map_ID <- sub(x = names(U_obs), pattern = "_BAMM_.*", replacement = "")
    BAMM_posterior_sample_ID <- paste0("BAMM_",sub(x = names(U_obs), pattern = ".*_BAMM_", replacement = ""))

    # Build perm_data_df
    perm_data_df <- data.frame(trait_map_ID = trait_map_ID, # Stochastic map ID
                               BAMM_posterior_sample_ID = BAMM_posterior_sample_ID, # As filtered in the master function
                               U_obs = as.numeric(U_obs),
                               U_perm = as.numeric(U_perm))

    if (two_tailed)
    { # For two-tailed test, distribution based on difference in absolute U-stats
      perm_data_df$abs_delta_U <- abs(U_obs) - abs(U_perm)
    } else { # For one-tailed test, distribution based on difference in U-stats
      perm_data_df$delta_U <- as.numeric(U_obs) - as.numeric(U_perm)
    }

    # Reorder per Trait_map_ID
    # perm_data_df <- perm_data_df[order(perm_data_df$trait_map_ID), ]
    # Store in output
    STRAPP_results$perm_data_df <- perm_data_df
  }

  ## Export output
  return(STRAPP_results)
}


### Sub-function to handle multinominal data ####

compute_STRAPP_test_for_multinominal_data <- function (
    BAMM_data, trait_data,
    trait_data_type,
    rate_type = "net_diversification",
    uncertainty_strategy = "paired",
    nb_permutations = NULL,
    alpha = 0.05,
    posthoc_pairwise_tests = FALSE,
    two_tailed = TRUE,
    p.adjust_method = "none",
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{

  ### Check input validity
  {

  }

  ## Extract rates data
  if (rate_type == "speciation")
  {
    rates_data <- BAMM_data$tipLambda
  }
  else if (rate_type == "extinction")
  {
    rates_data <- BAMM_data$tipMu
  }
  else if (rate_type == "net_diversification")
  {
    rates_data <- lapply(X = 1:length(BAMM_data$tipLambda), FUN = function (i) { BAMM_data$tipLambda[[i]] - BAMM_data$tipMu[[i]] })
  }

  ## Extract regime data
  regimes_data <- BAMM_data$tipStates

  ## Set number of permutations
  if (is.null(nb_permutations))
  {
    # If NULL, set to the number of posterior samples
    nb_permutations <- length(rates_data)
  }

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:nb_permutations)
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[l]],
                                                           regimes = regimes_data[[l]], stringsAsFactors = FALSE)
  }

  # Permute tip rates on tips using blocks defined by regime membership
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run permutations in parallel
    posterior_samples_permuted_rates_data <- parallel::parLapply(
      cl = cl, X = posterior_samples_random_rates_data, fun = block_permute_rates_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    posterior_samples_permuted_rates_data <- lapply(posterior_samples_random_rates_data, block_permute_rates_data)
  }
  names(posterior_samples_permuted_rates_data) <- BAMM_data$BAMM_posterior_sample_ID

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })
  names(posterior_samples_obs_rates_data) <- BAMM_data$BAMM_posterior_sample_ID

  ## Display hypothesis

  if (print_hypothesis)
  {
    cat(paste0("Selected Kruskal-Wallis's one-way ANOVA on ranks test:\n\n",
               "Null hypothesis: taxa have equal ",rate_type," rates independent from states.\n",
               "Alternative hypothesis: taxa have different ",rate_type," rates between states.\n\n",
               "'Estimate' stats is the ",alpha*100,"% quantile of differences in H-stats between observed and permuted data.\n",
               "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
  }

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

  ### Compute test for 'rates_only' strategy: run tests for each BAMM posterior with a unique trait dataset
  if (uncertainty_strategy == "rates_only")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed with a unique trait dataset of ML estimates matched across rates from ",nb_permutations," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",nb_permutations,".\n\n"))
    }

    ## Compute Kruskal-Wallis test on each permutation. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Compute test on observed data
      H_obs <- parallel::parLapply(cl = cl,
                                   X = posterior_samples_obs_rates_data,
                                   fun = kruskal_wallis_test,
                                   trait_data = trait_data)
      # Compute test on permuted data
      H_perm <- parallel::parLapply(cl = cl,
                                    X = posterior_samples_permuted_rates_data,
                                    fun = kruskal_wallis_test,
                                    trait_data = trait_data)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      H_obs <- lapply(X = posterior_samples_obs_rates_data,
                      FUN = kruskal_wallis_test,
                      trait_data = trait_data)
      H_perm <- lapply(X = posterior_samples_permuted_rates_data,
                       FUN = kruskal_wallis_test,
                       trait_data = trait_data)
    }

    ## Unlist outputs
    H_obs <- unlist(H_obs)
    H_perm <- unlist(H_perm)
    # Assign names to track origin of data
    names(H_obs) <- names(H_perm) <- paste0("Map_ML_",names(posterior_samples_obs_rates_data))
  }

  ### Compute test for 'paired' strategy: run tests for each (BAMM posterior, trait dataset) pair
  if (uncertainty_strategy == "paired")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed across trait data extracted from ",nb_permutations," stochastic maps paired with rates from ",nb_permutations," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",nb_permutations,".\n\n"))
    }

    ## Compute Kruskal-Wallis test on each pair of trait map X BAMM posterior. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Add the custom kruskal_wallis_test function to the environment of each worker
      parallel::clusterExport(cl, varlist = c("kruskal_wallis_test"), envir = environment())
      # Compute test on observed data
      H_obs <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_obs_rates_data), # Loop along pair indices
        # Function is pairing rates and traits
        fun = function(i, rates, traits)
        {
          kruskal_wallis_test(trait_data = traits[[i]], rates = rates[[i]])
        },
        rates = posterior_samples_obs_rates_data,
        traits = trait_data)
      # Compute test on permuted data
      H_perm <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_permuted_rates_data), # Loop along pair indices
        # Function is pairing rates and traits
        fun = function(i, rates, traits)
        {
          kruskal_wallis_test(trait_data = traits[[i]], rates = rates[[i]])
        },
        rates = posterior_samples_permuted_rates_data,
        traits = trait_data)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      H_obs <- purrr::map2(.x = posterior_samples_obs_rates_data,
                           .y = trait_data,
                           .f = kruskal_wallis_test)
      H_perm <- purrr::map2(.x = posterior_samples_permuted_rates_data,
                            .y = trait_data,
                            .f = kruskal_wallis_test)
    }

    ## Unlist outputs
    H_obs <- unlist(H_obs)
    H_perm <- unlist(H_perm)
    # Assign names to track origin of data
    names(H_obs) <- names(H_perm) <- paste0(names(trait_data),"_",names(posterior_samples_obs_rates_data))
  }

  ### Compute test for 'full' strategy: run tests for each combination of BAMM posterior × trait dataset
  if (uncertainty_strategy == "full")
  {
    if (print_hypothesis)
    {
      cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                 "Tests will be performed across trait data extracted from ",length(trait_data)," stochastic maps crossed with rates from ",length(rates_data)," BAMM posterior samples.\n",
                 "Total number of test stats in the null distribution = ",length(trait_data)*length(rates_data),".\n\n"))
    }

    ## Compute Kruskal-Wallis test on each nested combination. For observed data and permuted data.
    # Nested with BAMM samples > trait_data
    # This prevent to have to distribute the BAMM_data object across multiple workers to increase efficiency
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Add the custom kruskal_wallis_test function to the environment of each worker
      parallel::clusterExport(cl, varlist = c("kruskal_wallis_test"), envir = environment())
      # Compute test on observed data
      H_obs <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_obs_rates_data), # Loop along BAMM rates
        # Function is mapping rates X traits
        function(i, rates, traits)
        {
          lapply(X = traits,
                 # Function is matching traits with a given BAMM rates
                 FUN = function (tr)
                 {
                   kruskal_wallis_test(trait_data = tr, rates = rates[[i]])
                 })
        },
        rates = posterior_samples_obs_rates_data,
        traits = trait_data)
      # Compute test on permuted data
      H_perm <- parallel::parLapply(
        cl = cl,
        X = seq_along(posterior_samples_permuted_rates_data), # Loop along BAMM rates
        # Function is mapping rates X traits
        function(i, rates, traits)
        {
          lapply(X = traits,
                 # Function is matching traits with a given BAMM rates
                 FUN = function (tr)
                 {
                   kruskal_wallis_test(trait_data = tr, rates = rates[[i]])
                 })
        },
        rates = posterior_samples_permuted_rates_data,
        traits = trait_data)
      # Close cluster
      parallel::stopCluster(cl)

      # Add names of BAMM samples
      names(H_obs) <- names(H_perm) <- names(posterior_samples_obs_rates_data)
      # Unlist nested results
      H_obs <- unlist(H_obs)
      H_perm <- unlist(H_perm)

    } else { # In series
      H_obs <- purrr::map(
        .x = posterior_samples_obs_rates_data,
        .f = function (rate_sample)
        {
          purrr::map(
            .x = trait_data,
            .f = function (trait_sample)
            {
              kruskal_wallis_test(trait_data = trait_sample, rates = rate_sample)
            }
          )
        })
      H_obs <- unlist(H_obs)
      H_perm <- purrr::map(
        .x = posterior_samples_permuted_rates_data,
        .f = function (rate_sample)
        {
          purrr::map(
            .x = trait_data,
            .f = function (trait_sample)
            {
              kruskal_wallis_test(trait_data = trait_sample, rates = rate_sample)
            }
          )
        })
      H_perm <- unlist(H_perm)
    }

    # Assign names to track origin of data
    names(H_obs) <- sub(x = names(H_obs), pattern = "\\.Kruskal-Wallis chi-squared$", replacement = "")
    names_df <- do.call(rbind, strsplit(names(H_obs), split = "\\."))
    names(H_obs) <- paste0(names_df[,2],"_",names_df[,1])

    names(H_perm) <- sub(x = names(H_perm), pattern = "\\.Kruskal-Wallis chi-squared$", replacement = "")
    names_df <- do.call(rbind, strsplit(names(H_perm), split = "\\."))
    names(H_perm) <- paste0(names_df[,2],"_",names_df[,1])
  }

  ## Compute p-value

  # Ho: rate differences in ranks between states are equal in observed data and permuted data
  # Ha: rate differences in ranks between states are more extreme in observed data than in permuted data
  # Because the H-stat follows a Khi-squared distribution, thus is strictly positive and increases when
  # ranks between states are biased indifferently towards lower or higher ranks,
  # the comparison is made on H-stat differences.
  # P-value = frequency of cases where observed stats is higher than the permuted stats
  p_value <- sum(H_obs <= H_perm) / length(H_perm)

  ## Save test stats

  # Need to compare the delta_H with alpha % quantile to see if higher than zero.
  STRAPP_results <- list(
    estimate = stats::quantile(as.numeric(H_obs) - as.numeric(H_perm), p = alpha))

  ## Save test summary results
  STRAPP_results$stats_median <- stats::median(as.numeric(H_obs) - as.numeric(H_perm))
  STRAPP_results$nb_test_stats <- length(H_obs) # Number of test stats in the null distribution
  STRAPP_results$p_value <- p_value # P-value of the test
  STRAPP_results$method <- "Kruskal-Wallis" # Stats method
  STRAPP_results$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification
  STRAPP_results$trait_data_type <- trait_data_type # Type of trait data: continuous, categorical, or biogeographic
  STRAPP_results$trait_data_type_for_stats <- "multinominal" # Type of trait data used to select statistical method: continuous, binary, or multinominal
  STRAPP_results$uncertainty_strategy <- uncertainty_strategy # Type of strategy employed to account for uncertainty in estimates

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    # Extract Map_ID and BAMM_ID from names
    trait_map_ID <- sub(x = names(H_obs), pattern = "_BAMM_.*", replacement = "")
    BAMM_posterior_sample_ID <- paste0("BAMM_",sub(x = names(H_obs), pattern = ".*_BAMM_", replacement = ""))

    # Build perm_data_df
    perm_data_df <- data.frame(trait_map_ID = trait_map_ID, # Stochastic map ID
                               BAMM_posterior_sample_ID = BAMM_posterior_sample_ID, # As filtered in the master function
                               H_obs = as.numeric(H_obs),
                               H_perm = as.numeric(H_perm),
                               delta_H = as.numeric(H_obs) - as.numeric(H_perm))

    # Reorder per Trait_map_ID
    # perm_data_df <- perm_data_df[order(perm_data_df$trait_map_ID), ]
    # Store in output
    STRAPP_results$perm_data_df <- perm_data_df
  }

  #### For posthoc tests, need to adjust results to account only for map cases where the tested pair is present across all maps!!!! #####
  ## Careful of the post-hoc trick. Need to be run on maps that host each pair only
  # i.e, some maps will not have pair A/B despite the overall test being able to run
  # Think about consequences for downstream plots

  # See how example works. Maybe no need to worry, we can just not record any value for cases with missing pairs of states,
  # since by default tests for all pairs of present states are run.

  ## Deal with post hoc pairwise tests
  if (posthoc_pairwise_tests)
  {
    ## Initiate output elements for post hoc pairwise tests
    STRAPP_results$posthoc_pairwise_tests <- list()

    ## Print hypothesis if requested
    if (print_hypothesis)
    {
      if (two_tailed)
      {
        cat(paste0("# --------- Post hoc pairwise tests --------- #\n\n",
                   "Selected two-tailed Dunn's post hoc pairwise rank-sum test:\n",
                   "Tests will be ran across all possible unique pairs of states.\n\n",
                   "Null hypothesis: taxa have equal ",rate_type," rates independent from states.\n",
                   "Alternative hypothesis: taxa have different ",rate_type," rates between states.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of absolute differences in Z-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      } else {
        cat(paste0("# --------- Post hoc pairwise tests --------- #\n\n",
                   "Selected one-tailed Dunn's post hoc pairwise rank-sum test:\n",
                   "Tests will be ran across all possible asymmetric pairs of states.\n\n",
                   "Null hypothesis: taxa in the first state/range have lower or equal ",rate_type," rates than taxa in the second state/range in 'pairs'.\n",
                   "Alternative hypothesis: taxa in the first state/range have higher ",rate_type," rates than taxa in the second state/range in 'pairs'.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of differences in Z-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      }
    }

    ## Wrapped-up function to extract Z-stats from Dunn's post hoc pairwise rank-sum tests
    dunn_test <- function(rates, trait_data, two_tailed)
    {
      # If all rates are similar, dunn.test will throw an error
      # Avoid this by providing a dummy test_output_df
      if (length(unique(rates)) == 1)
      {
        # Run fake Dunn test to get pairs
        nb_taxa <- length(rates)
        fake_rates <- c(rep(x = 1, times = nb_taxa/2), rep(x = 2, times = nb_taxa/2))
        if (round(nb_taxa/2) != (nb_taxa/2)) { fake_rates <- c(fake_rates, 2) } # Add an extra rates for odd number of rates
        invisible(utils::capture.output(test_output <- dunn.test::dunn.test(x = fake_rates, g = trait_data)))
        Z_approximation <- stats::qnorm(p = 1 - 10^-9)
        # Create dummy test_output_df
        if (two_tailed) # For two-tailed tests
        {
        test_output_df <- data.frame(pairs = gsub(pattern = " - ", replacement = " != ", x = test_output$comparisons),
                                     Z_stats = Z_approximation)
        } else { # For one-tailed tests
          test_output_df <- data.frame(pairs = c(gsub(pattern = " - ", replacement = " > ", x = test_output$comparisons), gsub(pattern = " - ", replacement = " < ", x = test_output$comparisons)),
                                       Z_stats = c(rep(x = Z_approximation, times = length(test_output$comparisons)), rep(x = -Z_approximation, times = length(test_output$comparisons))))
        }
      } else {

        # Compute the Dunn test for all possible unique pairs of states
        invisible(utils::capture.output(test_output <- dunn.test::dunn.test(x = rates, g = trait_data)))

        # Reformat test output
        if (two_tailed) # For two-tailed tests
        {
          test_output_df <- data.frame(pairs = gsub(pattern = " - ", replacement = " != ", x = test_output$comparisons),
                                       Z_stats = test_output$Z)
        } else { # For one-tailed tests
          test_output_df <- data.frame(pairs = c(gsub(pattern = " - ", replacement = " > ", x = test_output$comparisons), gsub(pattern = " - ", replacement = " < ", x = test_output$comparisons)),
                                       Z_stats = c(test_output$Z, -test_output$Z))
        }

      }



      # If the test failed to provide a statistic because the value is reaching the ceiling for computation,
      # use the normal distribution and set an extremely high p-value to approximate a value
      if (any(is.na(test_output_df$Z_stats)))
      {
        Z_approximation <- stats::qnorm(p = 1 - 10^-9)
        test_output_df$Z_stats[is.na(test_output_df$Z_stats)] <- Z_approximation

      }

      # Export df
      return(test_output_df)
    }

    ### Compute test for 'rates_only' strategy: run tests for each BAMM posterior with a unique trait dataset
    if (uncertainty_strategy == "rates_only")
    {
      if (print_hypothesis)
      {
        cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                   "Pairwise post hoc tests will be performed with a unique trait dataset of ML estimates matched across rates from ",nb_permutations," BAMM posterior samples.\n",
                   "Total number of test stats in the null distribution of each pairwise post hoc test = ",nb_permutations,".\n\n"))
      }

      ## Compute Dunn test on each permutation. For observed data and permuted data.
      if (nthreads > 1) # In parallel
      {
        # Open cluster
        cl <- parallel::makePSOCKcluster(nthreads)
        # Compute test on observed data
        Dunn_obs <- parallel::parLapply(cl = cl,
                                        X = posterior_samples_obs_rates_data,
                                        fun = dunn_test,
                                        trait_data = trait_data,
                                        two_tailed = two_tailed)
        # Compute test on permuted data
        Dunn_perm <- parallel::parLapply(cl = cl,
                                         X = posterior_samples_permuted_rates_data,
                                         fun = dunn_test,
                                         trait_data = trait_data,
                                         two_tailed = two_tailed)
        # Close cluster
        parallel::stopCluster(cl)
      } else { # In series
        Dunn_obs <- lapply(X = posterior_samples_obs_rates_data,
                           FUN = dunn_test,
                           trait_data = trait_data,
                           two_tailed = two_tailed)
        Dunn_perm <- lapply(X = posterior_samples_permuted_rates_data,
                            FUN = dunn_test,
                            trait_data = trait_data,
                            two_tailed = two_tailed)
      }
      # Assign names to track origin of data
      names(Dunn_obs) <- names(Dunn_perm) <- paste0("Map_ML_",names(posterior_samples_obs_rates_data))

      ## Extract Z-scores from outputs

      # Get list of pairs (unique list as we use one ML trait estimate set)
      pairs_list <- Dunn_obs[[1]]$pairs
      # Initiate list for Z-scores
      Z_obs <- list()
      Z_perm <- list()

      # Loop per pairs of states tested
      for (i in seq_along(pairs_list))
      {
        # i <- 1
        Z_obs[[i]] <- unlist(lapply(X = Dunn_obs, FUN = function (x) { x$Z_stats[i] }))
        Z_perm[[i]] <- unlist(lapply(X = Dunn_perm, FUN = function (x) { x$Z_stats[i] }))
      }
      names(Z_obs) <- names(Z_perm) <- pairs_list
    }

    ### Compute test for 'paired' strategy: run tests for each (BAMM posterior, trait dataset) pair
    if (uncertainty_strategy == "paired")
    {
      if (print_hypothesis)
      {
        cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                   "Pairwise post hoc tests will be performed across trait data extracted from ",nb_permutations," stochastic maps paired with rates from ",nb_permutations," BAMM posterior samples.\n",
                   "Maximum number of test stats in the null distribution of each pairwise post hoc test = ",nb_permutations,".\n",
                   "The actual number of test stats depends on the availability of the tested pair of states across all stochastic maps.\n\n"))
      }

      ## Compute Dunn test on each pair of trait map X BAMM posterior. For observed data and permuted data.
      if (nthreads > 1) # In parallel
      {
        # Open cluster
        cl <- parallel::makePSOCKcluster(nthreads)
        # Add the custom kruskal_wallis_test function to the environment of each worker
        parallel::clusterExport(cl, varlist = c("kruskal_wallis_test"), envir = environment())
        # Compute test on observed data
        Dunn_obs <- parallel::parLapply(
          cl = cl,
          X = seq_along(posterior_samples_obs_rates_data), # Loop along pair indices
          # Function is pairing rates and traits
          fun = function(i, rates, traits, two_tailed)
          {
            dunn_test(trait_data = traits[[i]], rates = rates[[i]], two_tailed = two_tailed)
          },
          rates = posterior_samples_obs_rates_data,
          traits = trait_data,
          two_tailed = two_tailed)
        # Compute test on permuted data
        Dunn_perm <- parallel::parLapply(
          cl = cl,
          X = seq_along(posterior_samples_permuted_rates_data), # Loop along pair indices
          # Function is pairing rates and traits
          fun = function(i, rates, traits, two_tailed)
          {
            dunn_test(trait_data = traits[[i]], rates = rates[[i]], two_tailed = two_tailed)
          },
          rates = posterior_samples_permuted_rates_data,
          traits = trait_data,
          two_tailed = two_tailed)
        # Close cluster
        parallel::stopCluster(cl)
      } else { # In series
        Dunn_obs <- purrr::map2(.x = posterior_samples_obs_rates_data,
                             .y = trait_data,
                             .f = dunn_test,
                             two_tailed = two_tailed)
        Dunn_perm <- purrr::map2(.x = posterior_samples_permuted_rates_data,
                              .y = trait_data,
                              .f = dunn_test,
                              two_tailed = two_tailed)
      }
      # Assign names to track origin of data
      names(Dunn_obs) <- names(Dunn_perm) <- paste0(names(trait_data),"_",names(posterior_samples_obs_rates_data))

      ## Extract Z-scores from outputs

      # Get list of pairs by recording all unique pairs found across all stochastic maps
      pairs_list <- unique(unlist(lapply(X = Dunn_obs, FUN = function (x) { x$pairs })))
      pairs_list <- pairs_list[order(pairs_list)]

      # Initiate list for Z-scores
      Z_obs <- list()
      Z_perm <- list()

      # Loop per pairs of states tested
      for (i in seq_along(pairs_list))
      {
        # i <- 1
        Z_obs[[i]] <- unlist(lapply(X = Dunn_obs, FUN = function (x) { x$Z_stats[i] } ))
        Z_perm[[i]] <- unlist(lapply(X = Dunn_perm, FUN = function (x) { x$Z_stats[i] } ))
      }
      names(Z_obs) <- names(Z_perm) <- pairs_list

    }

    ### Compute test for 'full' strategy: run tests for each combination of BAMM posterior × trait dataset
    if (uncertainty_strategy == "full")
    {
      if (print_hypothesis)
      {
        cat(paste0("Uncertainty strategy = '",uncertainty_strategy,"'\n",
                   "Pairwise post hoc tests will be performed across trait data extracted from ",length(trait_data)," stochastic maps crossed with rates from ",length(rates_data)," BAMM posterior samples.\n",
                   "Maximum number of test stats in the null distribution = ",length(trait_data)*length(rates_data),".\n",
                   "The actual number of test stats depends on the availability of the tested pair of states across stochastic maps.\n"))
      }

      ## Compute Dunn test on each nested combination. For observed data and permuted data.
      # Nested with BAMM samples > trait_data
      # This prevent to have to distribute the BAMM_data object across multiple workers to increase efficiency
      if (nthreads > 1) # In parallel
      {
        # Open cluster
        cl <- parallel::makePSOCKcluster(nthreads)
        # Add the custom dunn_test function to the environment of each worker
        parallel::clusterExport(cl, varlist = c("dunn_test"), envir = environment())
        # Compute test on observed data
        Dunn_obs <- parallel::parLapply(
          cl = cl,
          X = seq_along(posterior_samples_obs_rates_data), # Loop along BAMM rates
          # Function is mapping rates X traits
          function(i, rates, traits, two_tailed)
          {
            lapply(X = traits,
                   # Function is matching traits with a given BAMM rates
                   FUN = function (tr)
                   {
                     dunn_test(trait_data = tr, rates = rates[[i]], two_tailed = two_tailed)
                   })
          },
          rates = posterior_samples_obs_rates_data,
          traits = trait_data,
          two_tailed = two_tailed)
        # Compute test on permuted data
        Dunn_perm <- parallel::parLapply(
          cl = cl,
          X = seq_along(posterior_samples_permuted_rates_data), # Loop along BAMM rates
          # Function is mapping rates X traits
          function(i, rates, traits, two_tailed)
          {
            lapply(X = traits,
                   # Function is matching traits with a given BAMM rates
                   FUN = function (tr)
                   {
                     dunn_test(trait_data = tr, rates = rates[[i]], two_tailed = two_tailed)
                   })
          },
          rates = posterior_samples_permuted_rates_data,
          traits = trait_data,
          two_tailed = two_tailed)
        # Close cluster
        parallel::stopCluster(cl)

        # Add names of BAMM samples
        names(Dunn_obs) <- names(Dunn_perm) <- names(posterior_samples_obs_rates_data)
        # Unlist nested results
        Dunn_obs <- unlist(Dunn_obs, recursive = FALSE, use.names = TRUE) # Unlist only the first level
        Dunn_perm <- unlist(Dunn_perm, recursive = FALSE, use.names = TRUE) # Unlist only the first level

      } else { # In series
        Dunn_obs <- purrr::map(
          .x = posterior_samples_obs_rates_data,
          .f = function (rate_sample)
          {
            purrr::map(
              .x = trait_data,
              .f = function (trait_sample)
              {
                dunn_test(trait_data = trait_sample, rates = rate_sample, two_tailed = two_tailed)
              }
            )
          })
        Dunn_obs <- unlist(Dunn_obs, recursive = FALSE, use.names = TRUE) # Unlist only the first level
        Dunn_perm <- purrr::map(
          .x = posterior_samples_permuted_rates_data,
          .f = function (rate_sample)
          {
            purrr::map(
              .x = trait_data,
              .f = function (trait_sample)
              {
                dunn_test(trait_data = trait_sample, rates = rate_sample, two_tailed = two_tailed)
              }
            )
          })
        Dunn_perm <- unlist(Dunn_perm, recursive = FALSE, use.names = TRUE) # Unlist only the first level
      }

      ## Assign names as 'Map_*_BAMM_*' to track origin of data
      names_df <- do.call(rbind, strsplit(names(Dunn_obs), split = "\\."))
      names(Dunn_obs) <- paste0(names_df[,2],"_",names_df[,1])

      names_df <- do.call(rbind, strsplit(names(Dunn_perm), split = "\\."))
      names(Dunn_perm) <- paste0(names_df[,2],"_",names_df[,1])

      ## Extract Z-scores from outputs

      # Get list of pairs by recording all unique pairs found across all stochastic maps
      pairs_list <- unique(unlist(lapply(X = Dunn_obs, FUN = function (x) { x$pairs })))
      pairs_list <- pairs_list[order(pairs_list)]

      # Initiate list for Z-scores
      Z_obs <- list()
      Z_perm <- list()

      # Loop per pairs of states tested
      for (i in seq_along(pairs_list))
      {
        # i <- 1
        Z_obs[[i]] <- unlist(lapply(X = Dunn_obs, FUN = function (x) { x$Z_stats[i] } ))
        Z_perm[[i]] <- unlist(lapply(X = Dunn_perm, FUN = function (x) { x$Z_stats[i] } ))
      }
      names(Z_obs) <- names(Z_perm) <- pairs_list

    }

    ## Filter NA for iterations without the tested state pairs
    Z_obs_filtered <- lapply(X = Z_obs, FUN = function (x) { x[!is.na(x)] } )
    Z_perm_filtered <- lapply(X = Z_perm, FUN = function (x) { x[!is.na(x)] } )

    ## Compute p-value for each pair of states
    p_values <- list()

    # Loop per pairs of states tested
    for (i in seq_along(Z_obs_filtered))
    {
      # i <- 1

      if (two_tailed)
      {
        # For two-tailed, compare differences of absolute values of Z-scores. alpha % should be higher than 0.
        # Ho: rate differences in ranks between states are equal in observed data and permuted data
        # Ha: rate differences in ranks between states are more extreme in observed data than in permuted data
        # P-value = frequency of cases where observed stats is less extreme
        # (closer from null hypothesis) than the permuted stats

        p_values[[i]] <- sum(abs(Z_obs_filtered[[i]]) <= abs(Z_perm_filtered[[i]])) / length(Z_perm_filtered[[i]])
      } else {
        # For one-tailed, compare differences of values of Z-scores. alpha % should be higher than 0
        # Ho: rate differences in ranks between states are lower or equal in observed data than in permuted data
        # Ha: rate differences in ranks between states are higher in observed data than in permuted data
        # P-value = frequency of cases where observed stats is lower than the permuted stats
        p_values[[i]] <- sum(Z_obs_filtered[[i]] <= Z_perm_filtered[[i]]) / length(Z_perm_filtered[[i]])
      }
    }
    p_values <- unlist(p_values)

    ## Adjust p-values for multiple comparisons
    if (two_tailed)
    {
      # For two-tailed, use the number of unique pairs for adjustment
      p_values_adjusted <- stats::p.adjust(p = p_values, method = p.adjust_method)
    } else {
      # For one-tailed, use the number of unique pairs by splitting p-values in two blocks
      # because reciprocal tests (A > B vs. B > A) are not independent tests.
      n_pairs <- length(p_values)
      p_values_forward <- p_values[1:(n_pairs/2)]
      p_values_backward <- p_values[((n_pairs/2)+1):n_pairs]
      p_values_forward <- stats::p.adjust(p = p_values_forward, method = p.adjust_method)
      p_values_backward <- stats::p.adjust(p = p_values_backward, method = p.adjust_method)
      p_values_adjusted <- c(p_values_forward, p_values_backward)
    }

    ## Save test stats
    estimates <- list()
    stats_median <- list()

    # Loop per pairs of states tested
    for (i in seq_along(Z_obs_filtered))
    {
      # i <- 1

      if (two_tailed)
      {
        # If two-tailed test, need to compare the abs_delta_Z with alpha % quantile to see if higher than zero.
        estimates[[i]] <- stats::quantile(abs(Z_obs_filtered[[i]]) - abs(Z_perm_filtered[[i]]), p = alpha)
        stats_median[[i]] <- stats::median(abs(Z_obs_filtered[[i]]) - abs(Z_perm_filtered[[i]]))
      } else {
        # If one-tailed test, need to compare the delta_Z with alpha % quantile to see if higher than zero.
        estimates[[i]] <- stats::quantile(as.numeric(Z_obs_filtered[[i]]) - as.numeric(Z_perm_filtered[[i]]), p = alpha)
        stats_median[[i]] <- stats::median(as.numeric(Z_obs_filtered[[i]]) - as.numeric(Z_perm_filtered[[i]]))
      }
    }

    ## Build summary df
    summary_df <- data.frame(pairs = pairs_list,
                             estimates = unlist(estimates),
                             stats_median = unlist(stats_median),
                             nb_test_stats = unlist(lapply(Z_obs_filtered, FUN = length)), # Number of test stats in each distribution. May differ depending on availability of trait states across stochastic maps
                             p_values = p_values,
                             p_values_adjusted = p_values_adjusted)
    row.names(summary_df) <- NULL

    ## Save test summary results
    STRAPP_results$posthoc_pairwise_tests$summary_df <- summary_df # Tests per pairs: estimates and p-values
    STRAPP_results$posthoc_pairwise_tests$method <- "Dunn" # Stats method
    STRAPP_results$posthoc_pairwise_tests$two_tailed <- two_tailed # Type of test: two-tailed or not

    ## Save permutation results in an array
    if (return_perm_data)
    {
      ## 3D-array to save permutation data
      # 1D = pairs
      # 2D = Stochastic map X BAMM posterior samples
      # 3D = Stats: trait_map_ID, BAMM_posterior_sample_ID, Z_obs, Z_perm, Z_delta/Z_abs_delta

      # Some combinations will be empty because the pairs was not found for a given stochastic map

      if (two_tailed)
      { # For two-tailed test, distribution based on difference in absolute Z-scores
        perm_data_array <- array(data = list(), # Need to use list to be able to store both character strings and numerical values
                                 dim = c(pairs = length(pairs_list), test_ID = length(Dunn_obs), stats = 5),
                                 dimnames = list(pairs = pairs_list, test_ID = 1:length(Dunn_obs), stats = c("trait_map_ID", "BAMM_posterior_sample_ID", "Z_obs", "Z_perm", "abs_delta_Z")))
      } else { # For one-tailed test, distribution based on difference in Z-scores
        perm_data_array <- array(data = list(), # Need to use list to be able to store both character strings and numerical values
                                 dim = c(pairs = length(pairs_list), test_ID = length(Dunn_obs), stats = 5),
                                 dimnames = list(pairs = pairs_list, test_ID = names(Dunn_obs), stats = c("trait_map_ID", "BAMM_posterior_sample_ID", "Z_obs", "Z_perm", "delta_Z")))
      }

      # Extract information on stochastic maps
      perm_data_array[,,"trait_map_ID"] <- matrix(data = sub(x = names(Dunn_obs), pattern = "_BAMM.*", replacement = ""),
                                                  nrow = dim(perm_data_array)[1],
                                                  ncol = dim(perm_data_array)[2],
                                                  byrow = TRUE)

      # Extract information on BAMM samples
      perm_data_array[,,"BAMM_posterior_sample_ID"] <- matrix(data = sub(x = names(Dunn_obs), pattern = ".*_BAMM", replacement = "BAMM"),
                                                              nrow = dim(perm_data_array)[1],
                                                              ncol = dim(perm_data_array)[2],
                                                              byrow = TRUE)

      # Extract Z_obs across pairs
      Z_obs_df <- do.call(rbind.data.frame, Z_obs)
      perm_data_array[,,"Z_obs"] <- as.matrix(Z_obs_df)

      # Extract Z_perm across pairs
      Z_perm_df <- do.call(rbind.data.frame, Z_perm)
      perm_data_array[,,"Z_perm"] <- as.matrix(Z_perm_df)

      # Extract delta-stats distribution
      if (two_tailed)
      { # For two-tailed test, distribution based on difference in absolute Z-scores
        perm_data_array[,,"abs_delta_Z"] <- as.matrix(abs(Z_obs_df) - abs(Z_perm_df))
      } else { # For one-tailed test, distribution based on difference in Z-scores
        perm_data_array[,,"delta_Z"] <- as.matrix((Z_obs_df - Z_perm_df))
      }

      # Store array
      STRAPP_results$posthoc_pairwise_tests$perm_data_array <- perm_data_array
    }
  }

  ## Export output
  return(STRAPP_results)
}



### Helper function to permute rates on tips using regime-blocks
# Input = data.frame with $regimes and $rates

#' @title Permutes rates on tips using regime-blocks
#'
#' @description Permutes rates on tips using regime membership to define
#'   blocks used for permutation as in a STRAPP test.
#'   Each block of tips assigned the a regime has the same rates at a given point in time.
#'   During permutation, each block of tips is assigned a unique rate from any regime drawn randomly.
#'
#' @param rates_regimes_df Data.frame with `$regimes` (integer) and `$rates` (numerical)
#'   found in each tip (one row per tip).
#'
#' @return The function returns a numerical vector with permuted rates across tips.
#'
#' @author Maël Doré, Dan Rabosky, Huateng Huang
#'
#' @references Rabosky, D. L. and Huang, H., 2015. A Robust Semi-Parametric Test for Detecting Trait-Dependent Diversification.
#'  Systematic Biology 65: 181-193. \doi{10.1093/sysbio/syv066}.
#'
#' @noRd
#'

block_permute_rates_data <- function (rates_regimes_df)
{
  # Extract regimes and rates data
  regimes <- rates_regimes_df$regimes
  rates <- rates_regimes_df$rates

  # Extract ID of unique regimes
  regimes_ID <- unique(regimes)

  # Extract rates for each regime
  # Rates are all equals in a given regime since they were extracted at the same focal_time
  rates_per_regimes <- numeric(length(regimes_ID))
  for (k in 1:length(regimes_ID))
  {
    rates_per_regimes[k] <- rates[regimes == regimes_ID[k]][1]
  }
  # Randomly rearrange regimes ID
  new_regimes_ID <- sample(regimes_ID, size = length(regimes_ID), replace = FALSE)

  # Replace rates with rates of newly assigned regimes = block-permutation based on regimes
  new_rates <- rep(0, length(regimes))
  for (k in 1:length(regimes_ID))
  {
    new_rates[which(regimes == regimes_ID[k])] <- rates_per_regimes[which(regimes_ID == new_regimes_ID[k])]
  }
  new_rates
}

