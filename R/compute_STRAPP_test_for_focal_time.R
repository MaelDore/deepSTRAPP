## Functions to compute STRAPP tests
# One master function to prepare data and select the proper test function according to data type
# There sub-functions carrying out tests according to data type

#' @title Computes STRAPP to test for a relationship between diversification rates and trait data
#'
#' @description Carries out the appropriate statistical method to test for a relationship between
#'   diversification rates and trait data for a given point in the past (i.e. the `focal_time`).
#'   Tests are based on block-permutations: rates data are randomized across tips following blocks
#'   defined by the diversification regimes identified on each tip (typically from a BAMM).
#'
#'   Such tests are called STructured RAte Permutations on Phylogenies (STRAPP) as described in
#'   Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification.
#'   Systematic biology, 65(2), 181-193. \url{https://doi.org/10.1093/sysbio/syv066}.
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
#'   that contains at least a `$trait_data` element, a `$focal_time` element, and a `$trait_data_type`.
#'   `$trait_data` is a named vector with the trait data found on the phylogeny at `focal_time`.
#'   `$focal_time` informs on the time in the past at which the trait and rates data will be tested.
#'   `$trait_data_type` informs on the type of trait data: continuous, categorical, or biogeographic.
#' @param rate_type A character string specifying the type of diversification rates to use. Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param nb_permutations Integer. To select the number of random permutations to perform during the tests.
#'   If NULL (default), all posterior samples will be used once.
#' @param replace_samples Logical. To specify whether to allow 'replacement' (i.e., multiple use) of a posterior sample when drawing samples used to carry out the test. Default is `FALSE`.
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
#'  and what significant level is used to rejected or not the null hypothesis. Default is `TRUE`.
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
#'   * Allow to choose if random sampling of posterior configurations must be done with replacement or not with `replace_samples`.
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
#' @return The function returns a list with at least eight elements.
#'
#'   Summary elements for the main test:
#'   * `$estimate` Named numerical. Value of the test statistic used to assess significance of the test
#'   according to the significance level provided (`alpha`). The test is significant if `$estimate` is higher than zero.
#'   * `$stats_median` Numerical. Median value of the distribution of test statistics across all selected posterior samples.
#'   * `$p-value` Numerical. P-value of the test. The test is considered significant if `$p-value` is lower than `alpha`.
#'   * `$method` Character string. The statistical method used to carry out the test.
#'   * `$rate_type` Character string. The type of diversification rates tested. One of 'speciation', 'extinction' or 'net_diversification'.
#'   * `$trait_data_type` Character string. The type of trait data as found in 'trait_data_list$trait_data_type'. One of 'continuous', 'categorical', or 'biogeographic'.
#'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. One of 'continuous', 'binary', or 'multinominal'.
#'   * `$focal_time` The time in the past at which the trait and rates data were tested.
#'
#'   If using continuous or binary data:
#'   * `$two-tailed` Logical. Record the type of test used: two-tailed if `TRUE`, one-tailed if `FALSE`.
#'   If `one_tailed_hypothesis` is provided (only for continuous and binary trait data):
#'   * `$one_tailed_hypothesis` Character string. Record of the alternative hypothesis used for the one-tailed tests.
#'
#'   If `posthoc_pairwise_tests = TRUE` (only for multinomial trait data):
#'   * `$posthoc_pairwise_tests` List of at least 3 sub-elements:
#'     + `$summary_df` Data.frame of five variables providing the summary results of post hoc pairwise tests
#'     + `$method` Character string. The statistical method used to carry out the test. Here, "Dunn".
#'     + `$two-tailed` Logical. Record the type of post hoc pairwise tests used: two-tailed if `TRUE`, one-tailed if `FALSE`.
#'
#'   If `return_perm_data = TRUE`, the stats data computed from the posterior samples for observed and permuted data are provided.
#'   This is needed to plot the histogram of the null distribution used to assess significance of the test with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'   * `$perm_data_df` A data.frame with four variables summarizing the data generated during the STRAPP test:
#'     + `$posterior_samples_random_ID` Integer. ID of the posterior samples randomly drawn and used for the STRAPP test.
#'     + `$*_obs` Numerical. Test stats computed from the observed data in the posterior samples. Name depends on the test used.
#'     + `$*_perm` Numerical. Test stats computed from the permuted data in the posterior samples. Name depends on the test used.
#'     + `$delta_*` OR `$abs_delta_*` Numerical. Test stats computed for the STRAPP test comparing observed stats and permuted stats.
#'       Name depends on the test used and the type of tests (two-tailed compare absolute values; one-tailed compare raw values).
#'   Combined with `posthoc_pairwise_tests = TRUE`, the stats data are also provided for the post hoc pairwise tests:
#'   * `$posthoc_pairwise_tests$perm_data_array` A 3D array containing stats data for all post hoc pairwise tests in a similar format that `$perm_data_df`.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#'
#'   Original function in BAMMtools: [BAMMtools::traitDependentBAMM()]
#'
#'   Statistical tests: [stats::cor.test()] [stats::wilcox.test()] [stats::kruskal.test()] [dunn.test::dunn.test()]
#'
#' @references For STRAPP: Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification.
#'   Systematic biology, 65(2), 181-193. \url{https://doi.org/10.1093/sysbio/syv066}.
#'
#'   For STRAPP in deep times: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'   Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'   \url{https://doi_of_Paper_to_provide.html}
#'
#' @examples
#' # ------ Prepare data ------ #
#'
#' ## Load the BAMM_object summarizing 1000 posterior samples of BAMM with diversification rates
#' # for ponerine ants extracted for 10My ago.
#' data(Ponerinae_BAMM_object_10My, package = "deepSTRAPP")
#'
#' # Plot the associated phylogeny with mapped rates
#' plot_BAMM_rates(Ponerinae_BAMM_object_10My)
#'
#' ## Load the object containing head width trait data for ponerine ants extracted for 10My ago.
#' data(Ponerinae_trait_data_10My, package = "deepSTRAPP")
#'
#' # Plot the associated contMap (continuous trait stochastic map)
#' plot(Ponerinae_trait_data_10My$contMap)
#'
#' # Check that objects are ordered in the same fashion
#' identical(names(Ponerinae_BAMM_object_10My$tipStates[[1]]),
#'           names(Ponerinae_trait_data_10My$trait_data))
#'
#' # Save continuous data
#' trait_data_continuous <- Ponerinae_trait_data_10My
#'
#' ## Transform trait data into binary and multinominal data
#'
#' # Binarize data into two states
#' trait_data_binary <- trait_data_continuous
#' trait_data_binary$trait_data[trait_data_continuous$trait_data < 0] <- "state_A"
#' trait_data_binary$trait_data[trait_data_continuous$trait_data >= 0] <- "state_B"
#' trait_data_binary$trait_data_type <- "categorical"
#'
#' table(trait_data_binary$trait_data)
#'
#' # Categorize data into three states
#' trait_data_multinominal <- trait_data_continuous
#' trait_data_multinominal$trait_data[trait_data_continuous$trait_data < 0] <- "state_B"
#' trait_data_multinominal$trait_data[trait_data_continuous$trait_data < -1] <- "state_A"
#' trait_data_multinominal$trait_data[trait_data_continuous$trait_data >= 0] <- "state_C"
#' trait_data_multinominal$trait_data_type <- "categorical"
#'
#' table(trait_data_multinominal$trait_data)
#'
#' \dontrun{  (May take a minute to run)
#' # ------ Compute STRAPP test for continuous data ------ #
#'
#' plot(x = trait_data_continuous$trait_data, y = Ponerinae_BAMM_object_10My$tipLambda[[1]])
#'
#' # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
#' # between "net_diversification" rates and trait data
#' STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_continuous,
#'    two_tailed = FALSE,
#'    one_tailed_hypothesis = "negative",
#'    return_perm_data = TRUE)
#' str(STRAPP_results, max.level = 2)
#' # Data from the posterior samples is available in STRAPP_results$perm_data_df
#' head(STRAPP_results$perm_data_df)
#'
#' # ------ Compute STRAPP test for binary data ------ #
#'
#' # Compute STRAPP test under the alternative hypothesis that "state_A" is associated
#' # with higher "net_diversification" that "state_B"
#' STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_binary,
#'    two_tailed = FALSE,
#'    one_tailed_hypothesis = c("state_A > state_B"))
#' str(STRAPP_results, max.level = 1)
#'
#' # Compute STRAPP test under the alternative hypothesis that "state_B" is associated
#' # with higher "net_diversification" that "state_A"
#' STRAPP_results <- compute_STRAPP_test_for_focal_time(BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_binary,
#'    two_tailed = FALSE,
#'    one_tailed_hypothesis = c("state_B > state_A"))
#' str(STRAPP_results, max.level = 1)
#'
#' # ------ Compute STRAPP test for multinominal data ------ #
#'
#' # Compute STRAPP test between all three states, and compute post hoc tests
#' # for differences in rates between all possible pairs of states
#' # with a p-value adjusted for multiple comparison using Bonferroni's correction
#' STRAPP_results <- compute_STRAPP_test_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    trait_data_list = trait_data_multinominal,
#'    posthoc_pairwise_tests = TRUE,
#'    two_tailed = TRUE,
#'    p.adjust_method = "bonferroni")
#' str(STRAPP_results, max.level = 3)
#' # All post hoc pairwise test summaries are available in $summary_df
#' STRAPP_results$posthoc_pairwise_tests$summary_df }
#'


### Master function to prepare data and select the proper test function according to data type ####

compute_STRAPP_test_for_focal_time <- function (BAMM_object, trait_data_list,
                                                rate_type = "net_diversification",
                                                nb_permutations = NULL,
                                                replace_samples = FALSE,
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
    ## BAMM_object
    # BAMM_object must be a 'bammdata' object
    if (!("bammdata" %in% class(BAMM_object)))
    {
      stop("'BAMM_object' must have the 'bammdata' class. See ?BAMMtools::getEventData() and ?deepSTRAPP::update_rates_and_regimes_for_focal_time() to learn how to generate those objects.")
    }
    # Number of posterior sample data must be equal between $tipStates, $tipLambda and $tipMu
    posterior_samples_length <- c(length(BAMM_object$tipStates), length(BAMM_object$tipLambda), length(BAMM_object$tipMu))
    if (length(unique(posterior_samples_length)) != 1)
    {
      stop("Number of posterior samples in 'BAMM_object' must be equal between $tipStates, $tipLambda and $tipMu.\nPlease check the structure of your 'BAMM_object' with str(BAMM_object, 1)")
    }
    # Number of branches in each posterior sample must be equal within $tipStates, $tipLambda and $tipMu
    tipStates_data_length <- unlist(lapply(X = BAMM_object$tipStates, FUN = length))
    if (length(unique(tipStates_data_length)) != 1)
    {
      stop("Number of branches in each posterior sample of 'BAMM_object$tipStates' must be equal.\nPlease check the structure of your 'BAMM_object' with str(BAMM_object$tipStates, 1)")
    }
    tipLambda_data_length <- unlist(lapply(X = BAMM_object$tipLambda, FUN = length))
    if (length(unique(tipLambda_data_length)) != 1)
    {
      stop("Number of branches in each posterior sample of 'BAMM_object$tipLambda' must be equal.\nPlease check the structure of your 'BAMM_object' with str(BAMM_object$tipLambda, 1)")
    }
    tipMu_data_length <- unlist(lapply(X = BAMM_object$tipMu, FUN = length))
    if (length(unique(tipMu_data_length)) != 1)
    {
      stop("Number of branches in each posterior sample of 'BAMM_object$tipMu' must be equal.\nPlease check the structure of your 'BAMM_object' with str(BAMM_object$tipMu, 1)")
    }
    # Number of branches in each posterior sample must be equal between $tipStates, $tipLambda and $tipMu
    posterior_samples_data_length <- c(unique(tipStates_data_length), unique(tipLambda_data_length), unique(tipMu_data_length))
    if (length(unique(posterior_samples_data_length)) != 1)
    {
      stop(paste0("Number of branches in posterior samples of 'BAMM_object$tipMu', 'BAMM_object$tipLambda', and 'BAMM_object$tipMu' must be equal.\nThere respective number of branches is: ",paste(posterior_samples_data_length, collapse = ", "),".\nPlease check the structure of your 'BAMM_object' with str(BAMM_object, 2)"))
    }

    ## trait_data_list
    # trait_data_list must be a list with $trait_data and $trait_data_type
    if (is.null(trait_data_list$trait_data) | is.null(trait_data_list$trait_data_type))
    {
      stop("'trait_data_list' must be a list with $trait_data and $trait_data_type elements.")
    }
    # $trait_data must be a named vector (can be numerical or character string)
    if (is.null(names(trait_data_list$trait_data)))
    {
      stop(paste0("'trait_data_list$trait_data' must be a named vector with names matching those found in BAMM_object$tipStates, BAMM_object$tipLambda, and BAMM_object$tipMu.\n",
                  "Names are either tip.label or tipward_node_ID of the branches cut at 'trait_data_list$focal_time' with deepSTRAPP::extract_most_likely_trait_values_for_focal_time."))
    }
    # $trait_data_type can only be "continuous", categorical" or "biogeographic"
    if (!(trait_data_list$trait_data_type %in% c("continuous", "categorical", "biogeographic")))
    {
      stop("'trait_data_list$trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
    }
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
      stop("Number of branches in 'trait_data_list$trait_data' must be equal to number of branches in posterior samples in 'BAMM_object$tipStates', 'BAMM_object$$tipLambda' and 'BAMM_object$tipMu'.\nPlease check the structure of your 'BAMM_object' with str(BAMM_object, 2)")
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

    ## nb_permutations
    # If nb_permutations is higher than number of posterior samples (length of $tipStates, $tipLambda and $tipMu) AND replace_samples = FALSE,
    # Send an error to say that replace_samples should be set to TRUE to allow multiple samplings of posterior in order to reach the requested number of permutations
    if (!is.null(nb_permutations))
    {
      if ((nb_permutations > unique(posterior_samples_length)) & !replace_samples)
      {
        stop(paste0("Number of permutations ('nb_permutations' = ",nb_permutations,") is higher than number of posterior samples in 'BAMM_object' (",unique(posterior_samples_length),").\n",
                    "'replace_samples' should be set to 'TRUE' so that multiple replicates of each posterior sample can be used in order to reach the requested number of permutations.\n",
                    "Alternatively, the 'nb_permutations' can be set to be equal or lower than the number of posterior samples. Default is to use the same number so each posterior is permutated once."))
      }
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

  ## Extract BAMM rates and regimes data
  BAMM_data <- list(tipStates = BAMM_object$tipStates, tipLambda = BAMM_object$tipLambda, tipMu = BAMM_object$tipMu)

  ## Extract trait data
  trait_data <- trait_data_list$trait_data

  ## Extract type of trait data
  trait_data_type <- trait_data_list$trait_data_type

  # If trait_data_type is "categorical" or "biogeographic", reclassify according to the number of states
  if (trait_data_type %in% c("categorical", "biogeographic"))
  {
    nb_levels <- nlevels(as.factor(trait_data))
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
  if (trait_data_type_for_stats != "multinominal" & posthoc_pairwise_tests == TRUE)
  {
    warning(paste0("There are only two states/ranges found at focal time = ",BAMM_object$focal_time,": ",paste(levels(as.factor(trait_data)), collapse = ", "),".\n",
                   "'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                   "If you want to test specific hypotheses with continuous or categorical binary data, use 'two_tailed = FALSE' and provide the 'one_tailed_hypothesis'.\n"))
  }

  ## Compute the appropriate internal function depending on the type of data

  switch(EXPR = trait_data_type_for_stats,
         continuous =   { # Case for continuous data
           # Stat test = Spearman's rank Rho test
           STRAPP_results <- compute_STRAPP_test_for_continuous_data(
             BAMM_data = BAMM_data,
             trait_data = trait_data,
             trait_data_type = trait_data_type,
             rate_type = rate_type,
             nb_permutations = nb_permutations,
             replace_samples = replace_samples,
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
             nb_permutations = nb_permutations,
             replace_samples = replace_samples,
             alpha = alpha,
             two_tailed = two_tailed,
             one_tailed_hypothesis = one_tailed_hypothesis,
             return_perm_data = return_perm_data,
             nthreads = nthreads,
             print_hypothesis = print_hypothesis)
         },
         multinominal = { # Case for multinominal data (Case of categorical/biogeographic data with more than two states)
           # Stat test = Kruskal-Wallis H test
           # Can define the post-hoc pairwise tests to compute
           STRAPP_results <- compute_STRAPP_test_for_multinominal_data(
             BAMM_data = BAMM_data,
             trait_data = trait_data,
             trait_data_type = trait_data_type,
             rate_type = rate_type,
             nb_permutations = nb_permutations,
             replace_samples = replace_samples,
             alpha = alpha,
             posthoc_pairwise_tests = posthoc_pairwise_tests, # See if I implement that for pairwise posthoc tests. Need to provide list of pairs with hypotheses
             two_tailed = two_tailed,
             p.adjust_method = p.adjust_method,
             return_perm_data = return_perm_data,
             nthreads = nthreads,
             print_hypothesis = print_hypothesis)
         }
  )

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
    nb_permutations = NULL,
    replace_samples = FALSE,
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

  # Randomly sample posteriors to use for each permutation
  posterior_samples_random_ID <- sample(x = 1:length(rates_data), size = nb_permutations, replace = replace_samples)

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:length(posterior_samples_random_ID))
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[posterior_samples_random_ID[l]]],
                                                           regimes = regimes_data[[posterior_samples_random_ID[l]]], stringsAsFactors = FALSE)
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

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })

  ## Print what is tested
  if (print_hypothesis)
  {
    if (two_tailed) # For two-tailed test
    {
      cat(paste0("Selected two-tailed Spearman's rank correlation test:\n\n",
                 "Null hypothesis: no correlation between trait data and diversification rates.\n\n",
                 "Alternative hypothesis: negative or positive correlation between trait data diversification rates.\n\n",
                 "'Estimate' stats is the ",alpha*100,"% quantile of abolute rho differences between observed and permuted data.\n",
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

  ## Unlist outputs
  rho_obs <- unlist(rho_obs)
  rho_perm <- unlist(rho_perm)

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
  STRAPP_results$p_value <- p_value # P-value of the test
  STRAPP_results$method <- "Spearman" # Stats method
  STRAPP_results$two_tailed <- two_tailed # Type of test: two-tailed or not
  STRAPP_results$one_tailed_hypothesis <- one_tailed_hypothesis # Type of hypothesis if one-tailed test
  STRAPP_results$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification
  STRAPP_results$trait_data_type <- trait_data_type # Type of trait data: continuous, categorical, or biogeographic
  STRAPP_results$trait_data_type_for_stats <- "continuous" # Type of trait data used to select statistical method: continuous, binary, or multinominal

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    perm_data_df <- data.frame(posterior_samples_random_ID = posterior_samples_random_ID,
                               rho_obs = as.numeric(rho_obs),
                               rho_perm = as.numeric(rho_perm))
    if (two_tailed)
    { # For two-tailed test, distribution based on difference in absolute correlations
      perm_data_df$abs_delta_rho <- abs(rho_obs) - abs(rho_perm)
    } else { # For one-tailed test, distribution based on difference in correlations
      perm_data_df$delta_rho <- as.numeric(rho_obs) - as.numeric(rho_perm)
    }
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
    nb_permutations = NULL,
    replace_samples = FALSE,
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

  # Randomly sample posteriors to use for each permutation
  posterior_samples_random_ID <- sample(x = 1:length(rates_data), size = nb_permutations, replace = replace_samples)

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:length(posterior_samples_random_ID))
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[posterior_samples_random_ID[l]]],
                                                           regimes = regimes_data[[posterior_samples_random_ID[l]]], stringsAsFactors = FALSE)
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

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })

  ## Prepare trait data

  # one_tailed_hypothesis <- c("state_A > state_B")

  obs_trait_states <- unique(trait_data)[order(unique(trait_data))]
  trait_states <- NA
  if (two_tailed)
  { # Case for two-tailed test

    if (print_hypothesis)
    {
      cat(paste0("Selected two-tailed Mann-Whitney-Wilcoxon rank-sum test:\n\n",
                 "Null hypothesis: taxa with trait '",
                 obs_trait_states[1], "' have equal ",rate_type," rates than those with trait '",
                 obs_trait_states[2], "'.\n",
                 "Alternative hypothesis: taxa with trait '",
                 obs_trait_states[1], "' have higher or lower ",rate_type," rates than those with trait '",
                 obs_trait_states[2], "'.\n\n",
                 "'Estimate' stats is the ",alpha*100,"% quantile of absolute differences in U-stats between observed and permuted data.\n",
                 "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
    }

  } else { # Case for one-tailed test

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
                   "Null hypothesis: taxa with trait '",
                   trait_states[1], "' have lower or equal ",rate_type," rates than those with trait '",
                   trait_states[2], "'.\n",
                   "Alternative hypothesis: taxa with trait '",
                   trait_states[1], "' have higher ",rate_type," rates than those with trait '",
                   trait_states[2],"'.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of differences in U-stats between observed and permuted data.\n",
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

  ## Center stats around location shift of the null hypothesis (mu)
  # Null hypothesis is that ranks of the values of the two groups are random
  # Compute location shift (mu) from state frequencies as average of the products of frequencies
  trait_data_counts <- table(trait_data)
  trait_data_counts <- trait_data_counts[!is.na(names(trait_data_counts))] # Remove NA
  stat_mu <- prod(trait_data_counts)/2
  # Center U-stats to get an estimate of how greater/lower (far away) than the null hypothesis (mu) are the calculated U-stats
  U_obs <- U_obs - stat_mu
  U_perm <- U_perm - stat_mu

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
  STRAPP_results$p_value <- p_value # P-value of the test
  STRAPP_results$method <- "Mann-Whitney-Wilcoxon" # Stats method
  STRAPP_results$two_tailed <- two_tailed # Type of test: two-tailed or not
  STRAPP_results$one_tailed_hypothesis <- one_tailed_hypothesis # Type of hypothesis if one-tailed test
  STRAPP_results$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification
  STRAPP_results$trait_data_type <- trait_data_type # Type of trait data: continuous, categorical, or biogeographic
  STRAPP_results$trait_data_type_for_stats <- "binary" # Type of trait data used to select statistical method: continuous, binary, or multinominal

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    perm_data_df <- data.frame(posterior_samples_random_ID = posterior_samples_random_ID,
                               U_obs = as.numeric(U_obs),
                               U_perm = as.numeric(U_perm))
    if (two_tailed)
    { # For two-tailed test, distribution based on difference in absolute U-stats
      perm_data_df$abs_delta_U <- abs(U_obs) - abs(U_perm)
    } else { # For one-tailed test, distribution based on difference in U-stats
      perm_data_df$delta_U <- as.numeric(U_obs) - as.numeric(U_perm)
    }
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
    nb_permutations = NULL,
    replace_samples = FALSE,
    alpha = 0.05,
    posthoc_pairwise_tests = FALSE,
    two_tailed = TRUE,
    p.adjust_method = "none",
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{

  ### Check input validity

  # Anything regarding posthoc_pairwise_tests, and two_tailed?

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

  # Randomly sample posteriors to use for each permutation
  posterior_samples_random_ID <- sample(x = 1:length(rates_data), size = nb_permutations, replace = replace_samples)

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:length(posterior_samples_random_ID))
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[posterior_samples_random_ID[l]]],
                                                           regimes = regimes_data[[posterior_samples_random_ID[l]]], stringsAsFactors = FALSE)
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

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })

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
      H_approximation <- stats::qchisq(p = 0.999, df = test_output$parameter)
      return(H_approximation)
    } else { # Otherwise, provide the computed H-stats
      return(test_output$statistic)
    }
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
  STRAPP_results$p_value <- p_value # P-value of the test
  STRAPP_results$method <- "Kruskal-Wallis" # Stats method
  STRAPP_results$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification
  STRAPP_results$trait_data_type <- trait_data_type # Type of trait data: continuous, categorical, or biogeographic
  STRAPP_results$trait_data_type_for_stats <- "multinominal" # Type of trait data used to select statistical method: continuous, binary, or multinominal

  ## Save permutation results
  if (return_perm_data)
  {
    perm_data_df <- data.frame(posterior_samples_random_ID = posterior_samples_random_ID,
                               H_obs = as.numeric(H_obs),
                               H_perm = as.numeric(H_perm),
                               delta_H = as.numeric(H_obs) - as.numeric(H_perm))
    STRAPP_results$perm_data_df <- perm_data_df
  }

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
                   "Null hypothesis: taxa with trait in the first state have lower or equal ",rate_type," rates than taxa with the second state in 'pairs'.\n",
                   "Alternative hypothesis: taxa with trait the first state have higher ",rate_type," rates than taxa with the second state in 'pairs'.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of differences in Z-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      }
    }

    ## Wrapped-up function to extract Z-stats from Dunn's post hoc pairwise rank-sum tests
    dunn_test <- function(rates, trait_data, two_tailed)
    {
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

      # If the test failed to provide a statistic because the value is reaching the ceiling for computation,
      # use the normal distribution and set an extremely high p-value to approximate a value
      if (any(is.na(test_output_df$Z_stats)))
      {
        Z_approximation <- stats::qnorm(p = 0.999)
        test_output_df$Z_stats[is.na(test_output_df$Z_stats)] <- Z_approximation

      }

      # Export df
      return(test_output_df)
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

    ## Extract Z-scores from outputs

    # Get list of pairs
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

    ## Compute p-value
    p_values <- list()

    # Loop per pairs of states tested
    for (i in seq_along(Z_obs))
    {
      # i <- 1

      if (two_tailed)
      {
        # For two-tailed, compare differences of absolute values of Z-scores. alpha % should be higher than 0.
        # Ho: rate differences in ranks between states are equal in observed data and permuted data
        # Ha: rate differences in ranks between states are more extreme in observed data than in permuted data
        # P-value = frequency of cases where observed stats is less extreme
        # (closer from null hypothesis) than the permuted stats

        p_values[[i]] <- sum(abs(Z_obs[[i]]) <= abs(Z_perm[[i]])) / length(Z_perm[[i]])
      } else {
        # For one-tailed, compare differences of values of Z-scores. alpha % should be higher than 0
        # Ho: rate differences in ranks between states are lower or equal in observed data than in permuted data
        # Ha: rate differences in ranks between states are higher in observed data than in permuted data
        # P-value = frequency of cases where observed stats is lower than the permuted stats
        p_values[[i]] <- sum(Z_obs[[i]] <= Z_perm[[i]]) / length(Z_perm[[i]])
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
    for (i in seq_along(Z_obs))
    {
      # i <- 1

      if (two_tailed)
      {
        # If two-tailed test, need to compare the abs_delta_Z with alpha % quantile to see if higher than zero.
        estimates[[i]] <- stats::quantile(abs(Z_obs[[i]]) - abs(Z_perm[[i]]), p = alpha)
        stats_median[[i]] <- stats::median(abs(Z_obs[[i]]) - abs(Z_perm[[i]]))
      } else {
        # If one-tailed test, need to compare the delta_Z with alpha % quantile to see if higher than zero.
        estimates[[i]] <- stats::quantile(as.numeric(Z_obs[[i]]) - as.numeric(Z_perm[[i]]), p = alpha)
        stats_median[[i]] <- stats::median(as.numeric(Z_obs[[i]]) - as.numeric(Z_perm[[i]]))
      }
    }

    ## Build summary df
    summary_df <- data.frame(pairs = pairs_list,
                             estimates = unlist(estimates),
                             stats_median = unlist(stats_median),
                             p_values = p_values,
                             p_values_adjusted = p_values_adjusted)

    ## Save test summary results
    STRAPP_results$posthoc_pairwise_tests$summary_df <- summary_df # Tests per pairs: estimates and p-values
    STRAPP_results$posthoc_pairwise_tests$method <- "Dunn" # Stats method
    STRAPP_results$posthoc_pairwise_tests$two_tailed <- two_tailed # Type of test: two-tailed or not

    ## Save permutation results in an array
    if (return_perm_data)
    {
      ## 3D-array to save permutation data
      # 1D = pairs
      # 2D = Posterior samples
      # 3D = Stats: Z_obs, Z_perm, Z_delta/Z_abs_delta

      if (two_tailed)
      { # For two-tailed test, distribution based on difference in absolute Z-scores
        perm_data_array <- array(data = NA,
                                 dim = c(pairs = length(pairs_list), posterior_samples = length(posterior_samples_random_ID), stats = 3),
                                 dimnames = list(pairs = pairs_list, posterior_samples = as.character(posterior_samples_random_ID), stats = c("Z_obs", "Z_perm", "abs_delta_Z")))
      } else { # For one-tailed test, distribution based on difference in Z-scores
        perm_data_array <- array(data = NA,
                                 dim = c(pairs = length(pairs_list), posterior_samples = length(posterior_samples_random_ID), stats = 3),
                                 dimnames = list(pairs = pairs_list, posterior_samples = as.character(posterior_samples_random_ID), stats = c("Z_obs", "Z_perm", "delta_Z")))
      }

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
#'  Systematic Biology 65: 181-193.  [https://doi.org/10.1093/sysbio/syv066].
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

