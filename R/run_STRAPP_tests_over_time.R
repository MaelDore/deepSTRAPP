
#' @title Runs a STRAPP to test for a relationship between diversification rates and trait data over multiple time steps
#'
#' @description Wrapper function to run STRAPP test workflows over multiple time steps in the past.
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
#'   Systematic biology, 65(2), 181-193. \url{https://doi.org/10.1093/sysbio/syv066}.
#'
#'   See the original [BAMMtools::traitDependentBAMM()] function used to
#'   carry out STRAPP test on extant time-calibrated phylogenies.
#'
#'   Tests can be carried out on speciation, extinction and net diversification rates.
#'
#' @param contMap Object of class `"contMap"`, typically generated with [phytools::contMap()],
#'   that contains a phylogenetic tree and associated continuous trait mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param ace Named numeric vector (Optional). Ancestral Character Estimates (ACE) of the internal nodes,
#'   typically generated with [phytools::fastAnc()], [phytools::anc.ML()], or [ape::ace()].
#'   Names are nodes_ID of the internal nodes. Values are ACE of the trait.
#'   Needed to provide accurate estimates of trait values.
#' @param tip_data Named numeric vector (Optional). Tip values of the trait.
#'   Names are nodes_ID of the internal nodes.
#'   Needed to provide accurate tip values.
#' @param trait_data_type Character string. Specify the type of trait data. Must be one of "continuous", "categorical", "biogeographic".
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [BAMMtools::getEventData()],
#'   that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#'   The phylogenetic tree must the same as the one associated with the `contMap`, `ace` and `tip_data`.
#' @param time_steps Numerical vector. Time steps at which the STRAPP tests should be carried out. If `NULL` (the default),
#'  `time_steps` will be generated from a combination of two arguments among `time_range`, `nb_time_steps`, and/or `time_step_duration`.
#' @param time_range Vector of two numerical values. Time boundaries within with the `time_steps` must be defined if not provided.
#'   If `NULL` (the default), and `time_range` is needed to generate the `time_steps`, the depth of the tree is used by default: `c(0, root_age)`.
#' @param nb_time_steps,time_step_duration Numerical. Number of time steps and duration of each time step used to generate `time_steps` if not provided.
#'   You must provide at least one of those two arguments to be able to generate `time_steps`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated phylogeny. Default is `TRUE`.
#' @param rate_type A character string specifying the type of diversification rates to use. Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param nb_permutations Integer. To select the number of random permutations to perform during the tests.
#'   If NULL (default), all posterior samples will be used once.
#' @param replace_samples Logical. To specify whether to allow 'replacement' (i.e., multiple use) of a posterior sample
#'   when drawing samples used to carry out the STRAPP test. Default is `FALSE`.
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
#'  and what significant level is used to rejected or not the null hypothesis. Default is `TRUE`.
#' @param extract_trait_data_melted_df Logical. Specify whether trait data must be extracted from the `updated_contMap` objects at each time step
#'  and returned in a melted data.frame. Default is `FALSE`.
#' @param extract_diversification_data_melted_df Logical. Specify whether diversification data (regimes ID and tip rates) must be extracted from the `updated_BAMM_object`
#'   at each time step and returned in a melted data.frame. Default is `FALSE`.
#' @param return_STRAPP_results Logical. Specify whether the `STRAPP_results` objects summarizing the results of the STRAPP tests carried out at each time step
#'   should be returned among the outputs in addition to the `$pvalues_summary_df` already providing test stat estimates and p-values obtained across all `time_steps`.
#' @param return_updated_trait_data_with_contMap Logical. Specify whether the `trait_data` extracted
#'   for the given `focal_time` and the updated version of mapped phylogeny (`contMap`) provided as input
#'   should be returned among the outputs. The updated `contMap` consists in cutting off branches and mapping
#'   that are younger than the `focal_time`. Default is `FALSE`.
#' @param return_updated_BAMM_object Logical. Specify whether the `updated_BAMM_object` with phylogeny and
#'   mapped diversification rates cut-off at the `focal_time` should be returned among the outputs.
#' @param verbose Logical. Should progression per `time_steps` be displayed? Default is `TRUE`.
#' @param verbose_extended Should progression per `time_steps` AND within each STRAPP workflow de displayed?
#'   In addition to printing progress along `time_steps`, a message will be printed at each step of the STRAPP workflow,
#'   and for every batch of 100 BAMM posterior samples whose rates are regimes are updated. If `extract_diversification_data_melted_df = TRUE`,
#'   a message for will also be printed when rates are extracted. Default is `FALSE`.
#'
#' @export
#' @importFrom phytools nodeHeights
#'
#' @details The function is a wrapper of [deepSTRAPP::run_STRAPP_test_for_focal_time()] that runs the
#'   STRAPP test workflow over multiple `time_steps`.
#'
#'   The STRAPP workflow is described step by step in the [deepSTRAPP::run_STRAPP_test_for_focal_time()] documentation.
#'
#'   Its main output is the `$pvalues_summary_df`: a data.frame providing test stat estimates and p-values obtained across all `time_steps`,
#'   that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
#'   If using multinominal data (with more than two states) and `posthoc_pairwise_tests = TRUE`, the output will also contain
#'   a data.frame providing test stat estimates and p-values for post hoc pairwise tests in `$pvalues_summary_df_for_posthoc_pairwise_tests`.
#'
#'   The function offers options to generate summary data.frames of the data extracted across `time_steps`:
#'   * If `extract_trait_data_melted_df = TRUE`, a data.frame of trait values found along branches at each time-step
#'     is provided in `$trait_data_df_over_time`.
#'   * If `extract_diversification_data_melted_df = TRUE`, a data.frame of diversification data (regimes ID and tip rates)
#'     found along branches at each time-step is provided in `$diversification_data_df_over_time`.
#'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
#'     the evolution diversification rates across trait values over time.
#'
#'   The function also allows to keep records of the intermediate objects generated during the STRAPP workflow:
#'   * If `return_STRAPP_results = TRUE`, a list of STRAPP test outputs is provided in `$STRAPP_results_over_time`.
#'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
#'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
#'   * If `return_updated_trait_data_with_contMap = TRUE`, a list of objects containing trait data and updated `contMap`
#'     is provided in `$updated_trait_data_with_contMap_over_time`. Updated `contMap` can be plotted with [phytools::plot.contMap()]
#'     to display a phylogeny mapped with trait values with branches cut at each `focal_time`.
#'   * If `return_updated_BAMM_object = TRUE`, a list of updated `BAMM_object` of class `"bammdata"` that contains rates and regimes ID
#'     found at each `focal_time`. Updated `BAMM_object` can be plotted with [BAMMtools::plot.bammdata()] to display
#'     a phylogeny mapped with diversification rates with branches cut at each `focal_time`.
#'
#' @return The function returns a list with at least five elements.
#'
#'   * `$pvalues_summary_df` Data.frame with three columns providing test stat `$estimate` and `$p_value` obtained for each time-step (i.e., `$focal_time`),
#'     that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
#'   * `$time_steps` Numerical vector. Time steps at which the STRAPP tests were carried out in the same order as the objects returned in the output lists.
#'   * `$trait_data_type` Character string. Specify the type of trait data. Possible values are: "continuous", "categorical", "biogeographic".
#'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. One of 'continuous', 'binary', or 'multinominal'.
#'   * `$rate_type` Character string. The type of diversification rates used in the tests: 'speciation', 'extinction' or 'net_diversification'.
#'
#'   Optional summary df for multinominal data, if `posthoc_pairwise_tests = TRUE`:
#'   * `$pvalues_summary_df_for_posthoc_pairwise_tests` Data.frame with four or five columns providing test stat `$estimate`, `$p_value`, and `$p_value_adjusted`
#'     (if `p.adjust_method` used is not "none") for each `$pair` of states involved in post hoc Dunn's tests obtained for each time-step (i.e., `$focal_time`).
#'      This data.frame can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the post hoc test results across time.
#'
#'   Optional melted data.frames:
#'   * `$trait_data_df_over_time` Data.frame with three columns providing `$trait_value` associated with each `$tip_ID` found along each time-step (i.e., `$focal_time`).
#'     Set `extract_trait_data_melted_df = TRUE` to include it in the output.
#'   * `$diversification_data_df_over_time` Data.frame with six columns providing diversification regimes (`$regime_ID`) and `$rates` sorted by `$rate_type` along tips (`$tip_ID`)
#'     found across all posterior samples (`$BAMM_sample_ID`) over each time-step (i.e., `$focal_time`).
#'     Set `extract_diversification_data_melted_df = TRUE` to include it in the output.
#'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
#'     the evolution diversification rates across trait values over time.
#'
#'   Optional objects generated for each time-step (i.e., `focal_time`) and ordered as in `$time_steps`:
#'   * `$STRAPP_results_over_time` List of objects summarizing the results of the STRAPP tests
#'     See [compute_STRAPP_test_for_focal_time()] for a detailed description of the elements in each object.
#'     Set `return_STRAPP_results = TRUE` to include it in the output.
#'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
#'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
#'   * `$updated_trait_data_with_contMap_over_time` List of objects containing trait data and updated `contMap`.
#'     Updated `contMap` can be plotted with [phytools::plot.contMap()] to display a phylogeny mapped with trait values
#'     with branches cut at each `focal_time`.
#'   * `$updated_BAMM_objects_over_time` List of objects containing rates and regimes ID mapped on phylogeny.
#'     Updated `BAMM_object` can be plotted with [BAMMtools::plot.bammdata()] to display a phylogeny mapped with
#'     diversification rates with branches cut at each `focal_time`.
#'
#' @author Maël Doré
#'
#' @seealso To run STRAPP workflow for a single `focal_time`: [deepSTRAPP::run_STRAPP_test_for_focal_time()]
#' [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#' [deepSTRAPP::extract_diversification_data_melted_df_for_focal_time()] [deepSTRAPP::compute_STRAPP_test_for_focal_time()]
#'
#' @examples
#' ## Load data
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
#' plot(Ponerinae_contMap)
#'
#' ## Set for three time steps of 5 My. Will generate STRAPP workflows for 0, 5 and 10 Mya.
#' nb_time_steps <- 3
#' time_step_duration <- 5
#'
#' \dontrun{  (May take several minutes to run)
#' ## Run STRAPP test on net diversification rates
#'
#' STRAPP_tests_over_time <- run_STRAPP_tests_over_time(
#'    contMap = Ponerinae_contMap,
#'    ace = Ponerinae_ACE, tip_data = Ponerinae_data_ln_HW,
#'    trait_data_type = "continuous",
#'    BAMM_object = Ponerinae_BAMM_object,
#'    nb_time_steps = nb_time_steps,
#'    time_step_duration = time_step_duration,
#'    return_perm_data = TRUE,
#'    extract_trait_data_melted_df = TRUE,
#'    extract_diversification_data_melted_df = TRUE,
#'    return_STRAPP_results = TRUE,
#'    return_updated_trait_data_with_contMap = TRUE,
#'    return_updated_BAMM_object = TRUE,
#'    verbose = TRUE,
#'    verbose_extended = TRUE)
#'
#' ## Explore output
#' str(STRAPP_tests_over_time, max.level = 1)
#'
#' # Display test summary
#' # Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
#' # showing the evolution of the test results across time.
#' STRAPP_tests_over_time$pvalues_summary_df
#'
#' # Access trait data in a melted data.frame
#' head(STRAPP_tests_over_time$trait_data_df_over_time)
#'
#' # Access the diversification data in a melted data.frame
#' head(STRAPP_tests_over_time$diversification_data_df_over_time)
#'
#' # Access STRAPP test results
#' str(STRAPP_tests_over_time$STRAPP_results, max.level = 2)
#'
#' # Plot updated contMap for time-step n°2
#' contMap_2 <- STRAPP_tests_over_time$updated_trait_data_with_contMap_over_time[[2]]
#' phytools::plot.contMap(contMap_2$contMap)
#' ape::nodelabels(text = contMap_2$contMap$tree$initial_nodes_ID)
#'
#' # Plot diversification rates on updated phylogeny for time-step n°2
#' BAMMtools::plot.bammdata(STRAPP_tests_over_time$updated_BAMM_objects_over_time[[2]], labels = TRUE)
#'
#' # Plot histogram of test stats for time-step n°2
#' plot_histogram_STRAPP_test_for_focal_time(
#'    STRAPP_results = STRAPP_tests_over_time$STRAPP_results_over_time[[2]]) }
#'


run_STRAPP_tests_over_time <- function (contMap, # Add densityMaps as alternative input and adjust doc
                                        ace = NULL,
                                        tip_data = NULL,
                                        trait_data_type,
                                        BAMM_object,
                                        time_steps = NULL,
                                        time_range = NULL,
                                        nb_time_steps = NULL,
                                        time_step_duration = NULL,
                                        keep_tip_labels = TRUE,
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
                                        print_hypothesis = TRUE,
                                        extract_trait_data_melted_df = FALSE,
                                        extract_diversification_data_melted_df = FALSE,
                                        return_STRAPP_results = FALSE,
                                        return_updated_trait_data_with_contMap = FALSE,
                                        return_updated_BAMM_object = FALSE,
                                        verbose = TRUE,
                                        verbose_extended = FALSE)
{
  ### Check input validity
  {
    # Should all already be included in the wrapped functions
  }

  ## Extract root_age
  root_age <- max(phytools::nodeHeights(contMap$tree)[,2])

  ## If provided, check that time_range is two positive numerical values, in increasing order, not above root_age.
  if (!is.null(time_range))
  {
    if (!all(time_range >= 0))
    {
      stop("'time_range' must be two positive numerical values.\n")
    } else {
      if (time_range[1] > time_range[2])
      {
        stop("'time_range' values must be provided in increasing order.\n")
      } else {
        if (time_range[2] > root_age)
        {
          stop("'time_range' older boundary must be equal or lower than the age of the root of the tree.\n",
               "Max time_range = ",max(time_range),". Root age = ",root_age, ".\n",)
        }
      }
    }
  }

  ## If provided, check that nb_time_steps is a positive integer (warning over 100 against time consuming option).
  if (!is.null(nb_time_steps))
  {
    if (!(round(nb_time_steps) == nb_time_steps) | !(nb_time_steps > 0))
    {
      stop("'nb_time_steps' must be a positive integer.\n")
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
      if (time_step_duration > root_age)
      {
        stop("'time_step_duration' must be a smaller than the root of the tree.\n",
             "'time_step_duration' = ",time_step_duration,". Root age = ",root_age,".\n",)
      }
    }
  }

  ## If time_step_duration AND time_range provided, check that time_step_duration is smaller than the range of time_range.
  if (!is.null(time_step_duration) && !is.null(time_range))
  {
    if (time_step_duration > diff(time_range))
    {
      stop("'time_step_duration' must be smaller than the range of 'time_range'.\n",
           "'time_step_duration' = ",time_step_duration,". 'time_range' = ",time_range,".\n",)
    }
  }

  ## If return_perm_data = TRUE, but return_STRAPP_results = FALSE, the distribution of the test statistics will not be returned in the output
  if (return_perm_data && !return_STRAPP_results)
  {
    stop("'return_perm_data' = TRUE but 'return_STRAPP_results' = FALSE.\n",
         "The distribution of the test statistics will not be returned in the output if 'return_STRAPP_results' = FALSE.\n",
         "Set 'return_STRAPP_results' = TRUE if you want to retrieve the test statistics in order to plot histograms of the STRAPP test results.\n")
  }
  # Only useful if return_STRAPP_results = TRUE

  ### Check validity of time-steps arguments while generating the time-steps vector

  if (is.null(time_steps))
  {
    warning("'time_steps' were not provided. Instead, they were tentatively generated from a combination of 'time_range', 'nb_time_steps', and/or 'time_step_duration'.\n")

    if (!is.null(time_range) && !is.null(nb_time_steps) && !is.null(time_step_duration))
    {
      stop("In order to generate the 'time_steps', please provide only two of the three arguments: 'time_range', 'nb_time_steps' and/or 'time_step_duration'.\n")
    }

    if (is.null(nb_time_steps) && is.null(time_step_duration))
    {
      stop("In order to generate the 'time_steps', you must provide at least one of those two arguments: 'nb_time_steps' or 'time_step_duration'.\n")
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
                    "Max time-step = ",max(time_steps),". Root age = ",root_age, ".\n",
                    "You must provide 'nb_time_steps' and 'time_step_duration' arguments that produce time-steps that are compatible with the root age.\n"))
      }

    } else {
      # Cases with 'nb_time_steps' OR 'time_step_duration'

      # Check time_range
      # If absent, use c(0, root_age) and send warning
      if (is.null(time_range))
      {
        warning("'time_range' was not provided. The depth of the tree was used by default: c(0, root_age).\n")
        time_range <- c(0, root_age)
      }

      # Case with 'nb_time_steps' and no 'time_step_duration'
      if (!is.null(nb_time_steps))
      {
        # Produce time_steps from time_range & nb_time_steps
        time_steps <- seq(from = time_range[1], to = time_range[2], length.out = nb_time_steps)
      } else {
        # Case with no 'nb_time_steps' and 'time_step_duration'
        # Produce time_steps from time_range & time_step_duration
        time_steps <- seq(from = time_range[1], to = time_range[2], by = time_step_duration)
      }
    }

    # Print time steps generated
    if (verbose | verbose_extended)
    {
      cat(paste0("Time-steps generated from a combination of 'time_range', 'nb_time_steps', and/or 'time_step_duration': ",paste(time_steps, collapse = ", "),".\n\n"))
    }
  }

  ## Check if there is more than one time-step
  if (length(time_steps) < 2)
  {
    stop(paste0("Only one time-step was provided/generated.\n",
                "Please use run_STRAPP_test_for_focal_time() to run the STRAPP workflow for a unique time-step/focal time\n\n"))
  }

  ## Detect if trait_data is needed because requested, or to extract trait data in a melted df
  need_trait_data <- (return_updated_trait_data_with_contMap | extract_trait_data_melted_df)

  ### Run STRAPP test workflow per time-steps
  STRAPP_test_outputs_over_time <- list()
  for (i in seq_along(time_steps))
  {
    # i <- 1

    # Extract current focal_time
    focal_time_i <- time_steps[i]

    # Print progress per time-steps
    if (verbose | verbose_extended)
    {
      cat(paste0(Sys.time(), " - STRAPP test running for focal_time = ",focal_time_i," - n\u00B0 ", i, "/", length(time_steps),"\n\n"))
    }

    ## Run STRAPP test workflow for the current focal_time
    STRAPP_test_output_i <- run_STRAPP_test_for_focal_time(
      contMap = contMap,
      ace = ace,
      tip_data = tip_data,
      trait_data_type = trait_data_type,
      BAMM_object = BAMM_object,
      focal_time = focal_time_i,
      keep_tip_labels = keep_tip_labels,
      rate_type = rate_type,
      nb_permutations = nb_permutations,
      replace_samples = replace_samples,
      alpha = alpha,
      two_tailed = two_tailed,
      one_tailed_hypothesis = one_tailed_hypothesis,
      posthoc_pairwise_tests = posthoc_pairwise_tests,
      p.adjust_method = p.adjust_method,
      return_perm_data = return_perm_data,
      nthreads = nthreads,
      print_hypothesis = print_hypothesis,
      extract_diversification_data_melted_df = extract_diversification_data_melted_df,
      return_updated_trait_data_with_contMap = need_trait_data,
      return_updated_BAMM_object = return_updated_BAMM_object,
      verbose = verbose_extended)

    ## Store output for current focal time
    STRAPP_test_outputs_over_time[[i]] <- STRAPP_test_output_i
  }


  ## Initiate list for final output
  final_ouput <- list()

  ## Extract p_values and estimates in summary_df
  p_values <- unlist(lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$STRAPP_results$p_value} ))
  estimates <- unlist(lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$STRAPP_results$estimate} ))

  pvalues_summary_df <- data.frame(focal_time = time_steps,
                                   estimate = estimates,
                                   p_value = p_values)

  # Store pvalues_summary_df in final output
  final_ouput$pvalues_summary_df <- pvalues_summary_df

  ## Store time steps in final output
  final_ouput$time_steps <- time_steps

  ## Store type of trait data tested
  final_ouput$trait_data_type <- trait_data_type

  ## Store the type of trait data used to select the appropriate statistical method
  trait_data_type_for_stats_list <- unlist(lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$STRAPP_results$trait_data_type_for_stats } ))
  final_ouput$trait_data_type_for_stats <- unique(trait_data_type_for_stats_list)

  ## Store type of diversification rates tested
  final_ouput$rate_type <- rate_type

  ## Extract p_values and estimates in summary_df for post hoc tests
  if (posthoc_pairwise_tests)
  {
    # Extract summary_df from $STRAPP_results
    pvalues_summary_df_for_posthoc_pairwise_tests <- lapply(
      X = STRAPP_test_outputs_over_time,
      FUN = function (x) { x$STRAPP_results$posthoc_pairwise_tests$summary_df } )

    # Add focal_time and rename columns to match with $pvalues_summary_df
    for (i in seq_along(pvalues_summary_df_for_posthoc_pairwise_tests))
    {
      pvalues_summary_df_i <- pvalues_summary_df_for_posthoc_pairwise_tests[[i]]
      pvalues_summary_df_i$focal_time <- STRAPP_test_outputs_over_time[[i]]$STRAPP_results$focal_time
      pvalues_summary_df_i <- pvalues_summary_df_i[, c("focal_time", "pairs", "estimates", "p_values", "p_values_adjusted")]
      names(pvalues_summary_df_i) <- c("focal_time", "pair", "estimate", "p_value", "p_value_adjusted")
      pvalues_summary_df_for_posthoc_pairwise_tests[[i]] <- pvalues_summary_df_i
    }
    # Bind across all time-steps
    pvalues_summary_df_for_posthoc_pairwise_tests <- do.call(what = rbind, args = pvalues_summary_df_for_posthoc_pairwise_tests)

    # Store pvalues_summary_df in final output
    final_ouput$pvalues_summary_df_for_posthoc_pairwise_tests <- pvalues_summary_df_for_posthoc_pairwise_tests
  }

  ## Extract melted data.frame with trait data if requested
  # Can be used to plot rates through time ~ trait values
  if (extract_trait_data_melted_df)
  {
    # Extract trait data in a list
    trait_data_over_time <- lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$updated_trait_data_with_contMap$trait_data } )
    # Convert to data.frame in each focal_time
    for (i in seq_along(time_steps))
    {
      trait_data_over_time_i <- trait_data_over_time[[i]]
      trait_data_df_i <- data.frame(focal_time = time_steps[i], tip_ID = names(trait_data_over_time_i), trait_value = trait_data_over_time_i)
      row.names(trait_data_df_i) <- NULL
      trait_data_over_time[[i]] <- trait_data_df_i
    }
    # Bind across all time-steps
    trait_data_df_over_time <- do.call(what = rbind, args = trait_data_over_time)

    # Store trait_data_df_over_time in final output
    final_ouput$trait_data_df_over_time <- trait_data_df_over_time
  }

  ## Extract melted data.frame with rates/regimes data if requested
  # Can be used to plot rates through time ~ trait values
  if (extract_diversification_data_melted_df)
  {
    diversification_data_over_time <- lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$diversification_data_df } )
    # Bind across all time-steps
    diversification_data_df_over_time <- do.call(what = rbind, args = diversification_data_over_time)

    # Store STRAPP_results in final output
    final_ouput$diversification_data_df_over_time <- diversification_data_df_over_time
  }

  ## Extract STRAPP_results if requested
  # Can be used to plot histograms of STRAPP test statistics for each focal_time
  if (return_STRAPP_results)
  {
    STRAPP_results_over_time <- lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$STRAPP_results } )

    # Store STRAPP_results in final output
    final_ouput$STRAPP_results_over_time <- STRAPP_results_over_time
  }

  ## Extract updated contMap/simmaps with associated trait_data if requested
  # Can be used to plot an updated contMap for each focal_time
  if (return_updated_trait_data_with_contMap)
  {
    updated_trait_data_with_contMap_over_time <- lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$updated_trait_data_with_contMap } )

    # Store STRAPP_results in final output
    final_ouput$updated_trait_data_with_contMap_over_time <- updated_trait_data_with_contMap_over_time
  }

  ## Extract updated BAMM_object with associated rates and regimes if requested
  # Can be used to plot diversification rates on updated an phylogeny for each focal_time
  if (return_updated_BAMM_object)
  {
    updated_BAMM_objects_over_time <- lapply(X = STRAPP_test_outputs_over_time, FUN = function (x) { x$updated_BAMM_object } )

    # Store STRAPP_results in final output
    final_ouput$updated_BAMM_objects_over_time <- updated_BAMM_objects_over_time
  }

  ## Export final output
  return(final_ouput)

}

