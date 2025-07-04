
#' @title Run a full BAMM (Bayesian Analysis of Macroevolutionary Mixtures) workflow
#'
#' @description Run a full BAMM (Bayesian Analysis of Macroevolutionary Mixtures) workflow
#'   to produce a `BAMM_object` that contains a phylogenetic tree and associated diversification rates
#'   mapped along branches, across selected posterior samples:
#'
#'   * Step 1: Set BAMM - Record BAMM settings and generate all input files needed for BAMM.
#'   * Step 2: Run BAMM - Run BAMM and move output files in dedicated directory.
#'   * Step 3: Evaluate BAMM - Produce evaluation plots and ESS data.
#'   * Step 4: Import BAMM outputs - Load `BAMM_object` in R and subset posterior samples.
#'   * Step 5: Clean BAMM files - Remove files generated during the BAMM run.
#'
#'   The `BAMM_object` output is typically used as input to run deepSTRAPP with [deepSTRAPP::run_deepSTRAPP_for_focal_time()]
#'   or [deepSTRAPP::run_deepSTRAPP_over_time()]. Diversification rates and regimes shift can be visualized with [deepSTRAPP::plot_BAMM_rates()].
#'
#'   BAMM is a model of diversification for time-calibrated phylogenies that explores complex diversification dynamics
#'   by allowing multiple regime shifts across clades without a priori hypotheses on the location of such shifts.
#'   It uses reversible jump Markov chain Monte Carlo (rjMCMC) to automatically explore a vast range of models with different
#'   speciation and extinction rates, and different number and location of regime shits.
#'
#'   This function will work only if you have the BAMM C++ program installed in your machine.
#'   See the BAMM website: \url{http://bamm-project.org/} and the companion R package [BAMMtools].
#'
#' @param BAMM_install_directory_path Character string. The path to the directory where BAMM is.
#'   Use '/' to separate directory and sub-directories. The path must end with '/'.
#' @param phylo Time-calibrated phylogeny. Object of class `"phylo"` as defined in R package [ape]. The phylogeny must be rooted and fully resolved.
#'   BAMM does not currently work with fossils, so the tree must also be ultrametric.
#' @param prefix_for_files Character string. Prefix to add to all BAMM files stored in the `BAMM_output_directory_path` if `keep_BAMM_outputs = TRUE`.
#'   Files will be exported such as 'prefix_*' with an underscore separating the prefix and the file name. Default is `NULL` (no prefix is added).
#' @param seed Integer. Set the seed to ensure reproducibility. Default is `NULL` (a random seed is used).
#' @param numberOfGenerations Integer. Number of steps in the MCMC run. It should be set high enough to reach the equilibrium distribution
#'  and allows posterior samples to be uncorrelated. Check the Effective Sample Size of parameters with coda::effectiveSize() in the Evaluation step.
#'  Default value is `10^7`.
#' @param globalSamplingFraction Numerical. Global sampling fraction representing the overall proportion of terminals in the phylogeny compared to
#'  the estimated overall richness in the clade. It acts as a multiplier on the rates needed to achieve such extant diversity.
#'  Default is `1.0` (assuming all taxa are in the phylogeny).
#' @param sampleProbsFilename Character string. The path to the `.txt` file used to provide clade-specific sampling fractions.
#'  See [BAMMtools::samplingProbs()] to generate such file. If provided, `globalSamplingFraction` is ignored.
#' @param expectedNumberOfShifts Integer. Set the expected number of regime shifts. It acts as an hyperparameter controlling the exponential prior distribution
#'  used to modulate reversible jumps across model configurations in the rjMCMC run.
#'  If set to `NULL` (default), an empirical rule will be used to define this value: 1 regime shift expected for every 100 tips in the phylogeny, with a minimum of 1.
#'  The best practice consists in trying several values and inspect the similarity of the prior and posterior distribution of the regime shift parameter.
#'  See [BAMMtools::plotPrior()] and the Evaluation step to produce such evaluation plot.
#' @param eventDataWriteFreq Integer. Set the frequency in which to write the event data to the output file = the sampling frequency of posterior samples.
#'  If set to `NULL` (default), will set frequency such as 2000 posterior samples are recorded such as `eventDataWriteFreq = numberOfGenerations / 2000`.
#' @param burn_in Numerical. Proportion of posterior samples removed from the BAMM output to ensure that the remaining samples where drawn once the equilibrium distribution was reached.
#'  This can be evaluated looking at the MCMC trace (see Evaluation step). Default is `0.25`.
#' @param nb_posterior_samples Numerical. Number of posterior samples to extract, after removing the burn-in, in the final `BAMM_object` to use for downstream analyses.
#'  Default = 1000.
#' @param additional_BAMM_settings List of named elements. Additional settings options for BAMM provided as a list of named arguments.
#'  Ex: `list(lambdaInit0 = 0.5, muInit0 = 0)`. See available settings in the template file provided within the deepSTRAPP package files as 'BAMM_template_diversification.txt'.
#'  The template can also be loaded directly in R with `utils::data(BAMM_template_diversification)` and displayed with `print(BAMM_template_diversification)`.
#' @param BAMM_output_directory_path Character string. The path to the directory used to store input/output files generated.
#' Use '/' to separate directory and subdirectories. It must end with '/'. Default is `./BAMM_outputs/`
#' @param keep_BAMM_outputs Logical. Whether the `BAMM_output_directory` should be kept after the run. Default = `TRUE`.
#' @param MAP_odd_ratio_threshold Numerical. Controls the definition of 'core-shifts' used to distinguish across configurations when fetching the MAP samples.
#'   Shifts that have an odd-ratio of marginal posterior probability / prior lower than `MAP_odd_ratio_threshold` are ignored. See [BAMMtools::getBestShiftConfiguration()]. Default = `5`.
#' @param skip_evaluations Logical. Whether to skip the Evaluation step including MCMC trace, ESS, and prior/posterior comparisons for expected number of shifts. Default = `FALSE`.
#' @param plot_evaluations Logical. Whether to display the plots generated during the Evaluation step: MCMC trace, and prior/posterior comparisons for expected number of shifts. Default = `TRUE`.
#' @param save_evaluations Logical. Whether to save the outputs of evaluations in a table (ESS), and PDFs (MCMC trace, and prior/posterior comparisons for expected number of shifts)
#'  in the `BAMM_output_directory`. Default = `TRUE`.
#'
#' @export
#' @importFrom ape write.tree
#' @importFrom stringr str_detect str_remove
#' @importFrom utils read.csv write.csv data
#' @importFrom ggplot2 ggplot geom_line geom_point geom_vline aes labs ggtitle theme element_line element_rect element_text unit margin
#' @importFrom grDevices pdf dev.off
#' @importFrom coda effectiveSize
#' @importFrom cowplot save_plot
#' @importFrom BAMMtools setBAMMpriors plotPrior plot.bammdata samplingProbs getEventData subsetEventData
#'
#' @details This function runs a full BAMM (Bayesian Analysis of Macroevolutionary Mixtures) workflow
#'   to produce a `BAMM_object` that contains a phylogenetic tree and associated diversification rates
#'   mapped along branches, across selected posterior samples.
#'
#'  Step 1: Set BAMM
#'   * Produces a tree file for the phylogeny. Default file: 'phylogeny.tree'.
#'   * Save configuration settings used for the BAMM run. Default file: 'config_file.txt'.
#'   * Save default priors generated by [BAMMtools::setBAMMpriors] based on the phylogeny. Default file: 'priors.txt'.
#'
#'  Step 2: Run BAMM
#'   * Run BAMM using the system console
#'   * Move output files in dedicated `BAMM_output_directory`. Default directory is `./BAMM_outputs/`.
#'     - 'run_info.txt' containing a summary of your parameters/settings.
#'     - 'mcmc_log.txt' containing raw MCMC information useful in diagnosing convergence.
#'     - 'event_data.txt' containing all evolutionary rate parameters and their topological mappings.
#'     - 'chain_swap.txt' containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).
#'     - 'acceptance_info.txt' containing the history of acceptance/proposal of MCMC steps (If additional setting `outputAcceptanceInfo` is set to 1).
#'
#'  Step 3: Evaluate BAMM
#'   * Plot the MCMC trace = evolution of logLik across MCMC generations. Output file = 'MCMC_trace_logLik.pdf'.
#'   * Compute the Effective Sample Size (ESS) across posterior samples (after removing burn-in) using [coda::effectiveSize()].
#'     This is a way to evaluate if your MCMC runs has enough generations to produce robust estimates. Ideally, ESS should be higher than 200.
#'     Output file = 'ESS_df.csv'.
#'   * Plot the comparison of prior and posterior distributions of the number of regime shifts with [BAMMtools::plotPrior].
#'     Output file = 'PP_nb_shifts_plot.pdf'.
#'     A good value for `expectedNumberOfShifts` is one with high similarities between the distributions
#'     hinting that the information in the data coincides with your expectations for the number of regime shifts.
#'     The best practice consists in trying several values to control if it affects or not the final output.
#'
#'  Step 4: Import BAMM outputs
#'   * Load BAMM outputs with [BAMMtools::getEventData].
#'   * Subset posterior samples to the requested `nb_posterior_samples` with [BAMMtools::subsetEventData].
#'   * Record the `$expectedNumberOfShifts` used to set the prior. This is useful for downstream analyses involving comparison of prior vs. posterior probabilities
#'     (See [BAMMtools::distinctShiftConfigurations()]).
#'   * Record the marginal posterior probability of regime shift along branches based on the proportion of samples harboring a regime shift along each branch.
#'     (See [BAMMtools::marginalShiftProbsTree()]). Result is stored in `$MSP_tree` as phylogenetic tree with `$edge.length` scaled to the marginal posterior probability.
#'   * Extract the Maximum A Posteriori probability (MAP) configuration = the configuration of regime shift location found the most frequently among the posterior samples.
#'     (See [BAMMtools::getBestShiftConfiguration()]). This ignores shifts that have an odd-ratio of marginal posterior probability / prior lower than `MAP_odd_ratio_threshold`
#'     to avoid noise from non-core shifts. MAP sample indices are stored in `$MAP_indices`. Diversification rates and shift locations on branches are then averaged across all MAP samples and
#'     recorded as an object of class `"bammdata"` in `$MAP_BAMM_object` with a single `$eventData` table used to plot regime shifts on the phylogeny with [deepSTRAPP::plot_BAMM_rates()].
#'   * Extract the Maximum Shift Credibility (MSC) configuration = the configuration of regime shift location with the highest product of marginal probabilities across branches.
#'     (See [BAMMtools::maximumShiftCredibility()]). MSC sample indices are stored in `$MSC_indices`. Diversification rates and shift locations on branches are then averaged across all MSC samples and
#'     recorded as an object of class `"bammdata"` in `$MSC_BAMM_object` with a single `$eventData` table used to plot regime shifts on the phylogeny with [deepSTRAPP::plot_BAMM_rates()].
#'
#'  Step 5: Clean BAMM files
#'   * Remove files generated in Steps 1 & 2 if `keep_BAMM_outputs = FALSE`.
#'   * Delete the `BAMM_output_directory` if empty after cleaning files.
#'
#'  The `BAMM_object` output:
#'   * is typically used as input to run deepSTRAPP with [deepSTRAPP::run_deepSTRAPP_for_focal_time()] or [deepSTRAPP::run_deepSTRAPP_over_time()].
#'   * can be used to extract rates and regimes for any `focal_time` in the past with [deepSTRAPP::update_rates_and_regimes_for_focal_time()].
#'   * can be used to map diversification rates and regime shifts on the phylogeny with [deepSTRAPP::plot_BAMM_rates()].
#'
#' # Note on diversification models for time-calibrated phylogenies
#'
#' This function relies on BAMM to provide a reliable solution to map diversification rates and regime shifts on a time-calibrated phylogeny
#' and obtain the `BAMM_object` object needed to run the deepSTRAPP workflow ([run_deepSTRAPP_for_focal_time], [run_deepSTRAPP_over_time]).
#' However, it is one option among others for modeling diversification on phylogenies.
#' You may wish to explore alternatives models such as LSBDS model in RevBayes (Höhna et al., 2016), the MTBD model (Barido-Sottani et al., 2020),
#' or the ClaDS2 model (Maliet et al., 2019) for our own data.
#' However, you will need Bayesian models that infer regime shifts to be able to perform STRAPP tests (Rabosky & Huang, 2016).
#' Additionally, you need to format the model output such as in `BAMM_object`, so it can be used in a deepSTRAPP workflow.
#'
#' This function perform a single BAMM run to infer diversification rates and regime shifts.
#' Due to the stochastic nature of the exploration of the parameter space with MCMC process,
#' best practice recommend to ran multiple runs and check for convergence of the MCMC traces,
#' ensuring that the region of high probability has been reached by your MCMC runs.
#'
#' @return The function returns a `BAMM_object` of class `"bammdata"` which is a list with at least 22 elements.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'   * `$node.label` Vector of character strings. Labels of all internal nodes. (Present only if present in the initial `BAMM_object`)
#'
#'   BAMM internal elements used for tree exploration:
#'   * `$begin` Vector of numerical. Absolute time since root of edge/branch start (rootward).
#'   * `$end` Vector of numerical.  Absolute time since root of edge/branch end (tipward).
#'   * `$downseq` Vector of integers. Order of node visits when using a pre-order tree traversal.
#'   * `$lastvisit` ID of the last node visited when starting from the node in the corresponding position in `$downseq`.
#'
#'   BAMM elements summarizing diversification data:
#'   * `$numberEvents` Vector of integer. Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts.
#'   * `$eventData` List of data.frames. One per posterior sample. Records shift events and macroevolutionary regimes parameters. 1st line = Background root regime.
#'   * `$eventVectors` List of integer vectors. One per posterior sample. Record regime ID per branches.
#'   * `$tipStates` List of named integer vectors. One per posterior sample. Record regime ID per tips.
#'   * `$tipLambda` List of named numerical vectors. One per posterior sample. Record speciation rates per tips.
#'   * `$tipMu` List of named numerical vectors. One per posterior sample. Record extinction rates per tips.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of named numerical. Mean tip speciation rates across all posterior configurations of tips.
#'   * `$meanTipMu` Vector of named numerical. Mean tip extinction rates across all posterior configurations of tips.
#'   * `$type` Character string. Set the type of data modeled with BAMM. Should be "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history.
#'
#'  The function also produces files listed in the Details section and stored in the the `BAMM_output_directory`.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::run_deepSTRAPP_for_focal_time()] [deepSTRAPP::run_deepSTRAPP_over_time()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()] [deepSTRAPP::prepare_trait_data()] [deepSTRAPP::plot_BAMM_rates()]
#'
#' @references For BAMM: Rabosky, D. L. (2014). Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees.
#'  PloS one, 9(2), e89543. DOI: \url{https://doi.org/10.1371/journal.pone.0089543}. Website: \url{http://bamm-project.org/}.
#'
#'  For BAMMtools: Rabosky, D. L., Grundler, M., Anderson, C., Title, P., Shi, J. J., Brown, J. W., ... & Larson, J. G. (2014).
#'   BAMM tools: an R package for the analysis of evolutionary dynamics on phylogenetic trees. Methods in Ecology and Evolution, 5(7), 701-707.
#'   DOI: \url{https://doi.org/10.1111/2041-210X.12199}
#'
#' @examples
#' # ----- Example 1: Whale phylogeny ----- #
#'
#' library(phytools)
#' data(whale.tree)
#'
#' \dontrun{
#' # Run BAMM workflow with deepSTRAPP
#' whale_BAMM_object <- prepare_diversification_data(
#'    BAMM_install_directory_path = "./software/bamm-2.5.0/",
#'    phylo = whale.tree,
#'    prefix_for_files = "whale",
#'    numberOfGenerations = 100000 # Set low for the example
#' )}
#'
#' # Load directly the result
#' data(whale_BAMM_object)
#'
#' # Explore output
#' str(whale_BAMM_object, 1)
#'
#' # Plot mean net diversification rates and regime shifts on the phylogeny
#' plot_BAMM_rates(whale_BAMM_object,
#'                 labels = TRUE, legend = TRUE)
#'
#' # ----- Example 2: Ponerinae phylogeny ----- #
#'
#' # Load phylogeny
#' data("Ponerinae_tree", package = "deepSTRAPP")
#' plot(Ponerinae_tree, show.tip.label = FALSE)
#'
#' \dontrun{
#' # Run BAMM workflow with deepSTRAPP
#' Ponerinae_BAMM_object <- prepare_diversification_data(
#'    BAMM_install_directory_path = "./software/bamm-2.5.0/",
#'    phylo = Ponerinae_tree,
#'    prefix_for_files = "Ponerinae",
#'    numberOfGenerations = 10^7 # Set high for optimal run, but will take ages
#' )}
#'
#' # Load directly the result
#' data(Ponerinae_BAMM_object)
#'
#' # Explore output
#' str(Ponerinae_BAMM_object, 1)
#'
#' # Plot mean net diversification rates and regime shifts on the phylogeny
#' plot_BAMM_rates(Ponerinae_BAMM_object,
#'                 labels = FALSE, legend = TRUE)
#'


prepare_diversification_data <- function (BAMM_install_directory_path,
                                          phylo,
                                          prefix_for_files = NULL,
                                          seed = NULL,
                                          numberOfGenerations = 10^7,
                                          globalSamplingFraction = 1.0,
                                          sampleProbsFilename = NULL,
                                          expectedNumberOfShifts = NULL,
                                          eventDataWriteFreq = NULL,
                                          burn_in = 0.25,
                                          nb_posterior_samples = 1000,
                                          additional_BAMM_settings = list(),
                                          BAMM_output_directory_path = "./BAMM_outputs/",
                                          keep_BAMM_outputs = TRUE,
                                          MAP_odd_ratio_threshold = 5,
                                          skip_evaluations = FALSE,
                                          plot_evaluations = TRUE,
                                          save_evaluations = TRUE)

{
  ### Check input validity
  {
    ## BAMM_install_directory_path
    # BAMM_install_directory_path must be a directory, so it must end with '/'
    if (!stringr::str_detect(BAMM_install_directory_path, pattern = "/$"))
    {
      stop(paste0("'BAMM_install_directory_path' must end with '/'"))
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

    ## seed
    if (!is.null(seed))
    {
      if (!is.numeric(seed))
      {
        stop(paste0("'seed' must be an interger."))
      }
    }

    ## numberOfGenerations
    # numberOfGenerations must be a positive integer
    if ((numberOfGenerations != abs(numberOfGenerations)) | (numberOfGenerations != round(numberOfGenerations)))
    {
      stop(paste0("'numberOfGenerations' must be a positive integer defining the number of steps in the MCMC run.\n",
                  "It should be set high enough to reach the equilibrium distribution, and allows posterior samples to be decorrelated.\n",
                  "Check the Effective Sample Size of parameters with coda::effectiveSize() in the Evaluation step.\n",
                  "Default is '10^7'."))
    }

    ## globalSamplingFraction
    # globalSamplingFraction should be a numerical between 0 and 1
    if ((globalSamplingFraction < 0) | (globalSamplingFraction > 1))
    {
      stop(paste0("'globalSamplingFraction' represents the overall proportion of terminals in the phylogeny compared to the estimated overall richness in the clade.\n",
                  "It acts as a multipliers on the rates needed to achieve such extant diversity. It must be between 0 and 1.\n",
                  "Current value of 'globalSamplingFraction' is ",globalSamplingFraction,"."))
    }
    if (globalSamplingFraction == 0)
    {
      stop(paste0("'globalSamplingFraction' represents the overall proportion of terminals in the phylogeny compared to the estimated overall richness in the clade.\n",
                  "It acts as a multipliers on the rates needed to achieve such extant diversity.\n",
                  "It must be between 0 and 1, but cannot be 0.\n"))
    }

    ## sampleProbsFilename
    # If provided, sampleProbsFilename must be a directory, so it must end with '/'
    if (!is.null(sampleProbsFilename))
    {
      if (!stringr::str_detect(sampleProbsFilename, pattern = "\\.txt$"))
      {
        stop(paste0("'sampleProbsFilename' must end with '.txt'"))
      }
    }

    ## expectedNumberOfShifts
    # If provided, expectedNumberOfShifts must be a positive integer
    if (!is.null(expectedNumberOfShifts))
    {
      if ((expectedNumberOfShifts != abs(expectedNumberOfShifts)) | (expectedNumberOfShifts != round(expectedNumberOfShifts)))
      {
        stop(paste0("'expectedNumberOfShifts' must be a positive integer defining the expected number of diversification regime shifts in the phylogeny.\n"),
             "This value is used to set the hyperprior from which the number of shifts is derived.\n",
             "Current value of 'expectedNumberOfShifts' is ",expectedNumberOfShifts,".")
      }
    }

    ## eventDataWriteFreq
    # If provided, eventDataWriteFreq must be a positive integer
    if (!is.null(eventDataWriteFreq))
    {
      if ((eventDataWriteFreq != abs(eventDataWriteFreq)) | (eventDataWriteFreq != round(eventDataWriteFreq)))
      {
        stop(paste0("'eventDataWriteFreq' must be a positive integer defining the frequency at which parameters are sampled in the MCMC run.\n",
                    "It is defined in number of steps between each sampling event.\n",
                    "Current value of 'eventDataWriteFreq' is ",eventDataWriteFreq,"."))
      }
    }

    ## burn_in
    # Burn-in between 0 and 1
    if ((burn_in < 0) | (burn_in > 1))
    {
      stop(paste0("'burn_in' represents the proportion of posterior samples removed to ensure the remaining samples where drawn once the equilibrium distribution was reached.\n",
                  "This can be evaluated looking at the MCMC trace (see Evaluation step). It must be between 0 and 1. Default is '0.25'.\n",
                  "Current value of 'burn_in' is ",burn_in,"."))
    }

    ## nb_posterior_samples
    # nb_posterior_samples must be a positive integer
    if ((nb_posterior_samples != abs(nb_posterior_samples)) | (nb_posterior_samples != round(nb_posterior_samples)))
    {
      stop(paste0("'nb_posterior_samples' must be a positive integer defining the number of posterior samples retainted in the 'BAMM_object' output used for downstream analyses.\n",
                  "Default is '1000'. Current value of 'nb_posterior_samples' is ",nb_posterior_samples,"."))
    }

    ## additional_BAMM_settings
    available_BAMM_settings <- c("runInfoFilename", "sampleFromPriorOnly", "runMCMC", "simulatePriorShifts", "loadEventData", "eventDataInfile",
                                 "initializeModel", "overwrite", "validateEventConfiguration", "lambdaInitPrior", "lambdaShiftPrior", "muInitPrior",
                                 "lambdaIsTimeVariablePrior", "mcmcOutfile", "mcmcWriteFreq", "eventDataOutfile", "printFreq", "outputAcceptanceInfo",
                                 "acceptanceInfoFileName", "acceptanceResetFreq", "updateLambdaInitScale", "updateLambdaShiftScale", "updateMuInitScale",
                                 "updateEventLocationScale", "updateEventRateScale", "updateRateEventNumber", "updateRateEventPosition", "updateRateEventRate",
                                 "updateRateLambda0", "updateRateLambdaShift", "updateRateMu0", "updateRateLambdaTimeMode", "localGlobalMoveRatio",
                                 "lambdaInit0", "lambdaShift0", "muInit0", "initialNumberEvents", "numberOfChains", "deltaT", "swapPeriod", "chainSwapFileName",
                                 "minCladeSizeForShift", "segLength")
    # Check if they match the list of available parameters.
    if (!all(names(additional_BAMM_settings) %in% available_BAMM_settings))
    {
      # Extract settings that do not match with available settings
      error_settings <- names(additional_BAMM_settings)[!(names(additional_BAMM_settings) %in% available_BAMM_settings)]
      initial_warning_length_options <- options()$warning.length
      options(warning.length = 2000L)
      stop(paste0("Names in 'additional_BAMM_settings' do not match with available settings: ",paste(error_settings, collapse = ", "),".\n",
                  "Available setting parameters are: ",paste(available_BAMM_settings, collapse = ", "),".\n",
                  "See details in the template file provided with the package files 'BAMM_template_diversification.txt'."))
      options(warning.length = initial_warning_length_options)
    }

    ## Incompatibility across parameters
    # Check combination of numberOfGenerations, eventDataWriteFreq, burn_in, and nb_posterior_samples to see if enough posterior samples remains after burn-in
    if (is.null(eventDataWriteFreq))
    {
      eventData_freq <- round(as.numeric(numberOfGenerations) / 2000) # Default = write 2000 posterior samples before burn-in
    } else {
      eventData_freq <- eventDataWriteFreq
    }
    nb_posterior_after_burn_in <- (as.numeric(numberOfGenerations)/eventData_freq) * (1-burn_in)
    if (nb_posterior_after_burn_in < nb_posterior_samples)
    {
      stop(paste0("'nb_posterior_samples' requested in not compatible with the number of posterior samples retainted after burn-in.\n",
                  "The number of remaining samples = 'numberOfGenerations'/'eventData_freq' x (1 - 'burn_in'). Here it is ",nb_posterior_after_burn_in,".\n",
                  "'nb_posterior_samples' requested is set to ", nb_posterior_samples,".\n",
                  "Please adjust 'numberOfGenerations', 'eventData_freq', 'burn_in', or 'nb_posterior_samples' to obtain compatible values."))

    }

    ## BAMM_output_directory_path
    # BAMM_output_directory_path must be a directory, so it must end with '/'
    if (!stringr::str_detect(BAMM_output_directory_path, pattern = "/$"))
    {
      stop(paste0("'BAMM_output_directory_path' must end with '/'"))
    }

    ## skip_evaluations & (plot_evaluations | save_evaluations)
    if (skip_evaluations & plot_evaluations)
    {
      cat(paste0("WARNING: 'plot_evaluations' is set to 'TRUE', but 'skip_evaluations' is set to TRUE.\n",
                 "Evaluations will not be plotted as they are skipped.\n\n"))
      warning(paste0("'plot_evaluations' was set to 'TRUE', but 'skip_evaluations' was set to TRUE.\n",
                     "Evaluations were not be plotted as they were skipped."))
    }
    if (skip_evaluations & save_evaluations)
    {
      cat(paste0("WARNING: 'save_evaluations' is set to 'TRUE', but 'skip_evaluations' is set to TRUE.\n",
                 "Evaluations will not be saved as they are skipped.\n\n"))
      warning(paste0("'save_evaluations' is set to 'TRUE', but 'skip_evaluations' is set to TRUE.\n",
                     "Evaluations were not be saved as they were skipped."))
    }

    ## keep_BAMM_outputs & save_evaluations
    if (!keep_BAMM_outputs & !skip_evaluations & save_evaluations)
    {
      warning(paste0("BAMM outputs files were removed ('keep_BAMM_outputs = FALSE'), but evaluation files were preserved in ",BAMM_output_directory_path," as 'save_evaluations' was set to 'TRUE'."))
    }
  }

  #### ----------- Step 1: Set BAMM ----------- ####

  cat(paste0("# ----------- Step 1: Set BAMM ----------- #\n\n"))

  ## Generate the phylo.tree, set the config_file and paths to all input files.

  {
    ## Get names of the additional BAMM settings
    add_settings_names <- names(additional_BAMM_settings)

    ## Create output directory if missing
    if (!dir.exists(paths = file.path(BAMM_output_directory_path)))
    {
      dir.create(path = file.path(BAMM_output_directory_path))
    }

    ## Export the phylogeny in a .tree file

    # Build path
    if (is.null(prefix_for_files))
    {
      phy_path <- file.path(paste0(BAMM_output_directory_path, "phylogeny.tree"))
    } else {
      phy_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files,"_phylogeny.tree"))
    }
    # Export tree
    ape::write.tree(phy = phylo, file = phy_path)

    ## Load control file template for diversification analyses

    # # Load it from the root as in binary/installed version of the package
    # BAMM_config_file <- tryCatch({
    #   readLines(con = file.path("./BAMM_template_diversification.txt"))
    # }, warning = function(w) { }, # Do nothing
    # error = function(e) { }, # Do nothing
    # finally = { } # Do nothing
    # )
    # # If failed, load it from the /inst/ directory as in source/bundled version of the package
    # if (is.null(BAMM_config_file))
    # {
    #   BAMM_config_file <- readLines(con = file.path("./inst/BAMM_template_diversification.txt"))
    # }

    # Load control file template from internal deepSTRAPP data
    # utils::data("BAMM_template_diversification", package = "deepSTRAPP")
    BAMM_config_file <- deepSTRAPP::BAMM_template_diversification

    # Initiate new config file for this analysis
    my_config_file <- BAMM_config_file

    ### 1.1/ Set general settings and data input ####

    # Set prefix to add to all output files (separated with "_")
    if (is.null(prefix_for_files))
    {
      # Do nothing = Add no prefix (line is commented out in the template)
    } else {
      outName_line <- which(stringr::str_detect(string = my_config_file, pattern = "outName = "))[1]
      my_config_file[outName_line] <- paste0("outName = ", prefix_for_files)
    }

    # Set path to the phylogenetic tree file
    phy_path_line <- which(stringr::str_detect(string = my_config_file, pattern = "treefile ="))
    my_config_file[phy_path_line] <- paste0("treefile = ", phy_path)

    # Set path to the file storing all information on this run
    if ("runInfoFilename" %in% add_settings_names)
    {
      runInfoFilename <- file.path(additional_BAMM_settings$runInfoFilename)
    } else {
      runInfoFilename <- "run_info.txt"
    }
    runInfoFilename_line <- which(stringr::str_detect(string = my_config_file, pattern = "runInfoFilename = "))
    my_config_file[runInfoFilename_line] <- paste0("runInfoFilename = ", runInfoFilename)

    # Should the run sample from prior only?
    # To check that the prior settings are fine
    # To compare with the posterior, to check the influence of the prior on the posterior
    if ("sampleFromPriorOnly" %in% add_settings_names)
    {
      sampleFromPriorOnly <- additional_BAMM_settings$sampleFromPriorOnly
    } else {
      sampleFromPriorOnly <- 0 # 0 = No. 1 = Yes
    }
    sampleFromPriorOnly_line <- which(stringr::str_detect(string = my_config_file, pattern = "sampleFromPriorOnly = "))
    my_config_file[sampleFromPriorOnly_line] <- paste0("sampleFromPriorOnly = ", sampleFromPriorOnly)

    # Should the MCMC simulation be performed?
    # If runMCMC = 0, the program will only check whether the data file can be read and the initial likelihood computed
    if ("runMCMC" %in% add_settings_names)
    {
      runMCMC <- additional_BAMM_settings$runMCMC
    } else {
      runMCMC <- 1  # 0 = No. 1 = Yes
    }
    runMCMC_line <- which(stringr::str_detect(string = my_config_file, pattern = "runMCMC = "))[1]
    my_config_file[runMCMC_line] <- paste0("runMCMC = ", runMCMC)

    # Should the prior distribution of the number of shift events, given the hyperprior on the Poisson rate parameter, be simulated?
    # This was necessary to compute Bayes factors
    # But is now disabled as the exact (analytical) prior is now implemented in BAMMtools.
    if ("simulatePriorShifts" %in% add_settings_names)
    {
      simulatePriorShifts <- additional_BAMM_settings$simulatePriorShifts
    } else {
      simulatePriorShifts <- 0  # 0 = No. 1 = Yes
    }
    simulatePriorShifts_line <- which(stringr::str_detect(string = my_config_file, pattern = "simulatePriorShifts = "))[1]
    my_config_file[simulatePriorShifts_line] <- paste0("simulatePriorShifts = ", simulatePriorShifts)

    # Whether to load a previous event data file
    if ("loadEventData" %in% add_settings_names)
    {
      loadEventData <- additional_BAMM_settings$loadEventData
    } else {
      loadEventData <- 0  # 0 = No. 1 = Yes
    }
    loadEventData_line <- which(stringr::str_detect(string = my_config_file, pattern = "loadEventData = "))[1]
    my_config_file[loadEventData_line] <- paste0("loadEventData = ", loadEventData)

    # Provides file name of the event data file to load, used only if loadEventData = 1
    if ("eventDataInfile" %in% add_settings_names)
    {
      eventDataInfile <- file.path(additional_BAMM_settings$eventDataInfile)
    } else {
      eventDataInfile <- "event_data_in.txt"
    }
    eventDataInfile_line <- which(stringr::str_detect(string = my_config_file, pattern = "eventDataInfile = "))[1]
    my_config_file[eventDataInfile_line] <- paste0("eventDataInfile = ", eventDataInfile)

    # Whether to initialize (but not run) the MCMC.
    # If initializeModel = 0, the program will only ensure that the data files (e.g., treefile) can be read
    if ("initializeModel" %in% add_settings_names)
    {
      initializeModel <- additional_BAMM_settings$initializeModel
    } else {
      initializeModel <- "event_data_in.txt"
    }
    initializeModel <- 1  # 0 = No. 1 = Yes
    initializeModel_line <- which(stringr::str_detect(string = my_config_file, pattern = "initializeModel = "))[1]
    my_config_file[initializeModel_line] <- paste0("initializeModel = ", initializeModel)

    # Whether to use a "global" sampling probability to assign the proportion of terminal represented in the tree.
    # If False (0), expects a file path for clade-specific sampling probabilities (see sampleProbsFilename)
    if (!is.null(sampleProbsFilename))
    {
      useGlobalSamplingProbability = 0  # 0 = No. 1 = Yes
    } else {
      useGlobalSamplingProbability = 1  # 0 = No. 1 = Yes
    }
    useGlobalSamplingProbability_line <- which(stringr::str_detect(string = my_config_file, pattern = "useGlobalSamplingProbability = "))[1]
    my_config_file[useGlobalSamplingProbability_line] <- paste0("useGlobalSamplingProbability = ", useGlobalSamplingProbability)

    # Provides the global sampling fraction
    # If useGlobalSamplingProbability = 0, this is ignored and BAMM looks for a file path to clade-specific sampling fractions
    globalSamplingFraction_line <- which(stringr::str_detect(string = my_config_file, pattern = "globalSamplingFraction = "))[1]
    my_config_file[globalSamplingFraction_line] <- paste0("globalSamplingFraction = ", globalSamplingFraction)

    # Provides path to the file containing clade-specific sampling fractions
    if (is.null(sampleProbsFilename))
    {
      sampleProbsFilename <- "taxa_sampling_probs.txt"
    } else {
      sampleProbsFilename <- file.path(sampleProbsFilename)
    }
    sampleProbsFilename_line <- which(stringr::str_detect(string = my_config_file, pattern = "sampleProbsFilename = "))[1]
    my_config_file[sampleProbsFilename_line] <- paste0("sampleProbsFilename = ", sampleProbsFilename)

    # Set the seed for the random number generator.
    # If not specified (or is -1), a seed is obtained from the system clock
    if (is.null(seed))
    {
      seed <- -1
    }
    seed_line <- which(stringr::str_detect(string = my_config_file, pattern = "seed = "))[1]
    my_config_file[seed_line] <- paste0("seed = ", seed)

    # Should the output files be overwritten?
    # If True (1), the program will overwrite any output files in the current directory (if present)
    if ("overwrite" %in% add_settings_names)
    {
      overwrite <- additional_BAMM_settings$overwrite
    } else {
      overwrite <- 1  # 0 = No. 1 = Yes
    }
    overwrite_line <- which(stringr::str_detect(string = my_config_file, pattern = "overwrite = "))[1]
    my_config_file[overwrite_line] <- paste0("overwrite = ", overwrite)

    # Set limits to valid configurations
    # If 1, rejects proposals that cause a branch and both of its direct descendants to have at least one event.
      # Such an event configuration may cause the parameters of the parent event to change to unrealistic values.
    # If 0, no such proposals are immediately rejected. The default value is 1.
    if ("validateEventConfiguration" %in% add_settings_names)
    {
      validateEventConfiguration <- additional_BAMM_settings$validateEventConfiguration
    } else {
      validateEventConfiguration <- 1  # 0 = Do not reject. 1 = Reject
    }
    validateEventConfiguration_line <- which(stringr::str_detect(string = my_config_file, pattern = "validateEventConfiguration = "))[1]
    my_config_file[validateEventConfiguration_line] <- paste0("validateEventConfiguration = ", validateEventConfiguration)

    ### 1.2/ Set (hyper)prior settings ####

    ## Can use this help function to automatically tune prior adapted to your data by scaling the prior distributions based on the age (root depth) of your tree
    # In practice, setBAMMpriors first estimates the rate of speciation for your full tree under a pure birth model of diversification.
    # Then assume, arbitrarily, that a reasonable prior distribution for the initial lambda0/mu0 rate parameters is an exponential distribution with a mean five times greater than this pure birth value.
    # Rationale = having a weakly informative prior that is still in the order of magnitude of the true rate
    # For the shift parameter (alpha), the sd of the normal prior is set such as mean +/- 2s gives an alpha parameter that results in
    # either a 90% decline in the evolutionary rate or a 190% increase in rate on the interval of time from the root to the tips of the tree.

    # Build path to priors file
    if (is.null(prefix_for_files))
    {
      priors_path <- file.path(paste0(BAMM_output_directory_path, "priors.txt"))
    } else {
      priors_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files,"_priors.txt"))
    }
    # Generate priors file
    invisible(capture.output(BAMMtools::setBAMMpriors(phy = phylo, outfile = priors_path, suppressWarning = TRUE)))
    # Read prior file
    default_tuned_priors <- readLines(con = priors_path)

    # Set the expected number of shifts used to set the exponential hyperprior for nb of rate shifts (from which the Λ is drawn)
    # Suggested values (if set to NULL): 1 regime shift for 100 tips.
    # Best practice is to try runs with several value and inspect the convergence of the posterior distribution of the regime shift parameter with the prior distribution defined by this hyperprior parameter.
    if (is.null(expectedNumberOfShifts))
    {
      expectedNumberOfShifts <- max(1, round(length(phylo$tip.label)/100)) # 1 expected regime shift for 100 tips
    }
    expectedNumberOfShifts_line <- which(stringr::str_detect(string = my_config_file, pattern = "expectedNumberOfShifts = "))[1]
    my_config_file[expectedNumberOfShifts_line] <- paste0("expectedNumberOfShifts = ", expectedNumberOfShifts)

    # Set the rate parameter of the exponential prior(s) of initial lambda parameters (lambda0) of speciation rate regimes
      # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
    if ("lambdaInitPrior" %in% add_settings_names)
    {
      # If provided as 'additional_BAMM_settings', use this value
      lambdaInitPrior <- additional_BAMM_settings$lambdaInitPrior
    } else {
      # If not provided, use the default value suggested by BAMMtools::setBAMMpriors
      lambdaInitPrior_default_line <- which(stringr::str_detect(string = default_tuned_priors, pattern = "lambdaInitPrior = "))[1]
      lambdaInitPrior <- as.numeric(stringr::str_remove(string = default_tuned_priors[lambdaInitPrior_default_line], pattern = "lambdaInitPrior = "))

    }
    lambdaInitPrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaInitPrior = "))[1]
    my_config_file[lambdaInitPrior_line] <- paste0("lambdaInitPrior = ", lambdaInitPrior)

    # Set the standard deviation of the normal distribution prior(s) of rate variation parameters (alpha) of speciation rate regimes
      # alpha in lambda(t) = lamba0 x exp(alpha*t)
    # Mean of this prior(s) are fixed to zero such as a constant rate diversification process is the most probable a priori
    if ("lambdaShiftPrior" %in% add_settings_names)
    {
      # If provided as 'additional_BAMM_settings', use this value
      lambdaShiftPrior <- additional_BAMM_settings$lambdaShiftPrior
    } else {
      # If not provided, use the default value suggested by BAMMtools::setBAMMpriors
      lambdaShiftPrior_default_line <- which(stringr::str_detect(string = default_tuned_priors, pattern = "lambdaShiftPrior = "))[1]
      lambdaShiftPrior <- as.numeric(stringr::str_remove(string = default_tuned_priors[lambdaShiftPrior_default_line], pattern = "lambdaShiftPrior = "))

    }
    lambdaShiftPrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaShiftPrior = "))[1]
    my_config_file[lambdaShiftPrior_line] <- paste0("lambdaShiftPrior = ", lambdaShiftPrior)

    # Set the rate parameter of the exponential prior(s) of initial lambda parameters (mu0) of extinction rate regimes
      # mu0 in mu(t) = mu0 x exp(alpha*t)
    # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
    if ("muInitPrior" %in% add_settings_names)
    {
      # If provided as 'additional_BAMM_settings', use this value
      muInitPrior <- additional_BAMM_settings$muInitPrior
    } else {
      # If not provided, use the default value suggested by BAMMtools::setBAMMpriors
      muInitPrior_default_line <- which(stringr::str_detect(string = default_tuned_priors, pattern = "muInitPrior = "))[1]
      muInitPrior <- as.numeric(stringr::str_remove(string = default_tuned_priors[muInitPrior_default_line], pattern = "muInitPrior = "))

    }
    muInitPrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "muInitPrior = "))[1]
    my_config_file[muInitPrior_line] <- paste0("muInitPrior = ", muInitPrior)

    # No prior for rate variation parameters (alpha) of extinction rate regimes as they are assumed to follow constant rates

    # Set the prior (probability) of the time mode (of speciation?) being time-variable (vs. time-constant)
    # By default, allows all regimes to be time-variable, as their rate can still be estimated as constant with alpha = 0
    if ("lambdaIsTimeVariablePrior" %in% add_settings_names)
    {
      lambdaIsTimeVariablePrior <- additional_BAMM_settings$lambdaIsTimeVariablePrior
    } else {
      lambdaIsTimeVariablePrior <- 1
    }
    lambdaIsTimeVariablePrior <- 1
    lambdaIsTimeVariablePrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaIsTimeVariablePrior = "))[1]
    my_config_file[lambdaIsTimeVariablePrior_line] <- paste0("lambdaIsTimeVariablePrior = ", lambdaIsTimeVariablePrior)

    ### 1.3/ Set the MCMC simulation settings, MCMC logs and output options ####

    # Set the number of generations to perform MCMC simulation
    numberOfGenerations <- format(numberOfGenerations, scientific = F) # 10^7
    numberOfGenerations_line <- which(stringr::str_detect(string = my_config_file, pattern = "numberOfGenerations = "))[1]
    my_config_file[numberOfGenerations_line] <- paste0("numberOfGenerations = ", numberOfGenerations)

    # Set the path to the MCMC output file
    # Includes only summary information about MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)
    if ("mcmcOutfile" %in% add_settings_names)
    {
      mcmcOutfile <- file.path(additional_BAMM_settings$mcmcOutfile)
    } else {
      mcmcOutfile <- "mcmc_log.txt"
    }
    mcmcOutfile_line <- which(stringr::str_detect(string = my_config_file, pattern = "mcmcOutfile = "))[1]
    my_config_file[mcmcOutfile_line] <- paste0("mcmcOutfile = ", mcmcOutfile)

    # Set the frequency in which to write the MCMC output to the log file
    # Aim for 500-5000 posterior samples ideally
    if ("mcmcWriteFreq" %in% add_settings_names)
    {
      mcmcWriteFreq <- additional_BAMM_settings$mcmcWriteFreq
    } else {
      mcmcWriteFreq <- round(as.numeric(numberOfGenerations) / 2000) # Default = record 2000 MCMC generations
    }
    mcmcWriteFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "mcmcWriteFreq = "))[1]
    my_config_file[mcmcWriteFreq_line] <- paste0("mcmcWriteFreq = ", mcmcWriteFreq)

    # Set the path to the main output file
    # The raw event data. ALL of the results are contained in this file,
    # All branch-specific speciation rates, shift positions, marginal distributions, etc. can be reconstructed from this output.
    # See ?BAMMtools::getEventData to import this output in R
    if ("eventDataOutfile" %in% add_settings_names)
    {
      eventDataOutfile <- file.path(additional_BAMM_settings$eventDataOutfile)
    } else {
      eventDataOutfile <- "event_data.txt"
    }
    eventDataOutfile_line <- which(stringr::str_detect(string = my_config_file, pattern = "eventDataOutfile = "))[1]
    my_config_file[eventDataOutfile_line] <- paste0("eventDataOutfile = ", eventDataOutfile)

    # Set frequency in which to write the event data to the output file = the sampling frequency of posterior samples
    # Aim for 500-5000 posterior samples ideally
    # Will need to remove some to account for the burn-in
    if (is.null(eventDataWriteFreq))
    {
      eventDataWriteFreq <- round(as.numeric(numberOfGenerations) / 2000) # Default = write 2000 posterior samples before burn-in
    }
    eventDataWriteFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "eventDataWriteFreq = "))[1]
    my_config_file[eventDataWriteFreq_line] <- paste0("eventDataWriteFreq = ", eventDataWriteFreq)

    # Set frequency in which to print MCMC status to the screen
    if ("printFreq" %in% add_settings_names)
    {
      printFreq <- additional_BAMM_settings$printFreq
    } else {
      if (numberOfGenerations >= 10^6)
      {
        printFreq <- 10000 # Print status every 10^4 generations for long runs
      } else {
        printFreq <- 1000 # Print status every 10^3 generations for short runs
      }
    }
    printFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "printFreq = "))[1]
    my_config_file[printFreq_line] <- paste0("printFreq = ", printFreq)

    # Whether acceptance/proposal history should be saved.
    # If 1, outputs whether each proposal was accepted. The number identifying the proposal matches the one in the code.
    # The default value is 0 (i.e., do not output this information).
    if ("outputAcceptanceInfo" %in% add_settings_names)
    {
      outputAcceptanceInfo <- additional_BAMM_settings$outputAcceptanceInfo
    } else {
      outputAcceptanceInfo <- 0 # Do not save acceptance/proposal history
    }
    outputAcceptanceInfo_line <- which(stringr::str_detect(string = my_config_file, pattern = "outputAcceptanceInfo = "))[1]
    my_config_file[outputAcceptanceInfo_line] <- paste0("outputAcceptanceInfo = ", outputAcceptanceInfo)

    # Set the path to the acceptance/proposal history file
    # The path of the file to which to write whether each proposal was accepted.
    # outputAcceptanceInfo must be set to 1 for this information to be written.
    if ("acceptanceInfoFileName" %in% add_settings_names)
    {
      acceptanceInfoFileName <- file.path(additional_BAMM_settings$acceptanceInfoFileName)
    } else {
      acceptanceInfoFileName <- "acceptance_info.txt"
    }
    acceptanceInfoFileName_line <- which(stringr::str_detect(string = my_config_file, pattern = "acceptanceInfoFileName = "))[1]
    my_config_file[acceptanceInfoFileName_line] <- paste0("acceptanceInfoFileName = ", acceptanceInfoFileName)

    # Set frequency in which to update the acceptance rate calculation
    # Acceptance rate = how often new proposal are accepted as the next step in the chain
    # Give information on how efficient is the movement of the chain in the parameter space
    # The acceptance rate is output to both the MCMC data file and print to the screen
    if ("acceptanceResetFreq" %in% add_settings_names)
    {
      acceptanceResetFreq <- additional_BAMM_settings$acceptanceResetFreq
    } else {
      if (numberOfGenerations >= 10^6)
      {
        acceptanceResetFreq <- 10000/2 # Update acceptance rate every 5000 generations for long runs
      } else {
        acceptanceResetFreq <- 1000/2 # Print status every 500 generations for short runs
      }
    }
    acceptanceResetFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "acceptanceResetFreq = "))[1]
    my_config_file[acceptanceResetFreq_line] <- paste0("acceptanceResetFreq = ", acceptanceResetFreq)

    ### 1.4/ Set the scaling operators = temperatures, to propose new values for sampled parameters ####

    # The highest scaling operators = temperatures = the bigger changes can be implemented
    # Advantages = allows to escape suboptimum
    # Cons = may be unstable / harder to reach convergence

    # Set scale parameter used for updating the initial speciation rate (lambda0) for each regime/process
      # Updated as lambda0_new ~ lambda0_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
    if ("updateLambdaInitScale" %in% add_settings_names)
    {
      updateLambdaInitScale <- additional_BAMM_settings$updateLambdaInitScale
    } else {
      updateLambdaInitScale <- 2.0
    }
    updateLambdaInitScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateLambdaInitScale = "))[1]
    my_config_file[updateLambdaInitScale_line] <- paste0("updateLambdaInitScale = ", updateLambdaInitScale)

    # Set window size parameter used for updating the rate variation parameter (alpha) of speciation rates for each regime/process
      # Updated as alpha_new ~ alpha_old + U with U a uniform distribution ranging between - window_size_par and + window_size_par
    if ("updateLambdaShiftScale" %in% add_settings_names)
    {
      updateLambdaShiftScale <- additional_BAMM_settings$updateLambdaShiftScale
    } else {
      updateLambdaShiftScale <- 0.1
    }
    updateLambdaShiftScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateLambdaShiftScale = "))[1]
    my_config_file[updateLambdaShiftScale_line] <- paste0("updateLambdaShiftScale = ", updateLambdaShiftScale)

    # Set scale parameter used for updating the initial extinction rate (mu0) for each regime/process
      # Updated as mu0_new ~ mu0_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
    if ("updateMuInitScale" %in% add_settings_names)
    {
      updateMuInitScale <- additional_BAMM_settings$updateMuInitScale
    } else {
      updateMuInitScale <- 2.0
    }
    updateMuInitScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateMuInitScale = "))[1]
    my_config_file[updateMuInitScale_line] <- paste0("updateMuInitScale = ", updateMuInitScale)

    # Set window size parameter used for updating LOCAL moves of the position of shifts on the tree
      # Updated as position_new ~ position_old + U with U a uniform distribution ranging between - window_size_par and + window_size_par
    # Unit = fraction of root_to_tip length. May lead to jump of the shift position to a new branch
    # Ex: For a tree of 100My, with parameter set to 0.05, the proposal window for local position change is +/- 5My around the previous value
    if ("updateEventLocationScale" %in% add_settings_names)
    {
      updateEventLocationScale <- additional_BAMM_settings$updateEventLocationScale
    } else {
      updateEventLocationScale <- 0.05
    }
    updateEventLocationScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateEventLocationScale = "))[1]
    my_config_file[updateEventLocationScale_line] <- paste0("updateEventLocationScale = ", updateEventLocationScale)

    # Set scale parameter used for updating the LAMBDA rate parameter of the Poisson process controlling the number of shifts in the submodels
      # Updated as LAMBDA_new ~ LAMBDA_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
    if ("updateEventRateScale" %in% add_settings_names)
    {
      updateEventRateScale <- additional_BAMM_settings$updateEventRateScale
    } else {
      updateEventRateScale <- 4.0
    }
    updateEventRateScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateEventRateScale = "))[1]
    my_config_file[updateEventRateScale_line] <- paste0("updateEventRateScale = ", updateEventRateScale)

    ### 1.5/ Set the relative frequencies of operator uses (frequency of parameter updates) at each generation ####

    # Set the relative frequency of MCMC moves that change the number of events in the submodel (shift from Mk to Mk+1 or Mk-1)
    if ("updateRateEventNumber" %in% add_settings_names)
    {
      updateRateEventNumber <- additional_BAMM_settings$updateRateEventNumber
    } else {
      updateRateEventNumber <- 0.1 # 1/10 generations
    }
    updateRateEventNumber_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateEventNumber = "))[1]
    my_config_file[updateRateEventNumber_line] <- paste0("updateRateEventNumber = ", updateRateEventNumber)

    # Set the relative frequency of MCMC moves that change the location of an event on the tree (update position parameters)
    if ("updateRateEventPosition" %in% add_settings_names)
    {
      updateRateEventPosition <- additional_BAMM_settings$updateRateEventPosition
    } else {
      updateRateEventPosition <- 1 # Every generation !
    }
    updateRateEventPosition_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateEventPosition = "))[1]
    my_config_file[updateRateEventPosition_line] <- paste0("updateRateEventPosition = ", updateRateEventPosition)

    # Set the relative frequency of MCMC moves that change the rate at which events occur (update LAMBDA parameter)
    if ("updateRateEventRate" %in% add_settings_names)
    {
      updateRateEventRate <- additional_BAMM_settings$updateRateEventRate
    } else {
      updateRateEventRate <- 1 # Every generation !
    }
    updateRateEventRate_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateEventRate = "))[1]
    my_config_file[updateRateEventRate_line] <- paste0("updateRateEventRate = ", updateRateEventRate)

    # Set the relative frequency of MCMC moves that change the initial speciation rates (lambda0) associated with a regime
      # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
    if ("updateRateLambda0" %in% add_settings_names)
    {
      updateRateLambda0 <- additional_BAMM_settings$updateRateLambda0
    } else {
      updateRateLambda0 <- 1 # Every generation !
    }
    updateRateLambda0_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateLambda0 = "))[1]
    my_config_file[updateRateLambda0_line] <- paste0("updateRateLambda0 = ", updateRateLambda0)

    # Set the relative frequency of MCMC moves that change the exponential shift parameter (alpha) of the speciation rate associated with a regime
      # alpha in lambda(t) = lamba0 x exp(alpha*t)
    if ("updateRateLambdaShift" %in% add_settings_names)
    {
      updateRateLambdaShift <- additional_BAMM_settings$updateRateLambdaShift
    } else {
      updateRateLambdaShift <- 1 # Every generation !
    }
    updateRateLambdaShift_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateLambdaShift = "))[1]
    my_config_file[updateRateLambdaShift_line] <- paste0("updateRateLambdaShift = ", updateRateLambdaShift)

    # Set the relative frequency of MCMC moves that change the (initial) extinction rate associated with a regime
      # mu0 in mu(t) = mu0 x exp(alpha*t)
    # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
    if ("updateRateMu0" %in% add_settings_names)
    {
      updateRateMu0 <- additional_BAMM_settings$updateRateMu0
    } else {
      updateRateMu0 <- 1 # Every generation !
    }
    updateRateMu0_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateMu0 = "))[1]
    my_config_file[updateRateMu0_line] <- paste0("updateRateMu0 = ", updateRateMu0)

    # Set the relative frequency of MCMC moves that flip the time mode (time-constant <=> time-variable)
    # By default, only use time-variable mode, so the frequency is set to 0.
    if ("updateRateLambdaTimeMode" %in% add_settings_names)
    {
      updateRateLambdaTimeMode <- additional_BAMM_settings$updateRateLambdaTimeMode
    } else {
      updateRateLambdaTimeMode <- 0 # Never shift to time-constant mode
    }
    updateRateLambdaTimeMode_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateLambdaTimeMode = "))[1]
    my_config_file[updateRateLambdaTimeMode_line] <- paste0("updateRateLambdaTimeMode = ", updateRateLambdaTimeMode)

    # Set the ratio of local to global moves used to propose new location of events on the tree (update position parameters)
    if ("localGlobalMoveRatio" %in% add_settings_names)
    {
      localGlobalMoveRatio <- additional_BAMM_settings$localGlobalMoveRatio
    } else {
      localGlobalMoveRatio <- 10.0 # Ten times more local changes than global changes
    }
    localGlobalMoveRatio_line <- which(stringr::str_detect(string = my_config_file, pattern = "localGlobalMoveRatio = "))[1]
    my_config_file[localGlobalMoveRatio_line] <- paste0("localGlobalMoveRatio = ", localGlobalMoveRatio)

    ### 1.6/ Set the initial parameter values to start the MCMC chain(s) ####

    # The MCMC chain start with a model with no shift (M0 submodel). So the initial parameter values are for this unique regime

    # Run a BD model to obtain credible starting values for the root process
    BD_fit <- suppressWarnings(phytools::fit.bd(tree = phylo))

    # Set the initial speciation rate (lambda0) for the first regime starting at the root of the tree (regime 0)
    # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
    if ("lambdaInit0" %in% add_settings_names)
    {
      lambdaInit0 <- additional_BAMM_settings$lambdaInit0
    } else {
      lambdaInit0 <- BD_fit$b # Provide estimates from BD model
    }
    lambdaInit0_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaInit0 = "))[1]
    my_config_file[lambdaInit0_line] <- paste0("lambdaInit0 = ", lambdaInit0)

    # Set the initial shift parameter (alpha) for the root process (regime 0)
      # alpha in lambda(t) = lamba0 x exp(alpha*t)
    if ("lambdaShift0" %in% add_settings_names)
    {
      lambdaShift0 <- additional_BAMM_settings$lambdaShift0
    } else {
      lambdaShift0 <- 0 # Initial value set to 0 such as the process is a constant rate
    }
    lambdaShift0_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaShift0 = "))[1]
    my_config_file[lambdaShift0_line] <- paste0("lambdaShift0 = ", lambdaShift0)

    # Set the intial extinction rate (mu0) for the first regime starting at the root of the tree (regime 0)
      # mu0 in mu(t) = mu0 x exp(alpha*t)
    # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
    if ("muInit0" %in% add_settings_names)
    {
      muInit0 <- additional_BAMM_settings$muInit0
    } else {
      muInit0 <- BD_fit$d # Provide estimates from BD model
    }
    muInit0_line <- which(stringr::str_detect(string = my_config_file, pattern = "muInit0 = "))[1]
    my_config_file[muInit0_line] <- paste0("muInit0 = ", muInit0)

    # Set the initial number of non-root processes/shifts = M0 submodel
    if ("initialNumberEvents" %in% add_settings_names)
    {
      initialNumberEvents <- additional_BAMM_settings$initialNumberEvents
    } else {
      initialNumberEvents <- 0 # Start with a single regime = M0 submodel
    }
    initialNumberEvents_line <- which(stringr::str_detect(string = my_config_file, pattern = "initialNumberEvents = "))[1]
    my_config_file[initialNumberEvents_line] <- paste0("initialNumberEvents = ", muInit0)

    ### 1.7/ Set the MCMC chain behavior ####

    # Set the number of Markov chains to run
    # Each chain will have a different temperature to favor different exploration behavior of the parameter space
    if ("numberOfChains" %in% add_settings_names)
    {
      numberOfChains <- additional_BAMM_settings$numberOfChains
    } else {
      numberOfChains <- 4 # Default = using 4 coupled-chains to explore pararameter space for a single run.
    }
    numberOfChains_line <- which(stringr::str_detect(string = my_config_file, pattern = "numberOfChains = "))[1]
    my_config_file[numberOfChains_line] <- paste0("numberOfChains = ", numberOfChains)

    # Set the temperature increment parameter that control the difference of temperatures between the chains.
    # This value should be > 0
    # The temperature for the i-th chain is computed as 1 / [1 + deltaT * (i - 1)]
    # Chain 1 is the coldest. Highest chain is the hottest
    if ("deltaT" %in% add_settings_names)
    {
      deltaT <- additional_BAMM_settings$deltaT
    } else {
      deltaT <- 0.01
    }
    deltaT_line <- which(stringr::str_detect(string = my_config_file, pattern = "deltaT = "))[1]
    my_config_file[deltaT_line] <- paste0("deltaT = ", deltaT)

    # Set the frequency of generations at which to propose a chain swap
    # The coupled-MCMC algorithm will check chain state and eventually swap for the one having reach the highest likelihood
    if ("swapPeriod" %in% add_settings_names)
    {
      swapPeriod <- additional_BAMM_settings$swapPeriod
    } else {
      swapPeriod <- 1000 # Check swapping every 10^3 generations
    }
    swapPeriod_line <- which(stringr::str_detect(string = my_config_file, pattern = "swapPeriod = "))[1]
    my_config_file[swapPeriod_line] <- paste0("swapPeriod = ", swapPeriod)

    # Set the path to the file where to store information about each chain swap proposal
    # The format of each line is [generation],[rank_1],[rank_2],[swap_accepted]
    # where [generation] is the generation in which the swap proposal was made,
    # [rank_1] and [rank_2] are the chains that were chosen, and [swap_accepted] is
    # whether the swap was made. The cold chain has a rank of 1.
    if ("chainSwapFileName" %in% add_settings_names)
    {
      chainSwapFileName <- file.path(additional_BAMM_settings$chainSwapFileName)
    } else {
      chainSwapFileName <- "chain_swap_log.txt"
    }
    chainSwapFileName_line <- which(stringr::str_detect(string = my_config_file, pattern = "chainSwapFileName = "))[1]
    my_config_file[chainSwapFileName_line] <- paste0("chainSwapFileName = ", chainSwapFileName)

    ### 1.8/ Set other parameters ####

    # Set the minimum size of a clade to allow a shift to occur
    # Constrain location of possible rate-change events to occur only on branches with at least this many descendant tips.
    # The default value of 1 allows shifts to occur on all branches.
    if ("minCladeSizeForShift" %in% add_settings_names)
    {
      minCladeSizeForShift <- additional_BAMM_settings$minCladeSizeForShift
    } else {
      minCladeSizeForShift <- 1 # Shift can occur on all branches, even the terminal ones.
    }
    minCladeSizeForShift_line <- which(stringr::str_detect(string = my_config_file, pattern = "minCladeSizeForShift = "))[1]
    my_config_file[minCladeSizeForShift_line] <- paste0("minCladeSizeForShift = ", minCladeSizeForShift)

    # Set the "grain" at which time-continuous calculations are discretized
    # The continuous-time change in diversification rates are approximated by breaking each branch into constant-rate diversification segments
    # with each segment given a length determined by the segLength parameter.
    # segLength is in fraction of the root-to-tip distance of the tree.
    # Ex: For an ultrametric tree of 100My, a segLength of 0.02 lead to a step size of 2My
    # If the value is greater than a given branch length BAMM will not break the branch into segments but use the mean rate across the entire branch.
    if ("segLength" %in% add_settings_names)
    {
      segLength <- additional_BAMM_settings$segLength
    } else {
      segLength = 0.01 # Use segments of length = 1% of tree depth to discretize calculations of time-continuous rates
    }
    segLength_line <- which(stringr::str_detect(string = my_config_file, pattern = "segLength = "))[1]
    my_config_file[segLength_line] <- paste0("segLength = ", segLength)

    ### 1.9/ Export the updated custom config file ####

    ## Export the config file with the prefix
    # Build path
    if (is.null(prefix_for_files))
    {
      config_file_path <- file.path(paste0(BAMM_output_directory_path, "config_file.txt"))
    } else {
      config_file_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files,"_config_file.txt"))
    }
    # Export my_config_file
    writeLines(text = my_config_file, con = config_file_path)

  }

  #### ----------- Step 2: Run BAMM ----------- ####

  cat(paste0("# ----------- Step 2: Run BAMM ----------- #\n\n"))

  ## Run BAMM and move output files in dedicated directory

  {
    ## run BAMM
    BAMM_path <- file.path(paste0(BAMM_install_directory_path, "bamm"))
    system(paste0(BAMM_path, " -c ", config_file_path))

    ### Outputs
    # run_info.txt file containing a summary of your parameters/settings
    # mcmc_log.txt file containing raw MCMC information useful in diagnosing convergence
    # event_data.txt file containing all evolutionary rate parameters and their topological mappings
    # chain_swap_log.txt file containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).
    # acceptance_info.txt containing the history of acceptance/proposal of MCMC steps (If additional parameter 'outputAcceptanceInfo' is set to 1)

    ## Clean outputs = move files to the dedicated directory
    # Detect output files
    if (is.null(prefix_for_files))
    {
      file.rename(from = runInfoFilename, to = file.path(paste0(BAMM_output_directory_path, runInfoFilename)))
      file.rename(from = mcmcOutfile, to = file.path(paste0(BAMM_output_directory_path, mcmcOutfile)))
      file.rename(from = eventDataOutfile, to = file.path(paste0(BAMM_output_directory_path, eventDataOutfile)))
      file.rename(from = chainSwapFileName, to = file.path(paste0(BAMM_output_directory_path, chainSwapFileName)))
      if (outputAcceptanceInfo == 1)
      {
        file.rename(from = acceptanceInfoFileName, to = file.path(paste0(BAMM_output_directory_path, acceptanceInfoFileName)))
      }
    } else {
      file.rename(from = file.path(paste0(prefix_for_files, "_", runInfoFilename)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", runInfoFilename)))
      file.rename(from = file.path(paste0(prefix_for_files, "_", mcmcOutfile)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", mcmcOutfile)))
      file.rename(from = file.path(paste0(prefix_for_files, "_", eventDataOutfile)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", eventDataOutfile)))
      file.rename(from = file.path(paste0(prefix_for_files, "_", chainSwapFileName)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", chainSwapFileName)))
      if (outputAcceptanceInfo == 1)
      {
        file.rename(from = acceptanceInfoFileName, to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", acceptanceInfoFileName)))
      }
    }
  }

  #### ----------- Step 3: Evaluate BAMM ----------- ####

  cat(paste0("# ----------- Step 3: Evaluate BAMM ----------- #\n\n"))

  ## Produce evaluation plots and ESS data

  {
    if (!skip_evaluations)
    {

      ### 3.1/ Plot MCMC trace for logLik ####

      # Load the MCMC log file
      if (is.null(prefix_for_files))
      {
        MCMC_log <- utils::read.csv(file.path(paste0(BAMM_output_directory_path, mcmcOutfile)), header = T)
      } else {
        MCMC_log <- utils::read.csv(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", mcmcOutfile)), header = T)
      }

      # Find generations used to cut-off the burn-in
      burn_in_threshold <- ceiling(burn_in * MCMC_log$generation[nrow(MCMC_log)])

      # Create binding of new variables to avoid Notes
      logLik <- generation <- NULL

      # Create MCMC trace plot for logLik
      MCMC_logLik_ggplot <- ggplot2::ggplot(data = MCMC_log,
                                            mapping = ggplot2::aes(y = logLik, x = generation)) +

        ggplot2::geom_line(linewidth = 1.0, alpha = 1.0) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +

        ggplot2::geom_vline(xintercept = burn_in_threshold, linewidth = 1.5, linetype = 2, color = "red") +

        ggplot2::labs(x = "Generations", y = "LogLikelihood") +

        ggplot2::ggtitle("MCMC trace for logLik") +

        ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
                       panel.background = ggplot2::element_rect(fill = NA, color = NA),
                       plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                                          margin = ggplot2::margin(b = 15, t = 5)),
                       axis.title = ggplot2::element_text(size = 20, color = "black"),
                       axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                       axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
                       axis.line = ggplot2::element_line(linewidth = 1.5),
                       axis.text = ggplot2::element_text(size = 18, color = "black"),
                       axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

      # Display plot if requested
      if (plot_evaluations)
      {
        print(MCMC_logLik_ggplot)
      }

      # Save plot if requested
      if (save_evaluations)
      {
        if (is.null(prefix_for_files))
        {
          MCMC_logLik_path <- file.path(paste0(BAMM_output_directory_path, "MCMC_trace_logLik.pdf"))
        } else {
          MCMC_logLik_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_MCMC_trace_logLik.pdf"))
        }
        cowplot::save_plot(plot = MCMC_logLik_ggplot,
                           filename = MCMC_logLik_path,
                           base_height = 8, base_width = 10)
      }

      ## Remove burn-in from MCMC log to explore ESS and prior/posterior of LAMBDA (parameter controlling the nb of rates)
      post_burn_MCMC_log <- MCMC_log[MCMC_log$generation >= burn_in_threshold, ]

      ### 3.2/ Explore effective sample sizes ####

      # Store ESS info
      ESS_df <- data.frame(ESS_logLik = coda::effectiveSize(post_burn_MCMC_log$logLik),
                           ESS_logPrior = coda::effectiveSize(post_burn_MCMC_log$logPrior),
                           ESS_N_shifts = coda::effectiveSize(post_burn_MCMC_log$N_shifts),
                           ESS_eventRate = coda::effectiveSize(post_burn_MCMC_log$eventRate),
                           ESS_acceptRate = coda::effectiveSize(post_burn_MCMC_log$acceptRate))
      row.names(ESS_df) <- NULL

      cat("\nEffective sample sizes recorded in the MCMC log after removing burn-in:\n")
      print(ESS_df)
      cat("Ideally, ESS should be higher than 200. Increase the 'numberOfGenerations' or 'eventDataWriteFreq' if needed.\n\n")

      if (save_evaluations)
      {
        if (is.null(prefix_for_files))
        {
          utils::write.csv(x = ESS_df, file = file.path(paste0(BAMM_output_directory_path, "ESS_df.csv")), row.names = FALSE)
        } else {
          utils::write.csv(x = ESS_df, file = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_ESS_df.csv")), row.names = FALSE)
        }
      }

      ### 3.3/ Compare prior and posterior distribution of LAMBDA (parameter controlling the nb of rates) ####

      # Display plot if requested
      if (plot_evaluations)
      {
        BAMMtools::plotPrior(mcmc = MCMC_log,
                             expectedNumberOfShifts = expectedNumberOfShifts,
                             burnin = burn_in,
                             main = paste0("Comparison prior/posterior distributions\n",
                                           "Number of shifts"))
      }
      # Save plot if requested
      if (save_evaluations)
      {
        if (is.null(prefix_for_files))
        {
          PP_nb_shifts_path <- file.path(paste0(BAMM_output_directory_path, "PP_nb_shifts_plot.pdf"))
        } else {
          PP_nb_shifts_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_PP_nb_shifts_plot.pdf"))
        }
        grDevices::pdf(file = file.path(PP_nb_shifts_path),
                       width = 10, height = 8)

        BAMMtools::plotPrior(mcmc = MCMC_log,
                             expectedNumberOfShifts = expectedNumberOfShifts,
                             burnin = burn_in,
                             main = paste0("Comparison prior/posterior distributions\n",
                                           "of the number of regime shifts"))

        grDevices::dev.off()
      }
    }
  }

  #### ----------- Step 4: Import BAMM outputs ----------- ####

  cat(paste0("# ----------- Step 4: Import BAMM outputs ----------- #\n\n"))

  ## Load BAMM object in R and subset posterior samples

  {
    ## Build path to eventData file
    if (is.null(prefix_for_files))
    {
      eventData_path <- file.path(paste0(BAMM_output_directory_path, eventDataOutfile))
    } else {
      eventData_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", eventDataOutfile))
    }

    ## Create the bammdata summarizing BAMM outputs
    BAMM_data_output <- BAMMtools::getEventData(phy = phylo,
                                                eventdata = eventData_path,
                                                burnin = burn_in,
                                                type = "diversification")

    ## Select the subset of posterior samples

    # Get a subset of a selected number of posterior samples
    if (!is.null(seed))
    {
      set.seed(seed = seed)
    }
    sample_indices <- sample(x = 1:length(BAMM_data_output$eventData), size = nb_posterior_samples)
    BAMM_posterior_samples_data <- BAMMtools::subsetEventData(BAMM_data_output, index = sample_indices)

    ## Add the expectedNumberOfShifts as information in the output
    BAMM_posterior_samples_data$expectedNumberOfShifts <- expectedNumberOfShifts

    ## Extract Marginal Shift Probability of each branch and scale branch length accordingly
    MSP_tree <- BAMMtools::marginalShiftProbsTree(BAMM_posterior_samples_data)
    BAMM_posterior_samples_data$MSP_tree <- MSP_tree

    ## Extract the Maximum A Posteriori probability (MAP) configuration = the configuration of shift location showing up the most in the posterior sample
    # Ignore shifts that have an odd-ratio of marginal posterior probability / prior < 'MAP_odd_ratio_threshold' to avoid noise from non-core shifts
    # Rates are then averaged across all samples with the most frequent shift configuration of core-shifts

    MAP_detection <- BAMMtools::credibleShiftSet(BAMM_posterior_samples_data,
                                                 expectedNumberOfShifts = expectedNumberOfShifts,
                                                 threshold = MAP_odd_ratio_threshold,
                                                 set.limit = 0.95)
    # Extract indices of MAP samples
    BAMM_posterior_samples_data$MAP_indices <- MAP_detection$indices[[1]]

    # Compute mean rates/regimes across MAP samples
    MAP_BAMM_object <- BAMMtools::getBestShiftConfiguration(BAMM_posterior_samples_data,
                                                            expectedNumberOfShifts = expectedNumberOfShifts,
                                                            threshold = MAP_odd_ratio_threshold) # Odd-ratio threshold used to select core-shifts used to compare configurations

    # Reorder elements to fit order in the main BAMM_object
    if ("node.label" %in% names(MAP_BAMM_object))
    {
      MAP_BAMM_object <- MAP_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length", "node.label",
                                           "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                           "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                           "meanTipLambda", "meanTipMu", "type")]
    } else {
      MAP_BAMM_object <- MAP_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length",
                                           "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                           "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                           "meanTipLambda", "meanTipMu", "type")]
    }
    class(MAP_BAMM_object) <- "bammdata"
    attr(x = MAP_BAMM_object, which = "order") <- "cladewise"

    BAMM_posterior_samples_data$MAP_BAMM_object <- MAP_BAMM_object

    ## Extract the Maximum Shift Credibility (MSC) configuration = the configuration of shift location with the highest product of marginal probability across branch-specific shifts

    MSC_detection <- BAMMtools::maximumShiftCredibility(BAMM_posterior_samples_data)
    # Extract indices of MSC samples
    BAMM_posterior_samples_data$MSC_indices <- MSC_detection$bestconfigs[[1]]

    # Compute mean rates/regimes across MSC samples
    MSC_BAMM_object <- get_mean_eventData(BAMM_object = BAMM_posterior_samples_data,
                                          sample_indices = MSC_detection$bestconfigs[[1]])

    # Reorder elements to fit order in the main BAMM_object
    if ("node.label" %in% names(MSC_BAMM_object))
    {
      MSC_BAMM_object <- MSC_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length", "node.label",
                                           "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                           "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                           "meanTipLambda", "meanTipMu", "type")]
    } else {
      MSC_BAMM_object <- MSC_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length",
                                           "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                           "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                           "meanTipLambda", "meanTipMu", "type")]
    }
    class(MSC_BAMM_object) <- "bammdata"
    attr(x = MSC_BAMM_object, which = "order") <- "cladewise"

    BAMM_posterior_samples_data$MSC_BAMM_object <- MSC_BAMM_object
  }

  #### ----------- Step 5: Clean BAMM files ----------- ####

  cat(paste0("# ----------- Step 5: Clean BAMM files ----------- #\n\n"))

  ## Remove files generated during the BAMM run

  {
    # If requested, remove all files generated
    if(!keep_BAMM_outputs)
    {
      # Remove phylo, config file and default priors file
      file.remove(phy_path)
      file.remove(config_file_path)
      file.remove(priors_path)

      # Remove files generated during the BAMM run
      if (is.null(prefix_for_files))
      {
        file.remove(file.path(paste0(BAMM_output_directory_path, runInfoFilename)))
        file.remove(file.path(paste0(BAMM_output_directory_path, mcmcOutfile)))
        file.remove(file.path(paste0(BAMM_output_directory_path, eventDataOutfile)))
        file.remove(file.path(paste0(BAMM_output_directory_path, chainSwapFileName)))

        if (outputAcceptanceInfo == 1)
        {
          file.remove(file.path(paste0(BAMM_output_directory_path, acceptanceInfoFileName)))
        }
      } else {
        file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", runInfoFilename)))
        file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", mcmcOutfile)))
        file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", eventDataOutfile)))
        file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", chainSwapFileName)))
        if (outputAcceptanceInfo == 1)
        {
          file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", acceptanceInfoFileName)))
        }
      }

      # If empty, remove the BAMM_output_directory
      if (length(list.files(path = file.path(BAMM_output_directory_path))) == 0)
      {
        unlink(x = file.path(BAMM_output_directory_path), force = TRUE)
      }
    }

  }

  cat(paste0("# ----------- End of BAMM workflow ----------- #\n\n"))

  ## Export BAMM object with posterior samples data
  return(invisible(BAMM_posterior_samples_data))
}


## Ignore BAMM directories during compilation
# usethis::use_build_ignore("software/*")
# usethis::use_build_ignore("BAMM_outputs/*")



### Internal function used to compute mean BAMM_object across a selected set of posterior samples

## Used to aggregate rates/regimes across posterior samples with the Maximum Shift Credibility (MSC) instead of Maximum A Posteriori probability (MAP)

# Source: BAMMtools::getBestShiftConfiguration()
# Author: Dan Rabosky

#' @importFrom ape as.phylo extract.clade
#' @importFrom BAMMtools maximumShiftCredibility subsetEventData getEventData

get_mean_eventData <- function (BAMM_object, sample_indices)
{
  # Extract BAMM_object only for the selected set of samples
  subb <- BAMMtools::subsetEventData(BAMM_object, index = sample_indices)

  # Aggregate all regime shifts events across samples in $eventData
  for (i in 1:length(subb$eventData))
  {
    if (i == 1)
    {
      ff <- subb$eventData[[i]]
    }
    ff <- rbind(ff, subb$eventData[[i]])
  }
  shifts_tipward_nodes_ID <- unique(ff$node)

  # Initiate summary df
  xn <- numeric(length(shifts_tipward_nodes_ID))
  xc <- character(length(shifts_tipward_nodes_ID))
  dff <- data.frame(generation = xn, leftchild = xc, rightchild = xc,
                    abstime = xn, lambdainit = xn, lambdashift = xn,
                    muinit = xn, mushift = xn, stringsAsFactors = F)
  # Extract mean information for each regime
  for (i in 1:length(shifts_tipward_nodes_ID))
  {
    # Extract most left/right descendant tips of the regime
    if (shifts_tipward_nodes_ID[i] <= length(BAMM_object$tip.label))
    {
      dset <- c(BAMM_object$tip.label[shifts_tipward_nodes_ID[i]], NA)
    }
    else {
      tmp <- ape::extract.clade(as.phylo(BAMM_object), node = shifts_tipward_nodes_ID[i])
      dset <- tmp$tip.label[c(1, length(tmp$tip.label))]
    }
    tmp2 <- ff[ff$node == shifts_tipward_nodes_ID[i], ]
    dff$leftchild[i] <- dset[1]
    dff$rightchild[i] <- dset[2]
    # Aggregate mean regime parameters
    dff$abstime[i] <- mean(tmp2$time)
    dff$lambdainit[i] <- mean(tmp2$lam1)
    dff$lambdashift[i] <- mean(tmp2$lam2)
    dff$muinit[i] <- mean(tmp2$mu1)
    dff$mushift[i] <- mean(tmp2$mu2)
  }
  # Extract under bammdata format
  mean_BAMM_object <- BAMMtools::getEventData(ape::as.phylo(BAMM_object), eventdata = dff)

  return(mean_BAMM_object)
}


