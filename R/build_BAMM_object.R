
## Function to build a BAMM object from a BAMM output eventdata.txt file ####

#' @title Build a BAMM object for a deepSTRAPP run
#'
#' @description Build a BAMM object of class `bammdata` based on the output file of a BAMM run
#'   that contains a phylogenetic tree and associated diversification rates mapped along branches across BAMM posterior samples.
#'
#'   The `BAMM_object` output is typically used as input to run deepSTRAPP with [deepSTRAPP::run_deepSTRAPP_for_focal_time()]
#'   or [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'   This is a wrapper of the original [BAMMtools::getEventData()] function that additionally provides information on:
#'    * the Marginal Shift Probability (MSP) = the probability of a regime shift to occur along each branch.
#'    * the Maximum A Posteriori probability (MAP) configurations among the posterior samples = the configurations of regimes shifts
#'      that was sampled most frequently (See [BAMMtools::getBestShiftConfiguration()]).
#'    * the Maximum Shift Credibility (MSC) configurations among the posterior samples = the configurations of regime shift location
#'      with the highest product of marginal probabilities across branches (See [BAMMtools::maximumShiftCredibility()]).
#'
#'   Those additional elements are used by [deepSTRAPP::plot_BAMM_rates()] to display regime shift probabilities and locations.
#'
#'   This function is meant to enable users to inject into the deepSTRAPP framework the results of their own BAMM analyses.
#'   Alternatively, a full BAMM analyses starting from a time-calibrated phylogeny alone can be carried within deepSTRAPP with [deepSTRAPP::prepare_diversification_data()].
#'
#' # Note on Bayesian Analysis of Macroevolutionary Mixtures (BAMM)
#'
#'  BAMM is a model of diversification for time-calibrated phylogenies that explores complex diversification dynamics
#'  by allowing multiple regime shifts across clades without a priori hypotheses on the location of such shifts.
#'  It uses reversible jump Markov chain Monte Carlo (rjMCMC) to automatically explore a vast range of models with different
#'  speciation and extinction rates, and different number and location of regime shits.
#'
#'  BAMM is one option among others for modeling diversification on phylogenies.
#'  You may wish to explore alternatives models such as LSBDS model in RevBayes (Höhna et al., 2016), the MTBD model (Barido-Sottani et al., 2020),
#'  or the ClaDS2 model (Maliet et al., 2019) for your own data.
#'  However, you will need Bayesian models that infer regime shifts to be able to perform STRAPP tests (Rabosky & Huang, 2016).
#'  Additionally, you need to format the model output such as in `BAMM_object`, so it can be used in a deepSTRAPP workflow.
#'
#' @param phylo Object of class `"phylo"` as defined in R package `{ape}`. Time-calibrated phylogeny that was used to produce the BAMM run.
#'   The phylogeny must be rooted and fully resolved.
#' @param eventdata Character string specifying the path to a BAMM event-data file. Alternatively, an object of class data.frame that includes the event data from a BAMM run.
#' @param burn_in Numerical. Proportion of posterior samples removed from the BAMM output to ensure that the remaining samples where drawn once the equilibrium distribution was reached. Default is `0.25`
#' @param nb_posterior_samples Numerical. Number of posterior samples to extract, after removing the burn-in, in the final `BAMM_object` to use for downstream analyses.
#'  If set to `NULL` (default), all samples remaining after removing the burn-in will be kept.
#' @param seed Integer. Set the seed to ensure reproducibility when drawing random posterior samples. Default is `NULL` (a random seed is used).
#' @param expectedNumberOfShifts Integer. The expected number of regime shifts sets during the BAMM run as an hyperparameter controlling the exponential prior distribution
#'  used to modulate reversible jumps across model configurations in the rjMCMC run. This is needed to compute priors for regime shift along branches.
#' @param MAP_odd_ratio_threshold Numerical. Controls the definition of 'core-shifts' used to distinguish across configurations when fetching the MAP samples.
#'   Shifts that have an odd-ratio of marginal posterior probability / prior lower than `MAP_odd_ratio_threshold` are ignored. See [BAMMtools::getBestShiftConfiguration()]. Default = `5`.
#' @param verbose Logical. Whether to display progress in the console. Default = `FALSE`.
#'
#' @export
#' @importFrom BAMMtools getEventData subsetEventData
#'
#' @return The function returns a `BAMM_object` of class `"bammdata"` which is a list with at least 23 elements.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'   * `$node.label` Vector of character strings. Labels of all internal nodes. (Present only if present in the initial `phylo`)
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
#' @author Maël Doré
#'
#' @seealso Initial functions in BAMMtools: [BAMMtools::getEventData()] [BAMMtools::getBestShiftConfiguration()] [BAMMtools::maximumShiftCredibility()]
#'
#' Associated functions in deepSTRAPP: [deepSTRAPP::prepare_diversification_data()] [deepSTRAPP::plot_BAMM_rates()]
#' [deepSTRAPP::run_deepSTRAPP_for_focal_time()] [deepSTRAPP::run_deepSTRAPP_over_time()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("model_diversification_dynamics", package = "deepSTRAPP")}
#'
#' @references For BAMM: Rabosky, D. L. (2014). Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees.
#'  PloS one, 9(2), e89543. \doi{10.1371/journal.pone.0089543}. Website: \url{http://bamm-project.org/}.
#'
#'  For `{BAMMtools}`: Rabosky, D. L., Grundler, M., Anderson, C., Title, P., Shi, J. J., Brown, J. W., ... & Larson, J. G. (2014).
#'   BAMM tools: an R package for the analysis of evolutionary dynamics on phylogenetic trees. Methods in Ecology and Evolution, 5(7), 701-707.
#'   \doi{10.1111/2041-210X.12199}
#'
#' @examples
#' # The key output from a BAMM is the 'event_data.txt' file
#' # It can be loaded in R to use as input for deepSTRAPP
#'
#' library(phytools)
#' data(whale.tree)
#'
#' \dontrun{
#' ## The 'whale_event_data.txt' file used for example here is not provided within deepSTRAPP
#' BAMM_object <- build_BAMM_object(
#'    phylo = whale.tree,
#'    eventdata = "./BAMM_outputs/whale_event_data.txt",
#'    burn_in = 0.25, # Remove 25% as burn-in
#'    nb_posterior_samples = 1000, # Retain 1000 samples
#'    expectedNumberOfShifts = 1,
#'    verbose = TRUE)
#' str(BAMM_object, 1)
#' }
#'

build_BAMM_object <- function (phylo,
                               eventdata,
                               burn_in = 0.25,
                               nb_posterior_samples = NULL,
                               seed = NULL,
                               expectedNumberOfShifts,
                               MAP_odd_ratio_threshold = 5,
                               verbose = FALSE)
{
  ### Check input validity
  {
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
    # phylo must be in "cladewise" order for BAMM to works
    phylo_order <- attr(x = phylo, which = "order")
    if (phylo_order != "cladewise")
    {
      phylo <- ape::reorder.phylo(phylo, order = "cladewise")
      cat(paste0("WARNING: the internal ordering of edges in 'phylo$edge' must follow the 'cladewise' order for BAMM to work.\n",
                 "Your 'phylo' object was in '", phylo_order,"'. It was modified to follow the 'cladewise' structure.\n\n"))
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
    # nb_posterior_samples must be a positive integer or NULL
    if (!is.null(nb_posterior_samples))
    {
      if ((nb_posterior_samples != abs(nb_posterior_samples)) | (nb_posterior_samples != round(nb_posterior_samples)))
      {
        stop(paste0("'nb_posterior_samples' must be a positive integer defining the number of posterior samples retained in the 'BAMM_object' output used for downstream analyses.\n",
                    "Alternatively, you can set 'nb_posterior_samples' to 'NULL' to retain all remaining posterior samples after removing burn-in.\n",
                    "You can also use [deepSTRAPP::subset_BAMM_object()] to subset the BAMM_object after having built it with this function.\n",
                    "Current value of 'nb_posterior_samples' is ",nb_posterior_samples,"."))
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

    ## MAP_odd_ratio_threshold
    # If provided, MAP_odd_ratio_threshold must be a positive numerical
    if (!is.numeric(MAP_odd_ratio_threshold) | (MAP_odd_ratio_threshold < 0))
    {
      stop(paste0("'expectedNumberOfShifts' must be a positive numerical value. It controls the definition of 'core-shifts' used to distinguish across configurations when fetching the MAP samples.\n",
                  "Shifts that have an odd-ratio of marginal posterior probability / prior lower than `MAP_odd_ratio_threshold` are ignored. See [BAMMtools::getBestShiftConfiguration()].\n"))
    }
  }

  ## Load BAMM object in R and subset posterior samples

  ## Build path to eventData file
  eventData_path <- file.path(paste0(eventdata))

  ## Create the bammdata summarizing BAMM outputs
  BAMM_data_output <- BAMMtools::getEventData(phy = phylo,
                                              eventdata = eventData_path,
                                              burnin = burn_in,
                                              type = "diversification",
                                              verbose = verbose)

  ## If not provided, set the number of posterior samples to equal the number of samples remaining after burn-in
  if (is.null(nb_posterior_samples))
  {
    nb_posterior_samples <- length(BAMM_data_output$eventData)
    cat(paste0("WARNING: you set 'nb_posterior_samples = NULL', thus all posterior samples remaining after burn-in have been retained.\n",
               "The final number of BAMM posterior samples is ",nb_posterior_samples,".\n\n"))
  }

  ## Check if the requested number of posterior samples is compatible with the loaded object and requested burnin.
  nb_samples_loaded <- length(BAMM_data_output$eventData)
  if (nb_posterior_samples > nb_samples_loaded)
  {
    stop(paste0("The 'nb_posterior_samples' exceeds the remaining number of BAMM samples retained after burn-in.\n",
                "Please adjust 'nb_posterior_samples' and/or 'burn_in' to compatible values.\n",
                "Currently, ",nb_samples_loaded," samples were loaded. 'nb_posterior_samples' = ",nb_posterior_samples,"; 'burn_in' = ",burn_in,"."))
  }

  ## Select the subset of posterior samples (if needed)

  if (nb_posterior_samples < nb_samples_loaded)
  {
    # Get a subset of a selected number of posterior samples
    if (!is.null(seed))
    {
      set.seed(seed = seed)
    }
    sample_indices <- sample(x = 1:length(BAMM_data_output$eventData), size = nb_posterior_samples)
    BAMM_posterior_samples_data <- BAMMtools::subsetEventData(BAMM_data_output, index = sample_indices)
  } else {
    # Case when the number or posterior samples requested = the number of samples loaded
    # No need to subsample the loaded BAMM object
    BAMM_posterior_samples_data <- BAMM_data_output
  }

  ## Name the tip-level elements with the tip labels
  # BAMMtools::getEventData() returns $tipStates, $tipLambda, $tipMu, $meanTipLambda and $meanTipMu
  # unnamed, ordered by tip ID. Naming them makes the mapping to tips explicit, and is required by the
  # deepSTRAPP functions that match tips by label, such as deepSTRAPP::compute_STRAPP_test_for_focal_time().
  BAMM_posterior_samples_data <- name_tip_elements(BAMM_posterior_samples_data)

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
  # MAP_BAMM_object <- BAMMtools::getBestShiftConfiguration(BAMM_posterior_samples_data,
  MAP_BAMM_object <- getBestShiftConfiguration_fixed(BAMM_posterior_samples_data,
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
  ## Name the tip-level elements with the tip labels, as for the main BAMM_object
  MAP_BAMM_object <- name_tip_elements(MAP_BAMM_object)

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
  ## Name the tip-level elements with the tip labels, as for the main BAMM_object
  MSC_BAMM_object <- name_tip_elements(MSC_BAMM_object)

  class(MSC_BAMM_object) <- "bammdata"
  attr(x = MSC_BAMM_object, which = "order") <- "cladewise"

  BAMM_posterior_samples_data$MSC_BAMM_object <- MSC_BAMM_object

  ## Export BAMM object with posterior samples data
  return(invisible(BAMM_posterior_samples_data))
}


## Function to subject a BAMM object while keeping deepSTRAPP info ####

#' @title Subset a BAMM object before a deepSTRAPP run
#'
#' @description Subset a BAMM object of class `bammdata` while keeping additional information for deepSTRAPP.
#'
#'   The `BAMM_object` is typically generated directly with [deepSTRAPP::prepare_diversification_data()]
#'   or from external BAMM output files with [deepSTRAPP::build_BAMM_object()].
#'
#'   This is a wrapper of the original [BAMMtools::subsetEventData()] function that additionally preserves information on:
#'    * the Marginal Shift Probability (MSP) = the probability of a regime shift to occur along each branch.
#'    * the Maximum A Posteriori probability (MAP) configurations among the posterior samples = the configurations of regimes shifts
#'      that was sampled most frequently (See [BAMMtools::getBestShiftConfiguration()]).
#'    * the Maximum Shift Credibility (MSC) configurations among the posterior samples = the configurations of regime shift location
#'      with the highest product of marginal probabilities across branches (See [BAMMtools::maximumShiftCredibility()])
#'    that are stored in a `BAMM_object` when build with deepSTRAPP functions.
#'
#'   Those additional elements are used by [deepSTRAPP::plot_BAMM_rates()] to display regime shift probabilities and locations.
#'
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::prepare_diversification_data()]
#'   or [deepSTRAPP::build_BAMM_object()], that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#' @param nb_posterior_samples Integer. Number of posterior samples to extract from `BAMM_object`. Default = `NULL`.
#' @param sample_indices Integer or vector of integers. Indices of the posterior samples to extract. Default = `NULL`.
#' @param seed Integer. Set the seed to ensure reproducibility when `sample_indices` is not provided and posterior samples are drawn randomly.
#'   Default is `NULL` (a random seed is used).
#'
#' @export
#' @importFrom BAMMtools subsetEventData
#'
#' @return The function returns a `BAMM_object` of class `"bammdata"` which is a list with at least 23 elements.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'   * `$node.label` Vector of character strings. Labels of all internal nodes. (Present only if present in the initial `phylo`)
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
#' @author Maël Doré
#'
#' @seealso Initial function in BAMMtools: [BAMMtools::subsetEventData()]
#'
#' Associated functions in deepSTRAPP: [deepSTRAPP::prepare_diversification_data()] [deepSTRAPP::build_BAMM_object()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("model_diversification_dynamics", package = "deepSTRAPP")}
#'
#' @references For BAMM: Rabosky, D. L. (2014). Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees.
#'  PloS one, 9(2), e89543. \doi{10.1371/journal.pone.0089543}. Website: \url{http://bamm-project.org/}.
#'
#'  For `{BAMMtools}`: Rabosky, D. L., Grundler, M., Anderson, C., Title, P., Shi, J. J., Brown, J. W., ... & Larson, J. G. (2014).
#'   BAMM tools: an R package for the analysis of evolutionary dynamics on phylogenetic trees. Methods in Ecology and Evolution, 5(7), 701-707.
#'   \doi{10.1111/2041-210X.12199}
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  ## Load BAMM object
#'  # data(Ponerinae_BAMM_object_old_calib)
#'  # This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  # Check structure of BAMM_object
#'  str(Ponerinae_BAMM_object_old_calib, 1)
#'  # Check current number of BAMM posterior samples
#'  length(Ponerinae_BAMM_object_old_calib$eventData)
#'  # We have initially 1000 posterior samples in the updated BAMM object
#'
#'  ## Subset BAMM_object
#'  BAMM_object_subset <- subset_BAMM_object(
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     nb_posterior_samples = 100,
#'     seed = 1234)
#'
#'  # Check structure of the updated BAMM_object
#'  str(BAMM_object_subset, 1)
#'  # Check updated number of BAMM posterior samples
#'  length(BAMM_object_subset$eventData)
#'  # We have now 100 posterior samples in the updated BAMM object
#' }
#'

subset_BAMM_object <- function (BAMM_object,
                                nb_posterior_samples = NULL,
                                sample_indices = NULL,
                                seed = NULL)
{
  ### Check input validity
  {
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

    ## nb_posterior_samples & sample_indices
    # One of the two must be provided, but not the two
    if (!is.null(nb_posterior_samples) & !is.null(sample_indices))
    {
      stop("You must provide one of 'nb_posterior_samples' and 'sample_indices', but not both.")
    }
    if (is.null(nb_posterior_samples) & is.null(sample_indices))
    {
      stop("You must provide at least one of 'nb_posterior_samples' and 'sample_indices'.")
    }

    ## nb_posterior_samples
    if (!is.null(nb_posterior_samples))
    {
      # nb_posterior_samples must be an integer
      if (!(nb_posterior_samples - round(nb_posterior_samples) == 0))
      {
        stop("'nb_posterior_samples' must be an integer.")
      }
      # nb_posterior_samples must be within the range of available BAMM samples
      if (!(nb_posterior_samples <= length(BAMM_object$eventData)))
      {
        stop("'nb_posterior_samples' must be lower than the number of available BAMM samples.\n",
             "Your BAMM_object has ",length(BAMM_object$eventData)," posterior samples available.\n",
             "'nb_posterior_samples' = ", nb_posterior_samples, ".")
      }
    }

    ## sample_indices
    if (!is.null(sample_indices))
    {
      # sample_indices must be an integer or vector of integers
      if (!all((sample_indices - round(sample_indices) == 0)))
      {
        stop("'sample_indices' must be an integer or vector of integers.")
      }
      # sample_indices must fit within the range of available BAMM samples
      if (!(all(sample_indices <= length(BAMM_object$eventData))))
      {
        stop("'sample_indices' must fit within the range of available BAMM samples.\n",
             "Your BAMM_object has ",length(BAMM_object$eventData)," posterior samples available.")
      }
    }

    ## seed
    if (!is.null(seed))
    {
      if (!(seed - round(seed) == 0))
      {
        stop(paste0("'seed' must be an integer."))
      }
    }
  }

  ## Draw random indices if needed
  if (!is.null(nb_posterior_samples))
  {
    sample_indices <- sample(x = 1:length(BAMM_object$eventData), size = nb_posterior_samples, replace = F)
  }

  ## Subset using BAMMtools::subsetEventData
  BAMM_object_subset <- BAMMtools::subsetEventData(ephy = BAMM_object, index = sample_indices)

  ## Concatenate the additional deepSTRAPP elements into the BAMM_object
  BAMM_object_subset[c("expectedNumberOfShifts", "MSP_tree", "MAP_indices", "MAP_BAMM_object", "MSC_indices", "MSC_BAMM_object")] <- BAMM_object[c("expectedNumberOfShifts", "MSP_tree", "MAP_indices", "MAP_BAMM_object", "MSC_indices", "MSC_BAMM_object")]

  ## Export output
  return(invisible(BAMM_object_subset))
}


## Helper function to name the tip-level elements of a BAMM object ####

#' @title Name the tip-level elements of a BAMM object
#'
#' @description Assign tip labels as names to every element of a `"bammdata"` object that holds
#'   one value per tip: `$tipStates`, `$tipLambda`, `$tipMu`, `$meanTipLambda` and `$meanTipMu`.
#'
#'   [BAMMtools::getEventData()] returns those elements unnamed, ordered by tip ID, so that the i-th
#'   value corresponds to the i-th label in `$tip.label`. Naming them makes the mapping explicit, and
#'   makes the object directly usable by the deepSTRAPP functions that match tips by label,
#'   such as [deepSTRAPP::compute_STRAPP_test_for_focal_time()], which requires the names of
#'   `trait_data_list$trait_data` to match those of `BAMM_object$tipStates`.
#'
#'   Elements that are already named are simply overwritten with the same labels, so the function
#'   is safe to apply to an object that went through
#'   [deepSTRAPP::update_rates_and_regimes_for_focal_time()].
#'
#' @param BAMM_object Object of class `"bammdata"`, or a list holding the same elements.
#'
#' @return The input object, with named tip-level elements.
#'
#' @author Maël Doré
#'
#' @noRd
#'

name_tip_elements <- function (BAMM_object)
{
  tip_labels <- BAMM_object$tip.label
  nb_tips <- length(tip_labels)

  ## Elements holding one value per tip, in each posterior sample
  for (element in c("tipStates", "tipLambda", "tipMu"))
  {
    if (!is.null(BAMM_object[[element]]))
    {
      # Safety check: the element must hold exactly one value per tip
      element_lengths <- unlist(lapply(X = BAMM_object[[element]], FUN = length))
      if (any(element_lengths != nb_tips))
      {
        stop(paste0("Internal error: 'BAMM_object$",element,"' must hold exactly one value per tip to be named with the tip labels.\n",
                    "The phylogeny holds ",nb_tips," tips, but the posterior samples hold ",paste(unique(element_lengths), collapse = ", ")," values."))
      }

      BAMM_object[[element]] <- lapply(X = BAMM_object[[element]],
                                       FUN = function (x) { names(x) <- tip_labels ; return(x) })
    }
  }

  ## Elements holding one value per tip, averaged across posterior samples
  for (element in c("meanTipLambda", "meanTipMu"))
  {
    if (!is.null(BAMM_object[[element]]))
    {
      # Safety check: the element must hold exactly one value per tip
      if (length(BAMM_object[[element]]) != nb_tips)
      {
        stop(paste0("Internal error: 'BAMM_object$",element,"' must hold exactly one value per tip to be named with the tip labels.\n",
                    "The phylogeny holds ",nb_tips," tips, but the element holds ",length(BAMM_object[[element]])," values."))
      }

      names(BAMM_object[[element]]) <- tip_labels
    }
  }

  return(BAMM_object)
}
