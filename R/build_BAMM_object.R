
#' @title Build a BAMM object for a deepSTRAPP run
#'
#' @description Build a BAMM object of class `bammdata` based on the output file of a BAMM run
#'   that contains a phylogenetic tree and associated diversification rates mapped along branches across BAMM posterior samples.
#'
#'   The `BAMM_object` output is typically used as input to run deepSTRAPP with [deepSTRAPP::run_deepSTRAPP_for_focal_time()]
#'   or [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'   This is a wrapper of the [BAMMtools::getEventData()] that additionally provides information on:
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
#' @param nb_posterior_samples Numerical. Number of posterior samples to extract, after removing the burn-in, in the final `BAMM_object` to use for downstream analyses. Default = `1000`.
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
                               nb_posterior_samples = 1000,
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
    # nb_posterior_samples must be a positive integer
    if ((nb_posterior_samples != abs(nb_posterior_samples)) | (nb_posterior_samples != round(nb_posterior_samples)))
    {
      stop(paste0("'nb_posterior_samples' must be a positive integer defining the number of posterior samples retained in the 'BAMM_object' output used for downstream analyses.\n",
                  "Default is '1000'. Current value of 'nb_posterior_samples' is ",nb_posterior_samples,"."))
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

  ## Check if the requested number of posterior samples is compatible with the loaded object and requested burnin.
  nb_samples_loaded <- length(BAMM_data_output$eventData)
  if (nb_posterior_samples > nb_samples_loaded)
  {
    stop(paste0("The 'nb_posterior_samples' exceeds the remaining number of BAMM samples retained after burn-in.\n",
                "Please adjust 'nb_posterior_samples' and/or 'burn_in' to compatible values.\n",
                "Currently, ",nb_samples_loaded," samples were loaded. 'nb_posterior_samples' = ",nb_posterior_samples,"; 'burn_in' = ",burn_in,"."))
  }

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

  ## Export BAMM object with posterior samples data
  return(invisible(BAMM_posterior_samples_data))
}
