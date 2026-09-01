
## Function to prune a BAMM object while keeping deepSTRAPP info ####

#' @title Prune a BAMM object to a subset of tips
#'
#' @description Prune a `BAMM_object` of class `bammdata` down to a subset of tips, and update all
#'   internal elements so that the pruned object still describes the BAMM diversification dynamics,
#'   including all deepSTRAPP elements, but only for the retained tips and branches.
#'
#'   The `BAMM_object` is typically generated directly with [deepSTRAPP::prepare_diversification_data()]
#'   or from external BAMM output files with [deepSTRAPP::build_BAMM_object()].
#'
#'   This function is an extension of the original function [BAMMtools::subtreeBAMM()], which is designed to extract a subclade.
#'   However, this new function also accepts any arbitrary set of tips, whether or not it forms a monophyletic group.
#'   When tips are removed, internal nodes that are left with no descendant are removed,
#'   and internal nodes left with a single descendant are suppressed, their parent and child branches being merged into a single branch.
#'
#'   All BAMM elements are updated accordingly:
#'    * regime shifts located on removed branches are dropped,
#'    * regime shifts located on merged internal branches are re-attached to the new merged branch that contains them,
#'    * macroevolutionary regimes are re-indexed, and branch segments are rebuilt so that each merged branch
#'      carries all the regimes it went through.
#'
#'  This function also preserves and updates the additional deepSTRAPP elements:
#'    * the Marginal Shift Probability (MSP) = the probability of a regime shift to occur along each branch.
#'    * the Maximum A Posteriori probability (MAP) configurations among the posterior samples = the configurations of regimes shifts
#'      that was sampled most frequently (See [BAMMtools::getBestShiftConfiguration()]).
#'    * the Maximum Shift Credibility (MSC) configurations among the posterior samples = the configurations of regime shift location
#'      with the highest product of marginal probabilities across branches (See [BAMMtools::maximumShiftCredibility()]).
#'
#'  Pruning is intended to restrict a downstream deepSTRAPP analysis (e.g., to the taxa for which trait or range data are available)
#'  while keeping the diversification dynamics inferred on the full phylogeny.
#'
#'  If you need diversification rates estimated for a specific clade or set of taxa,
#'  you should run a dedicated BAMM analysis on that subset with appropriate sampling fractions,
#'  see [deepSTRAPP::prepare_diversification_data()].
#'
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::prepare_diversification_data()]
#'   or [deepSTRAPP::build_BAMM_object()], that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#' @param tips_to_keep Vector of character strings. Tips to retain in the pruned `BAMM_object`,
#'   given as tip labels as found in `BAMM_object$tip.label`. Default = `NULL`.
#' @param tips_to_prune Vector of character strings. Tips to remove from the `BAMM_object`,
#'   given as tip labels as found in `BAMM_object$tip.label`. Default = `NULL`.
#' @param MRCA_node Integer. ID of a single internal node of `BAMM_object`, as found in `BAMM_object$edge`,
#'   whose descendant branches/tips must be retained. Use it to focus on the diversification dynamics of one subclade,
#'   as in [BAMMtools::subtreeBAMM()]. Default = `NULL`.
#'
#'   Exactly one of `tips_to_keep`, `tips_to_prune`, and `MRCA_node` must be provided.
#'
#' @param recompute_shift_configurations Logical. Whether the MAP and MSC configurations must be detected again
#'   from the pruned posterior samples.
#'   * If `FALSE` (default), `$MAP_indices` and `$MSC_indices` are left unchanged (pruning removes branches, not
#'     posterior samples, so those indices remain valid), and `$MAP_BAMM_object` and `$MSC_BAMM_object` are simply
#'     pruned along with the main object. Use this to keep the pruned object comparable with the analysis run on
#'     the full phylogeny.
#'   * If `TRUE`, the MAP and MSC configurations are detected again from the pruned posterior samples, as in
#'     [deepSTRAPP::build_BAMM_object()]. Shifts located on removed branches no longer contribute, so the
#'     configurations retained as MAP/MSC may differ from those of the full phylogeny.
#' @param MAP_odd_ratio_threshold Numerical. Controls the definition of 'core-shifts' used to distinguish across configurations when fetching the MAP samples.
#'   Shifts that have an odd-ratio of marginal posterior probability / prior lower than `MAP_odd_ratio_threshold` are ignored. See [BAMMtools::getBestShiftConfiguration()].
#'   Only used when `recompute_shift_configurations = TRUE`. Default = `5`.
#' @param verbose Logical. Whether to display progress in the console. Default = `FALSE`.
#'
#' @export
#' @importFrom ape as.phylo keep.tip reorder.phylo is.rooted is.binary node.depth.edgelength
#' @importFrom BAMMtools credibleShiftSet getmrca marginalShiftProbsTree maximumShiftCredibility
#'
#' @details Prune a `BAMM_object` down to a subset of tips, and update all
#'   internal elements so that the pruned object still describes the BAMM diversification dynamics,
#'   including all deepSTRAPP elements, but only for the retained tips and branches.
#'
#'   When the retained tips do not span the original root, the root of the pruned phylogeny becomes the
#'   Most Recent Common Ancestor (MRCA) of the retained tips. The macroevolutionary regime that was governing
#'   the branch subtending that MRCA becomes the new background/root regime, and its rate parameters are
#'   re-anchored on the new root age (the rates occurring at any given time along the retained branches are unchanged,
#'   only the time of reference at which the initial rates are recorded is shifted).
#'
#'   This function also preserves and updates the additional deepSTRAPP elements:
#'    * the Marginal Shift Probability (MSP) = the probability of a regime shift to occur along each branch.
#'    * the Maximum A Posteriori probability (MAP) configurations among the posterior samples = the configurations of regimes shifts
#'      that was sampled most frequently (See [BAMMtools::getBestShiftConfiguration()]).
#'    * the Maximum Shift Credibility (MSC) configurations among the posterior samples = the configurations of regime shift location
#'      with the highest product of marginal probabilities across branches (See [BAMMtools::maximumShiftCredibility()]).
#'
#'   Those additional elements are used by [deepSTRAPP::plot_BAMM_rates()] to display regime shift probabilities and locations.
#'
#' # Note on what is *not* recomputed
#'
#'   This function does **not** re-estimate diversification rates.
#'   Removing tips changes the incomplete taxon sampling of the phylogeny, but the sampling fractions used
#'   during the original BAMM run are not updated, and rates are not re-inferred.
#'   If you need diversification rates estimated for a specific clade or set of taxa, you should run a
#'   dedicated BAMM analysis on that subset with appropriate sampling fractions,
#'   see [deepSTRAPP::prepare_diversification_data()].
#'
#' # Note on the Marginal Shift Probability tree
#'
#'   The `$MSP_tree` is always recomputed from the pruned posterior samples, independently of
#'   `recompute_shift_configurations`. Marginal shift probabilities are a deterministic function of the
#'   posterior samples and of the topology, and are not a choice of configuration. When several branches are
#'   merged into a single one, the marginal shift probability of the merged branch is the proportion of
#'   posterior samples in which at least one shift occurred anywhere along that merged branch.
#'
#' @return The function returns a `BAMM_object` of class `"bammdata"` which is a list with at least 26 elements.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all retained tips.
#'   * `$edge.length` Vector of numerical. Length of edges/branches. Branches resulting from the merging of several
#'     initial branches have a length equal to the sum of the lengths of the initial branches.
#'   * `$node.label` Vector of character strings. Labels of all internal nodes. (Present only if present in the initial `BAMM_object`)
#'
#'   BAMM internal elements used for tree exploration updated for the new pruned tree:
#'   * `$begin` Vector of numerical. Absolute time since root of edge/branch start (rootward).
#'   * `$end` Vector of numerical.  Absolute time since root of edge/branch end (tipward).
#'   * `$downseq` Vector of integers. Order of node visits when using a pre-order tree traversal.
#'   * `$lastvisit` ID of the last node visited when starting from the node in the corresponding position in `$downseq`.
#'
#'   BAMM elements summarizing diversification data:
#'   * `$numberEvents` Vector of integer. Number of events/macroevolutionary regimes (k+1) retained in each posterior configuration. k = number of shifts.
#'   * `$eventData` List of data.frames. One per posterior sample. Records shift events and macroevolutionary regimes parameters. 1st line = Background root regime.
#'   * `$eventVectors` List of integer vectors. One per posterior sample. Record regime ID per branches.
#'   * `$tipStates` List of integer vectors. One per posterior sample. Record regime ID per tips.
#'     Tip vectors are named after the tips only when they were named in the initial `BAMM_object`:
#'     the naming convention of the input is preserved.
#'   * `$tipLambda` List of numerical vectors. One per posterior sample. Record speciation rates per tips.
#'   * `$tipMu` List of numerical vectors. One per posterior sample. Record extinction rates per tips.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of numerical. Mean tip speciation rates across all posterior configurations of tips.
#'   * `$meanTipMu` Vector of numerical. Mean tip extinction rates across all posterior configurations of tips.
#'   * `$type` Character string. Set the type of data modeled with BAMM. Should be "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch, recomputed on the pruned phylogeny.
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata"` recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configurations, pruned to the retained tips.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata"` recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations, pruned to the retained tips.
#'
#'   Elements tracking the pruning, that can be used to relate the pruned object to the initial one:
#'   * `$pruned_tip_labels` Vector of character strings. Labels of the tips that were removed.
#'   * `$pruning_root_shift` Numerical. Time elapsed between the root of the initial phylogeny and the root of the
#'      pruned phylogeny. Equals `0` when the retained tips span the initial root. Since pruning does not change the
#'      distance of any retained node to the present, this is also the difference between the initial and the pruned root ages.
#'   * `$pruning_nodes_ID_df` Data.frame with four columns providing the conversion table for node IDs:
#'      `$new_node_ID`, `$initial_node_ID`, `$node_type` (`"tip"`, `"root"`, or `"internal"`), and `$node_label`.
#'   * `$pruning_edges_ID_df` Data.frame with four columns providing the conversion table for edge IDs:
#'      `$new_edge_ID`, `$initial_edge_ID`, `$position_in_path` (`1` = the most rootward initial edge merged into the new edge),
#'      and `$nb_merged_edges`. A new edge resulting from the merging of several initial edges holds several rows.
#'
#' @author Maël Doré
#'
#' @seealso Related function in BAMMtools: [BAMMtools::subtreeBAMM()]
#'
#' Associated functions in deepSTRAPP: [deepSTRAPP::prepare_diversification_data()] [deepSTRAPP::build_BAMM_object()]
#' [deepSTRAPP::subset_BAMM_object()] [deepSTRAPP::plot_BAMM_rates()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("import_external_analyses", package = "deepSTRAPP")}
#'
#' @references For BAMM: Rabosky, D. L. (2014). Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees.
#'  PloS one, 9(2), e89543. \doi{10.1371/journal.pone.0089543}. Website: \url{http://bamm-project.org/}.
#'
#'  For `{BAMMtools}`: Rabosky, D. L., Grundler, M., Anderson, C., Title, P., Shi, J. J., Brown, J. W., ... & Larson, J. G. (2014).
#'   BAMM tools: an R package for the analysis of evolutionary dynamics on phylogenetic trees. Methods in Ecology and Evolution, 5(7), 701-707.
#'   \doi{10.1111/2041-210X.12199}
#'
#' @examples
#' # Load BAMM output
#' data(whale_BAMM_object, package = "deepSTRAPP")
#'
#' # Check structure of the initial BAMM_object
#' str(whale_BAMM_object, 1)
#' # Check the initial number of tips
#' length(whale_BAMM_object$tip.label)
#' # We have initially 87 tips in the phylogeny
#'
#' # ----- Example 1: Prune based on set of tips to remove ----- #
#'
#' ## Prune an arbitrary set of tips
#' # Typically, the tips for which no trait or range data is available
#' set.seed(seed = 1234)
#' tips_without_data <- sample(whale_BAMM_object$tip.label, size = 30)
#'
#' whale_BAMM_object_pruned <- prune_BAMM_object(
#'    BAMM_object = whale_BAMM_object,
#'    tips_to_prune = tips_without_data,
#'    verbose = TRUE)
#'
#' # Check structure of the pruned BAMM_object
#' str(whale_BAMM_object_pruned, 1)
#' # Check the updated number of tips
#' length(whale_BAMM_object_pruned$tip.label)
#' # We have now 87 - 30 = 57 tips in the pruned phylogeny
#'
#' # The tips that were removed are recorded in the pruned BAMM_object
#' head(whale_BAMM_object_pruned$pruned_tip_labels)
#'
#' # Branches leading to removed tips are dropped, and the remaining branches are merged
#' # The conversion table records which initial branches were merged into each new branch
#' head(whale_BAMM_object_pruned$pruning_edges_ID_df)
#' table(whale_BAMM_object_pruned$pruning_edges_ID_df$nb_merged_edges)
#'
#' # Diversification rates estimated at the retained tips are left untouched by the pruning
#' retained_tips <- match(whale_BAMM_object_pruned$tip.label, whale_BAMM_object$tip.label)
#' all.equal(whale_BAMM_object_pruned$meanTipLambda,
#'           whale_BAMM_object$meanTipLambda[retained_tips])
#'
#' ## Plot mean rates and MAP regime shifts on the initial vs. pruned phylogeny
#' old_par <- par()$mfrow
#' par(mfrow = c(1, 2))
#'
#' plot_BAMM_rates(whale_BAMM_object, regimes_size = 3)
#' plot_BAMM_rates(whale_BAMM_object_pruned, regimes_size = 3)
#'
#' par(mfrow = old_par)
#'
#' # ----- Example 2: Prune based on set of tips to keep ----- #
#'
#' tips_with_data <- setdiff(whale_BAMM_object$tip.label, tips_without_data)
#'
#' whale_BAMM_object_kept <- prune_BAMM_object(
#'    BAMM_object = whale_BAMM_object,
#'    tips_to_keep = tips_with_data)
#'
#' # Both ways of selecting the tips give the same pruned BAMM_object
#' all.equal(whale_BAMM_object_pruned, whale_BAMM_object_kept)
#'
#' # ----- Example 3: Prune to retain a single subclade ----- #
#'
#' # Plot the initial phylogeny to pick the MRCA node of the focal subclade
#' plot(as.phylo(whale_BAMM_object), cex = 0.5)
#' ape::nodelabels()
#'
#' # Subset BAMM object to focus on node 103 = Odontoceti Infra-order ("toothed whales")
#' whale_BAMM_object_subclade <- prune_BAMM_object(
#'    BAMM_object = whale_BAMM_object,
#'   MRCA_node = 103)
#'
#' # Check the number of tips retained in the subclade
#' length(whale_BAMM_object_subclade$tip.label)
#' # Only 72 odontocete species remain
#'
#' # The root of the pruned phylogeny is now the MRCA of the subclade,
#' # so the phylogeny is shallower than the initial one
#' whale_BAMM_object_subclade$pruning_root_shift # 2 My shift in root_age
#' max(whale_BAMM_object$end) ; max(whale_BAMM_object_subclade$end)
#' # Cetacae phylogeny is 35.4 My old; Odontoceti is phylogeny is 33.4 My old
#'
#' ## Plot mean rates and MAP regime shifts on the pruned phylogeny
#' plot_BAMM_rates(whale_BAMM_object_subclade,
#'                 add_regime_shifts = TRUE,
#'                 configuration_type = "MAP",
#'                 regimes_size = 3)
#'
#' # Since we subsetted diversification dynamics only for the Odontoceti subclade,
#' # we may wish to identify the MAP/MSC configurations based only on the events
#' # occurring along the remaining branches, and not across the full initial tree.
#'
#' # For this, we can set 'recompute_shift_configurations = TRUE'.
#'
#' \donttest{ # This may take several seconds to run
#' ## Detect the MAP/MSC configurations again, on the pruned phylogeny
#' identical(whale_BAMM_object_pruned$MAP_indices, whale_BAMM_object$MAP_indices)
#'
#' # Set 'recompute_shift_configurations = TRUE' to identify the configurations
#' # that are the most supported once the removed branches no longer contribute
#' whale_BAMM_object_recomputed <- prune_BAMM_object(
#'    BAMM_object = whale_BAMM_object,
#'    MRCA_node = 103,
#'    recompute_shift_configurations = TRUE,
#'    verbose = TRUE)
#'
#' # Compare the number of posterior samples supporting the MAP configuration
#' length(whale_BAMM_object$MAP_indices)
#' length(whale_BAMM_object_recomputed$MAP_indices)
#'
#' ## Plot mean rates and updated MAP regime shifts on the pruned phylogeny
#' plot_BAMM_rates(whale_BAMM_object_recomputed,
#'                 add_regime_shifts = TRUE,
#'                 configuration_type = "MAP",
#'                 regimes_size = 3)
#' }
#'

prune_BAMM_object <- function (BAMM_object,
                               tips_to_keep = NULL,
                               tips_to_prune = NULL,
                               MRCA_node = NULL,
                               recompute_shift_configurations = FALSE,
                               MAP_odd_ratio_threshold = 5,
                               verbose = FALSE)
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

    ## Ordering of $edge
    # $eventVectors, $eventBranchSegs, and every edge ID used during the pruning are tied to the ordering
    # of the rows of 'BAMM_object$edge'. That ordering must be "cladewise", as guaranteed by
    # deepSTRAPP::build_BAMM_object(). It is checked here because ape::as.phylo() re-orders the edges
    # when the "order" attribute is not "cladewise", which would silently break that correspondence.
    if (!identical(attr(x = BAMM_object, which = "order"), "cladewise"))
    {
      stop(paste0("The internal ordering of edges in 'BAMM_object$edge' must follow the 'cladewise' order, and be recorded in the 'order' attribute of 'BAMM_object'.\n",
                  "The current 'order' attribute is: ",ifelse(is.null(attr(x = BAMM_object, which = "order")), "NULL", attr(x = BAMM_object, which = "order")),".\n",
                  "Please rebuild your 'BAMM_object' with ?deepSTRAPP::build_BAMM_object() or ?deepSTRAPP::prepare_diversification_data()."))
    }

    ## Extract the phylogeny stored in the BAMM_object
    # ape::as.phylo() dispatches on BAMMtools:::as.phylo.bammdata(), which does not carry $node.label over.
    # Node labels are therefore retrieved directly from the BAMM_object.
    initial_phylo <- ape::as.phylo(BAMM_object)
    initial_node.label <- BAMM_object$node.label

    ## The phylogeny stored in the BAMM_object must be rooted, binary, and free of duplicated tip labels
    if (!(ape::is.rooted(initial_phylo)))
    {
      stop(paste0("The phylogeny stored in 'BAMM_object' must be rooted."))
    }
    if (!(ape::is.binary(initial_phylo)))
    {
      stop(paste0("The phylogeny stored in 'BAMM_object' must be a fully resolved/dichotomous/binary phylogeny."))
    }
    if (anyDuplicated(initial_phylo$tip.label))
    {
      stop(paste0("The phylogeny stored in 'BAMM_object' holds duplicated tip labels, which makes tips impossible to identify unambiguously.\n",
                  "Duplicated labels are: ",collapse_for_message(unique(initial_phylo$tip.label[duplicated(initial_phylo$tip.label)])),"."))
    }

    ## $downseq & $lastvisit are required to explore the initial phylogeny
    if (is.null(BAMM_object$downseq) | is.null(BAMM_object$lastvisit))
    {
      stop(paste0("'BAMM_object' must hold the '$downseq' and '$lastvisit' elements recording the pre-order tree traversal.\n",
                  "See ?BAMMtools::getEventData() and ?deepSTRAPP::build_BAMM_object() to learn how to generate those objects."))
    }

    ## type
    # Only the "diversification" type is supported, as regime parameters are read from $lam1/$lam2/$mu1/$mu2
    if (!identical(BAMM_object$type, "diversification"))
    {
      stop(paste0("'BAMM_object$type' must be 'diversification'. Trait models ('trait' type) are not supported by this function.\n",
                  "Current value of 'BAMM_object$type' is '",BAMM_object$type,"'."))
    }
  }

  initial_Ntip <- length(initial_phylo$tip.label)

  ## Detect the root node ID as the only rootward node that is never the tipward node of an edge
  # This is more robust than assuming the standard 'Ntip + 1' numbering
  initial_root_node_ID <- setdiff(unique(initial_phylo$edge[, 1]), initial_phylo$edge[, 2])
  if (length(initial_root_node_ID) != 1)
  {
    stop(paste0("The root of the phylogeny stored in 'BAMM_object' could not be identified unambiguously.\n",
                "Exactly one node must be a rootward node without ever being a tipward node, but ",length(initial_root_node_ID)," were found.\n",
                "Please check that the phylogeny stored in 'BAMM_object' is a valid rooted phylogeny."))
  }

  ### Check the validity of the arguments defining the tips to retain
  {
    ## tips_to_keep, tips_to_prune & MRCA_node
    # Exactly one of the three must be provided
    provided_inputs <- c("tips_to_keep", "tips_to_prune", "MRCA_node")[!c(is.null(tips_to_keep), is.null(tips_to_prune), is.null(MRCA_node))]
    if (length(provided_inputs) == 0)
    {
      stop(paste0("You must provide exactly one of 'tips_to_keep', 'tips_to_prune', and 'MRCA_node' to define which tips are retained in the pruned 'BAMM_object'."))
    }
    if (length(provided_inputs) > 1)
    {
      stop(paste0("You must provide only one of 'tips_to_keep', 'tips_to_prune', and 'MRCA_node', but not a combination of them.\n",
                  "You currently provided: ",paste(provided_inputs, collapse = " & "),"."))
    }

    ## tips_to_keep & tips_to_prune
    # Tips must be provided as valid tip labels
    for (focal_arg_name in c("tips_to_keep", "tips_to_prune"))
    {
      focal_tips <- get(focal_arg_name)
      if (!is.null(focal_tips))
      {
        if (!is.character(focal_tips))
        {
          stop(paste0("'",focal_arg_name,"' must be a vector of character strings providing tip labels, as found in 'BAMM_object$tip.label'.\n",
                      "Tip indices are not accepted, to avoid any ambiguity on which tips are selected."))
        }
        missing_tips <- setdiff(focal_tips, initial_phylo$tip.label)
        if (length(missing_tips) > 0)
        {
          stop(paste0("Some tips provided in '",focal_arg_name,"' were not found among the tip labels of the phylogeny stored in 'BAMM_object'.\n",
                      "Missing tips are: ",collapse_for_message(missing_tips),"."))
        }
      }
    }

    ## MRCA_node
    # MRCA_node must be a single internal node ID
    if (!is.null(MRCA_node))
    {
      if (!is.numeric(MRCA_node) | (length(MRCA_node) != 1))
      {
        stop(paste0("'MRCA_node' must be a single integer providing the ID of the internal node subtending the subclade to retain, as found in 'BAMM_object$edge'.\n",
                    "To retain several subclades at once, list their tips in 'tips_to_keep' instead."))
      }
      if (MRCA_node != round(MRCA_node))
      {
        stop(paste0("'MRCA_node' must be an integer providing the ID of the internal node subtending the subclade to retain, as found in 'BAMM_object$edge'."))
      }
      # Internal nodes are the nodes that subtend at least one branch
      if (!(MRCA_node %in% initial_phylo$edge[, 1]))
      {
        stop(paste0("'MRCA_node' is not an internal node of the phylogeny stored in 'BAMM_object'.\n",
                    "Internal node IDs are the ",length(unique(initial_phylo$edge[, 1]))," unique values found in 'BAMM_object$edge[, 1]'. Nodes with an ID lower than or equal to ",initial_Ntip," are tips.\n",
                    "Current value of 'MRCA_node' is ",MRCA_node,"."))
      }
    }

    ## recompute_shift_configurations
    if (!is.logical(recompute_shift_configurations) | (length(recompute_shift_configurations) != 1) | is.na(recompute_shift_configurations))
    {
      stop(paste0("'recompute_shift_configurations' must be a single logical value ('TRUE' or 'FALSE')."))
    }

    ## verbose
    if (!is.logical(verbose) | (length(verbose) != 1) | is.na(verbose))
    {
      stop(paste0("'verbose' must be a single logical value ('TRUE' or 'FALSE')."))
    }

    ## expectedNumberOfShifts & MAP_odd_ratio_threshold
    # Only needed when the MAP/MSC configurations must be detected again
    if (recompute_shift_configurations)
    {
      # The expected number of regime shifts is always read from the BAMM_object, so that the value used here
      # is guaranteed to be the one used for the BAMM run that produced it
      expectedNumberOfShifts <- BAMM_object$expectedNumberOfShifts
      if (is.null(expectedNumberOfShifts))
      {
        stop(paste0("'recompute_shift_configurations = TRUE' requires the expected number of regime shifts used to set the prior in BAMM,\n",
                    "but 'BAMM_object$expectedNumberOfShifts' was not found. This element is added by ?deepSTRAPP::build_BAMM_object() and ?deepSTRAPP::prepare_diversification_data().\n",
                    "Please rebuild your 'BAMM_object' with one of those functions, or keep 'recompute_shift_configurations = FALSE'."))
      }
      if ((expectedNumberOfShifts != abs(expectedNumberOfShifts)) | (expectedNumberOfShifts != round(expectedNumberOfShifts)))
      {
        stop(paste0("'BAMM_object$expectedNumberOfShifts' must be a positive integer defining the expected number of diversification regime shifts in the phylogeny.\n",
                    "This value is used to set the hyperprior from which the number of shifts is derived.\n",
                    "Current value of 'BAMM_object$expectedNumberOfShifts' is ",expectedNumberOfShifts,"."))
      }
      if (!is.numeric(MAP_odd_ratio_threshold) | (MAP_odd_ratio_threshold < 0))
      {
        stop(paste0("'MAP_odd_ratio_threshold' must be a positive numerical value. It controls the definition of 'core-shifts' used to distinguish across configurations when fetching the MAP samples.\n",
                    "Shifts that have an odd-ratio of marginal posterior probability / prior lower than 'MAP_odd_ratio_threshold' are ignored. See [BAMMtools::getBestShiftConfiguration()].\n"))
      }
    }
  }

  ### Resolve the set of tips to retain
  {
    ## Case 1: tips provided as tips to keep
    if (!is.null(tips_to_keep))
    {
      retained_tip_labels <- unique(tips_to_keep)
    }

    ## Case 2: tips provided as tips to prune
    if (!is.null(tips_to_prune))
    {
      retained_tip_labels <- setdiff(initial_phylo$tip.label, tips_to_prune)
    }

    ## Case 3: tips defined by the MRCA node of the subclade to keep
    if (!is.null(MRCA_node))
    {
      # Read the descendant tips from the pre-order tree traversal already stored in the BAMM_object
      visited_nodes <- BAMM_object$downseq[which(BAMM_object$downseq == MRCA_node):which(BAMM_object$downseq == BAMM_object$lastvisit[MRCA_node])]
      retained_tip_labels <- initial_phylo$tip.label[visited_nodes[visited_nodes <= initial_Ntip]]
    }

    ## Keep the initial tip ordering, for reproducibility
    retained_tip_labels <- initial_phylo$tip.label[initial_phylo$tip.label %in% retained_tip_labels]
    pruned_tip_labels <- setdiff(initial_phylo$tip.label, retained_tip_labels)

    ## At least 2 tips must be retained to keep a valid rooted phylogeny
    if (length(retained_tip_labels) < 2)
    {
      stop(paste0("At least 2 tips must be retained to build a valid pruned 'BAMM_object' with a root and two descending branches.\n",
                  "Your current selection retains ",length(retained_tip_labels)," tip(s)."))
    }
    # Warn the user when the retained phylogeny becomes too small for meaningful downstream STRAPP tests
    if (length(retained_tip_labels) < 10)
    {
      warning(paste0("Only ",length(retained_tip_labels)," tips are retained in the pruned 'BAMM_object'.\n",
                     "STRAPP tests have very little power on such small phylogenies. Interpret downstream results with caution."))
    }
  }

  ## Nothing to prune: return the initial object untouched
  if (length(pruned_tip_labels) == 0)
  {
    cat(paste0("WARNING: all tips of the phylogeny were retained. The 'BAMM_object' was returned unchanged.\n\n"))
    return(invisible(BAMM_object))
  }

  if (verbose)
  {
    cat(paste0(Sys.time(), " - Pruning the phylogeny: ",length(pruned_tip_labels)," tip(s) removed, ",length(retained_tip_labels)," tip(s) retained.\n"))
  }

  ### Build the pruned topology and the conversion tables between the initial and the pruned phylogeny
  pruning_map <- build_pruning_map(initial_phylo = initial_phylo,
                                   initial_node.label = initial_node.label,
                                   initial_root_node_ID = initial_root_node_ID,
                                   retained_tip_labels = retained_tip_labels)

  if (verbose)
  {
    cat(paste0(Sys.time(), " - Pruned phylogeny built. Root shifted by ",signif(pruning_map$root_shift, 4)," time units. ",
               sum(pruning_map$nb_merged_edges > 1)," branch(es) resulting from merging.\n"))
    cat(paste0(Sys.time(), " - Updating BAMM elements across ",length(BAMM_object$eventData)," posterior sample(s).\n"))
  }

  ### Prune all BAMM elements of the main object
  pruned_BAMM_object <- prune_bammdata_core(BAMM_object = BAMM_object, pruning_map = pruning_map, verbose = verbose)

  ## Set class and ordering attributes right away, as the BAMMtools functions used below require an
  ## object of class "bammdata" holding a valid "order" attribute
  class(pruned_BAMM_object) <- "bammdata"
  attr(x = pruned_BAMM_object, which = "order") <- "cladewise"

  ### Update the additional deepSTRAPP elements

  ## $expectedNumberOfShifts
  if ("expectedNumberOfShifts" %in% names(BAMM_object))
  {
    pruned_BAMM_object$expectedNumberOfShifts <- BAMM_object$expectedNumberOfShifts
  }

  ## $MSP_tree
  # Always recomputed from the pruned posterior samples: marginal shift probabilities of merged branches
  # cannot be recovered from the initial per-branch probabilities.
  if ("MSP_tree" %in% names(BAMM_object))
  {
    if (verbose)
    {
      cat(paste0(Sys.time(), " - Recomputing the Marginal Shift Probability (MSP) tree on the pruned phylogeny.\n"))
    }
    pruned_BAMM_object$MSP_tree <- BAMMtools::marginalShiftProbsTree(pruned_BAMM_object)
  }

  ## $MAP_* and $MSC_*
  if (!recompute_shift_configurations)
  {
    ## Keep the initial MAP/MSC posterior sample indices: pruning removes branches, not posterior samples
    if ("MAP_indices" %in% names(BAMM_object))
    {
      pruned_BAMM_object$MAP_indices <- BAMM_object$MAP_indices
    }
    if ("MAP_BAMM_object" %in% names(BAMM_object))
    {
      if (verbose)
      {
        cat(paste0(Sys.time(), " - Pruning the Maximum A Posteriori probability (MAP) BAMM object.\n"))
      }
      pruned_BAMM_object$MAP_BAMM_object <- format_sub_BAMM_object(prune_bammdata_core(BAMM_object = BAMM_object$MAP_BAMM_object,
                                                                                       pruning_map = pruning_map, verbose = FALSE))
    }
    if ("MSC_indices" %in% names(BAMM_object))
    {
      pruned_BAMM_object$MSC_indices <- BAMM_object$MSC_indices
    }
    if ("MSC_BAMM_object" %in% names(BAMM_object))
    {
      if (verbose)
      {
        cat(paste0(Sys.time(), " - Pruning the Maximum Shift Credibility (MSC) BAMM object.\n"))
      }
      pruned_BAMM_object$MSC_BAMM_object <- format_sub_BAMM_object(prune_bammdata_core(BAMM_object = BAMM_object$MSC_BAMM_object,
                                                                                       pruning_map = pruning_map, verbose = FALSE))
    }

  } else {

    ## Detect the MAP/MSC configurations again, from the pruned posterior samples

    if (verbose)
    {
      cat(paste0(Sys.time(), " - Detecting the Maximum A Posteriori probability (MAP) configurations on the pruned phylogeny.\n"))
    }

    MAP_detection <- BAMMtools::credibleShiftSet(pruned_BAMM_object,
                                                 expectedNumberOfShifts = expectedNumberOfShifts,
                                                 threshold = MAP_odd_ratio_threshold,
                                                 set.limit = 0.95)
    # Extract indices of MAP samples
    pruned_BAMM_object$MAP_indices <- MAP_detection$indices[[1]]

    # Compute mean rates/regimes across MAP samples
    MAP_BAMM_object <- getBestShiftConfiguration_fixed(pruned_BAMM_object,
                                                       expectedNumberOfShifts = expectedNumberOfShifts,
                                                       threshold = MAP_odd_ratio_threshold)
    pruned_BAMM_object$MAP_BAMM_object <- format_sub_BAMM_object(MAP_BAMM_object)

    if (verbose)
    {
      cat(paste0(Sys.time(), " - Detecting the Maximum Shift Credibility (MSC) configurations on the pruned phylogeny.\n"))
    }

    MSC_detection <- BAMMtools::maximumShiftCredibility(pruned_BAMM_object)
    # Extract indices of MSC samples
    pruned_BAMM_object$MSC_indices <- MSC_detection$bestconfigs[[1]]

    # Compute mean rates/regimes across MSC samples
    MSC_BAMM_object <- get_mean_eventData(BAMM_object = pruned_BAMM_object,
                                          sample_indices = MSC_detection$bestconfigs[[1]])
    pruned_BAMM_object$MSC_BAMM_object <- format_sub_BAMM_object(MSC_BAMM_object)
  }

  ### Store the elements tracking the pruning

  pruned_BAMM_object$pruned_tip_labels <- pruned_tip_labels
  pruned_BAMM_object$pruning_root_shift <- pruning_map$root_shift
  pruned_BAMM_object$pruning_nodes_ID_df <- pruning_map$nodes_ID_df
  pruned_BAMM_object$pruning_edges_ID_df <- pruning_map$edges_ID_df

  ### Finalize the output

  if (verbose)
  {
    cat(paste0(Sys.time(), " - Pruning of the BAMM object completed.\n"))
  }

  ## Export output
  return(invisible(pruned_BAMM_object))
}


## Core function to prune all BAMM elements of a "bammdata" object ####

#' @title Prune all BAMM elements of a "bammdata" object
#'
#' @description Transfer and update every element of a `"bammdata"` object onto a pruned phylogeny,
#'   following the conversion map built by `build_pruning_map()`.
#'
#'   This is the engine shared by the main `BAMM_object` and by the `$MAP_BAMM_object` and `$MSC_BAMM_object`
#'   sub-objects, which share the same topology and can therefore be pruned with the same map.
#'
#'   For each posterior sample, the function:
#'   * drops the macroevolutionary regimes whose shift is located on a removed branch.
#'   * re-attaches the regimes whose shift is located on a suppressed internal branch to the merged branch that contains them.
#'   * re-anchors the background/root regime on the new root when the root has been shifted tipward,
#'     re-computing the initial rates it records for the new root age.
#'   * re-indexes the remaining regimes, and rebuilds the branch segments so that merged branches carry all the regimes they went through.
#'
#' @param BAMM_object Object of class `"bammdata"` to prune.
#' @param pruning_map List. Conversion map returned by `build_pruning_map()`.
#' @param verbose Logical. Whether to display progress in the console.
#'
#' @return A list holding the pruned BAMM elements, in the same order as the output of [deepSTRAPP::build_BAMM_object()].
#'
#' @author Maël Doré
#'
#' @noRd
#'

prune_bammdata_core <- function (BAMM_object, pruning_map, verbose = FALSE)
{
  pruned_phylo <- pruning_map$pruned_phylo
  new_Ntip <- length(pruned_phylo$tip.label)
  new_root_node_ID <- pruning_map$new_root_node_ID
  new_edge <- pruned_phylo$edge

  initial_edge <- BAMM_object$edge
  root_shift <- pruning_map$root_shift
  initial_node_of_new_root <- pruning_map$initial_node_of_new_root

  nb_samples <- length(BAMM_object$eventData)

  ## Tolerance used to absorb floating-point errors when comparing times
  time_tolerance <- 10^-8

  ### Initiate the pruned object with the phylogeny-related elements

  pruned_BAMM_object <- list()
  pruned_BAMM_object$edge <- new_edge
  pruned_BAMM_object$Nnode <- pruned_phylo$Nnode
  pruned_BAMM_object$tip.label <- pruned_phylo$tip.label
  pruned_BAMM_object$edge.length <- pruned_phylo$edge.length
  if (!is.null(pruning_map$pruned_node.label))
  {
    pruned_BAMM_object$node.label <- pruning_map$pruned_node.label
  }

  ### BAMM internal elements used for tree exploration

  ## $begin & $end = Absolute time since root of the rootward/tipward end of each branch
  pruned_BAMM_object$begin <- pruning_map$begin
  pruned_BAMM_object$end <- pruning_map$end

  ## $downseq & $lastvisit = Pre-order tree traversal
  pruned_BAMM_object$downseq <- pruning_map$downseq
  pruned_BAMM_object$lastvisit <- pruning_map$lastvisit

  ### Initiate the BAMM elements summarizing diversification data

  pruned_BAMM_object$numberEvents <- integer(nb_samples)
  pruned_BAMM_object$eventData <- vector(mode = "list", length = nb_samples)
  pruned_BAMM_object$eventVectors <- vector(mode = "list", length = nb_samples)
  pruned_BAMM_object$tipStates <- vector(mode = "list", length = nb_samples)
  pruned_BAMM_object$tipLambda <- vector(mode = "list", length = nb_samples)
  pruned_BAMM_object$tipMu <- vector(mode = "list", length = nb_samples)
  pruned_BAMM_object$eventBranchSegs <- vector(mode = "list", length = nb_samples)

  ## Position, in the initial edge matrix, of the most tipward initial branch merged into each pruned branch.
  ## In BAMM, $eventVectors records the regime found at the tipward end of a branch, so the regime of a
  ## pruned branch is the one of the most tipward initial branch it merges.
  most_tipward_initial_edge <- vapply(X = pruning_map$edge_paths,
                                      FUN = function (x) { x[length(x)] },
                                      FUN.VALUE = integer(1))

  ## Position, in the initial edge matrix, of the branch subtending the new root (used to identify the new background regime)
  if (root_shift > 0)
  {
    initial_edge_of_new_root <- which(initial_edge[, 2] == initial_node_of_new_root)
  }

  ### Loop per posterior sample

  for (i in seq_len(nb_samples))
  {
    # i <- 1

    initial_eventData_i <- BAMM_object$eventData[[i]]
    initial_eventVectors_i <- BAMM_object$eventVectors[[i]]

    ## Identify the macroevolutionary regime that becomes the background/root regime of the pruned phylogeny
    if (root_shift > 0)
    {
      # The regime governing the tipward end of the branch subtending the new root
      root_regime_index_i <- initial_eventVectors_i[initial_edge_of_new_root]
    } else {
      # The root has not moved: the initial background regime is kept
      root_regime_index_i <- initial_eventData_i$index[1]
    }

    ## Locate each regime shift on the pruned phylogeny
    # Shifts located on removed branches get a '0' and are dropped, unless they define the new background regime
    new_node_of_event_i <- pruning_map$new_tipward_node_of_initial_node[initial_eventData_i$node]
    new_node_of_event_i[initial_eventData_i$index == root_regime_index_i] <- new_root_node_ID

    ## Keep only the regimes that are still present on the pruned phylogeny
    kept_events_i <- (new_node_of_event_i > 0)
    new_eventData_i <- initial_eventData_i[kept_events_i, , drop = FALSE]
    new_eventData_i$node <- new_node_of_event_i[kept_events_i]

    ## Shift all times into the time frame of the pruned phylogeny
    new_eventData_i$time <- new_eventData_i$time - root_shift

    ## Re-anchor the background/root regime on the new root
    root_regime_position_i <- which(new_eventData_i$index == root_regime_index_i)
    if (length(root_regime_position_i) != 1)
    {
      stop(paste0("Internal error while pruning: the background/root macroevolutionary regime could not be identified in posterior sample n\u00B0",i,".\n",
                  "Please report the issue at https://github.com/MaelDore/deepSTRAPP/issues."))
    }
    if (new_eventData_i$time[root_regime_position_i] < -time_tolerance)
    {
      ## The regime started before the new root: the initial rates it records are re-computed for the new root age.
      ## The rate reached at any given time along the retained branches is unchanged, only the time of reference is shifted.
      elapsed_time_i <- -new_eventData_i$time[root_regime_position_i]
      new_eventData_i$lam1[root_regime_position_i] <- compute_rate_at_time(rate_0 = new_eventData_i$lam1[root_regime_position_i],
                                                                           shape = new_eventData_i$lam2[root_regime_position_i],
                                                                           elapsed_time = elapsed_time_i)
      new_eventData_i$mu1[root_regime_position_i] <- compute_rate_at_time(rate_0 = new_eventData_i$mu1[root_regime_position_i],
                                                                          shape = new_eventData_i$mu2[root_regime_position_i],
                                                                          elapsed_time = elapsed_time_i)
    }
    new_eventData_i$time[root_regime_position_i] <- 0

    ## Absorb residual floating-point errors on the remaining times
    new_eventData_i$time[new_eventData_i$time < 0] <- 0

    ## Re-order the regimes: background/root regime first, then by increasing time, as expected by BAMM
    events_order_i <- order(new_eventData_i$index != root_regime_index_i, new_eventData_i$time)
    new_eventData_i <- new_eventData_i[events_order_i, , drop = FALSE]

    ## Re-index the regimes from 1 to the number of retained regimes
    initial_regime_indices_i <- new_eventData_i$index
    new_eventData_i$index <- seq_len(nrow(new_eventData_i))
    rownames(new_eventData_i) <- NULL

    ## $numberEvents = Number of events/macroevolutionary regimes retained
    pruned_BAMM_object$numberEvents[i] <- nrow(new_eventData_i)

    ## $eventData = Shift events and macroevolutionary regime parameters
    pruned_BAMM_object$eventData[[i]] <- new_eventData_i

    ## $eventVectors = Regime ID per branch
    new_eventVectors_i <- match(x = initial_eventVectors_i[most_tipward_initial_edge], table = initial_regime_indices_i)
    if (anyNA(new_eventVectors_i))
    {
      stop(paste0("Internal error while pruning: some branches of the pruned phylogeny were assigned a macroevolutionary regime that was dropped, in posterior sample n\u00B0",i,".\n",
                  "Please report the issue at https://github.com/MaelDore/deepSTRAPP/issues."))
    }
    pruned_BAMM_object$eventVectors[[i]] <- new_eventVectors_i

    ## $tipStates = Regime ID per tip
    new_tipStates_i <- match(x = BAMM_object$tipStates[[i]][pruning_map$initial_tip_of_new_tip], table = initial_regime_indices_i)
    if (anyNA(new_tipStates_i))
    {
      stop(paste0("Internal error while pruning: some tips of the pruned phylogeny were assigned a macroevolutionary regime that was dropped, in posterior sample n\u00B0",i,".\n",
                  "Please report the issue at https://github.com/MaelDore/deepSTRAPP/issues."))
    }
    # match() drops the names, which are restored from the initial object.
    # Whether tip vectors are named depends on how the 'BAMM_object' was produced:
    # BAMMtools::getEventData() returns unnamed tip vectors, while
    # deepSTRAPP::update_rates_and_regimes_for_focal_time() returns named ones.
    # The convention of the initial object is preserved, so that pruning stays a pure subsetting operation.
    names(new_tipStates_i) <- names(BAMM_object$tipStates[[i]])[pruning_map$initial_tip_of_new_tip]
    pruned_BAMM_object$tipStates[[i]] <- new_tipStates_i

    ## $tipLambda & $tipMu = Speciation and extinction rates per tip
    # Tip ages are unchanged by pruning, and re-anchoring the background regime leaves the rate reached
    # at any given time untouched. Tip rates are therefore only subset and re-ordered.
    # Subsetting carries the names of the initial object over, when present.
    pruned_BAMM_object$tipLambda[[i]] <- BAMM_object$tipLambda[[i]][pruning_map$initial_tip_of_new_tip]
    pruned_BAMM_object$tipMu[[i]] <- BAMM_object$tipMu[[i]][pruning_map$initial_tip_of_new_tip]

    ## $eventBranchSegs = Regime ID per segment of branch
    # Rebuilt from scratch: a merged branch can carry several regimes, and thus several segments.
    pruned_BAMM_object$eventBranchSegs[[i]] <- build_eventBranchSegs(new_edge = new_edge,
                                                                     new_begin = pruning_map$begin,
                                                                     new_end = pruning_map$end,
                                                                     new_eventData = new_eventData_i,
                                                                     new_eventVectors = new_eventVectors_i,
                                                                     new_root_node_ID = new_root_node_ID)

    ## Print progress
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - BAMM elements pruned for posterior sample n\u00B0", i, "/", nb_samples,"\n"))
    }
  }

  ### Mean tip rates across posterior samples

  ## $meanTipLambda & $meanTipMu = Mean tip speciation and extinction rates
  # Subsetting carries the names of the initial object over, when present
  pruned_BAMM_object$meanTipLambda <- BAMM_object$meanTipLambda[pruning_map$initial_tip_of_new_tip]
  pruned_BAMM_object$meanTipMu <- BAMM_object$meanTipMu[pruning_map$initial_tip_of_new_tip]

  ## $type = Type of data modeled with BAMM
  pruned_BAMM_object$type <- BAMM_object$type

  ### Export the pruned BAMM elements
  return(pruned_BAMM_object)
}


## Helper function to build the pruning map ####

#' @title Build the conversion map between an initial phylogeny and its pruned version
#'
#' @description Prune a phylogeny down to a set of tips, and build all the conversion tables needed
#'   to transfer BAMM elements from the initial phylogeny onto the pruned one.
#'
#'   The key output is `$edge_paths`: when tips are removed, internal nodes left with a single descendant
#'   are suppressed, so that a branch of the pruned phylogeny corresponds to an ordered *path* of one or
#'   several branches of the initial phylogeny. `$edge_paths` records that path for each pruned branch,
#'   ordered from the most rootward to the most tipward initial branch.
#'
#' @param initial_phylo Object of class `"phylo"`. The initial, rooted, binary phylogeny.
#'   The ordering of `$edge` must match the one of the initial `BAMM_object`, as edge IDs refer to it.
#' @param initial_node.label Vector of character strings, or `NULL`. Node labels of the initial phylogeny.
#' @param initial_root_node_ID Integer. ID of the root of the initial phylogeny.
#' @param retained_tip_labels Vector of character strings. Labels of the tips to retain.
#'
#' @return A list with the following elements:
#'   * `$pruned_phylo` Object of class `"phylo"`. The pruned phylogeny, in `"cladewise"` order.
#'   * `$pruned_node.label` Vector of character strings, or `NULL`. Node labels of the pruned phylogeny.
#'   * `$initial_node_of_new_node` Vector of integers, indexed by pruned node ID. ID of the matching node in the initial phylogeny.
#'   * `$initial_tip_of_new_tip` Vector of integers, indexed by pruned tip ID. ID of the matching tip in the initial phylogeny.
#'   * `$edge_paths` List of vectors of integers, one per pruned branch. IDs of the initial branches merged into that pruned branch,
#'      ordered from the most rootward to the most tipward.
#'   * `$nb_merged_edges` Vector of integers. Number of initial branches merged into each pruned branch.
#'   * `$new_tipward_node_of_initial_node` Vector of integers, indexed by initial node ID. ID of the tipward node of the pruned branch
#'      that contains the initial branch subtending that node. Holds `0` for initial nodes that are not covered by the pruned phylogeny.
#'   * `$root_shift` Numerical. Time elapsed between the root of the initial phylogeny and the root of the pruned phylogeny.
#'   * `$initial_node_of_new_root` Integer. ID, in the initial phylogeny, of the node that became the root of the pruned phylogeny.
#'   * `$new_root_node_ID` Integer. ID of the root of the pruned phylogeny.
#'   * `$begin`, `$end` Vectors of numerical. Start and stop times of each pruned branch, in the time frame of the pruned phylogeny.
#'   * `$downseq`, `$lastvisit` Vectors of integers. Pre-order tree traversal of the pruned phylogeny.
#'   * `$nodes_ID_df`, `$edges_ID_df` Data.frames. User-facing conversion tables, stored in the pruned `BAMM_object`.
#'
#' @author Maël Doré
#'
#' @noRd
#'

build_pruning_map <- function (initial_phylo, initial_node.label, initial_root_node_ID, retained_tip_labels)
{
  initial_Ntip <- length(initial_phylo$tip.label)

  ### Build the pruned phylogeny
  # ape::keep.tip() suppresses internal nodes left with a single descendant, and sums the lengths
  # of the merged branches
  pruned_phylo <- ape::keep.tip(phy = initial_phylo, tip = retained_tip_labels)
  pruned_phylo <- ape::reorder.phylo(pruned_phylo, order = "cladewise")

  new_Ntip <- length(pruned_phylo$tip.label)
  new_Nnode <- pruned_phylo$Nnode

  ## Detect the root of the pruned phylogeny the same robust way as for the initial phylogeny
  new_root_node_ID <- setdiff(unique(pruned_phylo$edge[, 1]), pruned_phylo$edge[, 2])
  if ((length(new_root_node_ID) != 1) | !identical(as.integer(new_root_node_ID), as.integer(new_Ntip + 1L)))
  {
    stop(paste0("Internal error while pruning: the pruned phylogeny does not follow the standard node numbering expected by BAMM,\n",
                "where tips hold IDs 1 to Ntip and the root holds the ID Ntip + 1.\n"))
  }

  ### Compute the pre-order tree traversal of the pruned phylogeny
  pruned_recursive_sequence <- getRecursiveSequence(phy = pruned_phylo)

  ### Match each node of the pruned phylogeny with its counterpart in the initial phylogeny

  initial_node_of_new_node <- integer(new_Ntip + new_Nnode)

  ## Tips are matched by label
  initial_node_of_new_node[1:new_Ntip] <- match(pruned_phylo$tip.label, initial_phylo$tip.label)

  ## Internal nodes are matched with the MRCA, in the initial phylogeny, of their descendant tips.
  # The MRCA of the first and the last descendant tips found in the pre-order traversal is enough:
  # those two tips descend from the two different children of the focal node, so their MRCA is the
  # MRCA of the whole clade, both in the pruned and in the initial phylogeny.
  new_internal_nodes_ID <- (new_Ntip + 1L):(new_Ntip + new_Nnode)
  first_descendant_tip <- last_descendant_tip <- integer(new_Nnode)

  for (j in seq_along(new_internal_nodes_ID))
  {
    focal_node_j <- new_internal_nodes_ID[j]
    visited_nodes_j <- pruned_recursive_sequence$downseq[which(pruned_recursive_sequence$downseq == focal_node_j):which(pruned_recursive_sequence$downseq == pruned_recursive_sequence$lastvisit[focal_node_j])]
    visited_tips_j <- visited_nodes_j[visited_nodes_j <= new_Ntip]
    first_descendant_tip[j] <- visited_tips_j[1]
    last_descendant_tip[j] <- visited_tips_j[length(visited_tips_j)]
  }

  initial_node_of_new_node[new_internal_nodes_ID] <- BAMMtools::getmrca(phy = initial_phylo,
                                                                        t1 = initial_node_of_new_node[first_descendant_tip],
                                                                        t2 = initial_node_of_new_node[last_descendant_tip])

  if (any(initial_node_of_new_node == 0))
  {
    stop(paste0("Internal error while pruning: some nodes of the pruned phylogeny could not be matched with the initial phylogeny.\n",
                "Please check that the phylogeny stored in 'BAMM_object' is rooted and fully dichotomous."))
  }

  ### Match each branch of the pruned phylogeny with the ordered path of initial branches it merges

  ## Describe the initial phylogeny: rootward node and subtending branch of each node
  # Both are indexed by node ID, and are therefore independent of the ordering of $edge
  parent_of_node <- integer(max(initial_phylo$edge))
  parent_of_node[initial_phylo$edge[, 2]] <- initial_phylo$edge[, 1]

  edge_of_node <- integer(max(initial_phylo$edge))
  edge_of_node[initial_phylo$edge[, 2]] <- seq_len(nrow(initial_phylo$edge))

  edge_paths <- vector(mode = "list", length = nrow(pruned_phylo$edge))

  for (k in seq_len(nrow(pruned_phylo$edge)))
  {
    rootward_node_k <- initial_node_of_new_node[pruned_phylo$edge[k, 1]]
    tipward_node_k <- initial_node_of_new_node[pruned_phylo$edge[k, 2]]

    ## Walk rootward from the tipward node, collecting branches until the rootward node is reached
    path_k <- integer(0)
    current_node <- tipward_node_k
    while (current_node != rootward_node_k)
    {
      # Prepend, so that the path ends up ordered from the most rootward to the most tipward branch
      path_k <- c(edge_of_node[current_node], path_k)
      current_node <- parent_of_node[current_node]
      if (current_node == 0)
      {
        stop(paste0("Internal error while pruning: the path between two nodes of the pruned phylogeny could not be traced back in the initial phylogeny.\n",
                    "Please report the issue at https://github.com/MaelDore/deepSTRAPP/issues."))
      }
    }
    edge_paths[[k]] <- path_k
  }

  nb_merged_edges <- lengths(edge_paths)

  ### Build the lookup giving, for each initial node, the tipward node of the pruned branch that covers it

  new_tipward_node_of_initial_node <- integer(max(initial_phylo$edge))
  for (k in seq_len(nrow(pruned_phylo$edge)))
  {
    covered_initial_nodes <- initial_phylo$edge[edge_paths[[k]], 2]
    new_tipward_node_of_initial_node[covered_initial_nodes] <- pruned_phylo$edge[k, 2]
  }

  ### Re-derive node labels from the conversion map

  # {ape} carries $node.label over when pruning, but its handling of suppressed internal nodes can be
  # ambiguous. Node labels are therefore re-derived from the node conversion map, which is unambiguous.
  if (!is.null(initial_node.label))
  {
    pruned_node.label <- initial_node.label[initial_node_of_new_node[new_internal_nodes_ID] - initial_Ntip]
  } else {
    pruned_node.label <- NULL
  }

  ### Compute the root shift and the new branch times

  ## Node heights = absolute time since root, which is the time frame used by BAMM
  initial_node_heights <- ape::node.depth.edgelength(initial_phylo)
  new_node_heights <- ape::node.depth.edgelength(pruned_phylo)

  initial_node_of_new_root <- initial_node_of_new_node[new_root_node_ID]
  root_shift <- initial_node_heights[initial_node_of_new_root]

  new_begin <- new_node_heights[pruned_phylo$edge[, 1]]
  new_end <- new_node_heights[pruned_phylo$edge[, 2]]

  ### Build the conversion tables

  ## Node labels: tip labels for tips, node labels for internal nodes when present
  if (!is.null(pruned_node.label))
  {
    new_node_labels <- c(pruned_phylo$tip.label, pruned_node.label)
  } else {
    new_node_labels <- c(pruned_phylo$tip.label, rep(NA_character_, new_Nnode))
  }

  nodes_ID_df <- data.frame(new_node_ID = seq_len(new_Ntip + new_Nnode),
                            initial_node_ID = initial_node_of_new_node,
                            node_type = c(rep("tip", new_Ntip), "root", rep("internal", new_Nnode - 1L)),
                            node_label = new_node_labels,
                            stringsAsFactors = FALSE)

  edges_ID_df <- data.frame(new_edge_ID = rep(seq_along(edge_paths), times = nb_merged_edges),
                            initial_edge_ID = unlist(edge_paths, use.names = FALSE),
                            position_in_path = unlist(lapply(X = edge_paths, FUN = seq_along), use.names = FALSE),
                            nb_merged_edges = rep(nb_merged_edges, times = nb_merged_edges),
                            stringsAsFactors = FALSE)

  ### Export the pruning map

  pruning_map <- list(pruned_phylo = pruned_phylo,
                      pruned_node.label = pruned_node.label,
                      initial_node_of_new_node = initial_node_of_new_node,
                      initial_tip_of_new_tip = initial_node_of_new_node[1:new_Ntip],
                      edge_paths = edge_paths,
                      nb_merged_edges = nb_merged_edges,
                      new_tipward_node_of_initial_node = new_tipward_node_of_initial_node,
                      root_shift = root_shift,
                      initial_node_of_new_root = initial_node_of_new_root,
                      initial_root_node_ID = initial_root_node_ID,
                      new_root_node_ID = as.integer(new_root_node_ID),
                      begin = new_begin,
                      end = new_end,
                      downseq = pruned_recursive_sequence$downseq,
                      lastvisit = pruned_recursive_sequence$lastvisit,
                      nodes_ID_df = nodes_ID_df,
                      edges_ID_df = edges_ID_df)

  return(pruning_map)
}



## Helper function to rebuild the matrix of branch segments ####

#' @title Rebuild the matrix of branch segments of a pruned phylogeny
#'
#' @description Rebuild the `$eventBranchSegs` matrix of a posterior sample, recording the macroevolutionary
#'   regime found along each segment of each branch.
#'
#'   Branches that carry no regime shift hold a single segment spanning the whole branch. Branches that carry
#'   `n` regime shifts are split into `n + 1` segments. The resulting matrix therefore holds
#'   `number of branches + number of shifts` rows.
#'
#'  Code was directly adapted from the original [BAMMtools::getEventData()] function.
#'
#' @param new_edge Matrix of integers. The `$edge` element of the pruned phylogeny.
#' @param new_begin,new_end Vectors of numerical. The `$begin` and `$end` elements of the pruned phylogeny.
#' @param new_eventData Data.frame. The pruned and re-indexed `$eventData` of the focal posterior sample.
#'   The background/root regime must be the first row, and hold `index = 1`.
#' @param new_eventVectors Vector of integers. The pruned `$eventVectors` of the focal posterior sample.
#' @param new_root_node_ID Integer. ID of the root of the pruned phylogeny.
#'
#' @return A matrix of numerical with four columns: tipward node ID of the branch, start time of the segment,
#'   stop time of the segment, and ID of the macroevolutionary regime found along the segment.
#'   Rows are ordered by tipward node ID, then from the most rootward to the most tipward segment.
#'
#' @author Maël Doré. Original [BAMMtools::getEventData()] function written by Dan Rabosky & Mike Grundler.
#'
#' @seealso Initial function in BAMMtools: [BAMMtools::getEventData()]
#'
#' @noRd
#'

build_eventBranchSegs <- function (new_edge, new_begin, new_end, new_eventData, new_eventVectors, new_root_node_ID)
{
  ## Initiate with one segment per branch, spanning the whole branch
  eventBranchSegs <- cbind(new_edge[, 2], new_begin, new_end, new_eventVectors)
  dimnames(eventBranchSegs) <- NULL

  ## Identify the branches carrying at least one regime shift
  # The background/root regime holds index 1 and is anchored on the root, which subtends no branch
  shift_events <- new_eventData[new_eventData$index != 1L, , drop = FALSE]
  branches_with_shifts <- unique(shift_events$node)

  ## Branches without shifts keep their single segment
  eventBranchSegs <- eventBranchSegs[!(eventBranchSegs[, 1] %in% branches_with_shifts), , drop = FALSE]

  ## Split the branches carrying shifts into as many segments as needed
  for (focal_node in branches_with_shifts)
  {
    branch_events <- shift_events[shift_events$node == focal_node, , drop = FALSE]
    branch_events <- branch_events[order(branch_events$time), , drop = FALSE]

    focal_edge <- which(new_edge[, 2] == focal_node)
    rootward_node <- new_edge[focal_edge, 1]

    ## Regime found at the rootward end of the branch
    if (rootward_node == new_root_node_ID)
    {
      # The branch descends directly from the root: it starts under the background/root regime
      current_regime <- 1L
    } else {
      current_regime <- new_eventVectors[which(new_edge[, 2] == rootward_node)]
    }
    current_time <- new_begin[focal_edge]

    ## Add one segment per shift
    branch_segments <- matrix(data = NA_real_, nrow = nrow(branch_events) + 1L, ncol = 4)
    for (j in seq_len(nrow(branch_events)))
    {
      branch_segments[j, ] <- c(focal_node, current_time, branch_events$time[j], current_regime)
      current_regime <- branch_events$index[j]
      current_time <- branch_events$time[j]
    }
    ## Add the last segment, running from the last shift to the tipward end of the branch
    branch_segments[nrow(branch_segments), ] <- c(focal_node, current_time, new_end[focal_edge], current_regime)

    eventBranchSegs <- rbind(eventBranchSegs, branch_segments)
  }

  ## Order by tipward node ID, then from the most rootward to the most tipward segment
  eventBranchSegs <- eventBranchSegs[order(eventBranchSegs[, 1], eventBranchSegs[, 2]), , drop = FALSE]
  dimnames(eventBranchSegs) <- NULL

  return(eventBranchSegs)
}



## Helper function to truncate long vectors displayed in error messages ####

#' @title Collapse a vector into a truncated character string
#'
#' @description Collapse a vector into a single character string to be displayed in an error message,
#'   truncating it when it holds more than `max_items` elements.
#'
#' @param x Vector. The elements to display.
#' @param max_items Integer. Maximum number of elements displayed. Default = `20`.
#'
#' @return A character string.
#'
#' @author Maël Doré
#'
#' @noRd
#'

collapse_for_message <- function (x, max_items = 20)
{
  if (length(x) > max_items)
  {
    return(paste0(paste(x[seq_len(max_items)], collapse = ", "), ", ... (", length(x), " elements in total)"))
  }

  return(paste(x, collapse = ", "))
}


## Helper function to compute a rate at a given time within a macroevolutionary regime ####

#' @title Compute the rate reached after a given time within a macroevolutionary regime
#'
#' @description Compute the value of a time-dependent BAMM rate after `elapsed_time` time units
#'   spent within a macroevolutionary regime defined by an initial rate and a shape parameter.
#'
#'   BAMM models time-dependence as an exponential decay when the shape parameter is negative,
#'   and as a saturating increase when it is positive:
#'   * if `shape <= 0`: `rate_t = rate_0 * exp(shape * t)`
#'   * if `shape > 0`: `rate_t = rate_0 * (2 - exp(-shape * t))`
#'
#'   This is the same parametrization as the one used in
#'   [deepSTRAPP::update_rates_and_regimes_for_focal_time()], and it reproduces the unexported
#'   `BAMMtools:::exponentialRate()`, which cannot be called from here.
#'
#' @param rate_0 Numerical. Initial rate of the macroevolutionary regime.
#' @param shape Numerical. Shape parameter controlling the time-dependence of the rate.
#' @param elapsed_time Numerical. Time elapsed since the start of the macroevolutionary regime.
#'
#' @return A numerical value. The rate reached after `elapsed_time`.
#'
#' @author Maël Doré
#'
#' @noRd
#'

compute_rate_at_time <- function (rate_0, shape, elapsed_time)
{
  if (is.na(shape) | is.na(rate_0) | is.na(elapsed_time))
  {
    return(NA_real_)
  }

  if (shape <= 0)
  {
    ## Decrease: rate_t = rate_0 * exp(shape * t)
    rate_t <- rate_0 * exp(shape * elapsed_time)
  } else {
    ## Increase: rate_t = rate_0 * (2 - exp(-shape * t))
    rate_t <- rate_0 * (2 - exp(-shape * elapsed_time))
  }

  return(rate_t)
}




## Helper function to format the MAP/MSC sub-objects ####

#' @title Format a MAP/MSC BAMM sub-object
#'
#' @description Re-order the elements of a `$MAP_BAMM_object` or `$MSC_BAMM_object` to match the structure
#'   used in [deepSTRAPP::build_BAMM_object()], and set its class and ordering attributes.
#'
#' @param sub_BAMM_object List. The MAP or MSC BAMM sub-object to format.
#'
#' @return An object of class `"bammdata"` holding 18 elements (19 with `$node.label`).
#'
#' @author Maël Doré
#'
#' @noRd
#'

format_sub_BAMM_object <- function (sub_BAMM_object)
{
  ## Reorder elements to fit the order used in the main BAMM_object
  if ("node.label" %in% names(sub_BAMM_object))
  {
    sub_BAMM_object <- sub_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length", "node.label",
                                         "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                         "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                         "meanTipLambda", "meanTipMu", "type")]
  } else {
    sub_BAMM_object <- sub_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length",
                                         "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                         "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                         "meanTipLambda", "meanTipMu", "type")]
  }

  class(sub_BAMM_object) <- "bammdata"
  attr(x = sub_BAMM_object, which = "order") <- "cladewise"

  return(sub_BAMM_object)
}


