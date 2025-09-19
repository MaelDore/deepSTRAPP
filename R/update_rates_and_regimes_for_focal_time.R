#' @title Update diversification rates/regimes mapped on a phylogeny up to a given time in the past
#'
#' @description Updates an object of class `"bammdata"` to obtain the diversification rates/regimes
#'   found along branches at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the object to display a mapped phylogeny such as
#'   branches overlapping the `focal_time` are shorten to the `focal_time`.
#'
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::prepare_diversification_data()],
#'   that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   at which the tree and rate mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_rates Logical. Specify whether diversification rates stored in
#'   `$tipLambda` (speciation) and `$tipMu` (extinction) must be updated to summarize
#'   rates found along branches at `focal_time`. Default is `TRUE`.
#' @param update_regimes Logical. Specify whether diversification regimes stored in
#'   `$tipStates` must be updated to summarize regimes found along branches at `focal_time`.
#'   Default is `TRUE`.
#' @param update_tree Logical. Specify whether to update the phylogeny such as
#'   all branches that are younger than the `focal_time` are cut-off. Default is `FALSE`.
#' @param update_plot Logical. Specify whether to update the phylogeny AND the elements
#'   used by [deepSTRAPP::plot_BAMM_rates()] to plot diversification rates on the phylogeny
#'   such as all branches that are younger than the `focal_time` are cut-off. Default is `FALSE`.
#'   If set to `TRUE`, it will override the `update_tree` parameter and update the phylogeny.
#' @param update_all_elements Logical. Specify whether to update all the elements in the object, including
#'   rates/regimes/phylogeny/elements for [deepSTRAPP::plot_BAMM_rates()]/all other elements. Default is `FALSE`.
#'   If set to `TRUE`, it will override other `update_*` parameters and update all elements.
#' @param keep_tip_labels Logical. Should terminal branches with a single descendant tip
#'   retain their initial `tip.label` on the updated phylogeny and diversification rate mapping?
#'   Default is `TRUE`. If set to `FALSE`, the tipward node ID will be used as label for all tips.
#' @param verbose Logical. Should progression be displayed? A message will be printed for every batch of
#'   100 BAMM posterior samples updated. Default is `TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights getDescendants
#' @importFrom BAMMtools plot.bammdata dtRates
#'
#' @details The object of class `"bammdata"` (`BAMM_object`) is cut-off at a specific time in the past
#'   (i.e. the `focal_time`) and the current diversification rate values of the overlapping edges/branches are extracted.
#'
#'   ----- Update diversification rate data -----
#'
#'   If `update_rates = TRUE`, diversification rates of branches overlapping with `focal_time`
#'   will be updated. Each cut-off branches form a new tip present at `focal_time` with updated associated
#'   diversification rate data. Fossils older than `focal_time` do not yield any data.
#'   * `$tipLambda` contains speciation rates.
#'   * `$tipMu` contains extinction rates.
#'
#'   If `update_regimes = TRUE`, diversification regimes of branches overlapping with `focal_time`
#'   will be updated. Each cut-off branches form a new tip present at `focal_time` with updated associated
#'   diversification regime ID found in `$tipStates`. Fossils older than `focal_time` do not yield any data.
#'
#'   ----- Update the phylogeny -----
#'
#'   If `update_tree = TRUE`, elements defining the phylogeny, as in an `"phylo"` object
#'   will be updated such as all branches that are younger than the `focal_time` are cut-off:
#'   * `$edge` defines the tree topology.
#'   * `$Nnode` defines the number of internal nodes.
#'   * `$tip.label` provides the labels of all tips, including fossils older than `focal_time` if present.
#'   * `$edge.length` provides length of all branches.
#'   * `$node.label` provides the labels of all internal nodes. (Optional)
#'
#'   ----- Update the plot from [deepSTRAPP::plot_BAMM_rates()] -----
#'
#'   If `update_plot = TRUE`, elements used to plot diversification rates on the phylogeny
#'   using [deepSTRAPP::plot_BAMM_rates()] will be updated such as all branches that are younger
#'   than the `focal_time` are cut-off:
#'   * `$begin` provides absolute time since root of edge/branch start (rootward).
#'   * `$end` provides absolute time since root of edge/branch end (tipward).
#'   * `$eventVectors` provides regime membership per branches in each posterior sample configuration.
#'   * `$eventBranchSegs` provides regime membership per segments of branches in each posterior sample configuration.
#'   * `$dtrates` provides mean speciation and extinction rates along segments of branches, and resolution fraction (tau) describing
#'   the fraction of each segment length compared to the full depth of the initial tree (i.e., the root_age).
#'
#' @return The function returns a list as an updated `BAMM_object` of class `"bammdata"`.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips, including fossils older than `focal_time` if present.
#'     + If `keep_tip_labels = TRUE`, cut-off branches with a single descendant tip retain their initial `tip.label`.
#'     + If `keep_tip_labels = FALSE`, all cut-off branches are labeled using their tipward node ID.
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
#'   * `$tipStates` List of named integer vectors. One per posterior sample. Record regime ID per tips present at `focal_time`. Updated if `update_regimes = TRUE`.
#'   * `$tipLambda` List of named numerical vectors. One per posterior sample. Record speciation rates per tips present at `focal_time`. Updated if `update_rates = TRUE`.
#'   * `$tipMu` List of named numerical vectors. One per posterior sample. Record extinction rates per tips present at `focal_time`. Updated if `update_rates = TRUE`.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of named numerical. Mean tip speciation rates across all posterior configurations of tips present at `focal_time` (does not includes older fossils).
#'   * `$meanTipMu` Vector of named numerical. Mean tip extinction rates across all posterior configurations of tips present at `focal_time` (does not includes older fossils).
#'   * `$type` Character string. Set the type of data modeled with BAMM. Should be "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'      whose origin is older that `focal_time`.
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history, updated for the `focal_time`.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history, updated for the `focal_time`.
#'
#'   New elements added to provide update information:
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'   * `$dtrates` List of three elements.
#'     + 1/ `$dtrates$tau` Numerical. Resolution factor describing the fraction of each segment length used in [deepSTRAPP::plot_BAMM_rates()]
#'       compared to the full depth of the initial tree (i.e., the root_age)
#'     + 2/ `$dtrates$rates` List of two numerical vectors. Speciation and extinction rates along segments used by [deepSTRAPP::plot_BAMM_rates()].
#'     + 3/ `$dtrates$tmat` Matrix of numerical. Start and end times of segments in term of distance to the root.
#'   * `$initial_colorbreaks` List of three vectors of numerical. Rate values of the percentiles delimiting the bins for mapping rates to colors with [BAMMtools::plot.bammdata()].
#'     Each element provides values for different type of rates (`$speciation`, `$extinction`, `$net_diversification`).
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the rates/regimes were extracted and the tree was eventually cut.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::plot_BAMM_rates()]
#'
#' @examples
#' # ----- Example 1: Extant whales (87 taxa) ----- #
#'
#' ## Load the BAMM_object summarizing 1000 posterior samples of BAMM
#' data(whale_BAMM_object, package = "deepSTRAPP")
#'
#' ## Set focal-time to 5 My
#' focal_time = 5
#'
#' ## Update the BAMM object
#' whale_BAMM_object_5My <- update_rates_and_regimes_for_focal_time(
#'    BAMM_object = whale_BAMM_object,
#'    focal_time = 5,
#'    update_rates = TRUE, update_regimes = TRUE,
#'    update_tree = TRUE, update_plot = TRUE,
#'    update_all_elements = TRUE,
#'    keep_tip_labels = TRUE,
#'    verbose = TRUE)
#'
#' # Add "phylo" class to be compatible with phytools::nodeHeights()
#' class(whale_BAMM_object) <- unique(c(class(whale_BAMM_object), "phylo"))
#' root_age <- max(phytools::nodeHeights(whale_BAMM_object)[,2])
#' # Remove temporary "phylo" class
#' class(whale_BAMM_object) <- setdiff(class(whale_BAMM_object), "phylo")
#'
#' ## Plot initial BAMM_object for t = 0 My
#' plot_BAMM_rates(whale_BAMM_object, add_regime_shifts = TRUE,
#'                 labels = TRUE, legend = TRUE,
#'                 par.reset = FALSE) # Keep plotting parameters in memory to use abline().
#' abline(v = root_age - focal_time,
#'       col = "red", lty = 2, lwd = 2)
#'
#' ## Plot updated BAMM_object for t = 5 My
#' plot_BAMM_rates(whale_BAMM_object_5My, add_regime_shifts = TRUE,
#'                 labels = TRUE, legend = TRUE)
#'
#' # ----- Example 2: Extant Ponerinae (1,534 taxa) ----- #
#'
#' ## Load the BAMM_object summarizing 1000 posterior samples of BAMM
#' data(Ponerinae_BAMM_object, package = "deepSTRAPP")
#'
#' ## Set focal-time to 10 My
#' focal_time = 10
#'
#' \dontrun{
#' ## Update the BAMM object (May take several minutes to run)
#' Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
#'   BAMM_object = Ponerinae_BAMM_object,
#'   focal_time = focal_time,
#'   update_rates = TRUE, update_regimes = TRUE,
#'   update_tree = TRUE, update_plot = TRUE,
#'   update_all_elements = TRUE,
#'   keep_tip_labels = TRUE,
#'   verbose = TRUE) }
#'
#' ## Load results to save time
#' data(Ponerinae_BAMM_object_10My, package = "deepSTRAPP")
#'
#' ## Extract root age
#' # Add temporarily the "phylo" class to be compatible with phytools::nodeHeights()
#' class(Ponerinae_BAMM_object) <- unique(c(class(Ponerinae_BAMM_object), "phylo"))
#' root_age <- max(phytools::nodeHeights(Ponerinae_BAMM_object)[,2])
#' # Remove temporary "phylo" class
#' class(Ponerinae_BAMM_object) <- setdiff(class(Ponerinae_BAMM_object), "phylo")
#'
#' ## Plot diversification rates on the initial tree
#' plot_BAMM_rates(Ponerinae_BAMM_object,
#'                 legend = TRUE, labels = FALSE)
#' abline(v = root_age - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' ## Plot diversification rates and regime shifts on the updated tree (cut-off for 10 My)
#' # Keep the initial color scheme
#' plot_BAMM_rates(Ponerinae_BAMM_object_10My, legend = TRUE, labels = FALSE,
#'                 colorbreaks = Ponerinae_BAMM_object_10My$initial_colorbreaks$net_diversification)
#'
#' # Use a new color scheme mapped on the new distribution of rates
#' plot_BAMM_rates(Ponerinae_BAMM_object_10My, legend = TRUE, labels = FALSE)
#'


# # ----- Example 3: Non-ultrametric tree including extinct mammal groups ----- #
#
# ### Ideally, run BAMM on motmot::mammals dataset so I have a real dataset with fossils
# # (but BAMM does not work with fossils yet)
# # Need to check if the BAMM_output object generated deal with fossils in a specific way...
#
# str(updated_BAMM_object$tipStates, max.level = 1)
# str(updated_BAMM_object, max.level = 1)
#
# pdf(file = "./test_BAMMplot_t0.pdf", width = 20, height = 150)
# plot_BAMM_rates(BAMM_object, labels = TRUE)
# dev.off()


update_rates_and_regimes_for_focal_time <- function (BAMM_object, focal_time,
                                                     update_rates = TRUE, update_regimes = TRUE,
                                                     update_tree = FALSE, update_plot = FALSE,
                                                     update_all_elements = FALSE,
                                                     keep_tip_labels = TRUE,
                                                     verbose = TRUE)
{
  ### Check input validity
  {
    ## Extract root age
    # Add "phylo" class to be compatible with phytools::nodeHeights()
    class(BAMM_object) <- unique(c(class(BAMM_object), "phylo"))
    root_age <- max(phytools::nodeHeights(BAMM_object)[,2])

    ## BAMM_object
    # BAMM_object must be of class "bammdata"
    if (!("bammdata" %in% class(BAMM_object)))
    {
      stop("'BAMM_object' must have the 'bammdata' class. See ?BAMMtools::getEventData() to learn how to generate those objects.")
    }
    # BAMM_object must contain at least the elements $eventData, $tipStates, $tipLambda and $tipMu.
    if (!all(c("eventData", "tipStates", "tipLambda", "tipMu") %in% names(BAMM_object)))
    {
      stop(paste0("'BAMM_object' must contain at least the elements $eventData, $tipStates, $tipLambda and $tipMu to extract diversification rates and regimes.\n",
                  "See ?BAMMtools::getEventData() to learn how to generate those objects."))
    }
    # Number of posterior sample data must be equal between $eventData, $tipStates, $tipLambda and $tipMu
    posterior_samples_length <- c(length(BAMM_object$eventData), length(BAMM_object$tipStates), length(BAMM_object$tipLambda), length(BAMM_object$tipMu))
    if (length(unique(posterior_samples_length)) != 1)
    {
      stop(paste0("Number of posterior samples in 'BAMM_object' must be equal between $eventData, $tipStates, $tipLambda and $tipMu.\n",
                  "Please check the structure of your 'BAMM_object' with str(BAMM_object, 1).\n",
                  "See ?BAMMtools::getEventData() to learn how to generate those objects."))
    }
    # BAMM_object must have an $tipStates that is a list of N posterior samples with integer vector of regime membership per tips
    # Check that all integer vectors have a length equal to $tip.label
    if (!all(unlist(lapply(BAMM_object$tipStates, FUN = length)) == length(BAMM_object$tip.label)))
    {
      stop(paste0("Number of values in posterior samples of 'BAMM_object$tipStates' must be equal to the number of tips in the phylogeny.\n",
                  "Please check the structure of your 'BAMM_object' with str(BAMM_object$tipStates, 1).\n",
                  "See ?BAMMtools::getEventData() to learn how to generate those objects."))
    }
    # BAMM_object must have an $tipLambda that is a list of N posterior samples with integer vector of final speciation rates at tips = current speciation rates
    if (!all(unlist(lapply(BAMM_object$tipLambda, FUN = length)) == length(BAMM_object$tip.label)))
    {
      stop(paste0("Number of values in posterior samples of 'BAMM_object$tipLambda' must be equal to the number of tips in the phylogeny.\n",
                  "Please check the structure of your 'BAMM_object' with str(BAMM_object$tipLambda, 1).\n",
                  "See ?BAMMtools::getEventData() to learn how to generate those objects."))
    }
    # BAMM_object must have an $tipMu that is a list of N posterior samples with integer vector of final extinction rates at tips = current extinction rates
    if (!all(unlist(lapply(BAMM_object$tipMu, FUN = length)) == length(BAMM_object$tip.label)))
    {
      stop(paste0("Number of values in posterior samples of 'BAMM_object$tipMu' must be equal to the number of tips in the phylogeny.\n",
                  "Please check the structure of your 'BAMM_object' with str(BAMM_object$tipMu, 1).\n",
                  "See ?BAMMtools::getEventData() to learn how to generate those objects."))
    }

    ## focal_time
    # focal_time must be positive and smaller than the root age
    if (focal_time < 0)
    {
      stop(paste0("'focal_time' must be a positive number. It represents the time as a distance from the present."))
    }
    if (focal_time >= root_age)
    {
      stop(paste0("'focal_time' must be smaller than the root age of the phylogeny.\n",
                  "'focal_time' = ",focal_time,"; root age = ",root_age,"."))
    }

    ## update_rates & update_regimes
    # At least one of "update_rates" and "update_regimes" must be TRUE
    if (!update_rates & !update_rates)
    {
      stop(paste0("At least one of 'update_rates' and 'update_regimes' must be 'TRUE' for the function to update either rates or regimes found at 'focal_time'."))
    }

    ## update_plot & update_tree
    # If update_tree = FALSE, but update_plot = TRUE, show a warning claiming that the tree will be updated anyway, alongside plotting elements.
    if (update_plot & !update_tree)
    {
      warning(paste0("'update_tree' is set to 'FALSE', but was ignore as 'update_plot = TRUE' requires the tree to be updated for the plotting elements to be updated too."))
    }

    ## update_all_elements
    # If update_all_elements = TRUE & any other update_* = FALSE; show a warning claiming that update_* = FALSE will be ignore and all elements including rates/regimes/tree/plotting elements/BAMM elements will all be updated
    if (update_all_elements & any(!update_rates, !update_regimes, !update_tree, !update_plot))
    {
      warning(paste0("'update_all_elements' is set to 'TRUE'. All components (rates/regimes/tree/plotting elements/BAMM elements) were updated, even if other 'update_*' arguments are set to 'FALSE'."))
    }
  }

  ## Initiate new BAMM_object to update
  updated_BAMM_object <- BAMM_object

  # Add "phylo" class temporarily
  class(updated_BAMM_object) <- unique(c(class(updated_BAMM_object), "phylo"))

  ## Identify edges present at focal time

  # Edge, rootward_node, tipward_node, length (once cut)

  # Define level of tolerance used to round ages
  tol <- root_age * 10^-5
  closest_power <- round(log10(tol))
  closest_power <- min(closest_power, 0) # Use 0 as the minimal power

  # Get node ages per branch (no root edge)
  all_edges_df <- phytools::nodeHeights(updated_BAMM_object)
  # all_edges_df <- as.data.frame(round(root_age - all_edges_df, 5)) # Used to ensure ultrametricity of extant tips, but may be an issue for trees with very short time span
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, -1*closest_power))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- updated_BAMM_object$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- updated_BAMM_object$edge[which.min(updated_BAMM_object$edge[, 1] %in% updated_BAMM_object$edge[, 2]), 1]

  # Merge tip.label to the edge df

  # If tipward node is a tip, use tip.label
  all_edges_df$tip.label <- updated_BAMM_object$tip.label[all_edges_df$tipward_node_ID]
  # If tipward node is an internal node, use node ID
  all_edges_df$tip.label[is.na(all_edges_df$tip.label)] <- all_edges_df$tipward_node_ID[is.na(all_edges_df$tip.label)]

  # Identify edges present at focal time
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > focal_time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= focal_time
  all_edges_df$time_test <- all_edges_df$rootward_test & all_edges_df$tipward_test
  all_edges_df$length <- all_edges_df$rootward_node_age - focal_time

  # Initiate regime ID
  all_edges_df$regime_ID <- NA

  # Extract only edges that are present at the focal time
  # present_edges_df <- all_edges_df[all_edges_df$time_test, ]

  ## Loop per Posterior sample
  for (i in seq_along(updated_BAMM_object$eventData))
  {
    # i <- 1

    # Extract eventData records = Macroevolutionary regime parameters
    eventData_i <- updated_BAMM_object$eventData[[i]]

    # Compute updated regime age and length
    eventData_i$age <- root_age - eventData_i$time
    eventData_i$updated_length <- eventData_i$age - focal_time

    if (update_regimes || update_all_elements)
    {
      ## Identify edge ID per regimes
      ## Loop per regime
      for (j in 1:nrow(eventData_i))
      {
        # j <- 2

        tipward_node_ID_j <- eventData_i$node[j] # Nodes are tipward nodes ID of the branch where the regime starts

        # Get descendant tipward nodes of regime j
        regime_nodes_j <- phytools::getDescendants(tree = updated_BAMM_object, node = tipward_node_ID_j)
        # Remove tipward node of the starting edge
        regime_nodes_j <- setdiff(regime_nodes_j, tipward_node_ID_j)

        # Assign regime ID
        all_edges_df$regime_ID[all_edges_df$tipward_node_ID %in% regime_nodes_j] <- j

        # Deal with special case of the edge where the process starts
        # Should the edge where the process starts be included in the regime at the focal time?
        if (j != 1) # No need for the root process
        {
          # Identify the starting edge
          starting_edge_j <- as.numeric(all_edges_df$edge_ID[all_edges_df$tipward_node_ID == tipward_node_ID_j])

          # Get relative position of the regime shift
          relative_position_shift_j <- all_edges_df$rootward_node_age[starting_edge_j] - eventData_i$age[j]
          # Assign starting edge to process only if the regime shift happen before the time cut
          if (relative_position_shift_j < all_edges_df$length[starting_edge_j])
          {
            all_edges_df$regime_ID[starting_edge_j] <- j
          }
        }
      }

      # Filter regimes for tips that are present at the focal time
      tipStates_i <- all_edges_df$regime_ID[all_edges_df$time_test]

      # Name tip regimes with tip.labels/tipward_edge_ID
      if (keep_tip_labels)
      {
        names(tipStates_i) <- all_edges_df$tip.label[all_edges_df$time_test]
      } else {
        names(tipStates_i) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
      }

      # Store updated tipStates
      updated_BAMM_object$tipStates[[i]] <- tipStates_i
    }

    ## If needed, also update tipRates
    if (update_rates || update_all_elements)
    {
      eventData_i$tip_speciation_rates <- NA
      eventData_i$tip_extinction_rates <- NA

      ## Loop per regime
      for (j in 1:nrow(eventData_i))
      {
        # Compute new tip speciation rates based on regime parameters
        lambda_0_j <- eventData_i$lam1[j]
        alpha_j <- eventData_i$lam2[j]
        time_j <- eventData_i$updated_length[j]

        if (alpha_j <= 0) # If alpha <= 0 (decrease): lambda_t = lambda_0 * exp(alpha*t)
        {
          eventData_i$tip_speciation_rates[j] <- lambda_0_j * exp(alpha_j*time_j)
        } else { # If alpha > 0 (increase): lambda_t = lambda_0 * (2 - exp(-alpha*t))
          eventData_i$tip_speciation_rates[j] <- lambda_0_j * (2 - exp(-alpha_j*time_j))
        }

        # Compute new tip extinction rates based on regime parameters
        # All extinction rates are constant within regime in the current BAMM settings
        eventData_i$tip_extinction_rates[j] <- eventData_i$mu1[j]

        if (time_j < 0)
        {
          eventData_i$tip_speciation_rates[j] <- NA
          eventData_i$tip_extinction_rates[j] <- NA
        }
      }

      # Assign rates to edge according to regime ID
      all_edges_df$tipLambda <- NA
      all_edges_df$tipLambda <- eventData_i$tip_speciation_rates[match(x = all_edges_df$regime_ID, table = eventData_i$index)]
      all_edges_df$tipMu <- NA
      all_edges_df$tipMu <- eventData_i$tip_extinction_rates[match(x = all_edges_df$regime_ID, table = eventData_i$index)]

      # Filter regimes for tips that are present at the focal time
      tipLambda_i <- all_edges_df$tipLambda[all_edges_df$time_test]
      tipMu_i <- all_edges_df$tipMu[all_edges_df$time_test]

      # Name tip regimes with tip.labels/tipward_edge_ID
      if (keep_tip_labels)
      {
        names(tipLambda_i) <- all_edges_df$tip.label[all_edges_df$time_test]
        names(tipMu_i) <- all_edges_df$tip.label[all_edges_df$time_test]
      } else {
        names(tipLambda_i) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
        names(tipMu_i) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
      }

      # Store updated tipLambda & tipMu
      updated_BAMM_object$tipLambda[[i]] <- tipLambda_i
      updated_BAMM_object$tipMu[[i]] <- tipMu_i
    }

    ## Print progress
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Tip states/rates updated for BAMM posterior sample n\u00B0", i, "/", length(updated_BAMM_object$eventData),"\n"))
    }
  }

  ## Updates elements of the phylo objects if needed
  if (update_tree || update_plot || update_all_elements)
  {
    # Update phylo elements
    updated_BAMM_object <- cut_phylo_for_focal_time(tree = updated_BAMM_object, focal_time = focal_time, keep_tip_labels = keep_tip_labels)

    ## Extract and plot updated "phylo" object
    # updated_tree <- list(edge = updated_BAMM_object$edge, Nnode = updated_BAMM_object$Nnode, tip.label = updated_BAMM_object$tip.label, edge.length = updated_BAMM_object$edge.length)
    # class(updated_tree) <- "phylo"
    # plot(updated_tree)
    # nodelabels()
  }

  ## Updates elements needed to plot a "bammdata" object with plot_BAMM_rates()
  if (update_plot || update_all_elements)
  {
    ## $begin = Absolute time since root of edge/branch start
    updated_BAMM_object$begin <- BAMM_object$begin[updated_BAMM_object$edges_ID_df$initial_edge_ID] # Extract values only for remaining edges
    ## $end = Absolute time since root of edge/branch end
    updated_BAMM_object$end <- BAMM_object$end[updated_BAMM_object$edges_ID_df$initial_edge_ID] # Extract values only for remaining edges
    updated_BAMM_object$end <- sapply(X = updated_BAMM_object$end, FUN = function (x) { min(x, root_age - focal_time) }) # Adjust end distance to cut-off at focal_time

    ## $eventVectors = List of integer vectors of regime membership per branches in each posterior configuration
    updated_BAMM_object$eventVectors <- lapply(X = BAMM_object$eventVectors, FUN = function (x) { x[updated_BAMM_object$edges_ID_df$initial_edge_ID] } )

    ## $eventBranchSegs = Same but for segments with matrix including tipward node ID (NOT the edge ID) and begin/end ages of the segments
    # An edge with a shift is split in multiple segments.
    # Number of rows = number of segments = nb of edges + nb of shifts (each shift adds a segment to an edge)
    # Filtered using focal_time to keep only segments that are not younger than the focal_time

    # Loop per BAMM posterior samples
    for (i in seq_along(updated_BAMM_object$eventData))
    {
      # i <- 1

      # Extract matrix of branch segments
      eventBranchSegs_i <- updated_BAMM_object$eventBranchSegs[[i]]
      # Remove segments that are younger than focal_time
      updated_eventBranchSegs_i <- eventBranchSegs_i[(root_age - eventBranchSegs_i[,2] > focal_time), ]

      # Update tipward nodes ID
      updated_eventBranchSegs_i[ ,1] <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(updated_eventBranchSegs_i[ ,1], updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # Reorder following tipward nodes ID, then older segments > younger segments
      updated_eventBranchSegs_i <- updated_eventBranchSegs_i[order(updated_eventBranchSegs_i[ ,1], updated_eventBranchSegs_i[ ,2]), ]

      # Store updated matrix of branch segments
      updated_BAMM_object$eventBranchSegs[[i]] <- updated_eventBranchSegs_i
    }

    # # Updated eventBranchSegs matrix should have a number of rows/segments = number of edges + number of shifts
    # updated_nb_regimes <- unlist(lapply(X = updated_BAMM_object$eventBranchSegs, FUN = function (x) { length(unique(x[, 4])) }))
    # updated_nb_segments <- unlist(lapply(X = updated_BAMM_object$eventBranchSegs, FUN = nrow))
    # all(updated_nb_segments == (nrow(updated_BAMM_object$edge) + (updated_nb_regimes - 1)))

    ## Create $dtRates and update it such as it contains only rates for segments older than the focal_time
    # Can be used to keep consistency with estimated rates and color scheme used in deepSTRAPP::plot_BAMM_rates()

    ## New version. Works but not sure why...
    updated_dtrates <- BAMMtools::dtRates(updated_BAMM_object, tau = 0.01, tmat = TRUE)$dtrates
    updated_BAMM_object$dtrates <- updated_dtrates

    # ## Former 'manual' version that ensure segments to be similar, but creates artifacts on the terminal segments (not sure why...)
    # # Get initial dtrates
    # dtrates_t0 <- BAMMtools::dtRates(BAMM_object, tau = 0.01, tmat = TRUE)$dtrates
    #
    # # Find segments to remove segments that are younger than the focal_time
    # dtrates_segments_to_remove <- dtrates_t0$tmat[ , 2] >= (root_age - focal_time)
    #
    # # Update dtrates to remove segments that are older than the focal_time
    # updated_dtrates <- dtrates_t0
    # updated_dtrates$rates <- lapply(X = updated_dtrates$rates, FUN = function (x) { y <-  x[!dtrates_segments_to_remove]} )
    # updated_dtrates$tmat <- updated_dtrates$tmat[!dtrates_segments_to_remove,]
    #
    # # Update dtrates$tmat to set distance to root of the terminal segments according to the focal_time
    # updated_dtrates$tmat[(updated_dtrates$tmat[ ,3] > (root_age - focal_time)), 3] <- (root_age - focal_time)
    #
    # # Update tipward nodes ID in dtrates
    # updated_dtrates$tmat[ ,1] <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(updated_dtrates$tmat[ ,1], updated_BAMM_object$nodes_ID_df$initial_node_ID)]
    # # Update dimnames
    # attr(updated_dtrates$tmat, which = "dimnames")[[1]] <- as.character(1:nrow(updated_dtrates$tmat))
    # # Update tau as the fraction of the total tree length represented by each segment
    # new_depth <- (root_age - focal_time)
    # depth_ratio <- new_depth/root_age
    # updated_dtrates$tau <- updated_dtrates$tau/depth_ratio
    #
    # # Store updated $dtrates
    # updated_BAMM_object$dtrates <- updated_dtrates

    ## Save initial_colorbreaks to use as colorbreaks in order to match color gradients from the initial full phylogeny
    initial_plot_speciation <- BAMMtools::plot.bammdata(BAMM_object, legend = TRUE, show = FALSE, spex = "s")
    initial_plot_extinction <- BAMMtools::plot.bammdata(BAMM_object, legend = TRUE, show = FALSE, spex = "e")
    initial_plot_net_div <- BAMMtools::plot.bammdata(BAMM_object, legend = TRUE, show = FALSE, spex = "netdiv")
    # updated_BAMM_object$initial_colorbreaks_range <- range(initial_plot$colorbreaks)
    updated_BAMM_object$initial_colorbreaks <- list(speciation = initial_plot_speciation$colorbreaks,
                                                    extinction = initial_plot_extinction$colorbreaks,
                                                    net_diversification = initial_plot_net_div$colorbreaks)

    # ## Updated BAMM_object can be plotted with deepSTRAPP::plot_BAMM_rates()
    # plot_BAMM_rates(BAMM_object, legend = TRUE)
    # plot_BAMM_rates(updated_BAMM_object, legend = TRUE)
    # plot_BAMM_rates(updated_BAMM_object, legend = TRUE, colorbreaks = updated_BAMM_object$initial_colorbreaks$net_diversification)

    # Update the Marginal Shift Probability tree if present (used to plot regime shifts)
    if ("MSP_tree" %in% names(updated_BAMM_object))
    {
      # Use the topology of the new updated tree
      updated_BAMM_object$MSP_tree$edge <- updated_BAMM_object$edge
      updated_BAMM_object$MSP_tree$Nnode <- updated_BAMM_object$Nnode
      updated_BAMM_object$MSP_tree$tip.label <- updated_BAMM_object$tip.label
      # Extract edge length (Marginal shift posterior probabilities) for the remaining edges
      updated_BAMM_object$MSP_tree$edge.length <- updated_BAMM_object$MSP_tree$edge.length[updated_BAMM_object$edges_ID_df$initial_edge_ID]
    }

    ## Update the BAMM_object for the Maximum A Posteriori probability (MAP) configuration if present (used to plot regime shifts)
    if ("MAP_BAMM_object" %in% names(updated_BAMM_object))
    {
      ## Add "phylo" class to be compatible with phytools::getDescendants()
      class(updated_BAMM_object$MAP_BAMM_object) <- unique(c(class(updated_BAMM_object$MAP_BAMM_object), "phylo"))

      ## Use the $MAP_BAMM_object following updates from the main BAMM_object
      updated_BAMM_object$MAP_BAMM_object$edge <- updated_BAMM_object$edge
      updated_BAMM_object$MAP_BAMM_object$Nnode <- updated_BAMM_object$Nnode
      updated_BAMM_object$MAP_BAMM_object$tip.label <- updated_BAMM_object$tip.label
      updated_BAMM_object$MAP_BAMM_object$edge.length <- updated_BAMM_object$edge.length
      updated_BAMM_object$MAP_BAMM_object$begin <- updated_BAMM_object$begin
      updated_BAMM_object$MAP_BAMM_object$end <- updated_BAMM_object$end

      ## Update information according to 'focal_type'

      ## $eventVectors
      ## $eventVectors = List of integer vectors of regime membership per branches in each posterior configuration
      updated_BAMM_object$MAP_BAMM_object$eventVectors[[1]] <- updated_BAMM_object$MAP_BAMM_object$eventVectors[[1]][updated_BAMM_object$edges_ID_df$initial_edge_ID]

      ## $tipStates
      # Extract eventData records = Macroevolutionary regime parameters
      MAP_eventData <- updated_BAMM_object$MAP_BAMM_object$eventData[[1]]
      # Compute updated regime age and length
      MAP_eventData$age <- root_age - MAP_eventData$time
      MAP_eventData$updated_length <- MAP_eventData$age - focal_time
      # Identify edge ID per regimes
      # Loop per regime
      for (j in 1:nrow(MAP_eventData))
      {
        # j <- 2

        tipward_node_ID_j <- MAP_eventData$node[j] # Nodes are tipward nodes ID of the branch where the regime starts

        # Get descendant tipward nodes of regime j
        regime_nodes_j <- phytools::getDescendants(tree = updated_BAMM_object$MAP_BAMM_object, node = tipward_node_ID_j)

        # Assign regime ID
        all_edges_df$regime_ID[all_edges_df$tipward_node_ID %in% regime_nodes_j] <- j

        # Deal with special case of the edge where the process starts
        # Should the edge where the process starts be included in the regime at the focal time?
        if (j != 1) # No need for the root process
        {
          # Identify the starting edge
          starting_edge_j <- as.numeric(all_edges_df$edge_ID[all_edges_df$tipward_node_ID == tipward_node_ID_j])

          # Get relative position of the regime shift
          relative_position_shift_j <- all_edges_df$rootward_node_age[starting_edge_j] - MAP_eventData$age[j]
          # Assign starting edge to process only if the regime shift happen before the time cut
          if (relative_position_shift_j < all_edges_df$length[starting_edge_j])
          {
            all_edges_df$regime_ID[starting_edge_j] <- j
          }
        }
      }
      # Filter regimes for tips that are present at the focal time
      MAP_tipStates <- all_edges_df$regime_ID[all_edges_df$time_test]
      # Name tip regimes with tip.labels/tipward_edge_ID
      if (keep_tip_labels)
      {
        names(MAP_tipStates) <- all_edges_df$tip.label[all_edges_df$time_test]
      } else {
        names(MAP_tipStates) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
      }
      # Store updated tipStates
      updated_BAMM_object$MAP_BAMM_object$tipStates[[1]] <- MAP_tipStates

      ## $tipLambda & $tipMu
      MAP_eventData$tip_speciation_rates <- NA
      MAP_eventData$tip_extinction_rates <- NA
      # Loop per regime
      for (j in 1:nrow(MAP_eventData))
      {
        # Compute new tip speciation rates based on regime parameters
        lambda_0_j <- MAP_eventData$lam1[j]
        alpha_j <- MAP_eventData$lam2[j]
        time_j <- MAP_eventData$updated_length[j]

        if (alpha_j <= 0) # If alpha <= 0 (decrease): lambda_t = lambda_0 * exp(alpha*t)
        {
          MAP_eventData$tip_speciation_rates[j] <- lambda_0_j * exp(alpha_j*time_j)
        } else { # If alpha > 0 (increase): lambda_t = lambda_0 * (2 - exp(-alpha*t))
          MAP_eventData$tip_speciation_rates[j] <- lambda_0_j * (2 - exp(-alpha_j*time_j))
        }

        # Compute new tip extinction rates based on regime parameters
        # All extinction rates are constant within regime in the current BAMM settings
        MAP_eventData$tip_extinction_rates[j] <- MAP_eventData$mu1[j]
        if (time_j < 0)
        {
          MAP_eventData$tip_speciation_rates[j] <- NA
          MAP_eventData$tip_extinction_rates[j] <- NA
        }
      }

      # Assign rates to edge according to regime ID
      all_edges_df$tipLambda <- NA
      all_edges_df$tipLambda <- MAP_eventData$tip_speciation_rates[match(x = all_edges_df$regime_ID, table = MAP_eventData$index)]
      all_edges_df$tipMu <- NA
      all_edges_df$tipMu <- MAP_eventData$tip_extinction_rates[match(x = all_edges_df$regime_ID, table = MAP_eventData$index)]

      # Filter regimes for tips that are present at the focal time
      MAP_tipLambda <- all_edges_df$tipLambda[all_edges_df$time_test]
      MAP_tipMu <- all_edges_df$tipMu[all_edges_df$time_test]

      # Name tip regimes with tip.labels/tipward_edge_ID
      if (keep_tip_labels)
      {
        names(MAP_tipLambda) <- all_edges_df$tip.label[all_edges_df$time_test]
        names(MAP_tipMu) <- all_edges_df$tip.label[all_edges_df$time_test]
      } else {
        names(MAP_tipLambda) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
        names(MAP_tipMu) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
      }

      # Store updated tipLambda & tipMu
      updated_BAMM_object$MAP_BAMM_object$tipLambda[[1]] <- MAP_tipLambda
      updated_BAMM_object$MAP_BAMM_object$tipMu[[1]] <- MAP_tipMu

      ## $eventData # Dataframe recording shift events and macroevolutionary regimes in the focal posterior configuration. 1st line = Background root regime
      # Filter to keep only events that happened before focal_time
      MAP_eventData <- updated_BAMM_object$MAP_BAMM_object$eventData[[1]]
      MAP_eventData <- MAP_eventData[((root_age - MAP_eventData$time) > focal_time), ]
      # Update tipward nodes ID
      MAP_eventData$node <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(MAP_eventData$node, updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # Store updated df of macroevolutionary regimes
      updated_BAMM_object$MAP_BAMM_object$eventData[[1]] <- MAP_eventData

      ## $eventBranchSegs
      # Extract matrix of branch segments
      MAP_eventBranchSegs <- updated_BAMM_object$MAP_BAMM_object$eventBranchSegs[[1]]
      # Remove segments that are younger than focal_time
      MAP_eventBranchSegs <- MAP_eventBranchSegs[(root_age - MAP_eventBranchSegs[,2] > focal_time), ]

      # Update tipward nodes ID
      MAP_eventBranchSegs[ ,1] <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(MAP_eventBranchSegs[ ,1], updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # Reorder following tipward nodes ID, then older segments > younger segments
      MAP_eventBranchSegs <- MAP_eventBranchSegs[order(MAP_eventBranchSegs[ ,1], MAP_eventBranchSegs[ ,2]), ]

      # Store updated matrix of branch segments
      updated_BAMM_object$MAP_BAMM_object$eventBranchSegs[[1]] <- MAP_eventBranchSegs

      ## Use the $MAP_BAMM_object following updates from the main BAMM_object
      updated_BAMM_object$MAP_BAMM_object$type <- updated_BAMM_object$type

      ## Create $dtRates and update it such as it contains only rates for segments older than the focal_time
      # Can be used to keep consistency with estimated rates and color scheme used in deepSTRAPP::plot_BAMM_rates()

      ## New version. Works but not sure why...
      ## New version. Does not work on updated tree! Rate estimates are crap. Issue with reordering of egdes?
      MAP_dtrates <- BAMMtools::dtRates(updated_BAMM_object$MAP_BAMM_object, tau = 0.01, tmat = TRUE)$dtrates
      updated_BAMM_object$MAP_BAMM_object$dtrates <- MAP_dtrates

      # # ## Former 'manual' version that ensure segments to be similar, but creates artifacts on the terminal segments (not sure why...)
      # # # Get initial dtrates
      # # dtrates_t0 <- BAMMtools::dtRates(BAMM_object$MAP_BAMM_object, tau = 0.01, tmat = TRUE)$dtrates
      # #
      # # # Find segments to remove segments that are younger than the focal_time
      # # dtrates_segments_to_remove <- dtrates_t0$tmat[ , 2] >= (root_age - focal_time)
      # #
      # # # Update dtrates to remove segments that are older than the focal_time
      # # MAP_dtrates <- dtrates_t0
      # # MAP_dtrates$rates <- lapply(X = MAP_dtrates$rates, FUN = function (x) { y <-  x[!dtrates_segments_to_remove]} )
      # # MAP_dtrates$tmat <- MAP_dtrates$tmat[!dtrates_segments_to_remove,]
      # #
      # # # Update tipward nodes ID in dtrates
      # # MAP_dtrates$tmat[ ,1] <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(MAP_dtrates$tmat[ ,1], updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # # # Update dimnames
      # # attr(MAP_dtrates$tmat, which = "dimnames")[[1]] <- as.character(1:nrow(MAP_dtrates$tmat))
      # # # Update tau as the fraction of the total tree length represented by each segment
      # # new_depth <- (root_age - focal_time)
      # # depth_ratio <- new_depth/root_age
      # # MAP_dtrates$tau <- MAP_dtrates$tau/depth_ratio
      #
      # # Store updated $dtrates
      # updated_BAMM_object$MAP_BAMM_object$dtrates <- MAP_dtrates

      ## Save initial_colorbreaks to use as colorbreaks in order to match color gradients from the initial full phylogeny
      initial_plot_speciation <- BAMMtools::plot.bammdata(BAMM_object$MAP_BAMM_object, legend = TRUE, show = FALSE, spex = "s")
      initial_plot_extinction <- BAMMtools::plot.bammdata(BAMM_object$MAP_BAMM_object, legend = TRUE, show = FALSE, spex = "e")
      initial_plot_net_div <- BAMMtools::plot.bammdata(BAMM_object$MAP_BAMM_object, legend = TRUE, show = FALSE, spex = "netdiv")
      # updated_BAMM_object$MAP_BAMM_object$initial_colorbreaks_range <- range(initial_plot$colorbreaks)
      updated_BAMM_object$MAP_BAMM_object$initial_colorbreaks <- list(speciation = initial_plot_speciation$colorbreaks,
                                                                      extinction = initial_plot_extinction$colorbreaks,
                                                                      net_diversification = initial_plot_net_div$colorbreaks)

      # Remove temporary "phylo" class
      class(updated_BAMM_object$MAP_BAMM_object) <- setdiff(class(updated_BAMM_object$MAP_BAMM_object), "phylo")
    }

    ## Update the BAMM_object for the Maximum Shift Credibility (MSC) configuration if present (used to plot regime shifts)
    if ("MSC_BAMM_object" %in% names(updated_BAMM_object))
    {
      ## Add "phylo" class to be compatible with phytools::getDescendants()
      class(updated_BAMM_object$MSC_BAMM_object) <- unique(c(class(updated_BAMM_object$MSC_BAMM_object), "phylo"))

      ## Use the $MSC_BAMM_object following updates from the main BAMM_object
      updated_BAMM_object$MSC_BAMM_object$edge <- updated_BAMM_object$edge
      updated_BAMM_object$MSC_BAMM_object$Nnode <- updated_BAMM_object$Nnode
      updated_BAMM_object$MSC_BAMM_object$tip.label <- updated_BAMM_object$tip.label
      updated_BAMM_object$MSC_BAMM_object$edge.length <- updated_BAMM_object$edge.length
      updated_BAMM_object$MSC_BAMM_object$begin <- updated_BAMM_object$begin
      updated_BAMM_object$MSC_BAMM_object$end <- updated_BAMM_object$end

      ## Update information according to 'focal_type'

      ## $eventVectors
      ## $eventVectors = List of integer vectors of regime membership per branches in each posterior configuration
      updated_BAMM_object$MSC_BAMM_object$eventVectors[[1]] <- updated_BAMM_object$MSC_BAMM_object$eventVectors[[1]][updated_BAMM_object$edges_ID_df$initial_edge_ID]

      ## $tipStates
      # Extract eventData records = Macroevolutionary regime parameters
      MSC_eventData <- updated_BAMM_object$MSC_BAMM_object$eventData[[1]]
      # Compute updated regime age and length
      MSC_eventData$age <- root_age - MSC_eventData$time
      MSC_eventData$updated_length <- MSC_eventData$age - focal_time
      # Identify edge ID per regimes
      # Loop per regime
      for (j in 1:nrow(MSC_eventData))
      {
        # j <- 2

        tipward_node_ID_j <- MSC_eventData$node[j] # Nodes are tipward nodes ID of the branch where the regime starts

        # Get descendant tipward nodes of regime j
        regime_nodes_j <- phytools::getDescendants(tree = updated_BAMM_object$MSC_BAMM_object, node = tipward_node_ID_j)

        # Assign regime ID
        all_edges_df$regime_ID[all_edges_df$tipward_node_ID %in% regime_nodes_j] <- j

        # Deal with special case of the edge where the process starts
        # Should the edge where the process starts be included in the regime at the focal time?
        if (j != 1) # No need for the root process
        {
          # Identify the starting edge
          starting_edge_j <- as.numeric(all_edges_df$edge_ID[all_edges_df$tipward_node_ID == tipward_node_ID_j])

          # Get relative position of the regime shift
          relative_position_shift_j <- all_edges_df$rootward_node_age[starting_edge_j] - MSC_eventData$age[j]
          # Assign starting edge to process only if the regime shift happen before the time cut
          if (relative_position_shift_j < all_edges_df$length[starting_edge_j])
          {
            all_edges_df$regime_ID[starting_edge_j] <- j
          }
        }
      }
      # Filter regimes for tips that are present at the focal time
      MSC_tipStates <- all_edges_df$regime_ID[all_edges_df$time_test]
      # Name tip regimes with tip.labels/tipward_edge_ID
      if (keep_tip_labels)
      {
        names(MSC_tipStates) <- all_edges_df$tip.label[all_edges_df$time_test]
      } else {
        names(MSC_tipStates) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
      }
      # Store updated tipStates
      updated_BAMM_object$MSC_BAMM_object$tipStates[[1]] <- MSC_tipStates

      ## $tipLambda & $tipMu
      MSC_eventData$tip_speciation_rates <- NA
      MSC_eventData$tip_extinction_rates <- NA
      # Loop per regime
      for (j in 1:nrow(MSC_eventData))
      {
        # Compute new tip speciation rates based on regime parameters
        lambda_0_j <- MSC_eventData$lam1[j]
        alpha_j <- MSC_eventData$lam2[j]
        time_j <- MSC_eventData$updated_length[j]

        if (alpha_j <= 0) # If alpha <= 0 (decrease): lambda_t = lambda_0 * exp(alpha*t)
        {
          MSC_eventData$tip_speciation_rates[j] <- lambda_0_j * exp(alpha_j*time_j)
        } else { # If alpha > 0 (increase): lambda_t = lambda_0 * (2 - exp(-alpha*t))
          MSC_eventData$tip_speciation_rates[j] <- lambda_0_j * (2 - exp(-alpha_j*time_j))
        }

        # Compute new tip extinction rates based on regime parameters
        # All extinction rates are constant within regime in the current BAMM settings
        MSC_eventData$tip_extinction_rates[j] <- MSC_eventData$mu1[j]
        if (time_j < 0)
        {
          MSC_eventData$tip_speciation_rates[j] <- NA
          MSC_eventData$tip_extinction_rates[j] <- NA
        }
      }

      # Assign rates to edge according to regime ID
      all_edges_df$tipLambda <- NA
      all_edges_df$tipLambda <- MSC_eventData$tip_speciation_rates[match(x = all_edges_df$regime_ID, table = MSC_eventData$index)]
      all_edges_df$tipMu <- NA
      all_edges_df$tipMu <- MSC_eventData$tip_extinction_rates[match(x = all_edges_df$regime_ID, table = MSC_eventData$index)]

      # Filter regimes for tips that are present at the focal time
      MSC_tipLambda <- all_edges_df$tipLambda[all_edges_df$time_test]
      MSC_tipMu <- all_edges_df$tipMu[all_edges_df$time_test]

      # Name tip regimes with tip.labels/tipward_edge_ID
      if (keep_tip_labels)
      {
        names(MSC_tipLambda) <- all_edges_df$tip.label[all_edges_df$time_test]
        names(MSC_tipMu) <- all_edges_df$tip.label[all_edges_df$time_test]
      } else {
        names(MSC_tipLambda) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
        names(MSC_tipMu) <- all_edges_df$tipward_node_ID[all_edges_df$time_test]
      }

      # Store updated tipLambda & tipMu
      updated_BAMM_object$MSC_BAMM_object$tipLambda[[1]] <- MSC_tipLambda
      updated_BAMM_object$MSC_BAMM_object$tipMu[[1]] <- MSC_tipMu

      ## $eventData # Dataframe recording shift events and macroevolutionary regimes in the focal posterior configuration. 1st line = Background root regime
      # Filter to keep only events that happened before focal_time
      MSC_eventData <- updated_BAMM_object$MSC_BAMM_object$eventData[[1]]
      MSC_eventData <- MSC_eventData[((root_age - MSC_eventData$time) > focal_time), ]
      # Update tipward nodes ID
      MSC_eventData$node <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(MSC_eventData$node, updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # Store updated df of macroevolutionary regimes
      updated_BAMM_object$MSC_BAMM_object$eventData[[1]] <- MSC_eventData

      ## $eventBranchSegs
      # Extract matrix of branch segments
      MSC_eventBranchSegs <- updated_BAMM_object$MSC_BAMM_object$eventBranchSegs[[1]]
      # Remove segments that are younger than focal_time
      MSC_eventBranchSegs <- MSC_eventBranchSegs[(root_age - MSC_eventBranchSegs[,2] > focal_time), ]

      # Update tipward nodes ID
      MSC_eventBranchSegs[ ,1] <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(MSC_eventBranchSegs[ ,1], updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # Reorder following tipward nodes ID, then older segments > younger segments
      MSC_eventBranchSegs <- MSC_eventBranchSegs[order(MSC_eventBranchSegs[ ,1], MSC_eventBranchSegs[ ,2]), ]

      # Store updated matrix of branch segments
      updated_BAMM_object$MSC_BAMM_object$eventBranchSegs[[1]] <- MSC_eventBranchSegs

      ## Use the $MSC_BAMM_object following updates from the main BAMM_object
      updated_BAMM_object$MSC_BAMM_object$type <- updated_BAMM_object$type

      ## Create $dtRates and update it such as it contains only rates for segments older than the focal_time
      # Can be used to keep consistency with estimated rates and color scheme used in deepSTRAPP::plot_BAMM_rates()

      ## New version. Works but not sure why...
      ## New version. Does not work on updated tree! Rate estimates are crap. Issue with reordering of egdes?
      MSC_dtrates <- BAMMtools::dtRates(updated_BAMM_object$MSC_BAMM_object, tau = 0.01, tmat = TRUE)$dtrates
      updated_BAMM_object$MSC_BAMM_object$dtrates <- MSC_dtrates

      # # ## Former 'manual' version that ensure segments to be similar, but creates artifacts on the terminal segments (not sure why...)
      # # # Get initial dtrates
      # # dtrates_t0 <- BAMMtools::dtRates(BAMM_object$MSC_BAMM_object, tau = 0.01, tmat = TRUE)$dtrates
      # #
      # # # Find segments to remove segments that are younger than the focal_time
      # # dtrates_segments_to_remove <- dtrates_t0$tmat[ , 2] >= (root_age - focal_time)
      # #
      # # # Update dtrates to remove segments that are older than the focal_time
      # # MSC_dtrates <- dtrates_t0
      # # MSC_dtrates$rates <- lapply(X = MSC_dtrates$rates, FUN = function (x) { y <-  x[!dtrates_segments_to_remove]} )
      # # MSC_dtrates$tmat <- MSC_dtrates$tmat[!dtrates_segments_to_remove,]
      # #
      # # # Update tipward nodes ID in dtrates
      # # MSC_dtrates$tmat[ ,1] <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(MSC_dtrates$tmat[ ,1], updated_BAMM_object$nodes_ID_df$initial_node_ID)]
      # # # Update dimnames
      # # attr(MSC_dtrates$tmat, which = "dimnames")[[1]] <- as.character(1:nrow(MSC_dtrates$tmat))
      # # # Update tau as the fraction of the total tree length represented by each segment
      # # new_depth <- (root_age - focal_time)
      # # depth_ratio <- new_depth/root_age
      # # MSC_dtrates$tau <- MSC_dtrates$tau/depth_ratio
      #
      # # Store updated $dtrates
      # updated_BAMM_object$MSC_BAMM_object$dtrates <- MSC_dtrates

      ## Save initial_colorbreaks to use as colorbreaks in order to match color gradients from the initial full phylogeny
      initial_plot_speciation <- BAMMtools::plot.bammdata(BAMM_object$MSC_BAMM_object, legend = TRUE, show = FALSE, spex = "s")
      initial_plot_extinction <- BAMMtools::plot.bammdata(BAMM_object$MSC_BAMM_object, legend = TRUE, show = FALSE, spex = "e")
      initial_plot_net_div <- BAMMtools::plot.bammdata(BAMM_object$MSC_BAMM_object, legend = TRUE, show = FALSE, spex = "netdiv")
      # updated_BAMM_object$MSC_BAMM_object$initial_colorbreaks_range <- range(initial_plot$colorbreaks)
      updated_BAMM_object$MSC_BAMM_object$initial_colorbreaks <- list(speciation = initial_plot_speciation$colorbreaks,
                                                                      extinction = initial_plot_extinction$colorbreaks,
                                                                      net_diversification = initial_plot_net_div$colorbreaks)


      # Remove temporary "phylo" class
      class(updated_BAMM_object$MSC_BAMM_object) <- setdiff(class(updated_BAMM_object$MSC_BAMM_object), "phylo")
    }

  }

  if (update_all_elements)
  {
    ## Info for tree exploration
    # $downseq # Order of node visits when using a pre-order tree traversal
    # $lastvisit # ID of the last node visited when starting from the node in the corresponding position in downseq.
    updated_BAMM_object <- getRecursiveSequence(updated_BAMM_object)

    ## $numberEvents # Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts
    # Extract number of regimes detected across the updated segments
    updated_BAMM_object$numberEvents <- unlist(lapply(X = updated_BAMM_object$eventBranchSegs, FUN = function (x) { length(unique(x[, 4])) }))

    ## $eventData # Dataframe recording shift events and macroevolutionary regimes in the focal posterior configuration. 1st line = Background root regime
    # Loop per BAMM posterior samples
    for (i in seq_along(updated_BAMM_object$eventData))
    {
      # i <- 1

      # Extract df of macroevolutionary regimes
      eventData_i <- updated_BAMM_object$eventData[[i]]
      # Filter to keep only events that happened before focal_time
      updated_eventData_i <- eventData_i[((root_age - eventData_i$time) > focal_time), ]
      # Update tipward nodes ID
      updated_eventData_i$node <- updated_BAMM_object$nodes_ID_df$new_node_ID[match(updated_eventData_i$node, updated_BAMM_object$nodes_ID_df$initial_node_ID)]

      # Store updated df of macroevolutionary regimes
      updated_BAMM_object$eventData[[i]] <- updated_eventData_i
    }

    ## Mean tip rates across all posterior configurations

    # Bind all tipLambda in a df
    tipLambda_df <- data.frame(do.call(what = rbind, args = updated_BAMM_object$tipLambda))
    # Compute mean
    updated_meanTipLambda <- apply(X = tipLambda_df, MARGIN = 2, FUN = mean)
    # Provides tip.labels/tipward_node_ID as names. Extract only current tip names (not older fossils !)
    if (keep_tip_labels)
    {
      names(updated_meanTipLambda) <- updated_BAMM_object$tip.label[all_edges_df$time_test[all_edges_df$tip.label %in% updated_BAMM_object$tip.label]]
    } else {
      names(updated_meanTipLambda) <- updated_BAMM_object$tip.label[all_edges_df$time_test[all_edges_df$tipward_node_ID %in% updated_BAMM_object$tip.label]]
    }
    # Store updated tipLambda
    updated_BAMM_object$meanTipLambda <- updated_meanTipLambda

    # Bind all tipMu in a df
    tipMu_df <- data.frame(do.call(what = rbind, args = updated_BAMM_object$tipMu))
    # Compute mean
    updated_meanTipMu <- apply(X = tipMu_df, MARGIN = 2, FUN = mean)
    # Provides tip.labels/tipward_node_ID as names
    if (keep_tip_labels)
    {
      names(updated_meanTipMu) <- updated_BAMM_object$tip.label[all_edges_df$time_test[all_edges_df$tip.label %in% updated_BAMM_object$tip.label]]
    } else {
      names(updated_meanTipMu) <- updated_BAMM_object$tip.label[all_edges_df$time_test[all_edges_df$tipward_node_ID %in% updated_BAMM_object$tip.label]]
    }
    # Store updated tipMu
    updated_BAMM_object$meanTipMu <- updated_meanTipMu

    ## Update the BAMM_object for the Maximum A Posteriori probability (MAP) configuration if present (used to plot regime shifts)
    if ("MAP_BAMM_object" %in% names(updated_BAMM_object))
    {
      # Update $downseq & $lastvisit from the main BAMM_object
      updated_BAMM_object$MAP_BAMM_object$downseq <- updated_BAMM_object$downseq
      updated_BAMM_object$MAP_BAMM_object$lastvisit <- updated_BAMM_object$lastvisit

      # Update $eventVectors by extracting information of remaning branches only
      updated_BAMM_object$MAP_BAMM_object$eventVectors[[1]] <- updated_BAMM_object$MAP_BAMM_object$eventVectors[[1]][updated_BAMM_object$edges_ID_df$initial_edge_ID]

      ## Update $numberEvents from $eventBranchSegs
      updated_BAMM_object$numberEvents <- length(unique(updated_BAMM_object$MAP_BAMM_object$eventBranchSegs[[1]][, 4]))

      ## Update $meanTipLambda and $meanTipMu as in $TipLambda and $TipMu
      updated_BAMM_object$MAP_BAMM_object$meanTipLambda <- updated_BAMM_object$MAP_BAMM_object$tipLambda[[1]]
      updated_BAMM_object$MAP_BAMM_object$meanTipMu <- updated_BAMM_object$MAP_BAMM_object$tipMu[[1]]

      ## Reorder elements to fit order in the main BAMM_object
      if ("node.label" %in% names(updated_BAMM_object$MAP_BAMM_object))
      {
        updated_BAMM_object$MAP_BAMM_object <- updated_BAMM_object$MAP_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length", "node.label",
            "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
            "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
            "meanTipLambda", "meanTipMu", "type", "dtrates", "initial_colorbreaks")]
      } else {
        updated_BAMM_object$MAP_BAMM_object <- updated_BAMM_object$MAP_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length",
            "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
            "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
            "meanTipLambda", "meanTipMu", "type", "dtrates", "initial_colorbreaks")]
      }
      class(updated_BAMM_object$MAP_BAMM_object) <- "bammdata"
      attr(x = updated_BAMM_object$MAP_BAMM_object, which = "order") <- "cladewise"
    }

    ## Update the BAMM_object for the Maximum Shift Credibility (MSC) configuration if present (used to plot regime shifts)
    if ("MSC_BAMM_object" %in% names(updated_BAMM_object))
    {
      # Update $downseq & $lastvisit from the main BAMM_object
      updated_BAMM_object$MSC_BAMM_object$downseq <- updated_BAMM_object$downseq
      updated_BAMM_object$MSC_BAMM_object$lastvisit <- updated_BAMM_object$lastvisit

      # Update $eventVectors by extracting information of remaning branches only
      updated_BAMM_object$MSC_BAMM_object$eventVectors[[1]] <- updated_BAMM_object$MSC_BAMM_object$eventVectors[[1]][updated_BAMM_object$edges_ID_df$initial_edge_ID]

      ## Update $numberEvents from $eventBranchSegs
      updated_BAMM_object$numberEvents <- length(unique(updated_BAMM_object$MSC_BAMM_object$eventBranchSegs[[1]][, 4]))

      ## Update $meanTipLambda and $meanTipMu as in $TipLambda and $TipMu
      updated_BAMM_object$MSC_BAMM_object$meanTipLambda <- updated_BAMM_object$MSC_BAMM_object$tipLambda[[1]]
      updated_BAMM_object$MSC_BAMM_object$meanTipMu <- updated_BAMM_object$MSC_BAMM_object$tipMu[[1]]

      ## Reorder elements to fit order in the main BAMM_object
      if ("node.label" %in% names(updated_BAMM_object$MSC_BAMM_object))
      {
        updated_BAMM_object$MSC_BAMM_object <- updated_BAMM_object$MSC_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length", "node.label",
                                                                                     "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                                                                     "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                                                                     "meanTipLambda", "meanTipMu", "type", "dtrates", "initial_colorbreaks")]
      } else {
        updated_BAMM_object$MSC_BAMM_object <- updated_BAMM_object$MSC_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length",
                                                                                     "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                                                                     "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                                                                     "meanTipLambda", "meanTipMu", "type", "dtrates", "initial_colorbreaks")]
      }
      class(updated_BAMM_object$MSC_BAMM_object) <- "bammdata"
      attr(x = updated_BAMM_object$MSC_BAMM_object, which = "order") <- "cladewise"
    }
  }

  # Inform focal time
  updated_BAMM_object$focal_time <- focal_time

  # Remove temporary "phylo" class
  class(updated_BAMM_object) <- setdiff(class(updated_BAMM_object), "phylo")

  # Export updated BAMM_object
  return(updated_BAMM_object)
}



## Helper function used to generate $mapped.edge from $edge and $maps

# Internal function copied from the R package BAMMtools
# Associated C scripts in /src/treetraverse.c
# Source: BAMMtools:::getRecursiveSequence()
# Authors: Dan Rabosky, Mike Grundler

# Declare the use of compiled C code in the package
# Add an import in NAMESPACE
#' @useDynLib deepSTRAPP

getRecursiveSequence <- function (phy)
{
  rootnd = as.integer(phy$Nnode + 2)
  anc = as.integer(phy$edge[, 1])
  desc = as.integer(phy$edge[, 2])
  ne = as.integer(dim(phy$edge)[1])
  L = .C("setrecursivesequence", anc, desc, rootnd, ne, integer(ne + 1), integer(ne + 1), PACKAGE = "deepSTRAPP")
  phy$downseq = as.integer(L[[5]])
  phy$lastvisit = as.integer(L[[6]])
  return(phy)
}
