#' @title Extracts continuous trait data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely trait values found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogeny (`contMap`) such as
#'   branches overlapping the `focal_time` are shorten to the `focal_time`, and
#'   the continuous trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
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
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_contMap Logical. Specify whether the mapped phylogeny (`contMap`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated contMap. Default is `TRUE`.
#'   Used only if `update_contMap = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.contMap
#' @importFrom ape nodelabels
#'
#' @details The mapped phylogeny (`contMap`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   If providing only the `contMap` trait values at tips and internal nodes will be extracted from
#'   the mapping of the `contMap` leading to a slight dependency with the actual tip data
#'   and estimated ancestral character values.
#'
#'   True ML estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'
#'   ----- Update the `contMap` -----
#'
#'   To obtain an updated `contMap` alongside the trait data, set `update_contMap = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The continuous trait mapping is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named numerical vector with ML trait values found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_contMap = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$contMap`.
#'   * `$contMap` An object of class that contains the updated `contMap` with  branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$contMap$tree` element.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMap_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Only extent taxa (Ultrametric tree) ----- #
#'
#' ## Prepare data
#'
#' # Load eel data from the R package phytools
#' # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract body size
#' eel_data <- setNames(eel.data$Max_TL_cm,
#'                      rownames(eel.data))
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' eel_ACE <- phytools::fastAnc(tree = eel.tree, x = eel_data)
#'
#' # Run a Stochastic Mapping based on a Brownian Motion model
#' # to interpolate values along branches and obtain a "contMap" object
#' eel_contMap <- phytools::contMap(eel.tree, x = eel_data,
#'                                  res = 100, # Number of time steps
#'                                  plot = FALSE)
#'
#' # Set focal time to 50 Mya
#' focal_time <- 50
#'
#' ## Extract trait data and update contMap for the given focal_time
#'
#' # Extract from the contMap (values are not exact ML estimates)
#' eel_test <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = eel_contMap,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#' # Extract from tip data and ML estimates of ancestral characters (values are true ML estimates)
#' eel_test <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = eel_contMap,
#'    ace = eel_ACE, tip_data = eel_data,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#'
#' ## Visualize outputs
#'
#' # Print trait data
#' eel_test$trait_data
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(eel_contMap)
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(eel_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' plot(eel_test$contMap)
#' ape::nodelabels(text = eel_test$contMap$tree$initial_nodes_ID)
#'
#' # ----- Example 2: Include fossils (Non-ultrametric tree) ----- #
#'
#' ## Test with non-ultrametric trees like mammals in motmot
#'
#' ## Prepare data
#'
#' # Load mammals phylogeny and data from the R package motmot
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#'
#' library(motmot)
#' data("mammals")
#' force(mammals)
#'
#' mammals_tree <- mammals$mammal.phy
#' mammals_data <- setNames(object = mammals$mammal.mass$mean,
#'                          nm = row.names(mammals$mammal.mass))[mammals_tree$tip.label]
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' mammals_ACE <- phytools::fastAnc(tree = mammals_tree, x = mammals_data)
#'
#' # Run a Stochastic Mapping based on a Brownian Motion model
#' # to interpolate values along branches and obtain a "contMap" object
#' mammals_contMap <- phytools::contMap(mammals_tree, x = mammals_data,
#'                                      res = 100, # Number of time steps
#'                                      plot = FALSE)
#'
#' # Set focal time to 80 Mya
#' focal_time <- 80
#'
#' ## Extract trait data and update contMap for the given focal_time
#'
#' # Extract from the contMap (values are not exact ML estimates)
#' mammals_test <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = mammals_contMap,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#' # Extract from tip data and ML estimates of ancestral characters (values are true ML)
#' mammals_test <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = mammals_contMap,
#'    ace = mammals_ACE, tip_data = mammals_data,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#'
#' ## Visualize outputs
#'
#' # Print trait data
#' mammals_test$trait_data
#'
#' # Plot node labels on initial stochastic map with cut-off
#' phytools::plot.contMap(mammals_contMap)
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(mammals_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' phytools::plot.contMap(mammals_test$contMap)
#' ape::nodelabels(text = mammals_test$contMap$tree$initial_nodes_ID)
#'

## Once datasets are included in my package, remove motmot from dependencies

## Currently, For continuous traits
# Input = contMap

## Make a different function for each type of data, then make a wrapper function for all types
# extract_most_likely_trait_values_for_focal_time() = wrapper
   # Detect type of data based on $trait_data_type, but check it match the input (contMap or densityMaps) and the type of data in ace and tip_data
# extract_most_likely_trait_values_from_contMap_for_focal_time() = For continuous traits
# extract_most_likely_trait_values_from_densityMaps_for_focal_time() = For categorical traits
# extract_most_likely_range_values_from_densityMaps_for_focal_time() = For biogeographic traits

## See if better to extract most likely range/state from densityMaps instead of simmaps
# Check my scripts in preparing trait data for the Macroevolution Course

## Update the reference to extract_most_likely_trait_values_for_focal_time() to link to the wrapper function!
# Ex: in the doc of compute_STRAPP_test_for_focal_time.R

# Note that the use of ace and tip_data is useful only for contMap! densityMaps should provide accurate ace and tip_data
# Adjust doc, to explain better what is mandatory (contMap or densityMaps), and what is optional (ace and tip_data, for continuous traits)

### Possible update: Make it work with non-dichotomous trees!!!


extract_most_likely_trait_values_for_focal_time <- function (contMap = NULL,
                                                             densityMaps = NULL,
                                                             ace = NULL,
                                                             tip_data = NULL,
                                                             trait_data_type = "continuous", # Remove entry from the wrapper. Remove argument from the sub-functions
                                                             focal_time,
                                                             update_contMap = FALSE, # Change for update_map in the wrapper. Adjust to update_contMap / update_densityMaps in the sub-functions
                                                             keep_tip_labels = TRUE)
{

  ### Check input validity

  {
    ## contMap / densityMaps
    # Must provide either contMap OR densityMaps depending on the type of traits
    if ((is.null(contMap) & is.null(densityMaps)))
    {
      stop(paste0("You must provide a 'contMap' (for continuous traits) OR a 'densityMaps' (for categorical and biogeographic traits).\n",
                  "See ?BAMMtools::prepare_trait_data(), ?phytools::contMap(), and ?phytools::densityMap() to learn how to generate those objects."))
    }
    if ((!is.null(contMap) & !is.null(densityMaps)))
    {
      stop(paste0("You must provide a 'contMap' (for continuous traits) OR a 'densityMaps' (for categorical and biogeographic traits) according to the type of trait data.\n",
                  "For continuous traits, 'ace' and tip_data' can also be provided in complement to use accurate ML estimates of trait values at nodes and tips."))
    }

    ## trait_data_type
    # trait_data_type must be "continuous", categorical" or "biogeographic"
    if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
    {
      stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
    }

    ## Check that what is provided contMap OR densityMaps match the trait_data_type
    if ((!is.null(contMap) & trait_data_type != "continuous"))
    {
      stop(paste0("You provided a 'contMap' but selected '",trait_data_type,"' as 'trait_data_type'. contMaps are used to map continuous traits.\n",
                  "If you wish to extract trait values for a continuous trait, provide 'trait_data_type = continuous'.\n",
                  "If you wish to extract trait states/ranges for ",trait_data_type," data, provide 'densityMaps' as input instead of 'contMap'"))
    }
    if ((!is.null(densityMaps) & !(trait_data_type %in% c("categorical", "biogeographic"))))
    {
      stop(paste0("You provided 'densityMaps' but selected '",trait_data_type,"' as 'trait_data_type'. densityMaps are used to map categorical or biogeographic data.\n",
                  "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                  "If you wish to extract trait values for a continuous trait, provide 'contMap' as input instead of 'densityMaps'.\n",
                  "\tYou can also provide 'ace' and tip_data' in complement to a 'contMap' to use accurate ML estimates of trait values at nodes and tips."))
    }

    # Extract root age
    if (trait_data_type == "continuous")
    {
      root_age <- max(phytools::nodeHeights(contMap$tree)[,2])
    } else {
      ### To adjust for densityMaps!!! ########
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

    ### Only for continuous traits (to move into the sub-function)

    ## contMap
    # contMap must be a "contMap" class object
    if (!("contMap" %in% class(contMap)))
    {
      stop("'contMap' must have the 'contMap' class. See ?phytools::contMap() and ?deepSTRAPP::prepare_data() to learn how to generate those objects.")
    }
    # contMap$tree must have a $maps element
    if (is.null(contMap$tree$maps))
    {
      stop(paste0("'contMap' must have a $tree$maps element that provides the mapping of the evolution of the continuous trait on the phylogeny.\n",
                  "See ?phytools::contMap() and ?deepSTRAPP::prepare_data() to learn how to generate those objects."))
    }

    ## ace
    if (!is.null(ace))
    {
      # ace must be a named numerical vector
      if (!is.numeric(ace))
      {
        stop(paste0("'ace' is used only for continuous traits. It must be a named numerical vector that provides trait values for internal nodes.\n",
                    "The object you provided is not a numerical vector."))
      }
      # ace should have many values as there are internal nodes in the contMap$tree
      if (length(ace) != contMap$tree$Nnode)
      {
        stop(paste0("'ace' should have as many values as there are internal nodes in the contMap$tree.\n",
                    "Number of values in 'ace' = ",length(ace),"; number of internal nodes in the contMap$tree = ",contMap$tree$Nnode,"."))
      }
      # names(ace) = internal node IDs
      if (!all(as.numeric(names(ace)) %in% (length(contMap$tree$tip.label) + 1):(length(contMap$tree$tip.label) + contMap$tree$Nnode)))
      {
        stop(paste0("'names(ace)' should match numerical ID of internal nodes in the contMap$tree."))
      }
      if (!all(as.numeric(names(ace)) == (length(contMap$tree$tip.label) + 1):(length(contMap$tree$tip.label) + contMap$tree$Nnode)))
      {
        warning(paste0("Values in 'ace' are not ordered in increasing numerical ID of internal nodes.\n",
                       "They were reordered to follow the numerical ID of internal nodes."))
      }
    }

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named numerical vector
      if (!is.numeric(tip_data))
      {
        stop(paste0("'tip_data' is used only for continuous traits. It must be a named numerical vector that provides trait values for tips.\n",
                    "The object you provided is not a numerical vector."))
      }
      # tip_data should have many values as there are tips in the contMap$tree
      if (length(tip_data) != length(contMap$tree$tip.label))
      {
        stop(paste0("'tip_data' should have as many values as there are tips in the contMap$tree.\n",
                    "Number of values in 'tip_data' = ",length(tip_data),"; number of tips in the contMap$tree = ",length(contMap$tree$tip.label),"."))
      }
      # names(tip_data) = contMap$tree$tip.label
      if (!all(names(tip_data) %in% contMap$tree$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the contMap$tree$tip.label."))
      }
      if (!all(names(tip_data) == contMap$tree$tip.label))
      {
        warning(paste0("Values in 'tip_data' are not ordered as tip labels in the contMap$tree.\n",
                       "They were reordered to follow tip labels."))
      }
    }


    ### Only for categorical traits (to move into the sub-function)

    # densityMap must be a list with "densityMap" class objects
    # names(densityMap) should be the states (just check that they are not NULL and are character strings as is.character())

    # Check that if present, the densityMaps$trait_data_type match the one provided in trait_data_type

    # ace and tip_data must be NULL as they are not used

    ### Only for biogeographic traits (to move into the sub-function)

    # densityMap must be a list with "densityMap" class objects
    # names(densityMap) should be the ranges (just check that they are not NULL and are character strings as is.character())

    # Check that if present, the densityMaps$trait_data_type match the one provided in trait_data_type

    # ace and tip_data must be NULL as they are not used
  }


  if (is.null(ace))
  {
    cat(paste0("WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using values interpolated in the contMap instead.\n"))
  }
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using values interpolated in the contMap instead.\n"))
  }


  ## Extract node values if provided with 'ace' and 'tip_data'
  if (!is.null(ace) & !is.null(tip_data))
  {
    node_data_is_provided <- T

    # Reorder and rename tips according to their node index
    tip_data <- tip_data[match(x = contMap$tree$tip.label, table = names(tip_data))]
    names(tip_data) <- 1:length(tip_data)

    # Reorder ace according to their node index
    ace <- ace[order(as.numeric(names(ace)))]

    # Concatenate ACE and current tip values
    node_data <- c(tip_data, ace)
  } else {
    node_data_is_provided <- F
  }

  ## Identify edges present at focal time

  # Edge, rootward_node, tipward_node, length (once cut)

  # Get node ages per edge (no root edge)
  all_edges_df <- phytools::nodeHeights(contMap$tree)
  root_age <- max(phytools::nodeHeights(contMap$tree)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 5)) # # May be an issue for trees with very short time span
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- contMap$tree$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  # Assign tip.labels. Use tipward_node_ID for internal edges
  all_edges_df$tip.label <- contMap$tree$tip.label[match(x = all_edges_df$tipward_node_ID, 1:length(contMap$tree$tip.label))]
  all_edges_df$tip.label[is.na(all_edges_df$tip.label)] <- all_edges_df$tipward_node_ID[is.na(all_edges_df$tip.label)]

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- contMap$tree$edge[which.min(contMap$tree$edge[, 1] %in% contMap$tree$edge[, 2]), 1]

  # Identify edges present at the focal time
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > focal_time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= focal_time
  all_edges_df$time_test <- all_edges_df$rootward_test & all_edges_df$tipward_test

  # If no edge present, send warning
  if (sum(all_edges_df$time_test) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_contMap)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "continuous"))
    } else {
      # Return a NULL object for contMap
      updated_contMap <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "continuous", contMap = updated_contMap))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$time_test, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

    # Compute node distances to focal time
    present_edges_df$rootward_node_dist <- abs(present_edges_df$rootward_node_age - focal_time)
    present_edges_df$tipward_node_dist <- abs(present_edges_df$tipward_node_age - focal_time)

    # Initiate fields for (scaled) ACE values at nodes
    present_edges_df$rootward_node_scaled_ACE <- NA
    present_edges_df$tipward_node_scaled_ACE <- NA
    present_edges_df$rootward_node_ACE <- NA
    present_edges_df$tipward_node_ACE <- NA

    # Create vector to convert scaled values from 0 to 1000 in the contMap to initial values
    if (!node_data_is_provided)
    {
      trans <- 0:1000/1000 * (contMap$lims[2] - contMap$lims[1]) + contMap$lims[1]
      names(trans) <- 0:1000
    }

    # Loop per edge
    for (i in 1:nrow(present_edges_df))
    {
      # i <- 1

      # If no ace provided, extract ace from contMap using the value interpolated for the closest edge segment to the nodes
      if (!node_data_is_provided)
      {
        ## Extract scaled ACE values at nodes from the contMap

        # Extract edge ID
        edge_ID_i <- as.numeric(present_edges_df$edge_ID[i])

        # Extract associated edge mapping
        edge_map_i <- contMap$tree$maps[[edge_ID_i]]

        # Extract rootward node scaled ACE values as the first mapped values on the edge
        # Discrepancy with actual rootward node ACE values as this is the expected value for the mean age of the first segment of the edge...
        present_edges_df$rootward_node_scaled_ACE[i] <- as.numeric(names(edge_map_i)[1])

        # Extract tipward node scaled ACE values as the first mapped values on the edge
        # Discrepancy with actual tipward node ACE values as this is the expected value for the mean age of the last segment of the edge...
        present_edges_df$tipward_node_scaled_ACE[i] <- as.numeric(names(edge_map_i)[length(edge_map_i)])

        ## Convert the scaled ACE values at nodes into initial values

        present_edges_df$rootward_node_ACE[i] <- trans[as.character(present_edges_df$rootward_node_scaled_ACE[i])]
        present_edges_df$tipward_node_ACE[i] <- trans[as.character(present_edges_df$tipward_node_scaled_ACE[i])]

      # If 'ace' and 'tip_data' provided, match node_data with appropriate node
      } else {

        # Extract rootward node scaled ACE values
        present_edges_df$rootward_node_ACE[i] <- as.numeric(node_data[as.character(present_edges_df$rootward_node_ID[i])])

        # Extract tipward node scaled ACE values
        present_edges_df$tipward_node_ACE[i] <- as.numeric(node_data[as.character(present_edges_df$tipward_node_ID[i])])
      }

      ## Interpolate trait value at focal time

      # Based on equations from Felsenstein, 1985
      # Estimate ACE along an edge at a specific time-step as a weighted mean of node values with weights being the inverse distance to the nodes
       # ACE = (Xr/Dr + Xt/Dt) / (1/Dr + 1/Dt)
         # Xr = Trait value at rootward node
         # Dr = Distance from focal time to rootward node
         # Xt = Trait value at tipward node
         # Dt = Distance from focal time to tipward node

      # Case when focal time is different from rootward/tipward time
      if (all(c(present_edges_df$rootward_node_dist[i], present_edges_df$tipward_node_dist[i]) != 0))
      {
        present_edges_df$ACE_at_focal_time[i] <- ((present_edges_df$rootward_node_ACE[i]/present_edges_df$rootward_node_dist[i]) + (present_edges_df$tipward_node_ACE[i]/present_edges_df$tipward_node_dist[i])) / ((1/present_edges_df$rootward_node_dist[i]) + (1/present_edges_df$tipward_node_dist[i]))

      }

      # Case when focal time is rootward
      if (present_edges_df$rootward_node_dist[i] == 0)
      {
        present_edges_df$ACE_at_focal_time[i] <- present_edges_df$rootward_node_ACE[i]

      }

      # Case when focal time is tipward
      if (present_edges_df$tipward_node_dist[i] == 0)
      {
        present_edges_df$ACE_at_focal_time[i] <- present_edges_df$tipward_node_ACE[i]

      }
    }

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ACE_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
    }


    ## Update contMap if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_contMap)
    {
      ## Cut contMap$tree at focal time and update trait mapping in contMap$tree$maps and contMap$tree$mapped.edge
      updated_contMap <- cut_contMap_for_focal_time(contMap = contMap, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_contMap)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "continuous"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "continuous", contMap = updated_contMap))
    }
  }
}


## Make unit tests for ultrametric (eel.tree / eel_contMap) and non-ultrametric trees (mammals$mammals.phy / mammals_contMap)

## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0


