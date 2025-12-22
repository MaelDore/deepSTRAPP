

### Functions to cut the phylogeny and posterior probability mapping of a categorical trait for a given focal time in the past
## cut_densityMap_for_focal_time() => cut a single densityMap object
## cut_densityMaps_for_focal_time() => cut all densityMaps in a list

## Function to cut a single densityMap object ####

#' @title Cut the phylogeny and posterior probability mapping of a categorical trait for a given focal time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'   Likewise, remove posterior probability mapping of the categorical trait for the cut off branches
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param densityMap Object of class `"densityMap"`, typically generated with [phytools::densityMap()],
#'   that contains a phylogenetic tree and associated posterior probability mapping of a categorical trait.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   for which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.densityMap
#'
#' @details The phylogenetic tree is cut for a specific time in the past (i.e. the `focal_time`).
#'
#'   When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'   the leaf left is labeled with the tip.label of the unique descendant tip.
#'
#'   When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'   the leaf left is labeled with the node ID of the unique descendant tip.
#'
#'   In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'   the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The posterior probability mapping of a categorical trait is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return The function returns the cut densityMap as an object of class `"densityMap"` with three elements.
#'
#'   It contains a `$tree` element of classes `"simmap"` and `"phylo"`. This function updates and adds multiple useful sub-elements to the `$tree` element.
#'     * `$maps` An updated list of named numerical vectors. Provides the mapping of posterior probability of the state along each remaining edge.
#'     * `$mapped.edge` An updated matrix. Provides the evolutionary time spent across posterior probabilities (columns) along the remaining edges (rows).
#'     * `$root_age` Integer. Stores the age of the root of the tree.
#'     * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#'   The `$col` element describes the colors used to map each possible posterior probability value from 0 to 1000.
#'
#'   The `$states` element provide the name of the states. Here, the first value is the absence of the state labeled as "Not X" with X being the state.
#'   The second value is the name of the state.
#'
#'   High posterior probability reflects high likelihood to harbor the state. Low probability reflects high likelihood to NOT harbor the state.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'  [deepSTRAPP::extract_most_likely_states_from_densityMaps_for_focal_time()]
#'
#'  For a guided tutorial, see this vignette: \code{vignette("cut_phylogenies", package = "deepSTRAPP")}
#'
#' @examples
#' # ----- Prepare data ----- #
#'
#' # Load mammals phylogeny and data from the R package motmot included within deepSTRAPP
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#' data("mammals", package = "deepSTRAPP")
#'
#' # Obtain mammal tree
#' mammals_tree <- mammals$mammal.phy
#' # Convert mass data into categories
#' mammals_mass <- setNames(object = mammals$mammal.mass$mean,
#'                          nm = row.names(mammals$mammal.mass))[mammals_tree$tip.label]
#' mammals_data <- mammals_mass
#' mammals_data[seq_along(mammals_data)] <- "small"
#' mammals_data[mammals_mass > 5] <- "medium"
#' mammals_data[mammals_mass > 10] <- "large"
#' table(mammals_data)
#'
#' \dontrun{  (May take several minutes to run)
#' # Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100,
#'                                        plot_map = FALSE)
#'
#' # Set focal time
#' focal_time <- 80
#'
#' # Extract the density map for small mammals (state 3 = "small")
#' mammals_densityMap_small <- mammals_cat_data$densityMaps[[3]]
#'
#' # ----- Example 1: keep_tip_labels = TRUE ----- #
#'
#' # Cut densityMap to 80 Mya while keeping tip.label
#' # on terminal branches with a unique descending tip.
#' updated_mammals_densityMap_small <- cut_densityMap_for_focal_time(
#'      densityMap = mammals_densityMap_small,
#'      focal_time = focal_time,
#'      keep_tip_labels = TRUE)
#'
#' # Plot node labels on initial stochastic map with cut-off
#' phytools::plot.densityMap(mammals_densityMap_small, fsize = 0.7, lwd = 2)
#' ape::nodelabels(cex = 0.7)
#' abline(v = max(phytools::nodeHeights(mammals_densityMap_small$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map
#' phytools::plot.densityMap(updated_mammals_densityMap_small, fsize = 0.8)
#' ape::nodelabels(cex = 0.8, text = updated_mammals_densityMap_small$tree$initial_nodes_ID)
#'
#' # ----- Example 2: keep_tip_labels = FALSE ----- #
#'
#' # Cut densityMap to 80 Mya while NOT keeping tip.label
#' updated_mammals_densityMap_small <- cut_densityMap_for_focal_time(
#'      densityMap = mammals_densityMap_small,
#'      focal_time = focal_time,
#'      keep_tip_labels = FALSE)
#'
#' # Plot node labels on initial stochastic map with cut-off
#' phytools::plot.densityMap(mammals_densityMap_small, fsize = 0.7, lwd = 2)
#' ape::nodelabels(cex = 0.7)
#' abline(v = max(phytools::nodeHeights(mammals_densityMap_small$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map
#' phytools::plot.densityMap(updated_mammals_densityMap_small, fsize = 0.8)
#' ape::nodelabels(cex = 0.8, text = updated_mammals_densityMap_small$tree$initial_nodes_ID) }
#'

cut_densityMap_for_focal_time <- function(densityMap, focal_time, keep_tip_labels = TRUE)
{
  ### Check input validity
  {
    ## densityMap
    # densityMap must be a "densityMap" class object
    if (!("densityMap" %in% class(densityMap)))
    {
      stop("'densityMap' must have the 'densityMap' class. See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }

    ## densityMap$tree
    # densityMap$tree must be an object with the two classes: "simmap" and "phylo"
    if (!(all(c("simmap", "phylo") %in% class(densityMap$tree))))
    {
      stop(paste0("'densityMap$tree' must have the 'simmap' and 'phylo' classes indicating a trait is mapped on the phylogeny.\n",
                  "See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # densityMap$tree must be rooted
    if (!(ape::is.rooted(densityMap$tree)))
    {
      stop(paste0("'densityMap$tree' must be a rooted phylogeny."))
    }
    # densityMap$tree must be fully resolved/dichotomous
    if (!(ape::is.binary(densityMap$tree)))
    {
      stop(paste0("'densityMap$tree' must be a fully resolved/dichotomous/binary phylogeny."))
    }

    ## Extract root age
    root_age <- max(phytools::nodeHeights(densityMap$tree)[,2])

    ## focal_time
    # focal_time must be positive and smaller to root age
    if (focal_time < 0)
    {
      stop(paste0("'focal_time' must be a positive number. It represents the time as a distance from the present."))
    }
    if (focal_time >= root_age)
    {
      stop(paste0("'focal_time' must be smaller than the root age of the phylogeny.\n",
                  "'focal_time' = ",focal_time,"; root age = ",root_age,"."))
    }
  }

  ## Cut densityMap$tree at focal time

  updated_densityMap_tree <- densityMap
  updated_densityMap_tree$tree <- cut_phylo_for_focal_time(tree = updated_densityMap_tree$tree, focal_time = focal_time, keep_tip_labels = keep_tip_labels)

  ## Update densityMap$tree$maps for focal time

  updated_densityMap_maps <- densityMap
  updated_densityMap_maps$tree <- update_maps_for_focal_time(tree_with_maps = updated_densityMap_maps$tree, focal_time = focal_time)

  ## Merge outputs
  updated_densityMap <- updated_densityMap_tree
  updated_densityMap$tree$maps <- updated_densityMap_maps$tree$maps

  ## Update $mapped.edge
  # Obtain mapped.edge
  updated_densityMap$tree$mapped.edge <- makeMappedEdge(edge = updated_densityMap$tree$edge, maps = updated_densityMap$tree$maps)
  # Order columns by increasing PP
  updated_densityMap$tree$mapped.edge <- updated_densityMap$tree$mapped.edge[, order(as.numeric(colnames(updated_densityMap$tree$mapped.edge))), drop = FALSE]

  # Export densityMap with updated tree and character mapping ($tree$maps and $tree$mapped.edge)
  return(updated_densityMap)
}

## Function to cut a list of densityMap objects ####

#' @title Cut phylogenies and posterior probability mapping of each state for a given focal time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'   Likewise, remove posterior probability mapping of the categorical trait for the cut off branches
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param densityMaps List of objects of class `"densityMap"`, typically generated with [deepSTRAPP::prepare_trait_data()].
#'   Each densityMap (see [phytools::densityMap()]) contains a phylogenetic tree and associated posterior probability mapping of a categorical trait.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   for which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.densityMap
#'
#' @details The phylogenetic tree is cut for a specific time in the past (i.e. the `focal_time`).
#'
#'   When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'   the leaf left is labeled with the tip.label of the unique descendant tip.
#'
#'   When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'   the leaf left is labeled with the node ID of the unique descendant tip.
#'
#'   In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'   the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The continuous trait mapping is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return The function returns an updated list of objects as cut densityMaps of class `"densityMap"`.
#'
#'   Each densityMap object contains three elements:
#'   * `$tree` element of classes `"simmap"` and `"phylo"`. This function updates and adds multiple useful sub-elements to the `$tree` element.
#'     - `$maps` An updated list of named numerical vectors. Provides the mapping of posterior probability of the state along each remaining edge.
#'     - `$mapped.edge` An updated matrix. Provides the evolutionary time spent across posterior probabilities (columns) along the remaining edges (rows).
#'     - `$root_age` Integer. Stores the age of the root of the tree.
#'     - `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     - `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     - `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     - `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#'   * `$col` element describes the colors used to map each possible posterior probability value from 0 to 1000.
#'
#'   * `$states` element provide the name of the states. Here, the first value is the absence of the state labeled as "Not X" with X being the state.
#'   The second value is the name of the state.
#'
#'   High posterior probability reflects high likelihood to harbor the state. Low probability reflects high likelihood to NOT harbor the state.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_densityMap_for_focal_time()]
#'  [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'  [deepSTRAPP::extract_most_likely_states_from_densityMaps_for_focal_time()]
#'
#' @examples
#' # ----- Prepare data ----- #
#'
#' # Load mammals phylogeny and data from the R package motmot, and implemented in deepSTRAPP
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#' data("mammals", package = "deepSTRAPP")
#'
#' # Obtain mammal tree
#' mammals_tree <- mammals$mammal.phy
#' # Convert mass data into categories
#' mammals_mass <- setNames(object = mammals$mammal.mass$mean,
#'                          nm = row.names(mammals$mammal.mass))[mammals_tree$tip.label]
#' mammals_data <- mammals_mass
#' mammals_data[seq_along(mammals_data)] <- "small"
#' mammals_data[mammals_mass > 5] <- "medium"
#' mammals_data[mammals_mass > 10] <- "large"
#' table(mammals_data)
#'
#' \donttest{ # (May take several minutes to run)
#' # Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100,
#'                                        plot_map = FALSE)
#'
#' # Set focal time
#' focal_time <- 80
#'
#' # Extract the density maps
#' mammals_densityMaps <- mammals_cat_data$densityMaps
#'
#' # ----- Example 1: keep_tip_labels = TRUE ----- #
#'
#' # Cut densityMaps to 80 Mya while keeping tip.label
#' # on terminal branches with a unique descending tip.
#' updated_mammals_densityMaps <- cut_densityMaps_for_focal_time(
#'     densityMaps = mammals_densityMaps,
#'     focal_time = focal_time,
#'     keep_tip_labels = TRUE)
#'
#' ## Plot density maps as overlay of all state posterior probabilities
#' # ?plot_densityMaps_overlay
#'
#' # Plot initial density maps
#' plot_densityMaps_overlay(densityMaps = mammals_densityMaps, fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated/cut density maps
#' plot_densityMaps_overlay(densityMaps = updated_mammals_densityMaps, fsize = 0.8)
#'
#' # ----- Example 2: keep_tip_labels = FALSE ----- #
#'
#' # Cut densityMap to 80 Mya while NOT keeping tip.label
#' updated_mammals_densityMaps <- cut_densityMaps_for_focal_time(
#'     densityMaps = mammals_densityMaps,
#'     focal_time = focal_time,
#'     keep_tip_labels = FALSE)
#'
#' # Plot initial density maps
#' plot_densityMaps_overlay(densityMaps = mammals_densityMaps, fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated/cut density maps
#' plot_densityMaps_overlay(densityMaps = updated_mammals_densityMaps, fsize = 0.8) }
#'


cut_densityMaps_for_focal_time <- function(densityMaps, focal_time, keep_tip_labels = TRUE)
{
  ## Loop across densityMaps
  updated_densityMaps <- lapply(X = densityMaps, FUN = cut_densityMap_for_focal_time,
                                focal_time = focal_time,
                                keep_tip_labels = keep_tip_labels)

  # Export densityMaps with updated trees and character mappings ($tree$maps and $tree$mapped.edge)
  return(updated_densityMaps)
}


### Possible update: Make it work with non-dichotomous trees!!!

## Make unit tests for ultrametric (eel.tree) and non-ultrametric trees (tortoise.tree)

## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0

