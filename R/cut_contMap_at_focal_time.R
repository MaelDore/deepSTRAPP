#' @title Cuts the phylogeny and continuous trait mapping at a given time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'   Likewise, remove continuous trait mapping for the cut off branches
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param contMap Object of class `"contMap"`, typically generated with [phytools::contMap()],
#'   that contains a phylogenetic tree and associated continuous trait mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights
#'
#' @details The phylogenetic tree is cut at a specific time in the past (i.e. the `focal_time`).
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
#' @return The function returns the cut contMap as an object of class `"contMap"`.
#'   It contains a `$tree` element of classes `"simmap"` and `"phylo"`. This function updates and adds multiple useful sub-elements to the `$tree` element.
#'   * `$maps` An updated list of named numerical vectors. Provides the mapping of trait values along each remaining edge.
#'   * `$mapped.edge` An updated matrix. Provides the evolutionary time spent across trait values (columns) along the remaining edges (rows).
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_at_focal_time()] [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
#' @examples
#' # ----- Prepare data ----- #
#'
#' # Load mammals phylogeny and data from the R package motmot
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#' library(motmot)
#'
#' data("mammals")
#' force(mammals)
#'
#' mammals_tree <- mammals$mammal.phy
#' mammals_data <- setNames(object = mammals$mammal.mass$mean,
#'                          nm = row.names(mammals$mammal.mass))[mammals_tree$tip.label]
#'
#' # Run a stochastic mapping based on a Brownian Motion model
#' # for Ancestral Trait Estimates to obtain a "contMap" object
#' mammals_contMap <- phytools::contMap(mammals_tree, x = mammals_data,
#'                                      res = 100, # Number of time steps
#'                                      plot = FALSE)
#'
#' # Set focal time
#' focal_time <- 80
#'
#' # ----- Example 1: keep_tip_labels = TRUE ----- #
#'
#' # Cut contMap to 80 Mya while keeping tip.label
#' # on terminal branches with a unique descending tip.
#' updated_contMap <- cut_contMap_at_focal_time(contMap = mammals_contMap,
#'                                              focal_time = focal_time,
#'                                              keep_tip_labels = TRUE)
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(mammals_contMap, lwd = 2)
#' nodelabels()
#' abline(v = max(phytools::nodeHeights(mammals_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)

#' # Plot initial node labels on cut stochastic map
#' plot(updated_contMap)
#' nodelabels(text = updated_contMap$tree$initial_nodes_ID)
#'
#' # ----- Example 2: keep_tip_labels = FALSE ----- #
#'
#' # Cut contMap to 80 Mya while NOT keeping tip.label.
#' updated_contMap <- cut_contMap_at_focal_time(contMap = mammals_contMap,
#'                                              focal_time = focal_time,
#'                                             keep_tip_labels = FALSE)
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(mammals_contMap)
#' nodelabels()
#' abline(v = max(phytools::nodeHeights(mammals_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map
#' plot(updated_contMap)
#' nodelabels(text = updated_contMap$tree$initial_nodes_ID)
#'

### Possible update: Make it work with non-dichotomous trees!!!


cut_contMap_at_focal_time <- function(contMap, focal_time, keep_tip_labels = TRUE)
{
  ### Check input validity

  # contMap must be a "contMap" class object

  # contMap$tree must be an object with the two classes: "simmap" and "phylo"
  # contMap$tree must be rooted
    # ape::is.rooted
  # contMap$tree must be fully resolved/dichotomous
    # ape::is.binary

  # focal_time must be positive and smaller to root age
    # If focal_time = root_age, throw a specific error
  # focal_time must be numerical

  # keep_tip_labels must be TRUE or FALSE

  ## Extract root age
  root_age <- max(phytools::nodeHeights(contMap$tree)[,2])

  # If focal_time equals root_age, send warning
  if (focal_time == root_age)
  {
    warning(paste0(focal_time, " equals root age = ",root_age," Mya. Return an empty object.\n"))

    # Return a NULL object
    updated_contMap <- NULL
    return(updated_contMap)

  } else {

    ## Cut contMap$tree at focal time

    updated_contMap_tree <- contMap
    updated_contMap_tree$tree <- cut_phylo_at_focal_time(tree = updated_contMap_tree$tree, focal_time = focal_time, keep_tip_labels = keep_tip_labels)

    ## Update contMap$tree$maps for focal time

    updated_contMap_maps <- contMap
    updated_contMap_maps$tree <- update_maps_at_focal_time(tree_with_maps = updated_contMap_maps$tree, focal_time = focal_time)

    ## Merge outputs
    updated_contMap <- updated_contMap_tree
    updated_contMap$tree$maps <- updated_contMap_maps$tree$maps

    ## Update $mapped.edge
    updated_contMap$tree$mapped.edge <- makeMappedEdge(edge = updated_contMap$tree$edge, maps = updated_contMap$tree$maps)
    updated_contMap$tree$mapped.edge <- updated_contMap$tree$mapped.edge[, order(as.numeric(colnames(updated_contMap$tree$mapped.edge)))]

    # Export contMap with updated tree and character mapping ($tree$maps and $tree$mapped.edge)
    return(updated_contMap)
  }
}


## Make unit tests for ultrametric (eel.tree) and non-ultrametric trees (tortoise.tree)


## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0


## Helper function used to generate $mapped.edge from $edge and $maps

# Internal function copied from the R package phytools
# Source: phytools:::makeMappedEdge()
# Author: Liam Revell: liam.revell@umb.edu

makeMappedEdge <- function (edge, maps)
{
  st <- sort(unique(unlist(sapply(maps, function(x) names(x)))))
  mapped.edge <- matrix(0, nrow(edge), length(st))
  rownames(mapped.edge) <- apply(edge, 1, function(x) paste(x, collapse = ","))
  colnames(mapped.edge) <- st
  for (i in 1:length(maps))
  {
    for (j in 1:length(maps[[i]]))
    {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }

  return(mapped.edge)
}
