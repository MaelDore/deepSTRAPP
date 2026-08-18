
### Functions to cut the phylogeny and ancestral trait mapping of a simmap for a given focal time in the past

## 1/ cut_simmap_for_focal_time() => cut a single simmap object
## 2/ cut_simmaps_for_focal_time() => cut all simmaps in a list


## 1/ Function to cut a single simmap object ####

#' @title Cut the phylogeny and categorical trait/range mapping for a given focal time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'   Likewise, remove mapping for the cut off branches
#'   by updating the `$maps` and `$mapped.edge` elements.
#'
#' @param simmap Object with the classes `"phylo"` and `"simmap"`, typically generated with [phytools::make.simmap()],
#'   that contains a phylogenetic tree and associated categorical trait/range mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   for which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#'
#' @details The mapped phylogenetic tree is cut for a specific time in the past (i.e. the `focal_time`).
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
#'   The stochastic trait mapping is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return The function returns the cut/updated simmap as a list of at least six elements with the classes `"phylo"` and `"simmap"`:
#'
#'  Initial updated elements:
#'    * `$edge` A numerical matrix with two columns listing the updated ID of rootward and tipward nodes of the remaining edges
#'    * `$edge.length` A numerical vector. Updated length of remaining edges.
#'    * `$Nnode` Integer. Number of remaining internal nodes.
#'    * `$tip.label` Vector of character strings. Labels of the branches cut at the `focal_time`.
#'      + If `keep_tip_labels = TRUE`, the cut terminal branches are labeled with the tip.label of their unique descendant tip.
#'      + If `keep_tip_labels = FALSE`, the cut terminal branches are labeled with the node ID of their unique descendant tip.
#'    * `$maps` An updated list of named numerical vectors. Provides the mapping of trait values along each remaining edge.
#'    * `$mapped.edge` An updated matrix. Provides the evolutionary time spent across trait values (columns) along the remaining edges (rows).
#'
#'  Additional elements:
#'    * `$root_age` Integer. Stores the age of the root of the tree.
#'    * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'    * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'    * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'    * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("cut_phylogenies", package = "deepSTRAPP")}
#'
#' @examples
#' # ----- Prepare data ----- #
#'
#' ## Load mammals phylogeny and data from the R package motmot, and implemented in deepSTRAPP
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

#' \donttest{ # (May take several seconds to run)
#' # Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100,
#'                                        return_simmaps = TRUE,
#'                                        plot_map = FALSE)
#'
#' # Extract the simulated evolutionary histories (simmaps)
#' mammals_simmaps <- mammals_cat_data$simmaps
#' # Extract simulation n°1
#' mammals_simmap_1 <- mammals_simmaps[[1]]
#'
#' # Set focal time to 80Mya
#' focal_time <- 80
#'
#' # ----- Example 1: keep_tip_labels = TRUE ----- #
#'
#' # Cut the simmap to 80 Mya while keeping tip.label
#' # on terminal branches with a unique descending tip.
#' updated_simmap_1 <- cut_simmap_for_focal_time(simmap = mammals_simmap_1,
#'                                               focal_time = focal_time,
#'                                               keep_tip_labels = TRUE)
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(mammals_simmap_1, lwd = 2, fsize = 0.5,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_simmap_1)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map
#' plot(updated_simmap_1, fsize = 0.8,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.8, text = updated_simmap_1$initial_nodes_ID)
#'
#' # ----- Example 2: keep_tip_labels = FALSE ----- #
#'
#' # Cut the simmap to 80 Mya while NOT keeping tip.label.
#' updated_simmap_1 <- cut_simmap_for_focal_time(simmap = mammals_simmap_1,
#'                                               focal_time = focal_time,
#'                                               keep_tip_labels = FALSE)
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(mammals_simmap_1, lwd = 2, fsize = 0.5,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_simmap_1)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map
#' plot(updated_simmap_1, fsize = 0.8,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.8, text = updated_simmap_1$initial_nodes_ID)
#' }
#'

cut_simmap_for_focal_time <- function(simmap, focal_time, keep_tip_labels = TRUE)
{
  ### Check input validity
  {
    ## simmap
    # simmap must be an object with the two classes: "simmap" and "phylo"
    if (!(all(c("simmap", "phylo") %in% class(simmap))))
    {
      stop(paste0("'simmap' must have the 'simmap' and 'phylo' classes indicating a trait is mapped on the phylogeny.\n",
                  "See ?phytools::make.simmap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # simmap must be rooted
    if (!(ape::is.rooted(simmap)))
    {
      stop(paste0("'simmap' must be a rooted phylogeny."))
    }
    # simmap must be fully resolved/dichotomous
    if (!(ape::is.binary(simmap)))
    {
      stop(paste0("'simmap' must be a fully resolved/dichotomous/binary phylogeny."))
    }

    ## Extract root age
    root_age <- max(phytools::nodeHeights(simmap)[,2])

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

  ## Cut phylogeny at focal time
  updated_simmap_tree <- cut_phylo_for_focal_time(tree = simmap, focal_time = focal_time, keep_tip_labels = keep_tip_labels)

  ## Update simmap $maps for focal time
  updated_maps <- update_maps_for_focal_time(tree_with_maps = simmap, focal_time = focal_time)

  ## Merge outputs
  updated_simmap <- updated_simmap_tree
  updated_simmap$maps <- updated_maps$maps

  ## Update $mapped.edge
  updated_simmap$mapped.edge <- makeMappedEdge(edge = updated_simmap$edge, maps = updated_simmap$maps)
  updated_simmap$mapped.edge <- updated_simmap$mapped.edge[, order(colnames(updated_simmap$mapped.edge))]

  # Export simmap with updated tree and character mapping ($maps and $mapped.edge)
  return(updated_simmap)
}

# Possible update: Make it work with non-dichotomous trees!!!


## 2/ Function to cut multiple simmaps stored in a list ####

#' @title Cut the phylogenies and categorical trait/range mappings of a list of simmaps for a given focal time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'   Likewise, remove continuous trait mapping for the cut off branches
#'   by updating the `$maps` and `$mapped.edge` elements.
#'
#' @param simmaps List of objects with the classes `"phylo"` and `"simmap"`, typically generated with [phytools::make.simmap()],
#'   that contains a phylogenetic tree and associated categorical trait/range mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   for which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#'
#' @details The phylogenetic trees are cut for a specific time in the past (i.e. the `focal_time`).
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
#'   The categorical trait/range mappings are updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return The function returns an updated list of objects as cut/updated simmaps of classes `"phylo"` and `"simmap"`.
#'
#'  Each simmap object represent an updated stochastic map cut/updated simmap as a list of at least six elements:
#'
#'  Initial updated elements:
#'    * `$edge` A numerical matrix with two columns listing the updated ID of rootward and tipward nodes of the remaining edges
#'    * `$edge.length` A numerical vector. Updated length of remaining edges.
#'    * `$Nnode` Integer. Number of remaining internal nodes.
#'    * `$tip.label` Vector of character strings. Labels of the branches cut at the `focal_time`.
#'      + If `keep_tip_labels = TRUE`, the cut terminal branches are labeled with the tip.label of their unique descendant tip.
#'      + If `keep_tip_labels = FALSE`, the cut terminal branches are labeled with the node ID of their unique descendant tip.
#'    * `$maps` An updated list of named numerical vectors. Provides the mapping of trait values along each remaining edge.
#'    * `$mapped.edge` An updated matrix. Provides the evolutionary time spent across trait values (columns) along the remaining edges (rows).
#'
#'  Additional elements:
#'    * `$root_age` Integer. Stores the age of the root of the tree.
#'    * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'    * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'    * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'    * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("cut_phylogenies", package = "deepSTRAPP")}
#'
#' @examples
#' # ----- Prepare data ----- #
#'
#' ## Load mammals phylogeny and data from the R package motmot, and implemented in deepSTRAPP
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
#' \donttest{ # (May take several seconds to run)
#' # Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100,
#'                                        return_simmaps = TRUE,
#'                                        plot_map = FALSE)
#'
#' # Extract the simulated evolutionary histories (simmaps)
#' mammals_simmaps <- mammals_cat_data$simmaps
#'
#' # Set focal time to 80Mya
#' focal_time <- 80
#'
#' # ----- Example 1: keep_tip_labels = TRUE ----- #
#'
#' # Cut the simmap to 80 Mya while keeping tip.label
#' # on terminal branches with a unique descending tip.
#' updated_simmaps <- cut_simmaps_for_focal_time(simmaps = mammals_simmaps,
#'                                               focal_time = focal_time,
#'                                               keep_tip_labels = TRUE)
#'
#' # Plot node labels on initial stochastic map n°1 with cut-off
#' plot(mammals_simmaps[[1]], lwd = 2, fsize = 0.5,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_simmaps[[1]])[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map n°1
#' plot(updated_simmaps[[1]], fsize = 0.8,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.8, text = updated_simmaps[[1]]$initial_nodes_ID)
#'
#' # ----- Example 2: keep_tip_labels = FALSE ----- #
#'
#' # Cut the simmap to 80 Mya while NOT keeping tip.label.
#' updated_simmaps <- cut_simmaps_for_focal_time(simmaps = mammals_simmaps,
#'                                               focal_time = focal_time,
#'                                               keep_tip_labels = FALSE)
#'
#' # Plot node labels on initial stochastic map n°1 with cut-off
#' plot(mammals_simmaps[[1]], lwd = 2, fsize = 0.5,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_simmaps[[1]])[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot initial node labels on cut stochastic map n°1
#' plot(updated_simmaps[[1]], fsize = 0.8,
#'      colors = c(large = "black", medium = "blue", small = "dodgerblue"))
#' ape::nodelabels(cex = 0.8, text = updated_simmaps[[1]]$initial_nodes_ID)
#' }
#'


cut_simmaps_for_focal_time <- function(simmaps, focal_time, keep_tip_labels = TRUE)
{
  ## Loop across simmaps
  updated_simmaps <- lapply(X = simmaps, FUN = cut_simmap_for_focal_time,
                            focal_time = focal_time,
                            keep_tip_labels = keep_tip_labels)

  # Export simmaps with updated trees and character mappings (e$maps and $mapped.edge)
  return(updated_simmaps)
}

