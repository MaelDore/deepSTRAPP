#' @title Cut off the character mapping on a phylogeny for a given focal time in the past
#'
#' @description Cuts off the mapping of a trait on a phylogeny over all branches
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   The phylogenetic tree is NOT updated. The mapping in `$mapped.edge` is NOT updated.
#'   Only the mapping in `$maps` is updated.
#'
#' @param tree_with_maps Object of class `"phylo"` with a character mapping stored in the `$maps` element.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   for which the character mapping must be cut. It must be smaller than the root age of the phylogeny.
#'
#' @importFrom phytools nodeHeights
#'
#' @return The function returns the mapped phylogeny with the mapping along edges/branches in `$maps` updated
#'  such as it ends at the `focal_time`.
#'
#' @author Maël Doré
#'
#' @noRd
#'

### Possible update: Make it work with non-dichotomous trees!!!


update_maps_for_focal_time <- function(tree_with_maps, focal_time)
{
  ### Check input validity
  {
    ## tree_with_maps
    # tree_with_maps must be an object with the two classes: "simmap" and "phylo".
    if (!all(c("phylo", "simmap") %in% class(tree_with_maps)))
    {
      stop("'tree_with_maps' must have the 'simmap' and 'phylo' classes. See ?phytools::make.simmap() to learn how to create simmap objects.")
    }
    # tree_with_maps must be rooted
    if (!(ape::is.rooted(tree_with_maps)))
    {
      stop(paste0("'tree_with_maps' must be a rooted phylogeny."))
    }
    # tree_with_maps must be fully resolved/dichotomous
    if (!(ape::is.binary(tree_with_maps)))
    {
      stop(paste0("'tree_with_maps' must be a fully resolved/dichotomous/binary phylogeny."))
    }
    # tree_with_maps must have a $maps element
    if (is.null(tree_with_maps$maps))
    {
      stop("'tree_with_maps' must be a list with a $maps element that contains the trait mapping.")
    }
    # tree_with_maps$maps must have be a list of N items. N = number of edges
    if (length(tree_with_maps$maps) != nrow(tree_with_maps$edge))
    {
      stop(paste0("'tree_with_maps$maps' does not match with the number of edges in the phylogeny.\n",
                  "'$maps' has mapping information for ",length(tree_with_maps$maps)," edges.\n",
                  "There are ",nrow(tree_with_maps$edge)," edges in the phylogeny."))
    }

    ## Extract root age
    root_age <- max(phytools::nodeHeights(tree_with_maps)[,2])

    ## focal_time
    # focal_time must be positive and smaller than root age
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

  ## Initiate new tree_with_maps
  tree_with_updated_maps <- tree_with_maps

  ## Identify edges present at focal time

  # Define level of tolerance used to round ages
  tol <- root_age * 10^-5
  closest_power <- round(log10(tol))
  closest_power <- min(closest_power, 0) # Use 0 as the minimal power

  # Get node ages per edge (no root edge)
  all_edges_df <- phytools::nodeHeights(tree_with_maps)
  # all_edges_df <- as.data.frame(round(root_age - all_edges_df, 5)) # May be an issue for trees with very short time span
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, -1*closest_power))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  all_edges_df$rootward_node_age

  # Get nodes ID per edge
  all_edges_ID_df <- tree_with_maps$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  root_node_ID <- tree_with_maps$edge[which.min(tree_with_maps$edge[, 1] %in% tree_with_maps$edge[, 2]), 1]

  ## Keep only edges present at focal time and before

  # Identify edges present at the focal time
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > focal_time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= focal_time

  # Remove edges that are younger than the focal_time
  focal_edges_df <- all_edges_df[all_edges_df$rootward_test, ]

  ## Adjust length of edge present at focal time
  focal_edges_df$initial_length <- tree_with_maps$edge.length[as.numeric(focal_edges_df$edge_ID)]
  focal_edges_df$updated_length <- focal_edges_df$initial_length
  focal_edges_df$updated_length[focal_edges_df$tipward_test] <- focal_edges_df$rootward_node_age[focal_edges_df$tipward_test] - focal_time

  ## Update maps edge per edge
  updated_maps <- list()

  # Loop per edge
  for (i in 1:nrow(focal_edges_df))
  {
    # i <- 1

    # Extract edge ID
    edge_ID_i <- as.numeric(focal_edges_df$edge_ID[i])
    # Extract time to cut
    edge_length_i <- focal_edges_df$updated_length[i]
    # Extract current edge map
    edge_map_i <- tree_with_maps$maps[[edge_ID_i]]

    # If updated length = current length, keep current maps
    if ((sum(edge_map_i) <= edge_length_i)) # Use <= instead of == to avoid issue with inaccuracy of rounded values from floating points
    {
      updated_maps[[i]] <- edge_map_i

    } else { # If updated length != current length (should be smaller), update current maps to remove the maps from time younger than the focal_time

      # Initiate updated map for edge i
      updated_edge_map_i <- edge_map_i

      # Identify segment matching the time cut
      cut_position_i <- which.min(cumsum(edge_map_i) < edge_length_i)
      # Remove younger segments
      updated_edge_map_i <- updated_edge_map_i[1:cut_position_i]
      # Compute new length of the last segment
      if (cut_position_i > 1)
      { # Case when time cut does not happen in first segment, need to remove the previous segment length
        last_segment_updated_length <- edge_length_i - cumsum(updated_edge_map_i)[cut_position_i - 1]
      } else { # Case when time cut happens in first segment, new edge length = new segment length
        last_segment_updated_length <- edge_length_i
      }
      names(last_segment_updated_length) <- names(updated_edge_map_i[cut_position_i])
      updated_edge_map_i[cut_position_i] <- last_segment_updated_length

      # Store updated map for edge i
      updated_maps[[i]] <- updated_edge_map_i
    }
  }

  ## Length of maps and updated length should be equal
  # cbind(focal_edges_df$updated_length, lapply(X = updated_maps, FUN = sum))

  # Store updated maps
  tree_with_updated_maps$maps <- updated_maps

  # Export phylogeny with updated $maps
  return(tree_with_updated_maps)
}


## Make unit tests for ultrametric (eel.tree / eel_contMap) and non-ultrametric trees (tortoise.tree / tortoise_contMap)

## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0
