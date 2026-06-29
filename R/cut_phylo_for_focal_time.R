#' @title Cut the phylogeny for a given time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'
#' @param tree Object of class `"phylo"`. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numerical. The time, in terms of time distance from the present,
#'   for which the tree must be cut. It must be smaller than the root age of the phylogeny.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights
#' @importFrom ape nodelabels edgelabels
#'
#' @details When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'   the leaf left is labeled with the tip.label of the unique descendant tip.
#'
#'   When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'   the leaf left is labeled with the node ID of the unique descendant tip.
#'
#'   In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'   the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#' @return The function returns the cut phylogeny as an object of class `"phylo"`. It adds multiple useful elements to the object.
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_contMap_for_focal_time()] [deepSTRAPP::cut_densityMaps_for_focal_time()]
#'
#' For a guided tutorial, see this vignette: \code{vignette("cut_phylogenies", package = "deepSTRAPP")}
#'
#' @examples
#' # Load eel phylogeny from the R package phytools
#' # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#' data("eel.tree", package = "phytools")
#'
#' # ----- Example 1: keep_tip_labels = TRUE ----- #
#'
#' # Cut tree to 30 Mya while keeping tip.label on terminal branches with a unique descending tip.
#' cut_eel.tree <- cut_phylo_for_focal_time(tree = eel.tree, focal_time = 30, keep_tip_labels = TRUE)
#'
#' # Plot internal node labels on initial tree with cut-off
#' plot(eel.tree, cex = 0.5)
#' abline(v = max(phytools::nodeHeights(eel.tree)[,2]) - 30, col = "red", lty = 2, lwd = 2)
#' nb_tips <- length(eel.tree$tip.label)
#' nodelabels_in_cut_tree <- (nb_tips + 1):(nb_tips + eel.tree$Nnode)
#' nodelabels_in_cut_tree[!(nodelabels_in_cut_tree %in% cut_eel.tree$initial_nodes_ID)] <- NA
#' ape::nodelabels(text = nodelabels_in_cut_tree)
#'
#' # Plot initial internal node labels on cut tree
#' plot(cut_eel.tree, cex = 0.8)
#' ape::nodelabels(text = cut_eel.tree$initial_nodes_ID)
#'
#' # Plot edge labels on initial tree with cut-off
#' plot(eel.tree, cex = 0.5)
#' abline(v = max(phytools::nodeHeights(eel.tree)[,2]) - 30, col = "red", lty = 2, lwd = 2)
#' edgelabels_in_cut_tree <- 1:nrow(eel.tree$edge)
#' edgelabels_in_cut_tree[!(1:nrow(eel.tree$edge) %in% cut_eel.tree$initial_edges_ID)] <- NA
#' ape::edgelabels(text = edgelabels_in_cut_tree)
#'
#' # Plot initial edge labels on cut tree
#' plot(cut_eel.tree, cex = 0.8)
#' ape::edgelabels(text = cut_eel.tree$initial_edges_ID)
#'
#' # ----- Example 2: keep_tip_labels = FALSE ----- #
#'
#' # Cut tree to 30 Mya without keeping tip.label on terminal branches with a unique descending tip.
#' # All tip.labels are converted to their descending/tipward node ID
#' cut_eel.tree <- cut_phylo_for_focal_time(tree = eel.tree, focal_time = 30, keep_tip_labels = FALSE)
#' plot(cut_eel.tree, cex = 0.8)
#'


### Possible update: Make it work with non-dichotomous trees!!!


cut_phylo_for_focal_time <- function(tree, focal_time, keep_tip_labels = TRUE)
{

  ### Check input validity
  {
    ## tree
    # tree must be a "phylo" class object
    if (!("phylo" %in% class(tree)))
    {
      stop("'tree' must have the 'phylo' class. See ?ape::read.tree() and ?ape::read.nexus() to learn how to import phylogenies in R.")
    }
    # tree must be rooted
    if (!(ape::is.rooted(tree)))
    {
      stop(paste0("'tree' must be a rooted phylogeny."))
    }
    # tree must be fully resolved/dichotomous
    if (!(ape::is.binary(tree)))
    {
      stop(paste0("'tree' must be a fully resolved/dichotomous/binary phylogeny."))
    }

    ## Extract root age
    root_age <- max(phytools::nodeHeights(tree)[,2])

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

  ## Initiate new tree
  cut_tree <- tree

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = cut_tree, focal_time = focal_time, tolerance = 10^-5)

  # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  root_node_ID <- tree$edge[which.min(tree$edge[, 1] %in% tree$edge[, 2]), 1]

  ## Merge tip.label to the edge df

  # If tipward node is a tip, use tip.label
  all_edges_df$tip.label <- tree$tip.label[all_edges_df$tipward_node_ID]
  # If tipward node is an internal node, use node ID
  all_edges_df$tip.label[is.na(all_edges_df$tip.label)] <- all_edges_df$tipward_node_ID[is.na(all_edges_df$tip.label)]

  # Identify edges present at the focal time
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > focal_time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= focal_time

  ## Keep only edges present at focal time and before
  new_edges_df <- all_edges_df[all_edges_df$rootward_test, ]

  ## Adjust length of edge present at focal time
  new_edges_df$initial_length <- tree$edge.length[as.numeric(new_edges_df$edge_ID)]
  new_edges_df$updated_length <- new_edges_df$initial_length
  new_edges_df$updated_length[new_edges_df$tipward_test] <- new_edges_df$rootward_node_age[new_edges_df$tipward_test] - focal_time

  # Record updated edge length
  cut_tree$edge.length <- new_edges_df$updated_length

  ## Update Nnode
  cut_tree$Nnode <- as.integer(nrow(new_edges_df)/2) # Works only for dichotomic tree

  ## Update node ID to form a valid phylo object

  # Nodes = 1:N for tips. N+1 for root. N+2:2N for internal nodes

  # Detect tips node ID as the tipward nodes that are not also the rootward node of any edge
  tips_node_ID <- new_edges_df$tipward_node_ID[!(new_edges_df$tipward_node_ID %in% new_edges_df$rootward_node_ID)]
  # Detect internal nodes ID as rootward nodes that are not the root node
  internal_nodes_ID <- setdiff(new_edges_df$rootward_node_ID, root_node_ID)

  # Get new ID for the tips
  # nodes_ID_df <- data.frame(new_node_ID = 1:(cut_tree$Nnode+1), initial_node_ID = new_edges_df$tipward_node_ID[new_edges_df$tipward_test]) # Only for ultrametric trees.
  nodes_ID_df <- data.frame(new_node_ID = 1:(cut_tree$Nnode+1), initial_node_ID = tips_node_ID)
  # Get new ID for the root
  nodes_ID_df <- rbind(nodes_ID_df, c(nrow(nodes_ID_df)+1, root_node_ID))
  # Get new ID for the internal nodes if any
  if (length(internal_nodes_ID) > 0)
  {
    # internal_nodes_ID_df <- data.frame(new_node_ID = (nrow(nodes_ID_df)+1):(nrow(new_edges_df)+1), initial_node_ID = new_edges_df$tipward_node_ID[!new_edges_df$tipward_test]) # Only for ultrametric trees.
    internal_nodes_ID_df <- data.frame(new_node_ID = (nrow(nodes_ID_df)+1):(nrow(new_edges_df)+1), initial_node_ID = new_edges_df$tipward_node_ID[new_edges_df$tipward_node_ID %in% internal_nodes_ID])
    nodes_ID_df <- rbind(nodes_ID_df, internal_nodes_ID_df)
  }
  # Convert to numeric
  nodes_ID_df$new_node_ID <- as.numeric(nodes_ID_df$new_node_ID)
  nodes_ID_df$initial_node_ID <- as.numeric(nodes_ID_df$initial_node_ID)

  # Store conversion table for nodes ID
  cut_tree$nodes_ID_df <- nodes_ID_df

  ## Add previous internal node ID as initial_nodes_ID
  cut_tree$initial_nodes_ID <- as.character(nodes_ID_df$initial_node_ID[nodes_ID_df$new_node_ID > (cut_tree$Nnode+1)])

  ## Update node.label if present
  if (!is.null(cut_tree$node.label))
  {
    initial_node.label <- cut_tree$node.label
    names(initial_node.label) <- as.character((length(tree$tip.label)+1):(nrow(tree$edge)+1))
    new_node.label <- initial_node.label[cut_tree$initial_nodes_ID]
    cut_tree$node.label <- unname(new_node.label)
  }

  ## Update edge ID to form a valid phylo object
  edges_ID_df <- data.frame(new_edge_ID = 1:nrow(new_edges_df), initial_edge_ID = new_edges_df$edge_ID)
  # Convert to numeric
  edges_ID_df$new_edge_ID <- as.numeric(edges_ID_df$new_edge_ID)
  edges_ID_df$initial_edge_ID <- as.numeric(edges_ID_df$initial_edge_ID)

  # Store conversion table for edges ID
  cut_tree$edges_ID_df <- edges_ID_df

  ## Add previous edge ID as initial_edges_ID
  cut_tree$initial_edges_ID <- new_edges_df$edge_ID

  ## Update tip.label (use the tipward node ID and/or tip.label depending on "keep_tip_labels")

  if (keep_tip_labels) # If keep_tip_labels == T: Use tip.label for edge leading to a unique leaf/tip, and tipward node ID for terminal edges leading to a clade.
  {
    # cut_tree$tip.label <- new_edges_df$tip.label[new_edges_df$tipward_test] # Only for ultrametric trees.
    cut_tree$tip.label <- new_edges_df$tip.label[new_edges_df$tipward_node_ID %in% tips_node_ID]
  } else { # If keep_tip_labels == F: Use tipward node ID for all terminal edges.
    # cut_tree$tip.label <- new_edges_df$tipward_node_ID[new_edges_df$tipward_test] # Only for ultrametric trees.
    cut_tree$tip.label <- new_edges_df$tipward_node_ID[new_edges_df$tipward_node_ID %in% tips_node_ID]
  }

  ## Update $edge

  # Extract edge table using former nodes ID
  cut_tree$edge <- as.matrix(new_edges_df[, c("rootward_node_ID", "tipward_node_ID")])
  dimnames(cut_tree$edge) <- NULL
  # Update nodes ID with the matching table
  cut_tree$edge[ , 1] <- as.integer(nodes_ID_df$new_node_ID[match(x = new_edges_df$rootward_node_ID, table = nodes_ID_df$initial_node_ID)])
  cut_tree$edge[ , 2] <- as.integer(nodes_ID_df$new_node_ID[match(x = new_edges_df$tipward_node_ID, table = nodes_ID_df$initial_node_ID)])

  ## Export output
  return(cut_tree)
}


## Helper function to identify edges present at a given focal time ####

#' @title Identify edges present at a given focal time
#'
#' @description Identify edges of a phylogeny that are present at a given `focal time`.
#'   A branch/edge is consider present if it overlaps with the given `focal time`, or if its tipward node.
#'
#'   For instance, on an ultrametric phylogeny, all terminal tips will be considered present at `focal time = 0`.
#'   At the crown age of the phylogeny, the two descending branches/edges will be considered present.
#'
#' @param phylo Object of class `"phylo"`. The phylogenetic tree must be rooted,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Numeric. The time, in terms of time distance from the present,
#'   at which present edges must be identified.
#' @param tolerance Numeric. Fraction of the total phylogeny depth used to round ages.
#'   This helps avoid floating-point precision errors by treating nearly identical node ages as equal.
#'   Default = `10^-5`.
#'
#' @importFrom phytools nodeHeights
#'
#' @return Returns a data.frame with nine columns.
#'
#'   * `$edge_ID` Integer. ID of the egdes as listed in `phylo$edge`.
#'   * `$rootward_node_ID` Integer. ID of the node found on the rootward end of the edge.
#'   * `$tipward_node_ID` Integer. ID of the node found on the tipward end of the edge.
#'   * `$rootward_node_age` Numeric. Age of the rootward node of the edge, in terms of time distance from the present.
#'   * `$tipward_node_age` Numeric. Age of the tipward node of the edge, in terms of time distance from the present.
#'   * `$tip.label` Character string. Label of the tipward node of the edge.
#'     + For terminal tips, this is the tip label as in `phylo$tip.label`.
#'     + For internal branches, this is the tipward node label as in `phylo$node.label`, or if absent, the `tipward_node_ID.`
#'   * `$rootward_test` Logical. Whether the rootward node of the edge existed before `focal time`.
#'   * `$tipward_test` Logical. Whether the tipward node of the edge existed after or at `focal time`.
#'   * `$edge_present` Logical. Whether the edge was present at `focal time`.
#'
#' @author Maël Doré
#'
#' @noRd
#'

identify_edges_at_focal_time <- function (phylo, focal_time,
                                          tolerance = 10^-5)
{

  # Get node ages per edge (no root edge)
  all_edges_df <- phytools::nodeHeights(phylo)
  root_age <- max(phytools::nodeHeights(phylo)[,2])

  # Define level of tolerance used to round ages
  scaled_tol <- root_age * tolerance
  closest_power <- round(log10(scaled_tol))
  closest_power <- min(closest_power, 0) # Use 0 as the minimal power

  # Round edge ages and assign names
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, -1*closest_power))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- phylo$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  ## Assign labels

  # If tipward node is a tip, use tip.label
  all_edges_df$tip.label <- phylo$tip.label[match(x = all_edges_df$tipward_node_ID, 1:length(phylo$tip.label))]
  # If tipward node is an internal node, use node ID
  all_edges_df$tip.label[is.na(all_edges_df$tip.label)] <- all_edges_df$tipward_node_ID[is.na(all_edges_df$tip.label)]

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- phylo$edge[which.min(phylo$edge[, 1] %in% phylo$edge[, 2]), 1]

  # Identify edges present at the focal time
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > focal_time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= focal_time
  all_edges_df$edge_present <- all_edges_df$rootward_test & all_edges_df$tipward_test

  # Export df
  return(all_edges_df)
}


## Make unit tests for ultrametric (eel.tree / eel_contMap) and non-ultrametric trees (mammals$mammals.phy / mammals_contMap)


## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0


