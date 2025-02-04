#' @title Cut the phylogeny at a given time in the past
#'
#' @description Cuts off all the branches of the phylogeny which are
#'   younger than a specific time in the past (i.e. the `focal_time`).
#'   Branches overlapping the `focal_time` are shorten to the `focal_time`.
#'
#' @param tree Object of class `"phylo"`. The phylogenetic tree does not need
#'   to be ultrametric and fully dichotomous.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree must be cut.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip must retained their initial `tip.label`. Default is `TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights
#' @importFrom ape nodelabels
#' @importFrom ape edgelabels
#'
#' @details When an entire lineage is cut (i.e. one or more nodes along a path) and \code{keep.lineages = TRUE},
#'   the leaves left are labeled as "l" followed by a number.
#'
#' @return The function returns the cut phylogeny as an object of class `"phylo"`. It adds multiple useful items to the object.
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @examples
#' # Load eel phylogeny from the R package phytools
#' library(phytools)
#' data(eel.tree)
#'
#' # ----- Example 1: keep_tip_labels = T ----- #
#'
#' # Cut tree to 30 Mya while keeping tip.label on terminal branches with a unique descending tip.
#' cut_eel.tree <- cut_phylo_at_focal_time(tree = eel.tree, focal_time = 30, keep_tip_labels = T)
#'
#' # Plot internal node labels on initial tree with cut-off
#' plot(eel.tree)
#' abline(v = max(phytools::nodeHeights(eel.tree)[,2]) - 30, col = "red", lty = 2, lwd = 2)
#' nodelabels_in_cut_tree <- (length(eel.tree$tip.label)+1):(length(eel.tree$tip.label) + eel.tree$Nnode)
#' nodelabels_in_cut_tree[!(nodelabels_in_cut_tree %in% cut_eel.tree$initial_nodes_ID)] <- NA
#' ape::nodelabels(text = nodelabels_in_cut_tree)
#'
#' # Plot initial internal node labels on cut tree
#' plot(cut_tree)
#' ape::nodelabels(text = cut_tree$initial_nodes_ID)
#'
#' # Plot edge labels on initial tree with cut-off
#' plot(eel.tree)
#' abline(v = max(phytools::nodeHeights(eel.tree)[,2]) - 30, col = "red", lty = 2, lwd = 2)
#' edgelabels_in_cut_tree <- 1:nrow(eel.tree$edge)
#' edgelabels_in_cut_tree[!(1:nrow(eel.tree$edge) %in% cut_eel.tree$initial_edges_ID)] <- NA
#' ape::edgelabels(text = edgelabels_in_cut_tree)
#'
#' # Plot initial edge labels on cut tree
#' plot(cut_eel.tree)
#' ape::edgelabels(text = cut_tree$initial_edges_ID)
#'
#' # ----- Example 2: keep_tip_labels = F ----- #
#'
#' # Cut tree to 30 Mya without keeping tip.label on terminal branches with a unique descending tip.
#' # All tip.labels are converted to their descending/tipward node ID
#' cut_eel.tree <- cut_phylo_at_focal_time(tree = eel.tree, focal_time = 30, keep_tip_labels = F)
#' plot(cut_eel.tree)
#'


### Try with non-dichotomous trees!!!


library(phytools)
library(ape)
data(package = "phytools")
data(bat.tree)
plot(bat.tree)

is.binary(bat.tree)

cut_bat.tree <- cut_phylo_at_focal_time(tree = bat.tree, focal_time = 0.3)

cut_bat.tree$edge
cut_bat.tree$tip.label

# Plot initial internal node labels on cut treeplot(cut_bat.tree)
plot(cut_bat.tree)
ape::nodelabels(text = cut_bat.tree$initial_nodes_ID)

# Plot edge labels on initial tree with cut-off
plot(bat.tree)
abline(v = max(phytools::nodeHeights(bat.tree)[,2]) - 0.05, col = "red", lty = 2, lwd = 2)
edgelabels_in_cut_bat.tree <- 1:nrow(bat.tree$edge)
edgelabels_in_cut_bat.tree[!(1:nrow(bat.tree$edge) %in% cut_bat.tree$initial_edges_ID)] <- NA
ape::edgelabels(text = edgelabels_in_cut_bat.tree)

# Plot initial edge labels on cut tree
plot(cut_bat.tree)
ape::edgelabels(text = cut_bat.tree$initial_edges_ID)

# ----- Example 2: keep_tip_labels = F ----- #

# Cut tree to 30 Mya without keeping tip.label on terminal branches with a unique descending tip.
# All tip.labels are converted to their descending/tipward node ID
cut_bat.tree <- cut_phylo_at_focal_time(tree = bat.tree, focal_time = 0.05, keep_tip_labels = F)
plot(cut_bat.tree)



cut_phylo_at_focal_time <- function(tree, focal_time, keep_tip_labels = T)
{
  ### Check input validity

  # tree must be a "phylo" class object
  # tree must be rooted
  # tree must be fully resolved/dichotomic
    # ape::is.binary

  # focal_time must be positive and smaller to root age
  # If focal_time = root_age, throw a specific error
  # focal_time must be numerical

  # keep_tip_labels must be TRUE or FALSE


  ## Initiate new tree
  cut_tree <- tree

  ## Identify edges present at focal time

  # Get node ages per edge (no root edge)
  all_edges_df <- phytools::nodeHeights(tree)
  root_age <- max(phytools::nodeHeights(tree)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 2))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Inform root_age
  cut_tree$root_age <- root_age

  # Get nodes ID per edge
  all_edges_ID_df <- tree$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

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
  cut_tree$Nnode <- nrow(new_edges_df)/2 # Works only for dichotomic tree

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
  cut_tree$edge[ , 1] <- nodes_ID_df$new_node_ID[match(x = new_edges_df$rootward_node_ID, table = nodes_ID_df$initial_node_ID)]
  cut_tree$edge[ , 2] <- nodes_ID_df$new_node_ID[match(x = new_edges_df$tipward_node_ID, table = nodes_ID_df$initial_node_ID)]

  # Export output
  return(cut_tree)

}

cut_tree <- cut_phylo_at_focal_time(tree = tree, focal_time = 0, keep_tip_labels = F)


