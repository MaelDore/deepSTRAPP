

## Dependencies: phytools, dplyr

## For continuous traits

# Input = contMap

# contMap <- eel_contMap
# ace = ACE_ln_TL$ace
# tip_data = ln_TL
# focal_time <- 50


extract_most_likely_trait_values_for_focal_time <- function (contMap, ace = NULL, tip_data = NULL, focal_time, update_contMap = F)
{

  ### Check input validity

  # focal_time must be positive and smaller or equal to root age
  # focal_time must be numerical

  # contMap must be a "contMap" class object and have a $maps slot

  # ace must be a named num vector with as many values as their is internal nodes
  # names = node ID
  # names must match node ID and be ordered (add a reordering step and send a warning if not)

  # Same for tip data, but names are tip labels
  # Ensure that all names match with tip labels (order can differ, but send a warning)

  if (is.null(ace))
  {
    warning("No ancestral character estimates (ace) for internal nodes have been provided. Using values interpolated in the contMap instead.\n")
    # message?
  }
  if (is.null(tip_data))
  {
    warning("No tip data have been provided. Using values interpolated in the contMap instead.\n")
    # message?
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
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 2))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- contMap$tree$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

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

    # Return a NULL object
    trait_data <- NULL
    return(trait_data)

    if (update_contMap)
    {
      updated_contMap <- NULL
      return(list(trait_data = trait_data, contMap = updated_contMap))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$time_test, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

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

    ## Update contMap if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_contMap)
    {
      ## Cut tree

      # Use public function cut_phylo_at_focal_time()
      # Need to include validity checks
      # Make it such as it will keep other items in the list, not just create a phylo object with the mandatory items
      # Make it such as it keeps all classes (and add "phylo" if needed)

      ## Update maps

      # Use public function update_maps_at_focal_time()
      # Need to include validity checks

      ## Make all of this another public function: cut_contMap_at_focal_time()
      # Need to include validity checks
      # Make it such as it will keep other items in the list, not just create a contMap object with the mandatory items
      # Make it such as it keeps all classes (and add "contMap" if needed)

    }

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ACE_at_focal_time
    names(trait_data) <- present_edges_df$edge_ID

    ## Export outputs
    if (!update_contMap)
    {
      return(trait_data)

    } else {
      return(list(trait_data = trait_data, contMap = updated_contMap))
    }
  }
}


## Does it make sense to ask only for the ancestral state values and the tree in the workflow? (Technically, no need for the contMap, but needed for plotting)
# Could use the output of	prepare_trait_data() when not asking for stochastic mapping
  # Output = df of trait value x nodes including internal nodes and tips
  # But the rationale of the workflow is to have a contMap for visualization, paralleling the need for stochastic mapping for categorical/biogeographic traits

## Add the option to update the ContMap (not needed for STRAPP test. Useful only for visualization)
  # Provide a cut contMap as output, aside the trait_data vector, with updated maps and tree

