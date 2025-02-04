

# BAMM_posterior_samples_data <- readRDS(file = "../Ponerinae_Historical_Biogeography/outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")
# BAMM_object <- BAMM_posterior_samples_data
# focal_time = 30


update_rates_and_regimes_for_focal_time <- function (BAMM_object, focal_time, update_rates = T, update_regimes = T, verbose = T)
{

  ### Check input validity

  # focal_time must be positive and smaller or equal to root age
  # focal_time must be numerical

  # BAMM_object must be of class "bammdata"
  # BAMM_object must have an $eventData that is a list of N posterior samples
  # BAMM_object must have an $eventData that is a list of N posterior samples with df of regimes
    # Check that all df have 7 columns
  # BAMM_object must have an $tipStates that is a list of N posterior samples with integer vector of regime membership per tips
    # Check that all integer vectors have a length equal to $tip.label
  # BAMM_object must have an $tipLambda that is a list of N posterior samples with integer vector of final speciation rates at tips = current speciation rates
  # BAMM_object must have an $tipMu that is a list of N posterior samples with integer vector of final extinction rates at tips = current extinction rates

  # At least one of "update_rates" and "update_regimes" must be T



  # Add "phylo" class to the BAMM_object
  class(BAMM_object) <- c("bammdata", "phylo")

  ## Identify edges present at focal time

  # Edge, rootward_node, tipward_node, length (once cut)

  # Get node ages per branch (no root edge)
  all_edges_df <- phytools::nodeHeights(BAMM_object)
  root_age <- max(phytools::nodeHeights(BAMM_object)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 2))
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- BAMM_object$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- BAMM_object$edge[which.min(BAMM_object$edge[, 1] %in% BAMM_object$edge[, 2]), 1]

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
  for (i in seq_along(BAMM_object$eventData))
  {
    # i <- 1

    # Extract eventData records = Macroevolutionary regime parameters
    eventData_i <- BAMM_object$eventData[[i]]

    # Compute updated regime age and length
    eventData_i$age <- root_age - eventData_i$time
    eventData_i$updated_length <- eventData_i$age - focal_time

    if (update_regimes)
    {
      ## Identify edge ID per regimes
      ## Loop per regime
      for (j in 1:nrow(eventData_i))
      {
        # j <- 2

        tipward_node_ID_j <- eventData_i$node[j] # Nodes are tipward nodes ID of the branch where the regime starts

        # Get descendant tipward nodes of regime j
        regime_nodes_j <- phytools::getDescendants(tree = BAMM_object, node = tipward_node_ID_j)

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

      # Update tipStates by providing only regimes for tips that are present at the focal time
      tipStates_i <- all_edges_df$regime_ID[all_edges_df$time_test]
      names(tipStates_i) <- all_edges_df$edge_ID[all_edges_df$time_test]
      BAMM_object$tipStates[[i]] <- tipStates_i
    }

    ## If needed, also update tipRates
    if (update_rates)
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

      # Update tipLambda and tipMu by providing only regimes for tips that are present at the focal time
      tipLambda_i <- all_edges_df$tipLambda[all_edges_df$time_test]
      names(tipLambda_i) <- all_edges_df$edge_ID[all_edges_df$time_test]
      BAMM_object$tipLambda[[i]] <- tipLambda_i

      tipMu_i <- all_edges_df$tipMu[all_edges_df$time_test]
      names(tipMu_i) <- all_edges_df$edge_ID[all_edges_df$time_test]
      BAMM_object$tipMu[[i]] <- tipMu_i
    }

    ## Print progress
    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Tip states/rates updated for BAMM posterior sample n\u00B0", i, "/", length(BAMM_object$eventData),"\n"))
    }
  }

  # Update tip labels
  BAMM_object$tip.label <- all_edges_df$edge_ID[all_edges_df$time_test]
  # Inform focal time
  BAMM_object$focal_time <- focal_time

  # Export updated BAMM_object
  return(BAMM_object)
}


## Add the option to update the tree in the BAMM_object (not needed for STRAPP test. Useful only for visualization)
# Update elements of a phylogeny: edge, Nnode, tip.label, edge.length
  # + BAMM stuff (Probably not needed???): begin (Absolute time since root of edge/branch start), end (Absolute time since root of edge/branch end), downseq (?), lastvisit (?), numberEvents (Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts)
  # Update only what is needed to be able to use BAMM plotting function so it become possible to vizualize rates on the cut tree!
  # Use cut_phylo_at_focal_time() to update the phylogeny elements
