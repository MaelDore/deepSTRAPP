### Helper functions to convert BioGeoBEARS BSM outputs into phytools.simmaps ####

#' @title Convert Biogeographic Stochastic Map (BSM) to phytools SIMMAP stochastic map (SM) format
#'
#' @description These functions converts a Biogeographic Stochastic Map (BSM) output from BioGeoBEARS into
#'  a `simmap` object from R package `{phytools}` (See [phytools::make.simmap()]).
#'
#'  They require a model fit with `BioGeoBEARS::bears_optim_run()` and the output of a Biogeographic Stochastic Mapping
#'  performed with `BioGeoBEARS::runBSM()` to produce `simmap` objects as phylogenies with the associated
#'  mapping of range evolution along branches across simulations.
#'
#'   * [deepSTRAPP::convert_BSM_to_simmap()]: Produce one `simmap` for the required simulation (index of the simulation provided with `sim_index`).
#'   * [deepSTRAPP::convert_BSMs_to_simmaps()]: Produce all `simmap` objects for all simulations stored in a unique `multiSimmap` object.
#'
#'  Initial functions in R package BioGeoBEARS by Nicholas J. Matzke:
#'  * `BioGeoBEARS::BSM_to_phytools_SM()`
#'  * `BioGeoBEARS::BSMs_to_phytools_SMs()`
#'
#' # Notes on BioGeoBEARS
#'
#'  The R package `BioGeoBEARS` is needed for this function to work with biogeographic data.
#'  Please install it manually from: \href{https://github.com/nmatzke/BioGeoBEARS}{https://github.com/nmatzke/BioGeoBEARS}.
#'
#' # Notes on using the resulting simmap object in phytools (adapted from Nicholas J. Matzke)
#'
#'  The phytools functions, like [phytools::countSimmap()], will only count the anagenetic events
#'  (range transitions occurring along branches) as it was written assuming purely anagenetic models.
#'
#'  It remains possible to extract cladogenetic events (range transitions occurring at speciation)
#'  by comparing the last-state-below-a-node with the descendant-pairs-above-a-node.
#'  However, it is recommended to use the built-in functions from BioGeoBEARS to summarize
#'  the biogeographic history based on the tables of cladogenetic and anagenetic events obtained
#'  from `BioGeoBEARS::runBSM()`. `simmap` objects should primarily be considered as a tool for visualization.
#'
#'  Associated functions in R package `BioGeoBEARS`:
#'
#'  * `BioGeoBEARS::simulate_source_areas_ana_clado()`: To select randomly a unique area source for transition from a multi-area state to a single area.
#'  * `BioGeoBEARS::get_dmat_times_from_res()`: To generate matrices of range expansion from source area to destination area.
#'  * `BioGeoBEARS::count_ana_clado_events()`: To count the number and type of events from BSM tables.
#'  * `BioGeoBEARS::hist_event_counts()`: To plot histograms of event counts across BSM tables.
#'
#'  Please note carefully that area-to-area dispersal events are not identical with the state transitions.
#'  For example, a state can be a geographic range with multiple areas, but under the logic of DEC-type models,
#'  a range-expansion event like  ABC->ABCD actually means that a dispersal happened from some specific area (A, B, or C)
#'  to the new area. BSMs track this area-to-area sourcing in its cladogenetic and anagenetic event tables,
#'  at least if `BioGeoBEARS::simulate_source_areas_ana_clado()` has been run on the output of `BioGeoBEARS::runBSM()`.
#'
#' @param model_fit A BioGeoBEARS results object, produced by ML inference via `BioGeoBEARS::bears_optim_run()`.
#' @param phylo Time-calibrated phylogeny used in the BioGeoBEARS analyses to produce the historical biogeographic inference
#'  and run the Biogeographic Stochastic Mapping. Object of class `"phylo"` as defined in `{ape}`.
#' @param BSM_output A list with two objects, a cladogenetic events table and an anagenetic events table, as the result of
#'   Biogeographic Stochastic Mapping conducted with `BioGeoBEARS::runBSM()`.
#' @param sim_index Integer. Index of the biogeographic simulation targeted to produce the `simmap` with [deepSTRAPP::convert_BSM_to_simmap()].
#'
#' @export
#' @importFrom plyr ldply
#'
#' @details These functions are slight adaptations of original functions from the R Package `BioGeoBEARS` by N. Matzke.
#'
#'  Initial functions: `BioGeoBEARS::BSM_to_phytools_SM()` `BioGeoBEARS::BSMs_to_phytools_SMs()`
#'
#'  Changes:
#'  * Solves issue with differences in ranges allowed across time-strata.
#'  * Requires directly the output of `BioGeoBEARS::runBSM()` instead of separated cladogenetic and anagenetic event tables.
#'  * Update the documentation.
#'
#' @return The [deepSTRAPP::convert_BSM_to_simmap()] function returns a list with two elements:
#'   * `$simmap` A unique `simmap` for a given biogeographic simulation as an object of classes `c("simmap", "phylo")`.
#'     This is a modified `{ape}` tree with additional elements to report range mapping, model parameters and likelihood.
#'     - `$maps` A list of named numerical vectors. Provides the mapping of ranges along each remaining edge.
#'       Names are the ranges. Values are residence times in each state across segments
#'     - `$mapped.edge` A numerical matrix. Provides the evolutionary time spent across ranges (columns) along the edges (rows).
#'       row.names() are the node ID at the rootward and tipward ends of each edge.
#'     - `$Q` Numerical matrix. The transition rates across ranges calculated from the ML parameter estimates of the model.
#'     - `$logL` Numeric. The log-likelihood of the data under the ML model.
#'   * `$residence_times` Data.frame with two rows. Summarizes the residence time spent in each range along all branches,
#'    in (raw) evolutionary time (i.e., branch lengths), and in percentage (perc).
#'
#' The [deepSTRAPP::convert_BSMs_to_simmaps()] function loop around the [deepSTRAPP::convert_BSM_to_simmap()] function to aggregate all `simmaps`
#'   from all biogeographic simulations in a unique list of classes `c("multiSimmap", "multiPhylo")`.
#'   * Each element in the `$simmap` of a biogeographic simulation obtained with [deepSTRAPP::convert_BSM_to_simmap()].
#'   * `$residence_times` summary data.frames are not preserved.
#'
#' @seealso [phytools::countSimmap()] [phytools::make.simmap()]
#' `BioGeoBEARS::simulate_source_areas_ana_clado()` `BioGeoBEARS::get_dmat_times_from_res()`
#' `BioGeoBEARS::count_ana_clado_events()` `BioGeoBEARS::hist_event_counts()`
#'
#' @author Nicholas J. Matzke. Contact: \email{matzke@@berkeley.edu}
#' @author Changes by Maël Doré (see Details)
#'
#' @references For BioGeoBEARS: Matzke, Nicholas J. (2018). BioGeoBEARS: BioGeography with Bayesian (and likelihood) Evolutionary Analysis with R Scripts.
#'    version 1.1.1, published on GitHub on November 6, 2018. \doi{10.5281/zenodo.1478250}. Website: \url{http://phylo.wikidot.com/biogeobears}.
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  ## Run only if you have R package 'BioGeoBEARS' installed.
#'  # Please install it manually from: https://github.com/nmatzke/BioGeoBEARS")
#'
#'  ## Load phylogeny and tip data
#'  library(phytools)
#'  data(eel.tree)
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Load directly output of prepare_trait_data() run on biogeographic data
#'  data(eel_biogeo_data, package = "deepSTRAPP")
#'
#'  ## Convert BSM output into a unique simmap, including residence times
#'  simmap_1 <- convert_BSM_to_simmap(model_fit = eel_biogeo_data$best_model_fit,
#'                                     phylo = eel.tree,
#'                                     BSM_output = eel_biogeo_data$BSM_output,
#'                                     sim_index = 1)
#'  # Explore output
#'  str(simmap_1, max.level = 1)
#'  # Print residence times in each range
#'  simmap_1$residence_times
#'  # Plot simmap
#'  plot(simmap_1$simmap)
#'
#'  ## Convert BSM output into all simmaps in a multiSimmap/multiPhylo object
#'  all_simmaps <- convert_BSMs_to_simmaps(model_fit = eel_biogeo_data$best_model_fit,
#'                                         phylo = eel.tree,
#'                                         BSM_output = eel_biogeo_data$BSM_output)
#'  # Explore output
#'  str(all_simmaps, max.level = 1)
#'  # Plot simmap n°1
#'  plot(all_simmaps[[1]]) }
#' }
#'

## importFrom for BioGeoBEARS
# importFrom BioGeoBEARS get_Qmat_COOmat_from_BioGeoBEARS_run_object get_sum_statetime_on_branch prt


convert_BSM_to_simmap <- function(model_fit, phylo, BSM_output, sim_index)
{
  ## Control for BioGeoBEARS install
  if (!requireNamespace("BioGeoBEARS", quietly = TRUE))
  {
    stop("Package 'BioGeoBEARS' is needed to work with biogeographic data.
       Please install it manually from: https://github.com/nmatzke/BioGeoBEARS;
       or from the alterative deepSTRAPP repository (https://maeldore.github.io/drat).
       For instructions, see the Dependencies section on deepSTRAPP homepage (https://github.com/MaelDore/deepSTRAPP).")
  }

  # Extract the tables of cladogenetic and anagenetic events
  clado_events_table <- BSM_output$RES_clado_events_tables[[sim_index]]
  ana_events_table <- BSM_output$RES_ana_events_tables[[sim_index]]

  run_parent_brs_TF = TRUE
  tr <- phylo

  # Error check
  if ((inherits(ana_events_table, 'data.frame')) && (dim(ana_events_table)[1] > 0))
  {
    run_parent_brs_TF = TRUE
  } else {
    ana_events_table = NA
    run_parent_brs_TF = FALSE
  }

  if (is.null(ana_events_table) == TRUE)
  {
    ana_events_table = NA
    run_parent_brs_TF = FALSE
  }

  # Is it time-stratified?
  stratTF = (length(model_fit$inputs$timeperiods) > 0)

  returned_mats = BioGeoBEARS::get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object = model_fit$inputs, include_null_range = model_fit$inputs$include_null_range)
  returned_mats
  areanames = returned_mats$areanames
  ranges_list = returned_mats$ranges_list

  # Fix for issue with using only reduced range list
  reduced_ranges_list = returned_mats$ranges_list
  max_range_size <- max(nchar(reduced_ranges_list))
  ranges_list <- generate_list_ranges(areas_list = areanames, max_range_size = max_range_size, include_null_range = model_fit$inputs$include_null_range)

  # Get list of edges:
  trtable = BioGeoBEARS::prt(tr, printflag=FALSE)
  trtable

  # View(trtable)
  # View(clado_events_table)

  # Convert pruningwise edge numbers to cladewise edgenums
  # The clado_events_table (non-stratified) has the edge numbers in
  # pruningwise order
  if (stratTF == FALSE)
  {
    pruningwise_edgenums = clado_events_table$parent_br
    cladewise_edgenums = trtable$parent_br
    translation_pruning_to_clade_edgenums = as.data.frame(cbind(pruningwise_edgenums, cladewise_edgenums), stringsAsFactors=FALSE)
    translation_pruning_to_clade_edgenums

    # Error trap for when there are no anagenetic events
    if (run_parent_brs_TF == TRUE)
    {
      ana_events_edgenums_indexes_in_clado_events_table = match(x = ana_events_table$parent_br, table = clado_events_table$parent_br)
      ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[ana_events_edgenums_indexes_in_clado_events_table]
    }
    clado_events_table$parent_br = trtable$parent_br
    # } else {
    #   # If same issue with time-stratified clado_events_table, also have duplicates due to branches split across time-strata!
    #   pruningwise_edgenums = unique(clado_events_table$parent_br)
    #   cladewise_edgenums = trtable$parent_br
    #   translation_pruning_to_clade_edgenums = as.data.frame(cbind(pruningwise_edgenums, cladewise_edgenums), stringsAsFactors=FALSE)
    #   translation_pruning_to_clade_edgenums
    #   # Error trap for when there are no anagenetic events
    #   if (run_parent_brs_TF == TRUE)
    #   {
    #     # ana_events_edgenums_indexes_in_clado_events_table = match(x = ana_events_table$parent_br, table = clado_events_table$parent_br)
    #     # ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[ana_events_edgenums_indexes_in_clado_events_table]
    #     ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[match(x = ana_events_table$parent_br, table = translation_pruning_to_clade_edgenums$pruningwise_edgenums)]
    #   }
    #   clado_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[match(x = clado_events_table$parent_br, table = translation_pruning_to_clade_edgenums$pruningwise_edgenums)]
  }


  # Get the edgenums, exclude the "NA" for the root branch
  # Order edgenums from smallest to largest
  nonroot_TF = trtable$node.type != "root"
  edgenums = trtable$parent_br[nonroot_TF]
  edgenums_order = order(edgenums)
  edgenums = edgenums[edgenums_order]

  # Get the "ancestor node, descendant node"
  ancnodenums = trtable$ancestor[nonroot_TF][edgenums_order]
  decnodenums = trtable$node[nonroot_TF][edgenums_order]

  # rownames for mapped.edge
  rownames_for_mapped_edge = paste0(ancnodenums, ",", decnodenums)
  rownames_for_mapped_edge

  # instantiate "maps" for phytools (a list, with array of state residence times
  maps = list()

  if (stratTF == TRUE)
  {
    time_tops = sort(unique(clado_events_table$time_top))
    time_bots = sort(unique(clado_events_table$time_bot))
  }


  # Loop through the edges, record any anagenetic events on the branches
  # fill in "maps"
  for (i in 1:length(edgenums))
  {
    # i <- 410
    edgenum = edgenums[i]

    # Trap for if ana_events_table is NA (common, if there are no anagenetic events in
    # the tree at all)
    if ( (length(ana_events_table) == 1) && (is.na(ana_events_table)) )
    {
      edgefound_TF = FALSE
    } else {
      edgefound_TF = ana_events_table$parent_br == edgenum
    }

    # If no anagenetic events are found, the whole branchlength is in the
    # starting state, as specified in "clado_events_table"
    if (sum(edgefound_TF) == 0)
    {
      # The states should be the same at the branch bottom and top:
      clado_row_TF = clado_events_table$parent_br == edgenum
      clado_row_TF[is.na(clado_row_TF)] = FALSE

      # 			if (stratTF == TRUE)
      # 				{
      # 				match_time_tops_TF = as.numeric(tmptable_rows$abs_event_time[nnr]) >= as.numeric(trtable$time_top)
      # 				match_time_bots_TF = as.numeric(tmptable_rows$abs_event_time[nnr]) < as.numeric(trtable$time_bot)
      # 				}


      ## NJM 2019-03-12_ fix: doubles can be found in time-strat, FIX
      clado_events_table[clado_row_TF,]


      # Error check
      if (sum(clado_row_TF) < 1)
      {
        txt = paste0("STOP ERROR #1 in BSM_to_phytools_SM(): 0 rows in clado_events_table match edge number/branch number (parent_br==", edgenum, ").")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }

      if (sum(clado_row_TF) == 1)
      {
        bottom_state_num_1based = clado_events_table$sampled_states_AT_brbots[clado_row_TF]
        top_state_num_1based = clado_events_table$sampled_states_AT_nodes[clado_row_TF]

        # Error handle for the root
        if (clado_events_table$node.type[clado_row_TF] == "root")
        {
          bottom_state_num_1based <- top_state_num_1based
        }
      }

      if (sum(clado_row_TF) > 1)
      {
        bottom_state_num_1based = unique(clado_events_table$sampled_states_AT_brbots[clado_row_TF])
        top_state_num_1based = unique(clado_events_table$sampled_states_AT_nodes[clado_row_TF])

        # Error check: there should be only 1 unique state corresponding to this
        # node (because we are in the section where no anagenetic histories were
        # found on the branch).
        if (length(top_state_num_1based) != 1)
        {
          txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'top_state_num_1based' corresponding to the node specified. Printing the matching rows of 'clado_events_table'."
          cat("\n\n")
          cat(txt)
          cat("\n\n")

          print("clado_events_table[clado_row_TF,]")
          print(clado_events_table[clado_row_TF,])

          print("top_state_num_1based")
          print(top_state_num_1based)

          stop(txt)
        } # END if (length(top_state_num_1based) != 1)

        if (length(bottom_state_num_1based) != 1)
        {
          txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'bottom_state_num_1based' corresponding to the node specified. Printing the matching rows of 'clado_events_table'."
          cat("\n\n")
          cat(txt)
          cat("\n\n")

          print("clado_events_table[clado_row_TF,]")
          print(clado_events_table[clado_row_TF,])

          print("bottom_state_num_1based")
          print(bottom_state_num_1based)

          stop(txt)
        } # END if (length(bottom_state_num_1based) != 1)
      } # END if (sum(clado_row_TF) > 1)

      # Error check
      if (bottom_state_num_1based != top_state_num_1based)
      {
        txt = paste0("STOP ERROR #2 in BSM_to_phytools_SM(): the top_state_num_1based (", top_state_num_1based, "), and bottom_state_num_1based (", bottom_state_num_1based, ") have to match at edge number/branch number (parent_br==", edgenum, "), because no anagenetic events were recorded on this branch.")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }

      # No events detected, so put in the states_array just one state
      # But, since there are NO events on this branch, don't add to it if
      # there is already something there (e.g. don't add the same branch for
      # multiple branch segments)
      names_of_states_array = c(ranges_list[bottom_state_num_1based])
      times_in_each_state_array = unique(c(clado_events_table$edge.length[clado_row_TF]))

      if (length(times_in_each_state_array) != 1)
      {
        txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'times_in_each_state_array' corresponding to the node specified. There should only be one time, because no anagenetic events were detected on this branch. Printing the matching rows of 'clado_events_table'."
        cat("\n\n")
        cat(txt)
        cat("\n\n")

        print("clado_events_table[clado_row_TF,]")
        print(clado_events_table[clado_row_TF,])

        print("times_in_each_state_array")
        print(times_in_each_state_array)

        stop(txt)
      } # END if (length(times_in_each_state_array) != 1)


      names(times_in_each_state_array) = names_of_states_array
      times_in_each_state_array
    } # END if (sum(edgefound_TF) == 0)

    # If some anagenetic events are found, the whole branchlength is in the
    # starting state, as specified in "clado_events_table"
    if (sum(edgefound_TF) > 0)
    {
      # The states will be listed in the ana_events_table
      rows_matching_edgenum_TF = ana_events_table$parent_br == edgenum
      tmp_ana_events_table = ana_events_table[rows_matching_edgenum_TF,]

      # Make sure the tmp_ana_events_table is sorted by REVERSE event_time (along branch)
      # (Do REVERSE, because older events have larger event ages)
      tmp_ana_events_table = tmp_ana_events_table[rev(order(tmp_ana_events_table$abs_event_time)),]
      tmp_ana_events_table

      numevents = sum(rows_matching_edgenum_TF)

      # 1 event, 2 states on branch
      if (numevents == 1)
      {
        first_state_name = c(tmp_ana_events_table$current_rangetxt[1])
        abs_time_bp_at_branch_top = tmp_ana_events_table$time_bp[1]
        abs_time_bp_at_branch_bot = abs_time_bp_at_branch_top + tmp_ana_events_table$edge.length[1]
        first_state_timelength = c(abs_time_bp_at_branch_bot - tmp_ana_events_table$abs_event_time[1])

        # The rest of the branch is the 2nd state
        further_state_name = c(tmp_ana_events_table$new_rangetxt[1])
        further_state_time = c(tmp_ana_events_table$edge.length[1] - first_state_timelength)

        times_in_each_state_array = c(first_state_timelength, further_state_time)
        names_of_states_array = c(first_state_name, further_state_name)
      } # END if (numevents == 1)

      # 2+ events, 3+ states on branch
      if (numevents >= 2)
      {
        first_state_name = c(tmp_ana_events_table$current_rangetxt[1])
        #first_state_time = c(tmp_ana_events_table$event_time[1])
        abs_time_bp_at_branch_top = tmp_ana_events_table$time_bp[1]
        abs_time_bp_at_branch_bot = abs_time_bp_at_branch_top + tmp_ana_events_table$edge.length[1]
        first_state_timelength = c(abs_time_bp_at_branch_bot - tmp_ana_events_table$abs_event_time[1])

        nonfirst_rows = 2:numevents
        nonlast_rows = 1:(numevents-1)

        further_state_names = c(tmp_ana_events_table$new_rangetxt)
        further_state_times = tmp_ana_events_table$abs_event_time[nonlast_rows] - tmp_ana_events_table$abs_event_time[nonfirst_rows]

        # How much of the branch is left
        last_time = c(tmp_ana_events_table$edge.length[numevents] - sum(c(first_state_timelength,further_state_times)))

        times_in_each_state_array = c(first_state_timelength, further_state_times, last_time)
        names_of_states_array = c(first_state_name, further_state_names)
      } # END if (numevents >= 2)

      names(times_in_each_state_array) = names_of_states_array
    } # END if (sum(edgefound_TF) > 0)

    # Store
    maps[[i]] = times_in_each_state_array
  } # END for (i in 1:length(edgenums))

  maps


  # Check that the sums of state residence times add up to the branch lengths
  simmap_times <- unlist(lapply(X=maps, FUN=sum))
  tree_times <- tr$edge.length
  all(simmap_times == tree_times)
  cbind(simmap_times, tree_times)

  # Check if issue with order of edges
  simmap_times <- simmap_times[order(simmap_times)]
  tree_times <- tree_times[order(tree_times)]
  cbind(simmap_times, tree_times)

  # Compute residence times per states
  edge_states_matrix <- plyr::ldply(.data = maps, .fun = base::rbind)
  state_times <- apply(X = edge_states_matrix, MARGIN = 2, FUN = sum, na.rm = T)
  state_times_reordered <- state_times[ranges_list]
  names(state_times_reordered) <- ranges_list
  state_times_reordered[is.na(state_times_reordered)] <- 0

  # Add total
  state_times_df <- as.data.frame(t(state_times_reordered))
  state_times_df$total <- sum(state_times_reordered)

  # Add percentage of times
  state_times_perc <- state_times_df / state_times_df$total * 100
  state_times_df <- rbind(state_times_df, state_times_perc)
  row.names(state_times_df) <- c("raw", "perc")

  # # Remove empty states
  # state_times_df <- state_times_df[, state_times_perc != 0]

  # Make the mapped.edge output
  mapped.edge = matrix(data=0.0, nrow=length(edgenums), ncol=length(ranges_list))
  row.names(mapped.edge) = rownames_for_mapped_edge
  colnames(mapped.edge) = ranges_list
  mapped.edge

  # For each branch,
  # 1. Get the list of observed states
  i = 3
  observed_states = sort(unique(names(maps[[i]])))
  observed_states


  # sapply to get the sum of each
  sapply(X = observed_states, FUN = BioGeoBEARS::get_sum_statetime_on_branch, branch_history_map = maps[[i]])

  # Fill in the table for each branch
  for (i in 1:nrow(mapped.edge))
  {
    observed_states = sort(unique(names(maps[[i]])))
    total_residence_times = sapply(X = observed_states, FUN = BioGeoBEARS::get_sum_statetime_on_branch, branch_history_map = maps[[i]])
    names_observed_states = names(total_residence_times)
    mapped.edge[i,names_observed_states] = unname(total_residence_times)
  }
  mapped.edge


  # Get the transition matrix and logL
  Q = returned_mats$Qmat
  # row.names(Q) = ranges_list
  # colnames(Q) = ranges_list
  row.names(Q) = reduced_ranges_list
  colnames(Q) = reduced_ranges_list
  logL = model_fit$total_loglikelihood

  tr_wSimmap = tr
  tr_wSimmap$maps = maps
  tr_wSimmap$mapped.edge = mapped.edge
  tr_wSimmap$Q = Q
  tr_wSimmap$logL = logL
  class(tr_wSimmap) = c("simmap", "phylo")
  tr_wSimmap

  tr_wSimmap$maps
  tr_wSimmap$mapped.edge

  return(list(simmap = tr_wSimmap, residence_times = state_times_df))
}

### Same but for list of multiple BSM maps

#'
#' @rdname convert_BSM_to_simmap
#'

convert_BSMs_to_simmaps <- function(model_fit, phylo, BSM_output)
{
  # Initiate final output
  simmaps_list = list()

  for (i in 1:length(BSM_output$RES_clado_events_tables))
  {
    # i <- 1
    simmaps_list[[i]] = convert_BSM_to_simmap(model_fit = model_fit, phylo = phylo, BSM_output = BSM_output, sim_index = i)$simmap
  }

  class(simmaps_list) = c("multiSimmap", "multiPhylo")
  return(simmaps_list)
}
