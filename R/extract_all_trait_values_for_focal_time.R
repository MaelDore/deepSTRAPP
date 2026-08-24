## Functions to extract all trait data from multiple stochastic maps at a given focal time
# One master function to select the proper pipeline according to data type
# Three sub-functions extracting trait evolution according to data type

### Master function to select the proper sub-function according to data type ####

#' @title Extract all trait data from stochastic maps at a given time in the past
#'
#' @description Extracts all trait values/states/ranges found along branches
#'   of multiple stochastic maps at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogenies (`contMaps` or `densityMaps`)
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$maps` and `$mapped.edge` elements.
#'
#' @param contMaps For continuous trait data. List of objects of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   each containing a phylogenetic tree and associated continuous trait mapping that
#'   represents an independent evolutionary history of ancestral trait evolution,
#'   conditioned to observed trait data and model fit (i.e., stochastic maps).
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param densityMaps For categorical trait or biogeographic data. List of objects of class `"densityMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state/range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state/range. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param simmaps For categorical trait or biogeographic data.
#'   List of objects of classes `"phylo"` and `"simmap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()] or [phytools::make.simmap()],
#'   that represent discrete character/geographic evolutionary history
#'   (i.e., transitions in character states/geographic ranges) mapped along branches.
#'   This is needed to be able to track which simulated history provided which trait data in downstream analyses
#'   employing the "paired" or "full" strategies to account for uncertainty in ancestral trait estimates.
#' @param nb_simulations For categorical trait or biogeographic data. Integer. The number of stochastic maps used to simulate
#'   trait evolution. This is needed for the "paired" and "full" strategies to account for trait estimate uncertainty,
#'   if only densityMaps summarizing posterior state/range density are provided, but not the simmaps representing all evolutionary histories.
#' @param tip_data (Optional) Named vector of tip values of the trait.
#'   * For continuous trait data: Named numerical vector of trait values.
#'   * For categorical trait or biogeographic data: Character string vector of states/ranges
#'   Names are nodes_ID of the internal nodes. Ensure accurate tip values are used.
#' @param trait_data_type Character string. Specify the type of trait data. Must be one of "continuous", "categorical", "biogeographic".
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_Map Logical. Specify whether the mapped phylogeny (`contMap`, `densityMaps`, and/or `simmaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated contMap. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.contMap
#' @importFrom ape nodelabels
#'
#' @details The mapped phylogenies (`contMaps`, `densityMaps`, or `simmaps`) are cut at a specific time in the past
#'   (i.e. the `focal_time`) and the associated trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   For continuous trait data:
#'
#'   Simulated trait data are extracted from `contMaps` (i.e. continuous stochastic maps).
#'
#'   For categorical trait and biogeographic data:
#'
#'   If `simmaps` are provided as input, all states/ranges are extracted directly from the stochastic mapes.
#'
#'   If only `densityMaps` are provided, the posterior probabilities of states/ranges are extracted.
#'   Posterior probabilities are multiplied by the `nb_simulations` to record the distribution of states/ranges
#'   across dummy stochastic maps and assigned to each tip and cut branches at `focal_time`.
#'
#'   True tip states/ranges will be used if `tip_data` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `contMap`/`densityMaps`/`simmaps` -----
#'
#'   To obtain an updated `contMap`/`densityMaps`/`simmaps` alongside the trait data, set `update_Map = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The mapping in `contMap`/`densityMaps`/`simmaps` (`$maps` and `$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#'   A specific sub-function (that can be used independently) is called according to the type of trait data and inputs provided:
#'   * For continuous traits: [deepSTRAPP::extract_all_trait_values_from_contMaps_for_focal_time()]
#'   * For categorical traits: [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()] and [deepSTRAPP::extract_all_states_from_simmaps_for_focal_time()]
#'   * For biogeographic ranges: [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()] and [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#'   To extract only the most likely trait value/state/range, see this associated function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A list of named numerical vector with simulated trait values found along branches overlapping the `focal_time` across each of the stochastic maps.
#'     Names are the tip.label/tipward node ID. ID of the associated stochastic maps associated with each item of the list can only be provided if simmaps are provided as inputs.
#'     If only densityMaps are provided, the distribution of states/ranges is associated with dummy maps_ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous", "categorical", or "biogeographic". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_Map = TRUE`, the output is a list with four to five elements: `$trait_data`, `$focal_time`, `$trait_data_type`, `$contMaps` or `$densityMaps`, and/or `$simmaps`.
#'
#'   For continuous trait data:
#'
#'   * `$contMaps` A list of objects of class `"contMap"` that contains the updated `contMaps` with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$contMap$tree` elements.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#'  For categorical trait and biogeographic data:
#'
#'   * `$densityMaps` A list of objects of class `"densityMap"` that contains the updated `densityMap` of each state/range,
#'      with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$densityMaps$tree` elements.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#'   * `$simmaps` A list of objects with the class `"simmap"` that contains the updated `simmap` = stochastic maps
#'      representing discrete character/geographic evolutionary history,
#'      with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to each `simmap`.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMaps_for_focal_time()] [deepSTRAPP::cut_densityMaps_for_focal_time()] [deepSTRAPP::cut_simmaps_for_focal_time()]
#'
#' Equivalent function to extract the most likely trait value/state/range from mapped phylogenies:
#' [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
#' Associated sub-functions per type of trait data:
#'
#' [deepSTRAPP::extract_all_trait_values_from_contMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_states_from_simmaps_for_focal_time()]
#' [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Continuous trait ----- #
#'
#' if (deepSTRAPP::is_dev_version())
#' {
#'  ## The R package 'contsimmap' is needed for this example to work.
#'  # Please install it manually from: https://github.com/bstaggmartin/contsimmap.
#'
#'  ## Prepare data
#'
#'  # Load eel data from the R package phytools
#'  # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#'  library(phytools)
#'  data(eel.tree)
#'  data(eel.data)
#'
#'  # Extract body size
#'  eel_data <- setNames(eel.data$Max_TL_cm,
#'                       rownames(eel.data))
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Map continuous trait evolution on the phylogeny
#'  eel_cont_data  <- prepare_trait_data(
#'     tip_data = eel_data,
#'     trait_data_type = "continuous",
#'     phylo = eel.tree,
#'     # evolutionary_models = c("BM", "OU", "lambda", "kappa"),
#'     evolutionary_models = c("BM"),
#'     # Perform stochastic mapping to obtain multiple evolutionary histories
#'     run_stochastic_maps = TRUE,
#'     nb_simulations = 100,
#'     verbose = TRUE)
#'
#'  # Extract continuous stochastic maps (contMaps)
#'  eel_contMaps <- eel_cont_data$contMaps
#'
#'  # Plot the interpolated map of ML estimates
#'  plot_contMap(contMap = eel_cont_data$contMap)
#'
#'  # Plot several continuous stochastic maps (contMaps)
#'  plot_contMap(contMap = eel_contMaps[[1]])
#'  plot_contMap(contMap = eel_contMaps[[10]])
#'  plot_contMap(contMap = eel_contMaps[[100]])
#'
#'  # Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  ## Extract all trait values for focal time = 10 Mya
#'  extract_trait_data_10My <- extract_all_trait_values_from_contMaps_for_focal_time(
#'     contMaps = eel_contMaps, focal_time = 10, update_contMaps = T)
#'
#'  # Convert in data.frame
#'   # Rows = Stochastic maps
#'   # Columns = Cut branches at 10 Mya
#'  trait_data_df <- as.data.frame(do.call(rbind, extract_trait_data_10My$trait_data))
#'  head(trait_data_df)
#'
#'  ## Plot updated contMaps
#'
#'  updated_contMaps_10My <- extract_trait_data_10My$contMaps
#'  plot_contMap(contMap = updated_contMaps_10My[[1]])
#'  plot_contMap(contMap = updated_contMaps_10My[[10]])
#'  plot_contMap(contMap = updated_contMaps_10My[[100]])
#'  }
#' }
#'
#'
#' # ----- Example 2: Categorical trait ----- #
#'
#' ## Load categorical trait data mapped on a phylogeny
#' data(eel_cat_3lvl_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_cat_3lvl_data, 1)
#' eel_cat_3lvl_data$densityMaps # Three density maps: one per state
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_cat_3lvl_data_10My <- extract_all_states_from_densityMaps_for_focal_time(
#'    densityMaps = eel_cat_3lvl_data$densityMaps,
#'    nb_simulations = 100,
#'    focal_time = focal_time,
#'    update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_3lvl_data_10My, 1)
#' eel_cat_3lvl_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Dummy stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_cat_3lvl_data_10My$trait_data))
#' trait_data_df[50:60, ]
#'
#' # Distributions of states across dummy stochastic maps reflect posterior probabilities
#' # recorded in densityMaps, but they are not true evolutionary histories.
#' # If you want to relate trait states with a simulated evolutionary history,
#' # you need to provide simmaps as input.
#' }
#'
#'
#'# ----- Example 3: Biogeographic ranges ----- #
#'
#' ## Load biogeographic range data mapped on a phylogeny
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_biogeo_data, 1)
#' length(eel_biogeo_data$simmaps) # 100 stochastic maps: one per simulated biogeographic history
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update simmaps for the given focal_time
#'
#' # Extract from the simmaps
#' eel_biogeo_data_10My <- extract_all_ranges_from_simmaps_for_focal_time(
#'    simmaps = eel_biogeo_data$simmaps,
#'    focal_time = focal_time,
#'    update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_biogeo_data_10My$trait_data))
#' trait_data_df[1:10, ]
#'
#' # Distributions of ranges are recorded across true stochastic maps
#' # as you provided simmaps as input.
#'
#' ## Plot updated stochastic maps
#'
#' # Plot initial stochastic map n°1
#' plot(eel_biogeo_data$simmaps[[1]], fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(eel_biogeo_data$simmaps[[1]])[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated stochastic map n°1, cut at 10 Mya
#' plot(eel_biogeo_data_10My$simmaps[[1]], fsize = 0.7)
#' }
#'


extract_all_trait_values_for_focal_time <- function (contMaps = NULL,
                                                     densityMaps = NULL,
                                                     simmaps = NULL,
                                                     nb_simulations = NULL,
                                                     tip_data = NULL,
                                                     trait_data_type,
                                                     focal_time,
                                                     update_Map = FALSE, # Change for update_Map in the wrapper. Adjust to update_contMap / update_densityMaps / update_simmap in the sub-functions
                                                     keep_tip_labels = TRUE)
{
  ### Check input validity

  ## contMap OR (densityMaps OR simmap)
  if (!is.null(contMaps) & (!is.null(densityMaps) | !is.null(simmaps)))
  {
    stop(paste0("You must provide 'contMaps' (for continuous traits) OR 'densityMaps'/'simmaps' (for categorical and biogeographic traits) according to the type of trait data.\n",
                "'tip_data' can also be provided in complement to use accurate trait values at tips."))
  }

  ## trait_data_type
  # trait_data_type must be "continuous", categorical" or "biogeographic"
  if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
  {
    stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
  }

  ## Check that what is provided contMap OR densityMaps/simmaps match the trait_data_type
  if ((!is.null(contMaps) & trait_data_type != "continuous"))
  {
    stop(paste0("You provided 'contMaps' but selected '",trait_data_type,"' as 'trait_data_type'. contMaps are used to map continuous traits.\n",
                "If you wish to extract trait values for a continuous trait, provide 'trait_data_type = continuous'.\n",
                "If you wish to extract trait states/ranges for ",trait_data_type," data, provide 'densityMaps' or 'simmaps' as input instead of 'contMaps'"))
  }
  if ((is.null(contMaps) & trait_data_type == "continuous"))
  {
    stop(paste0("You selected 'trait_data_type = continuous' but did not provide 'contMaps' as input.\n",
                "contMaps are needed to extract continuous trait data.\n",
                "See ?deepSTRAPP::prepare_trait_data(), ?contsimmap::make.contsimmap(), and ?deepSTRAPP::convert_contsimmap_to_contMaps_list() to learn how to generate those objects."))
  }

  if ((!is.null(densityMaps) & !(trait_data_type %in% c("categorical", "biogeographic"))))
  {
    stop(paste0("You provided 'densityMaps' but selected '",trait_data_type,"' as 'trait_data_type'. densityMaps are used to map categorical or biogeographic data.\n",
                "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                "If you wish to extract trait values for a continuous trait, provide 'contMaps' as input instead of 'densityMaps'.\n"))
  }
  if ((!is.null(simmaps) & !(trait_data_type %in% c("categorical", "biogeographic"))))
  {
    stop(paste0("You provided 'simmaps' but selected '",trait_data_type,"' as 'trait_data_type'. simmaps are used to map categorical or biogeographic data.\n",
                "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                "If you wish to extract trait values for a continuous trait, provide 'contMaps' as input instead of 'simmaps'.\n"))
  }
  if (is.null(densityMaps) & is.null(simmaps) & (trait_data_type %in% c("categorical", "biogeographic")))
  {
    stop(paste0("You selected 'trait_data_type = ",trait_data_type,"' but did not provide 'densityMaps' or 'simmaps' as input.\n",
                "'simmaps' are needed to extract state/range data from stochastic maps.\n",
                "Alternatively, 'densityMaps' and 'nb_simulations' can be provided to reconstruct states/ranges from posterior probabilities across dummy stochastic maps."))
  }

  ## nb_simulations
  if ((trait_data_type %in% c("categorical", "biogeographic")) & is.null(simmaps) & is.null(nb_simulations))
  {
    stop(paste0("'nb_simulations' must be provided alongside 'densityMaps' to reconstruct states/ranges from posterior probabilities.\n",
                "Alternatively, you can provide 'simmaps' to extract states/ranges directly from stochastic maps."))
  }

  ## Compute the appropriate sub-function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous = { # Case for continuous data
           # Input = contMap. Trait values are interpolated along branches.
           trait_data_extract <- extract_all_trait_values_from_contMaps_for_focal_time(
             contMaps = contMaps,
             focal_time = focal_time,
             update_contMaps = update_Map,
             keep_tip_labels = keep_tip_labels
           )
         },
         categorical = { # Case for categorical data

           # Input = simmaps. Based on stochastic simulations of trait states.
           if (!is.null(simmaps))
           {
             trait_data_extract <- extract_all_states_from_simmaps_for_focal_time(
               simmaps = simmaps,
               tip_data = tip_data,
               focal_time = focal_time,
               update_simmaps = update_Map,
               keep_tip_labels = keep_tip_labels)

             # Special case when both densityMaps and simmaps are provided and update_Map = TRUE
             # Need to update densityMaps too and add them to the simmap output
             if (!is.null(densityMaps) & update_Map)
             {
               ## Cut densityMap$tree at focal time and update trait mapping in densityMap$tree$maps and densityMap$tree$mapped.edge for all densityMaps in the list
               updated_densityMaps <- cut_densityMaps_for_focal_time(densityMaps = densityMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
               trait_data_extract$densityMaps <- updated_densityMaps
             }

           } else {
             # Input = densityMaps. Based on posterior distribution of trait states.
             trait_data_extract <- extract_all_states_from_densityMaps_for_focal_time(
               densityMaps = densityMaps,
               nb_simulations = nb_simulations,
               tip_data = tip_data,
               focal_time = focal_time,
               update_densityMaps = update_Map,
               keep_tip_labels = keep_tip_labels)
           }
         },
         biogeographic = { # Case for biogeographic data

           # Input = simmaps. Based on stochastic simulations of ranges.
           if (!is.null(simmaps))
           {
             trait_data_extract <- extract_all_ranges_from_simmaps_for_focal_time(
               simmaps = simmaps,
               tip_data = tip_data,
               focal_time = focal_time,
               update_simmaps = update_Map,
               keep_tip_labels = keep_tip_labels)

             # Special case when both densityMaps and simmaps are provided and update_Map = TRUE
             # Need to update densityMaps to and add them to the simmap output
             if (!is.null(densityMaps) & update_Map)
             {
               ## Cut densityMap$tree at focal time and update trait mapping in density$tree$maps and density$tree$mapped.edge for all densityMaps in the list
               updated_densityMaps <- cut_densityMaps_for_focal_time(densityMaps = densityMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
               trait_data_extract$densityMaps <- updated_densityMaps
             }

            } else {
             # Input = densityMaps. Based on posterior distribution of trait states.
             trait_data_extract <- extract_all_ranges_from_densityMaps_for_focal_time(
               densityMaps = densityMaps,
               nb_simulations = nb_simulations,
               tip_data = tip_data,
               focal_time = focal_time,
               update_densityMaps = update_Map,
               keep_tip_labels = keep_tip_labels)
           }
         }
  )

  ## Export the output
  return(invisible(trait_data_extract))
}


### Sub-function for continuous trait data ####

#' @title Extract trait data mapped across continuous stochastic maps at a given time in the past
#'
#' @description Extracts trait data across multiple continuous stochastic maps
#'   found along branches at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the continuous stochastic maps (`contMaps`) such as
#'   branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the continuous trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param contMaps List of object of class `"contMap"`, typically generated with [deepSTRAPP::prepare_trait_data()]
#'   or [phytools::contMap()], that contains a phylogenetic tree and associated continuous trait mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_contMaps Logical. Specify whether the continuous stochqstic maps (`contMaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated contMaps. Default is `TRUE`.
#'   Used only if `update_contMaps = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.contMap
#' @importFrom ape nodelabels
#'
#' @details The continuous stochastic maps (`contMaps`) are cut at a specific time in the past
#'   (i.e. the `focal_time`) and the recorded trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Trait data are extracted from each `contMap` and stored in the `$trait_data` list.
#'
#'   If providing only `contMaps`, trait values at tips will be extracted from
#'   the mapping of the `contMaps` leading to a slight dependency with the actual tip data.
#'
#'   True tip data will be used if `tip_data` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `contMaps` -----
#'
#'   To obtain updated `contMaps` alongside the trait data, set `update_contMaps = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The continuous trait mapping in `contMaps` (`$tree$maps` and `$tree$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A list of named numerical vector. Each item corresponds to trait data extracted from a stochastic map while found along branches overlapping the `focal_time`.
#'     Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_contMaps = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$contMaps`.
#'   * `$contMaps` A list that contains the updated `contMaps` with  branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$contMap$tree` element.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMaps_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()] [deepSTRAPP::extract_all_states_from_simmaps_for_focal_time()]
#' [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()] [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  ## The R package 'contsimmap' is needed for this example to work.
#'  # Please install it manually from: https://github.com/bstaggmartin/contsimmap.
#'
#'  ## Prepare data
#'
#'  # Load eel data from the R package phytools
#'  # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#'  library(phytools)
#'  data(eel.tree)
#'  data(eel.data)
#'
#'  # Extract body size
#'  eel_data <- setNames(eel.data$Max_TL_cm,
#'                       rownames(eel.data))
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Map continuous trait evolution on the phylogeny
#'  eel_cont_data  <- prepare_trait_data(
#'     tip_data = eel_data,
#'     trait_data_type = "continuous",
#'     phylo = eel.tree,
#'     # evolutionary_models = c("BM", "OU", "lambda", "kappa"),
#'     evolutionary_models = c("BM"),
#'     # Perform stochastic mapping to obtain multiple evolutionary histories
#'     run_stochastic_maps = TRUE,
#'     nb_simulations = 100,
#'     verbose = TRUE)
#'
#'  # Extract continuous stochastic maps (contMaps)
#'  eel_contMaps <- eel_cont_data$contMaps
#'
#'  # Plot the interpolated map of ML estimates
#'  plot_contMap(contMap = eel_cont_data$contMap)
#'
#'  # Plot several continuous stochastic maps (contMaps)
#'  plot_contMap(contMap = eel_contMaps[[1]])
#'  plot_contMap(contMap = eel_contMaps[[10]])
#'  plot_contMap(contMap = eel_contMaps[[100]])
#'
#'  # Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  ## Extract all trait values for focal time = 10 Mya
#'  extract_trait_data_10My <- extract_all_trait_values_from_contMaps_for_focal_time(
#'     contMaps = eel_contMaps, focal_time = 10, update_contMaps = T)
#'
#'  # Convert in data.frame
#'   # Rows = Stochastic maps
#'   # Columns = Cut branches at 10 Mya
#'  trait_data_df <- as.data.frame(do.call(rbind, extract_trait_data_10My$trait_data))
#'  head(trait_data_df)
#'
#'  ## Plot updated contMaps
#'
#'  updated_contMaps_10My <- extract_trait_data_10My$contMaps
#'  plot_contMap(contMap = updated_contMaps_10My[[1]])
#'  plot_contMap(contMap = updated_contMaps_10My[[10]])
#'  plot_contMap(contMap = updated_contMaps_10My[[100]])
#'  }
#' }
#'

extract_all_trait_values_from_contMaps_for_focal_time <- function (
    contMaps,
    focal_time,
    update_contMaps = FALSE,
    keep_tip_labels = TRUE)
{

  ### Check input validity
  {
    ## contMaps
    # Must provide contMaps for continuous traits
    if (is.null(contMaps))
    {
      stop(paste0("You must provide 'contMaps' for continuous traits.\n",
                  "See ?deepSTRAPP::prepare_trait_data(), ?contsimmap::make.contsimmap(), and ?deepSTRAPP::convert_contsimmap_to_contMaps_list() to learn how to generate those objects."))
    }
    # contMaps must be a list of "contMap" class object
    if (!all(unlist(lapply(X = contMaps, FUN = inherits, what = 'contMap'))))
    {
      stop("'contMaps' must be a list of objects with the 'contMap' class. See ?phytools::contMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    # contMap[[i]]$tree must have a $maps element
    if (!all(unlist(lapply(X = contMaps, FUN = function (x) { !is.null(x$tree$maps) } ))))
    {
      stop(paste0("Each 'contMap' must have a $tree$maps element that provides the mapping of the evolution of the continuous trait on the phylogeny.\n",
                  "See ?phytools::contMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # All contMap objects must be built on the same phylogeny
    all_trees <- lapply(X = contMaps, FUN = function (x) { x$tree[c("edge", "edge.length", "Nnode", "tip.label")] })
    if (!all(unlist(lapply(X = all_trees, FUN = identical, all_trees[[1]]))))
    {
      stop(paste0("All 'contMaps' must be built on the same phylogeny.\n",
                  "Check '?deepSTRAPP::prepare_trait_data()' to learn how to produce this kind of object."))
    }
    # All contMap objects must be built with the same time points
    all_maps <- lapply(X = contMaps, FUN = function (x) { unname(unlist(x$tree$maps)) })
    if (!all(unlist(lapply(X = all_maps, FUN = identical, all_maps[[1]]))))
    {
      stop(paste0("All 'contMaps' must be built with the same time points.\n",
                  "Check '?deepSTRAPP::prepare_trait_data()' to learn how to produce this kind of object."))
    }

    ## focal_time

    # Extract root age
    root_age <- max(phytools::nodeHeights(contMaps[[1]]$tree)[,2])

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
  }

  ## Identify edges present at focal time

  # Use contMaps[[1]] as template to identify edges
  contMap <- contMaps[[1]]

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = contMap$tree, focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- contMap$tree$edge[which.min(contMap$tree$edge[, 1] %in% contMap$tree$edge[, 2]), 1]

  # If no edge present, send warning
  if (sum(all_edges_df$edge_present) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_contMaps)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "continuous"))
    } else {
      # Return a NULL object for contMaps
      updated_contMaps <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "continuous", contMaps = updated_contMaps))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$edge_present, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

    # Compute node distances to focal time
    present_edges_df$rootward_node_dist <- abs(present_edges_df$rootward_node_age - focal_time)
    present_edges_df$tipward_node_dist <- abs(present_edges_df$tipward_node_age - focal_time)

    # Initiate fields for (scaled) ACE values at nodes
    present_edges_df$rootward_node_scaled_ACE <- NA
    present_edges_df$tipward_node_scaled_ACE <- NA
    present_edges_df$rootward_node_ACE <- NA
    present_edges_df$tipward_node_ACE <- NA

    # Initiate list of trait_data per contMaps
    trait_data_list <- list()

    ## Loop per contMaps
    for (i in seq_along(contMaps))
    {
      # i <- 1

      # Extract the contMap for simulation n°i
      contMap_i <- contMaps[[i]]

      # Initiate present_edges_df for simulation n°i
      present_edges_df_i <- present_edges_df

      # Loop per edge
      for (k in 1:nrow(present_edges_df_i))
      {
        # k <- 1

        ## Extract scaled ACE values at nodes from the contMap

        # Extract edge ID
        edge_ID_k <- as.numeric(present_edges_df$edge_ID[k])

        # Extract associated edge mapping
        edge_map_k <- contMap_i$tree$maps[[edge_ID_k]]

        # Extract rootward node scaled ACE values as the first mapped values on the edge
        # With contMap, discrepancy with actual rootward node ACE values as this is the expected value for the mean age of the first segment of the edge...
        # But no discrepancy when using contsimmap which compute values at edge ends
        present_edges_df_i$rootward_node_scaled_ACE[k] <- as.numeric(names(edge_map_k)[1])

        # Extract tipward node scaled ACE values as the first mapped values on the edge
        # With contMap, discrepancy with actual tipward node ACE values as this is the expected value for the mean age of the last segment of the edge...
        # But no discrepancy when using contsimmap which compute values at edge ends
        present_edges_df_i$tipward_node_scaled_ACE[k] <- as.numeric(names(edge_map_k)[length(edge_map_k)])

        ## Convert the scaled ACE values at nodes into initial values

        present_edges_df_i$rootward_node_ACE[k] <- unscale_0_1000(x_scaled = present_edges_df_i$rootward_node_scaled_ACE[k], max_val = contMap_i$lims[2], min_val = contMap_i$lims[1])
        present_edges_df_i$tipward_node_ACE[k] <- unscale_0_1000(x_scaled = present_edges_df_i$tipward_node_scaled_ACE[k], max_val = contMap_i$lims[2], min_val = contMap_i$lims[1])

        ## Interpolate trait value at focal time

        # Based on equations from Felsenstein, 1985
        # Estimate ACE along an edge at a specific time-step as a weighted mean of node values with weights being the inverse distance to the nodes
        # ACE = (Xr/Dr + Xt/Dt) / (1/Dr + 1/Dt)
        # Xr = Trait value at rootward node
        # Dr = Distance from focal time to rootward node
        # Xt = Trait value at tipward node
        # Dt = Distance from focal time to tipward node

        # Case when focal time is different from rootward/tipward time
        if (all(c(present_edges_df_i$rootward_node_dist[k], present_edges_df_i$tipward_node_dist[k]) != 0))
        {
          present_edges_df_i$ACE_at_focal_time[k] <- ((present_edges_df_i$rootward_node_ACE[k]/present_edges_df_i$rootward_node_dist[k]) + (present_edges_df_i$tipward_node_ACE[k]/present_edges_df_i$tipward_node_dist[k])) / ((1/present_edges_df_i$rootward_node_dist[k]) + (1/present_edges_df_i$tipward_node_dist[k]))

        }

        # Case when focal time is rootward
        if (present_edges_df_i$rootward_node_dist[k] == 0)
        {
          present_edges_df_i$ACE_at_focal_time[k] <- present_edges_df_i$rootward_node_ACE[k]

        }

        # Case when focal time is tipward
        if (present_edges_df_i$tipward_node_dist[k] == 0)
        {
          present_edges_df_i$ACE_at_focal_time[k] <- present_edges_df_i$tipward_node_ACE[k]

        }
      }

      ## Format "trait_data" output = named vector of most likely values at focal time
      trait_data_i <- present_edges_df_i$ACE_at_focal_time
      # names(trait_data) <- present_edges_df$edge_ID
      if (keep_tip_labels) # Names = tip.labels of tipward nodes
      {
        names(trait_data_i) <- present_edges_df_i$tip.label
      } else { # Names = tipward nodes ID
        names(trait_data_i) <- present_edges_df_i$tipward_node_ID
      }

      ## Store output in trait_data_list
      trait_data_list <- append(x = trait_data_list, values = list(trait_data_i))
      names(trait_data_list) <- paste0("Map_", 1:i)
    }

    ## Update contMap if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_contMaps)
    {
      ## Cut contMap$tree at focal time and update trait mapping in contMap$tree$maps and contMap$tree$mapped.edge
      updated_contMaps <- cut_contMaps_for_focal_time(contMaps = contMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_contMaps)
    {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "continuous"))

    } else {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "continuous", contMaps = updated_contMaps))
    }
  }
}


### Sub-functions for categorical trait data ####

#' @title Extract categorical trait data mapped across densityMaps at a given time in the past
#'
#' @description Extracts all states found along branches across densityMaps
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, `densityMaps` summarize the posterior probabilities of each state, as observed across multiple stochastic maps.
#'   Optionally, the function can update the `densityMaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the posterior state density mapped on the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements in each `densityMap`.
#' @param densityMaps List of objects of class `"densityMap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param nb_simulations Integer. The number of stochastic maps used to simulate trait evolution.
#'   This is needed to reconstruct the distribution of states across dummy stochastic maps based on
#'   the posterior state density summarized in `densityMaps`.
#' @param tip_data (Optional) Named character string vector of tip states.
#'   Names are nodes_ID of the internal nodes. Needed to provide accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_densityMaps Logical. Specify whether the mapped phylogeny (`densityMaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated densityMaps. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.densityMap
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The mapped phylogeny (`densityMaps`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the recorded states of the overlapping edges/branches are extracted.
#'
#'   Since `densityMaps` record only the posterior probabilities of states observed along branches.
#'   the function reconstructs the distribution of states across dummy stochastic maps based on
#'   the posterior state densities and the `nb_simulations` provided.
#'   If you wish to assign real stochastic maps ID to the state values extracted,
#'   you need to provide the full simulated evolutionary histories as `simmaps`
#'   and use [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()] instead.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Posterior probabilities states are extracted from the `densityMaps`.
#'   A distribution of states across dummy stochastic maps is generated based on
#'   the posterior state densities and the `nb_simulations`.
#'   For each dummy stochastic map, states are assigned to each tip and cut branches at `focal_time`.
#'
#'   True tip states will be used if `tip_data` are provided as optional inputs.
#'   Otherwise, state probabilities as recorded in the `densityMaps` will be used.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `densityMaps` -----
#'
#'   To obtain updated `densityMaps` alongside the trait data, set `update_densityMaps = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The posterior state densities of each state recorded in `densityMaps` (`$tree$maps` and `$tree$mapped.edge`)
#'   are updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A list of named character string vector. Each item corresponds to the states extracted from a dummy stochastic map as found along branches overlapping the `focal_time`.
#'     Names are the tip.label/tipward node ID. Names of each item in the list are dummy stochastic maps (`Dummy_map_X`) and do not represent true simulated evolutionary histories.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "categorical". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_densityMaps = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$densityMaps`.
#'
#'   * `$densityMaps` A list of objects of class `"densityMap"` that contains the updated `densityMap` of each state,
#'      with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$densityMaps$tree` elements.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_densityMaps_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' Sub-function to extract states directly from simmaps: [deepSTRAPP::extract_all_states_from_simmaps_for_focal_time()]
#'
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_all_trait_values_from_contMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()] [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Only extent taxa (Ultrametric tree) ----- #
#'
#' ## Load categorical trait data mapped on a phylogeny
#' data(eel_cat_3lvl_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_cat_3lvl_data, 1)
#' eel_cat_3lvl_data$densityMaps # Three density maps: one per state
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_cat_3lvl_data_10My <- extract_all_states_from_densityMaps_for_focal_time(
#'    densityMaps = eel_cat_3lvl_data$densityMaps,
#'    nb_simulations = 100,
#'    focal_time = focal_time,
#'    update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_3lvl_data_10My, 1)
#' eel_cat_3lvl_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Dummy stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_cat_3lvl_data_10My$trait_data))
#' trait_data_df[50:60, ]
#'
#' # Distributions of states across dummy stochastic maps reflect posterior probabilities
#' # recorded in densityMaps, but they are not true evolutionary histories.
#' # If you want to relate trait states with a simulated evolutionary history,
#' # you need to provide simmaps as input.
#' }
#'
#'
#' # ----- Example 2: Include fossils (Non-ultrametric tree) ----- #
#'
#' ## Prepare data
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
#' \donttest{ # (May take several minutes to run)
#' ## Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100,
#'                                        plot_map = FALSE)
#'
#' # Set focal time
#' focal_time <- 80
#'
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' mammals_cat_data_80My <- extract_all_states_from_densityMaps_for_focal_time(
#'     densityMaps = mammals_cat_data$densityMaps,
#'     nb_simulations = 100,
#'     focal_time = focal_time,
#'     update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(mammals_cat_data_80My, 1)
#' mammals_cat_data_80My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Dummy stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, mammals_cat_data_80My$trait_data))
#' trait_data_df[50:60, ]
#'
#' # Distributions of states across dummy stochastic maps reflect posterior probabilities
#' # recorded in densityMaps, but they are not true evolutionary histories.
#' # If you want to relate trait states with a simulated evolutionary history,
#' # you need to provide simmaps as input.
#' }
#'

extract_all_states_from_densityMaps_for_focal_time <- function (
    densityMaps,
    nb_simulations,
    tip_data = NULL,
    focal_time,
    update_densityMaps = FALSE,
    keep_tip_labels = TRUE)
{
  ### Check input validity
  {
    ## densityMaps
    # Must provide densityMaps for categorical traits
    if (is.null(densityMaps))
    {
      stop(paste0("You must provide 'densityMaps' for categorical traits).\n",
                  "See ?deepSTRAPP::prepare_trait_data(), ?phytools::make.simmap(), and ?phytools::densityMap() to learn how to generate those objects."))
    }
    # densityMaps must be a list of "densityMap" class objects
    if (!is.list(densityMaps))
    {
      stop("'densityMaps' must be a list that contains only objects of the 'densityMap' class. See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    all_classes <- unlist(lapply(X = densityMaps, FUN = class))
    if (!all("densityMap" == all_classes))
    {
      stop("'densityMaps' must be a list that contains only objects of the 'densityMap' class. See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    # densityMaps[[i]]$tree must have a $maps element
    maps_check <- unlist(lapply(X = densityMaps, FUN = function (x) { is.null(x$tree$maps) }))
    if (any(maps_check))
    {
      stop(paste0("'densityMaps' objects must have a $tree$maps element that provides the mapping of the evolution of the categorical trait on the phylogeny
                  as posterior probabilty for each edge to harbour a given state.\n",
                  "See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # names(densityMap) should be the states
    if (is.null(names(densityMaps)))
    {
      stop(paste0("'densityMaps' objects must be named after the associated states in this format: 'Density_map_X' where X is the state name.\n",
                  "See ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }

    # Extract state list
    states_list <- names(densityMaps)
    states_list <- str_remove(string = states_list, pattern = "Density_map_")

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named character string vector
      if (!is.character(tip_data))
      {
        if (is.factor(tip_data))
        {
          cat("WARNING: 'tip_data' was provided as factors. It is converted to a vector of character strings.\n")

          tip_data_names <- names(tip_data)
          tip_data <- as.character(tip_data)
          names(tip_data) <- tip_data_names

        } else {
          stop(paste0("For categorical traits, 'tip_data' must be a character string vector that provides states for tips.\n",
                      "The object you provided is not a character string vector."))
        }
      }
      # tip_data should have many states as there are tips in the densityMaps[[i]]$tree
      if (length(tip_data) != length(densityMaps[[1]]$tree$tip.label))
      {
        stop(paste0("'tip_data' should have as many states as there are tips in the densityMaps[[i]]$tree.\n",
                    "Number of states in 'tip_data' = ",length(tip_data),"; number of tips in the densityMaps[[i]]$tree = ",length(densityMaps[[1]]$tree$tip.label),"."))
      }
      # names(tip_data) = densityMaps[[i]]$tree$tip.label
      if (!all(names(tip_data) %in% densityMaps[[1]]$tree$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the densityMaps[[i]]$tree$tip.label."))
      }
      if (!all(names(tip_data) == densityMaps[[1]]$tree$tip.label))
      {
        warning(paste0("States in 'tip_data' are not ordered as tip labels in the densityMaps[[i]]$tree.\n",
                       "They were reordered to follow tip labels."))
      }
    }

    ## focal_time

    # Extract root age
    root_age <- max(phytools::nodeHeights(densityMaps[[1]]$tree)[,2])

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
  }

  ## Warn against not providing tip_data
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using states extracted from the densityMaps instead.\n"))
  }

  ## Extract tip states if provided in tip_data
  if (!is.null(tip_data))
  {
    # Reorder states in tip_data to match tip.label
    tip_data <- tip_data[densityMaps[[1]]$tree$tip.label]

    # Use them only for focal_time = 0
    tip_data_is_provided <- T
  } else {
    tip_data_is_provided <- F
  }

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = densityMaps[[1]]$tree, focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- densityMaps[[1]]$tree$edge[which.min(densityMaps[[1]]$tree$edge[, 1] %in% densityMaps[[1]]$tree$edge[, 2]), 1]

  # If no edge present, send warning
  if (sum(all_edges_df$edge_present) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_densityMaps)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "categorical"))
    } else {
      # Return a NULL object for densityMaps
      updated_densityMaps <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "categorical", densityMaps = updated_densityMaps))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$edge_present, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

    # Compute node distances to focal time
    present_edges_df$rootward_node_dist <- abs(present_edges_df$rootward_node_age - focal_time)
    present_edges_df$tipward_node_dist <- abs(present_edges_df$tipward_node_age - focal_time)

    ## Initiate list of trait_data per dummy stochastic maps

    # Format "trait_data" output = named vector of states recorded at focal time
    dummy_trait_data <- rep(x = NA, times = nrow(present_edges_df))
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(dummy_trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(dummy_trait_data) <- present_edges_df$tipward_node_ID
    }

    # Repeat across dummy stochastic maps
    trait_data_list <- rep(x = list(dummy_trait_data), times = nb_simulations)
    names(trait_data_list) <- paste0("Dummy_map_", 1:nb_simulations)

    ## Loop per edge
    for (k in 1:nrow(present_edges_df))
    {
      # k <- 4

      ## Extract posterior probabilities at focal time from densityMaps

      # Extract edge ID
      edge_ID_k <- as.numeric(present_edges_df$edge_ID[k])

      # Extract associated edge mappings across states
      edge_maps_k <- lapply(X = densityMaps, FUN = function (x) { x$tree$maps[[edge_ID_k]] } )

      # Compute rootward ages of segments
      segment_rootward_ages_k <- rev(cumsum(rev(edge_maps_k[[1]])) + present_edges_df$tipward_node_age[k])
      # Identify segment matching the given focal time
      if (all(!(segment_rootward_ages_k < focal_time)))
      {
        # Case where all rootward ages are lower than focal_time, then focal_segment is the last one
        focal_segment_ID <- length(segment_rootward_ages_k)
      } else {
        # Otherwise focal_segment is the last to have a rootward age > to focal_time
        # focal_segment_ID <- which.max(segment_rootward_ages_k < focal_time) - 1
        focal_segment_ID <- which.min(segment_rootward_ages_k >= focal_time) - 1
      }

      # Extract posterior probability for focal segments
      edge_PP_k <- as.numeric(unlist(lapply(X = edge_maps_k, FUN = function (x) { names(x)[focal_segment_ID] } )))
      edge_PP_k <- edge_PP_k / 1000 # Rescale to true proportion ranging from 0 to 1
      names(edge_PP_k) <- lapply(X = densityMaps, function (x) { x$states[2] } )

      # Convert into frequencies of state observations across dummy stochastic maps
      edge_frequencies_k <- round(edge_PP_k * nb_simulations, 0)

      # Generate state distribution across dummy stochastic maps for edge n°k
      dummy_states_edge_k <- c()
      for (i in seq_along(edge_frequencies_k))
      {
        # i <- 1

        # Extract state value
        state_i <- names(edge_frequencies_k)[i]

        # Repeat according to recorded frequencies
        dummy_state_i_edge_k <- rep(state_i, times = edge_frequencies_k[i])
        dummy_states_edge_k <- c(dummy_states_edge_k, dummy_state_i_edge_k)
      }

      # Record states distributed across dummy stochastic maps for edge n°k
      for (j in 1:nb_simulations)
      {
        trait_data_list[[j]][k] <- dummy_states_edge_k[j]
      }

      # Check if properly recorded
      # unlist(lapply(X = trait_data_list, FUN = function(x) { x[k] } ))
    }

    ## Match states from tip_data if needed to correct for possible discrepancy from the densityMaps

    # Build df for tips to adjust
    if (tip_data_is_provided)
    {
      tip_data_df <- as.data.frame(tip_data)
      tip_data_df$node_label <- row.names(tip_data_df)
      names(tip_data_df) <- c("state", "node_label")
      accurate_states_df <- tip_data_df[, c("node_label", "state")]
      row.names(accurate_states_df) <- NULL

      # Retrieve node ages
      accurate_states_df <- dplyr::left_join(x = accurate_states_df,
                                             y = all_edges_df[, c("tip.label", "tipward_node_age")],
                                             by = dplyr::join_by("node_label" == "tip.label"))
      # Remove root to avoid issue with NA
      accurate_states_df <- accurate_states_df[!is.na(accurate_states_df$tipward_node_age), ]

      # Detect matches based on focal time (apply a 10^-5 tolerance)
      if (any(abs(accurate_states_df$tipward_node_age - focal_time) < 1e-05))
      {
        # Extract only matched tips
        accurate_states_df_to_patch <- accurate_states_df[(abs(accurate_states_df$tipward_node_age - focal_time) < 1e-05), ]

        # Replace recorded states distributed across dummy stochastic maps with provided tip states
        for (j in 1:nb_simulations)
        {
          trait_data_list[[j]][match(x = accurate_states_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_states_df_to_patch$state
        }
      }
    }

    ## Update densityMaps if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_densityMaps)
    {
      ## Cut densityMap$tree at focal time and update trait mapping in density$tree$maps and density$tree$mapped.edge for all densityMaps in the list
      updated_densityMaps <- cut_densityMaps_for_focal_time(densityMaps = densityMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_densityMaps)
    {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "categorical"))

    } else {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "categorical", densityMaps = updated_densityMaps))
    }
  }
}

#' @title Extract categorical trait data mapped across stochastic maps at a given time in the past
#'
#' @description Extracts all states found along branches across multiple stochastic maps
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, `simmaps` represent independent evolutionary histories simulated conditioned to observed data and model fit (i.e, stochastic maps)
#'   Optionally, the function can update the `simmaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the states mapped on the cut off branches are removed
#'   by updating the `$maps` and `$mapped.edge` in each `simmap`.
#' @param simmaps List of objects of classes `"phylo"` and `"simmap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains multiple evolutionary histories mapping state evolution along branches.
#'   Each object (i.e., `simmap`) corresponds to an independent evolutionary histories simulated conditioned to observed data and model fit.
#' @param tip_data (Optional) Named character string vector of tip states.
#'   Names are nodes_ID of the internal nodes. Needed to provide accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_simmaps Logical. Specify whether the stochastic maps (`simmaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated simmaps. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.densityMap
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The stochastic maps (`simmaps`) are cut at a specific time in the past
#'   (i.e. the `focal_time`) and the recorded states of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Ancestral states are extracted from the `simmaps`.
#'   For each stochastic map, states are assigned to each tip and cut branches at `focal_time`.
#'
#'   True tip states will be used if `tip_data` are provided as optional inputs.
#'   Otherwise, states as recorded in the `simmaps` will be used.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `simmaps` -----
#'
#'   To obtain updated `simmaps` alongside the trait data, set `update_simmaps = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The states recorded in `simmaps` (`$maps` and `$mapped.edge`)
#'   are updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A list of named character string vector. Each item corresponds to the states extracted from a stochastic map as found along branches overlapping the `focal_time`.
#'     Names are the tip.label/tipward node ID. Names of each item in the list are the stochastic maps ID (`Map_X`)
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "categorical". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_simmaps = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$simmaps`.
#'
#'   * `$simmaps` A list of objects with the class `"simmap"` that contains the updated stochastic map from each trait evolution simulation,
#'      with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements each `simmap`.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_simmap_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' Sub-function to extract states from densityMaps: [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()]
#'
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_all_trait_values_from_contMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()] [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Only extent taxa (Ultrametric tree) ----- #
#'
#' ## Load categorical trait data mapped on a phylogeny
#' data(eel_cat_3lvl_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_cat_3lvl_data, 1)
#' length(eel_cat_3lvl_data$simmaps) # 100 stochastic maps: one per simulated evolutionary history
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the simmaps
#' eel_cat_3lvl_data_10My <- extract_all_states_from_simmaps_for_focal_time(
#'    simmaps = eel_cat_3lvl_data$simmaps,
#'    focal_time = focal_time,
#'    update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_3lvl_data_10My, 1)
#' eel_cat_3lvl_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_cat_3lvl_data_10My$trait_data))
#' trait_data_df[1:10, ]
#'
#' ## Plot updated stochastic maps
#'
#' # Plot initial stochastic map n°1
#' plot(eel_cat_3lvl_data$simmaps[[1]], fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(eel_cat_3lvl_data$simmaps[[1]])[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated stochastic map n°1, cut at 10 Mya
#' plot(eel_cat_3lvl_data_10My$simmaps[[1]], fsize = 0.7)
#' }
#'
#'
#' # ----- Example 2: Include fossils (Non-ultrametric tree) ----- #
#' ## Test with non-ultrametric trees like mammals in motmot
#'
#' ## Prepare data
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
#' \donttest{ # (May take several minutes to run)
#' ## Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100,
#'                                        # Keep the simmaps in the output
#'                                        return_simmaps = TRUE,
#'                                        plot_map = FALSE)
#'
#' # Set focal time
#' focal_time <- 80
#'
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' mammals_cat_data_80My <- extract_all_states_from_simmaps_for_focal_time(
#'     simmaps = mammals_cat_data$simmaps,
#'     tip_data = mammals_data,
#'     focal_time = focal_time,
#'     update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(mammals_cat_data_80My, 1)
#' mammals_cat_data_80My$trait_data[1:2]
#'
#' # Convert in data.frame
#' # Rows = Stochastic maps
#' # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, mammals_cat_data_80My$trait_data))
#' trait_data_df[1:10, ]
#'
#' ## Plot updated stochastic maps
#'
#' # Plot initial stochastic map n°1
#' plot(mammals_cat_data$simmaps[[1]], fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(mammals_cat_data$simmaps[[1]])[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated stochastic map n°1, cut at 80 Mya
#' plot(mammals_cat_data_80My$simmaps[[1]], fsize = 0.7)
#' }
#'

extract_all_states_from_simmaps_for_focal_time <- function (
    simmaps,
    tip_data = NULL,
    focal_time,
    update_simmaps = FALSE,
    keep_tip_labels = TRUE)
{
  ### Check input validity
  {
    ## simmaps
    # Must provide simmaps for categorical traits
    if (is.null(simmaps))
    {
      stop(paste0("You must provide 'simmaps' for categorical traits).\n",
                  "See ?deepSTRAPP::prepare_trait_data() and ?phytools::make.simmap() to learn how to generate those objects."))
    }
    # simmaps must be a list of "simmaps" class objects
    if (!is.list(simmaps))
    {
      stop("'simmaps' must be a list that contains only objects with the 'simmap' class. See ?deepSTRAPP::prepare_trait_data() and ?phytools::make.simmap() to learn how to generate those objects.")
    }
    all_classes_test <- unlist(lapply(X = simmaps, FUN = function (x) { inherits(x = x, what = "simmap") } ))
    if (!all(all_classes_test))
    {
      stop("'simmaps' must be a list that contains only objects with the 'simmap' class. See ?deepSTRAPP::prepare_trait_data() and ?phytools::make.simmap() to learn how to generate those objects.")
    }
    # simmaps[[i]] must have a $maps element
    maps_check <- unlist(lapply(X = simmaps, FUN = function (x) { is.null(x$maps) }))
    if (any(maps_check))
    {
      stop(paste0("'simmaps' objects must have a $maps element that provides the mapping of the evolution of the categorical trait on the phylogeny.\n",
                  "?deepSTRAPP::prepare_trait_data() and ?phytools::make.simmap() to learn how to generate those objects."))
    }

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named character string vector
      if (!is.character(tip_data))
      {
        if (is.factor(tip_data))
        {
          cat("WARNING: 'tip_data' was provided as factors. It is converted to a vector of character strings.\n")

          tip_data_names <- names(tip_data)
          tip_data <- as.character(tip_data)
          names(tip_data) <- tip_data_names

        } else {
          stop(paste0("For categorical traits, 'tip_data' must be a character string vector that provides states for tips.\n",
                      "The object you provided is not a character string vector."))
        }
      }
      # tip_data should have many states as there are tips in the simmaps[[i]]
      if (length(tip_data) != length(simmaps[[1]]$tip.label))
      {
        stop(paste0("'tip_data' should have as many states as there are tips in the simmaps.\n",
                    "Number of states in 'tip_data' = ",length(tip_data),"; number of tips in the simmaps[[1]] = ",length(simmaps[[1]]$tip.label),"."))
      }
      # names(tip_data) = simmaps[[i]]$tip.label
      if (!all(names(tip_data) %in% simmaps[[1]]$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the simmaps."))
      }
      if (!all(names(tip_data) == simmaps[[1]]$tip.label))
      {
        warning(paste0("States in 'tip_data' are not ordered as tip labels in the simmaps[[1]].\n",
                       "They were reordered to follow tip labels."))
      }
    }

    ## focal_time

    # Extract root age
    root_age <- max(phytools::nodeHeights(simmaps[[1]])[,2])

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
  }

  ## Warn against not providing tip_data
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using states extracted from the simmaps instead.\n"))
  }

  ## Extract tip states if provided in tip_data
  if (!is.null(tip_data))
  {
    # Reorder states in tip_data to match tip.label
    tip_data <- tip_data[simmaps[[1]]$tip.label]

    # Use them only for focal_time = 0
    tip_data_is_provided <- T
  } else {
    tip_data_is_provided <- F
  }

  ## Compute nb of simulations
  nb_simulations <- length(simmaps)

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = simmaps[[1]], focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- simmaps[[1]]$edge[which.min(simmaps[[1]]$edge[, 1] %in% simmaps[[1]]$edge[, 2]), 1]

  # If no edge present, send warning
  if (sum(all_edges_df$edge_present) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_simmaps)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "categorical"))
    } else {
      # Return a NULL object for simmaps
      updated_simmaps <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "categorical", simmaps = updated_simmaps))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$edge_present, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

    ## Initiate list of trait_data per stochastic maps

    # Format "trait_data" output = named vector of states recorded at focal time
    dummy_trait_data <- rep(x = NA, times = nrow(present_edges_df))
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(dummy_trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(dummy_trait_data) <- present_edges_df$tipward_node_ID
    }

    # Repeat across stochastic maps
    trait_data_list <- rep(x = list(dummy_trait_data), times = nb_simulations)
    names(trait_data_list) <- paste0("Map_", 1:nb_simulations)
    # trait_data_list[[1]]

    ## Loop per edge
    for (k in 1:nrow(present_edges_df))
    {
      # k <- 4

      ## Extract states of all edge segments from simmaps

      # Extract edge ID
      edge_ID_k <- as.numeric(present_edges_df$edge_ID[k])

      # Extract associated edge mappings across states
      edge_maps_k <- lapply(X = simmaps, FUN = function (x) { x$maps[[edge_ID_k]] } )

      ## Identify the edge segment associated with focal time in each simmap
      states_edge_k <- c()
      for (i in seq_along(edge_maps_k))
      {
        # i <- 1

        # Extract mapping for simulation n°i
        edge_maps_ki <- edge_maps_k[[i]]

        # Compute rootward ages of segments
        segment_rootward_ages_ki <- rev(cumsum(rev(edge_maps_ki)) + present_edges_df$tipward_node_age[k])
        # Identify segment matching the given focal time
        if (all(!(segment_rootward_ages_ki < focal_time)))
        {
          # Case where all rootward ages are lower than focal_time, then focal_segment is the last one
          focal_segment_ID_i <- length(segment_rootward_ages_ki)
        } else {
          # Otherwise focal_segment is the last to have a rootward age > to focal_time
          # focal_segment_ID_i <- which.max(segment_rootward_ages_ki < focal_time) - 1
          focal_segment_ID_i <- which.min(segment_rootward_ages_ki >= focal_time) - 1
        }

        ## Extract states for focal segments only
        states_edge_ki <- names(edge_maps_k[[i]][focal_segment_ID_i])
        states_edge_k <- c(states_edge_k, states_edge_ki)
      }

      # Record states distributed across stochastic maps for edge n°k
      for (j in 1:nb_simulations)
      {
        trait_data_list[[j]][k] <- states_edge_k[j]
      }

      # Check if properly recorded
      # unlist(lapply(X = trait_data_list, FUN = function(x) { x[k] } ))
    }

    ## Match states from tip_data if needed to correct for possible discrepancy from the simmaps

    # Build df for tips to adjust
    if (tip_data_is_provided)
    {
      tip_data_df <- as.data.frame(tip_data)
      tip_data_df$node_label <- row.names(tip_data_df)
      names(tip_data_df) <- c("state", "node_label")
      accurate_states_df <- tip_data_df[, c("node_label", "state")]
      row.names(accurate_states_df) <- NULL

      # Retrieve node ages
      accurate_states_df <- dplyr::left_join(x = accurate_states_df,
                                             y = all_edges_df[, c("tip.label", "tipward_node_age")],
                                             by = dplyr::join_by("node_label" == "tip.label"))
      # Remove root to avoid issue with NA
      accurate_states_df <- accurate_states_df[!is.na(accurate_states_df$tipward_node_age), ]

      # Detect matches based on focal time (apply a 10^-5 tolerance)
      if (any(abs(accurate_states_df$tipward_node_age - focal_time) < 1e-05))
      {
        # Extract only matched tips
        accurate_states_df_to_patch <- accurate_states_df[(abs(accurate_states_df$tipward_node_age - focal_time) < 1e-05), ]

        # Replace recorded states distributed across stochastic maps with provided tip states
        for (j in 1:nb_simulations)
        {
          trait_data_list[[j]][match(x = accurate_states_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_states_df_to_patch$state
        }
      }
    }

    ## Update simmaps if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_simmaps)
    {
      ## Cut simmaps at focal time and update trait mapping in simmaps$maps and simmaps$mapped.edge for all simmaps in the list
      updated_simmaps <- cut_simmaps_for_focal_time(simmaps = simmaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_simmaps)
    {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "categorical"))

    } else {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "categorical", simmaps = updated_simmaps))
    }
  }
}


### Sub-function for biogeographic range data ####

#' @title Extract biogeographic range data mapped across densityMaps at a given time in the past
#'
#' @description Extracts all ranges found along branches across densityMaps
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, `densityMaps` summarize the posterior probabilities of each range, as observed across multiple biogeographic stochastic maps.
#'   Optionally, the function can update the `densityMaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the posterior range density mapped on the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements in each `densityMap`.
#' @param densityMaps List of objects of class `"densityMap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being associated with a given range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a biogeographic range. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param nb_simulations Integer. The number of stochastic maps used to simulate biogeographic range evolution.
#'   This is needed to reconstruct the distribution of ranges across dummy stochastic maps based on
#'   the posterior range density summarized in `densityMaps`.
#' @param tip_data (Optional) Named character string vector of tip ranges.
#'   Names are nodes_ID of the internal nodes. Needed to provide accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_densityMaps Logical. Specify whether the mapped phylogeny (`densityMaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and range mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated densityMaps. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.densityMap
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The mapped phylogeny (`densityMaps`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the recorded biogeographic ranges of the overlapping edges/branches are extracted.
#'
#'   Since `densityMaps` record only the posterior probabilities of ranges observed along branches.
#'   the function reconstructs the distribution of states across dummy stochastic maps based on
#'   the posterior range densities and the `nb_simulations` provided.
#'   If you wish to assign real stochastic maps ID to the range values extracted,
#'   you need to provide the full simulated biogeographic histories as `simmaps`
#'   and use [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()] instead.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Posterior probabilities ranges are extracted from the `densityMaps`.
#'   A distribution of ranges across dummy stochastic maps is generated based on
#'   the posterior range densities and the `nb_simulations`.
#'   For each dummy stochastic map, ranges are assigned to each tip and cut branches at `focal_time`.
#'
#'   True tip ranges will be used if `tip_data` are provided as optional inputs.
#'   Otherwise, range probabilities as recorded in the `densityMaps` will be used.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `densityMaps` -----
#'
#'   To obtain updated `densityMaps` alongside the trait data, set `update_densityMaps = TRUE`.
#'   The update consists in cutting off branches and range mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The posterior densities of each range recorded in `densityMaps` (`$tree$maps` and `$tree$mapped.edge`)
#'   are updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A list of named character string vector. Each item corresponds to the ranges extracted from a dummy stochastic map as found along branches overlapping the `focal_time`.
#'     Names are the tip.label/tipward node ID. Names of each item in the list are dummy stochastic maps (`Dummy_map_X`) and do not represent true simulated biogeogaphic histories.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "biogeographic". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_densityMaps = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$densityMaps`.
#'
#'   * `$densityMaps` A list of objects of class `"densityMap"` that contains the updated `densityMap` of each range,
#'      with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$densityMaps$tree` elements.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_densityMaps_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' Sub-function to extract ranges directly from simmaps: [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_all_trait_values_from_contMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()] [deepSTRAPP::extract_all_states_from_simmaps_for_focal_time()]
#'
#' @examples
#' ## Load biogeographic range data mapped on a phylogeny
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_biogeo_data, 1)
#' eel_biogeo_data$densityMaps # Two density maps: one per unique area: A, B.
#' eel_biogeo_data$densityMaps_all_ranges # Three density maps: one per range: A, B, and AB.
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' # ----- Example 1: Using only unique areas ----- #
#'
#' \donttest{ # (May take several seconds to run)
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_biogeo_data_10My <- extract_all_ranges_from_densityMaps_for_focal_time(
#'    densityMaps = eel_biogeo_data$densityMaps,
#'    nb_simulations = 100,
#'    focal_time = focal_time,
#'    update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Dummy stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_biogeo_data_10My$trait_data))
#' trait_data_df[35:45, ]
#'
#' # Distributions of ranges across dummy stochastic maps reflect posterior probabilities
#' # recorded in densityMaps, but they are not true evolutionary histories.
#' # If you want to relate ranges with a simulated evolutionary history,
#' # you need to provide simmaps as input.
#'
#' # ----- Example 2: Using all ranges ----- #
#'
#' ## Extract trait data and update densityMaps_all_ranges for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_biogeo_data_10My <- extract_all_ranges_from_densityMaps_for_focal_time(
#'   densityMaps = eel_biogeo_data$densityMaps_all_ranges,
#'   nb_simulations = 100,
#'   focal_time = focal_time,
#'   update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'   # Rows = Dummy stochastic maps
#'   # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_biogeo_data_10My$trait_data))
#' trait_data_df[35:45, ]
#'
#' # Includes AB among possible ranges
#'
#' # Distributions of ranges across dummy stochastic maps reflect posterior probabilities
#' # recorded in densityMaps, but they are not true evolutionary histories.
#' # If you want to relate ranges with a simulated evolutionary history,
#' # you need to provide simmaps as input.
#' }
#'

extract_all_ranges_from_densityMaps_for_focal_time <- function (
    densityMaps,
    nb_simulations,
    tip_data = NULL,
    focal_time,
    update_densityMaps = FALSE,
    keep_tip_labels = TRUE)
{
  ### Check input validity
  {
    ## densityMaps
    # Must provide densityMaps for biogeographic traits
    if (is.null(densityMaps))
    {
      stop(paste0("You must provide 'densityMaps' for biogeographic traits).\n",
                  "See ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # densityMaps must be a list of "densityMap" class objects
    if (!is.list(densityMaps))
    {
      stop("'densityMaps' must be a list that contains only objects of the 'densityMap' class. See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    all_classes <- unlist(lapply(X = densityMaps, FUN = class))
    if (!all("densityMap" == all_classes))
    {
      stop("'densityMaps' must be a list that contains only objects of the 'densityMap' class. See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    # densityMaps[[i]]$tree must have a $maps element
    maps_check <- unlist(lapply(X = densityMaps, FUN = function (x) { is.null(x$tree$maps) }))
    if (any(maps_check))
    {
      stop(paste0("'densityMaps' objects must have a $tree$maps element that provides the mapping of the evolution of the biogeographic ranges on the phylogeny
                  as posterior probabilty for each edge to harbour a given range.\n",
                  "See ?phytools::densityMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # names(densityMap) should be the ranges
    if (is.null(names(densityMaps)))
    {
      stop(paste0("'densityMaps' objects must be named after the associated ranges in this format: 'Density_map_X' where X is the range name.\n",
                  "See ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }

    # Extract range list
    ranges_list <- names(densityMaps)
    ranges_list <- str_remove(string = ranges_list, pattern = "Density_map_")

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named character string vector
      if (!is.character(tip_data))
      {
        if (is.factor(tip_data))
        {
          cat("WARNING: 'tip_data' was provided as factors. It is converted to a vector of character strings.\n")

          tip_data_names <- names(tip_data)
          tip_data <- as.character(tip_data)
          names(tip_data) <- tip_data_names

        } else {
          stop(paste0("For biogeographic ranges, 'tip_data' must be a character string vector that provides ranges for tips.\n",
                      "The object you provided is not a character string vector."))
        }
      }
      # tip_data should have as many ranges as there are tips in the densityMaps[[i]]$tree
      if (length(tip_data) != length(densityMaps[[1]]$tree$tip.label))
      {
        stop(paste0("'tip_data' should have as many ranges as there are tips in the densityMaps[[i]]$tree.\n",
                    "Number of ranges in 'tip_data' = ",length(tip_data),"; number of tips in the densityMaps[[i]]$tree = ",length(densityMaps[[1]]$tree$tip.label),"."))
      }
      # names(tip_data) = densityMaps[[i]]$tree$tip.label
      if (!all(names(tip_data) %in% densityMaps[[1]]$tree$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the densityMaps[[i]]$tree$tip.label."))
      }
      if (!all(names(tip_data) == densityMaps[[1]]$tree$tip.label))
      {
        warning(paste0("Ranges in 'tip_data' are not ordered as tip labels in the densityMaps[[i]]$tree.\n",
                       "They were reordered to follow tip labels."))
      }
    }

    ## focal_time

    # Extract root age
    root_age <- max(phytools::nodeHeights(densityMaps[[1]]$tree)[,2])

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
  }

  ## Warn against not providing tip_data
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using ranges extracted from the densityMaps instead.\n"))
  }

  ## Extract tip ranges if provided in tip_data
  if (!is.null(tip_data))
  {
    # Reorder ranges in tip_data to match tip.label
    tip_data <- tip_data[densityMaps[[1]]$tree$tip.label]

    # Use them only for focal_time = 0
    tip_data_is_provided <- T
  } else {
    tip_data_is_provided <- F
  }

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = densityMaps[[1]]$tree, focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- densityMaps[[1]]$tree$edge[which.min(densityMaps[[1]]$tree$edge[, 1] %in% densityMaps[[1]]$tree$edge[, 2]), 1]

  # If no edge present, send warning
  if (sum(all_edges_df$edge_present) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_densityMaps)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "biogeographic"))
    } else {
      # Return a NULL object for densityMaps
      updated_densityMaps <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "biogeographic", densityMaps = updated_densityMaps))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$edge_present, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

    # Compute node distances to focal time
    present_edges_df$rootward_node_dist <- abs(present_edges_df$rootward_node_age - focal_time)
    present_edges_df$tipward_node_dist <- abs(present_edges_df$tipward_node_age - focal_time)

    ## Initiate list of trait_data per dummy stochastic maps

    # Format "trait_data" output = named vector of ranges recorded at focal time
    dummy_trait_data <- rep(x = NA, times = nrow(present_edges_df))
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(dummy_trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(dummy_trait_data) <- present_edges_df$tipward_node_ID
    }

    # Repeat across dummy stochastic maps
    trait_data_list <- rep(x = list(dummy_trait_data), times = nb_simulations)
    names(trait_data_list) <- paste0("Dummy_map_", 1:nb_simulations)

    ## Loop per edge
    for (k in 1:nrow(present_edges_df))
    {
      # k <- 4

      ## Extract posterior probabilities at focal time from densityMaps

      # Extract edge ID
      edge_ID_k <- as.numeric(present_edges_df$edge_ID[k])

      # Extract associated edge mappings across ranges
      edge_maps_k <- lapply(X = densityMaps, FUN = function (x) { x$tree$maps[[edge_ID_k]] } )

      # Compute rootward ages of segments
      segment_rootward_ages_k <- rev(cumsum(rev(edge_maps_k[[1]])) + present_edges_df$tipward_node_age[k])
      # Identify segment matching the given focal time
      if (all(!(segment_rootward_ages_k < focal_time)))
      {
        # Case where all rootward ages are lower than focal_time, then focal_segment is the last one
        focal_segment_ID <- length(segment_rootward_ages_k)
      } else {
        # Otherwise focal_segment is the last to have a rootward age > to focal_time
        # focal_segment_ID <- which.max(segment_rootward_ages_k < focal_time) - 1
        focal_segment_ID <- which.min(segment_rootward_ages_k >= focal_time) - 1
      }

      # Extract posterior probability for focal segments
      edge_PP_k <- as.numeric(unlist(lapply(X = edge_maps_k, FUN = function (x) { names(x)[focal_segment_ID] } )))
      edge_PP_k <- edge_PP_k / 1000 # Rescale to true proportion ranging from 0 to 1
      names(edge_PP_k) <- lapply(X = densityMaps, function (x) { x$states[2] } )

      # Convert into frequencies of range observations across dummy stochastic maps
      edge_frequencies_k <- round(edge_PP_k * nb_simulations, 0)

      # Generate range distribution across dummy stochastic maps for edge n°k
      dummy_ranges_edge_k <- c()
      for (i in seq_along(edge_frequencies_k))
      {
        # i <- 1

        # Extract range value
        range_i <- names(edge_frequencies_k)[i]

        # Repeat according to recorded frequencies
        dummy_range_i_edge_k <- rep(range_i, times = edge_frequencies_k[i])
        dummy_ranges_edge_k <- c(dummy_ranges_edge_k, dummy_range_i_edge_k)
      }

      # Record ranges distributed across dummy stochastic maps for edge n°k
      for (j in 1:nb_simulations)
      {
        trait_data_list[[j]][k] <- dummy_ranges_edge_k[j]
      }

      # Check if properly recorded
      # unlist(lapply(X = trait_data_list, FUN = function(x) { x[k] } ))
    }

    ## Match ranges from tip_data if needed to correct for possible discrepancy from the densityMaps

    # Build df for tips to adjust
    if (tip_data_is_provided)
    {
      tip_data_df <- as.data.frame(tip_data)
      tip_data_df$node_label <- row.names(tip_data_df)
      names(tip_data_df) <- c("range", "node_label")
      accurate_ranges_df <- tip_data_df[, c("node_label", "range")]
      row.names(accurate_ranges_df) <- NULL

      # Retrieve node ages
      accurate_ranges_df <- dplyr::left_join(x = accurate_ranges_df,
                                             y = all_edges_df[, c("tip.label", "tipward_node_age")],
                                             by = dplyr::join_by("node_label" == "tip.label"))
      # Remove root to avoid issue with NA
      accurate_ranges_df <- accurate_ranges_df[!is.na(accurate_ranges_df$tipward_node_age), ]

      # Detect matches based on focal time (apply a 10^-5 tolerance)
      if (any(abs(accurate_ranges_df$tipward_node_age - focal_time) < 1e-05))
      {
        # Extract only matched tips
        accurate_ranges_df_to_patch <- accurate_ranges_df[(abs(accurate_ranges_df$tipward_node_age - focal_time) < 1e-05), ]

        # Replace recorded ranges distributed across dummy stochastic maps with provided tip ranges
        for (j in 1:nb_simulations)
        {
          trait_data_list[[j]][match(x = accurate_ranges_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_ranges_df_to_patch$range
        }
      }
    }

    ## Update densityMaps if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_densityMaps)
    {
      ## Cut densityMap$tree at focal time and update trait mapping in density$tree$maps and density$tree$mapped.edge for all densityMaps in the list
      updated_densityMaps <- cut_densityMaps_for_focal_time(densityMaps = densityMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_densityMaps)
    {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "biogeographic"))

    } else {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "biogeographic", densityMaps = updated_densityMaps))
    }
  }
}


#' @title Extract range data mapped across biogeographic stochastic maps at a given time in the past
#'
#' @description Extracts all ranges found along branches across multiple biogeographic stochastic maps
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, `simmaps` represent independent biogeographic histories simulated conditioned to observed data and model fit (i.e, stochastic maps)
#'   Optionally, the function can update the `simmaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the states mapped on the cut off branches are removed
#'   by updating the `$maps` and `$mapped.edge` in each `simmap`.
#' @param simmaps List of objects of class `"simmap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains multiple biogeographic histories mapping range evolution along branches.
#'   Each object (i.e., `simmap`) corresponds to an independent biogeographic histories simulated conditioned to observed data and model fit.
#' @param tip_data (Optional) Named character string vector of tip states.
#'   Names are nodes_ID of the internal nodes. Needed to provide accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_simmaps Logical. Specify whether the stochastic maps (`simmaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated simmaps. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @export
#' @importFrom phytools nodeHeights plot.densityMap
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The stochastic maps (`simmaps`) are cut at a specific time in the past
#'   (i.e. the `focal_time`) and the recorded states of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Ancestral states are extracted from the `simmaps`.
#'   For each stochastic map, states are assigned to each tip and cut branches at `focal_time`.
#'
#'   True tip states will be used if `tip_data` are provided as optional inputs.
#'   Otherwise, states as recorded in the `simmaps` will be used.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `simmaps` -----
#'
#'   To obtain updated `simmaps` alongside the trait data, set `update_simmaps = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The states recorded in `simmaps` (`$maps` and `$mapped.edge`)
#'   are updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A list of named character string vector. Each item corresponds to the states extracted from a stochastic map as found along branches overlapping the `focal_time`.
#'     Names are the tip.label/tipward node ID. Names of each item in the list are the stochastic maps ID (`Map_X`)
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "categorical". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_simmaps = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$simmaps`.
#'
#'   * `$simmaps` A list of objects with the class `"simmap"` that contains the updated stochastic map from each biogeographic range evolution simulation,
#'      with branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements each `simmap`.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_simmap_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' Sub-function to extract states from densityMaps: [deepSTRAPP::extract_all_states_from_densityMaps_for_focal_time()]
#'
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_all_trait_values_from_contMaps_for_focal_time()]
#' [deepSTRAPP::extract_all_ranges_from_densityMaps_for_focal_time()] [deepSTRAPP::extract_all_ranges_from_simmaps_for_focal_time()]
#'
#' @examples
#' ## Load biogeographic range data mapped on a phylogeny
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_biogeo_data, 1)
#' length(eel_biogeo_data$simmaps) # 100 stochastic maps: one per simulated biogeographic history
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update simmaps for the given focal_time
#'
#' # Extract from the simmaps
#' eel_biogeo_data_10My <- extract_all_ranges_from_simmaps_for_focal_time(
#'    simmaps = eel_biogeo_data$simmaps,
#'    focal_time = focal_time,
#'    update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data[1:2]
#'
#' # Convert in data.frame
#'  # Rows = Stochastic maps
#'  # Columns = Cut branches at 10 Mya
#' trait_data_df <- as.data.frame(do.call(rbind, eel_biogeo_data_10My$trait_data))
#' trait_data_df[1:10, ]
#'
#' # Distributions of ranges are recorded across true stochastic maps
#' # as you provided simmaps as input.
#'
#' ## Plot updated stochastic maps
#'
#' # Plot initial stochastic map n°1
#' plot(eel_biogeo_data$simmaps[[1]], fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(eel_biogeo_data$simmaps[[1]])[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated stochastic map n°1, cut at 10 Mya
#' plot(eel_biogeo_data_10My$simmaps[[1]], fsize = 0.7)
#' }
#'

extract_all_ranges_from_simmaps_for_focal_time <- function (
    simmaps,
    tip_data = NULL,
    focal_time,
    update_simmaps = FALSE,
    keep_tip_labels = TRUE)
{
  ### Check input validity
  {
    ## simmaps
    # Must provide simmaps for biogeographic ranges
    if (is.null(simmaps))
    {
      stop(paste0("You must provide 'simmaps' for biogeographic ranges).\n",
                  "See ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }
    # simmaps must be a list of "simmaps" class objects
    if (!is.list(simmaps))
    {
      stop("'simmaps' must be a list that contains only objects with the 'simmap' class. See ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    all_classes_test <- unlist(lapply(X = simmaps, FUN = function (x) { inherits(x = x, what = "simmap") } ))
    if (!all(all_classes_test))
    {
      stop("'simmaps' must be a list that contains only objects with the 'simmap' class. See ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    # simmaps[[i]] must have a $maps element
    maps_check <- unlist(lapply(X = simmaps, FUN = function (x) { is.null(x$maps) }))
    if (any(maps_check))
    {
      stop(paste0("'simmaps' objects must have a $maps element that provides the mapping of the evolution of biogeographic ranges on the phylogeny.\n",
                  "?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named character string vector
      if (!is.character(tip_data))
      {
        if (is.factor(tip_data))
        {
          cat("WARNING: 'tip_data' was provided as factors. It is converted to a vector of character strings.\n")

          tip_data_names <- names(tip_data)
          tip_data <- as.character(tip_data)
          names(tip_data) <- tip_data_names

        } else {
          stop(paste0("For biogeographic ranges, 'tip_data' must be a character string vector that provides ranges for tips.\n",
                      "The object you provided is not a character string vector."))
        }
      }
      # tip_data should have many ranges as there are tips in the simmaps[[i]]
      if (length(tip_data) != length(simmaps[[1]]$tip.label))
      {
        stop(paste0("'tip_data' should have as many ranges as there are tips in the simmaps.\n",
                    "Number of ranges in 'tip_data' = ",length(tip_data),"; number of tips in the simmaps[[1]] = ",length(simmaps[[1]]$tip.label),"."))
      }
      # names(tip_data) = simmaps[[i]]$tip.label
      if (!all(names(tip_data) %in% simmaps[[1]]$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the simmaps."))
      }
      if (!all(names(tip_data) == simmaps[[1]]$tip.label))
      {
        warning(paste0("Ranges in 'tip_data' are not ordered as tip labels in the simmaps[[1]].\n",
                       "They were reordered to follow tip labels."))
      }
    }

    ## focal_time

    # Extract root age
    root_age <- max(phytools::nodeHeights(simmaps[[1]])[,2])

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
  }

  ## Warn against not providing tip_data
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using ranges extracted from the simmaps instead.\n"))
  }

  ## Extract tip ranges if provided in tip_data
  if (!is.null(tip_data))
  {
    # Reorder ranges in tip_data to match tip.label
    tip_data <- tip_data[simmaps[[1]]$tip.label]

    # Use them only for focal_time = 0
    tip_data_is_provided <- T
  } else {
    tip_data_is_provided <- F
  }

  ## Compute nb of simulations
  nb_simulations <- length(simmaps)

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = simmaps[[1]], focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- simmaps[[1]]$edge[which.min(simmaps[[1]]$edge[, 1] %in% simmaps[[1]]$edge[, 2]), 1]

  # If no edge present, send warning
  if (sum(all_edges_df$edge_present) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_simmaps)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "biogeographic"))
    } else {
      # Return a NULL object for simmaps
      updated_simmaps <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "biogeographic", simmaps = updated_simmaps))
    }

  } else {

    # Extract only edges that are present at the focal time
    present_edges_df <- all_edges_df[all_edges_df$edge_present, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

    ## Initiate list of trait_data per stochastic maps

    # Format "trait_data" output = named vector of ranges recorded at focal time
    dummy_trait_data <- rep(x = NA, times = nrow(present_edges_df))
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(dummy_trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(dummy_trait_data) <- present_edges_df$tipward_node_ID
    }

    # Repeat across stochastic maps
    trait_data_list <- rep(x = list(dummy_trait_data), times = nb_simulations)
    names(trait_data_list) <- paste0("Map_", 1:nb_simulations)
    # trait_data_list[[1]]

    ## Loop per edge
    for (k in 1:nrow(present_edges_df))
    {
      # k <- 4

      ## Extract ranges of all edge segments from simmaps

      # Extract edge ID
      edge_ID_k <- as.numeric(present_edges_df$edge_ID[k])

      # Extract associated edge mappings across ranges
      edge_maps_k <- lapply(X = simmaps, FUN = function (x) { x$maps[[edge_ID_k]] } )

      ## Identify the edge segment associated with focal time in each simmap
      ranges_edge_k <- c()
      for (i in seq_along(edge_maps_k))
      {
        # i <- 1

        # Extract mapping for simulation n°i
        edge_maps_ki <- edge_maps_k[[i]]

        # Compute rootward ages of segments
        segment_rootward_ages_ki <- rev(cumsum(rev(edge_maps_ki)) + present_edges_df$tipward_node_age[k])
        # Identify segment matching the given focal time
        if (all(!(segment_rootward_ages_ki < focal_time)))
        {
          # Case where all rootward ages are lower than focal_time, then focal_segment is the last one
          focal_segment_ID_i <- length(segment_rootward_ages_ki)
        } else {
          # Otherwise focal_segment is the last to have a rootward age > to focal_time
          # focal_segment_ID_i <- which.max(segment_rootward_ages_ki < focal_time) - 1
          focal_segment_ID_i <- which.min(segment_rootward_ages_ki >= focal_time) - 1
        }

        ## Extract ranges for focal segments only
        ranges_edge_ki <- names(edge_maps_k[[i]][focal_segment_ID_i])
        ranges_edge_k <- c(ranges_edge_k, ranges_edge_ki)
      }

      # Record ranges distributed across stochastic maps for edge n°k
      for (j in 1:nb_simulations)
      {
        trait_data_list[[j]][k] <- ranges_edge_k[j]
      }

      # Check if properly recorded
      # unlist(lapply(X = trait_data_list, FUN = function(x) { x[k] } ))
    }

    ## Match ranges from tip_data if needed to correct for possible discrepancy from the simmaps

    # Build df for tips to adjust
    if (tip_data_is_provided)
    {
      tip_data_df <- as.data.frame(tip_data)
      tip_data_df$node_label <- row.names(tip_data_df)
      names(tip_data_df) <- c("range", "node_label")
      accurate_ranges_df <- tip_data_df[, c("node_label", "range")]
      row.names(accurate_ranges_df) <- NULL

      # Retrieve node ages
      accurate_ranges_df <- dplyr::left_join(x = accurate_ranges_df,
                                             y = all_edges_df[, c("tip.label", "tipward_node_age")],
                                             by = dplyr::join_by("node_label" == "tip.label"))
      # Remove root to avoid issue with NA
      accurate_ranges_df <- accurate_ranges_df[!is.na(accurate_ranges_df$tipward_node_age), ]

      # Detect matches based on focal time (apply a 10^-5 tolerance)
      if (any(abs(accurate_ranges_df$tipward_node_age - focal_time) < 1e-05))
      {
        # Extract only matched tips
        accurate_ranges_df_to_patch <- accurate_ranges_df[(abs(accurate_ranges_df$tipward_node_age - focal_time) < 1e-05), ]

        # Replace recorded ranges distributed across stochastic maps with provided tip ranges
        for (j in 1:nb_simulations)
        {
          trait_data_list[[j]][match(x = accurate_ranges_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_ranges_df_to_patch$range
        }
      }
    }

    ## Update simmaps if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_simmaps)
    {
      ## Cut simmaps at focal time and update trait mapping in simmaps$maps and simmaps$mapped.edge for all simmaps in the list
      updated_simmaps <- cut_simmaps_for_focal_time(simmaps = simmaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_simmaps)
    {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "biogeographic"))

    } else {
      return(list(trait_data = trait_data_list, focal_time = focal_time, trait_data_type = "biogeographic", simmaps = updated_simmaps))
    }
  }
}


### Possible update: Make it work with non-dichotomous trees!!!

## Make unit tests for ultrametric (eel.tree / eel_contMap) and non-ultrametric trees (mammals$mammals.phy / mammals_contMap)

## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0


