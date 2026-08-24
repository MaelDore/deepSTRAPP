## Functions to extract the most likely trait data from a mapped phylogeny at a given focal time
# One master function to select the proper pipeline according to data type
# Three internal sub-functions extracting trait evolution according to data type

### Master function to select the proper sub-function according to data type ####

#' @title Extract most likely trait data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely trait values/states/ranges found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogeny (`contMap` or `densityMaps`)
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param contMap For continuous trait data. Object of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()] or [phytools::contMap()],
#'   that contains a phylogenetic tree and associated continuous trait mapping,
#'   representing interpolated ML estimates of ancestral trait values.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param contMaps For continuous trait data. List of objects of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   each containing a phylogenetic tree and associated continuous trait mapping that
#'   represents an independent evolutionary history of ancestral trait evolution,
#'   conditioned to observed trait data and model fit (i.e., stochastic maps).
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
#' @param ace (Optional) Ancestral Character Estimates (ACE) at the internal nodes.
#'   Obtained with [deepSTRAPP::prepare_trait_data()] as output in the `$ace` slot.
#'   * For continuous trait data: Named numerical vector typically generated with [phytools::fastAnc()], [phytools::anc.ML()], or [ape::ace()].
#'     Names are nodes_ID of the internal nodes. Values are ACE of the trait.
#'   * For categorical trait or biogeographic data: Matrix that record the posterior probabilities of ancestral states/ranges.
#'     Rows are internal nodes_ID. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   Needed in all cases to provide accurate estimates of trait values.
#' @param tip_data (Optional) Named vector of tip values of the trait.
#'   * For continuous trait data: Named numerical vector of trait values.
#'   * For categorical trait or biogeographic data: Character string vector of states/ranges
#'   Names are nodes_ID of the internal nodes. Needed to provide accurate tip values.
#' @param trait_data_type Character string. Specify the type of trait data. Must be one of "continuous", "categorical", "biogeographic".
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_Map Logical. Specify whether the mapped phylogeny (`contMap`, `densityMaps`, `simmaps`)
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
#' @details The mapped phylogeny (`contMap` or `densityMaps`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   For continuous trait data:
#'
#'   `contMap` is the default input representing the interpolated ML estimates of trait values along
#'   the phylogeny.
#'
#'   If providing only the `contMap` trait values at tips and internal nodes will be extracted from
#'   the mapping of the `contMap` leading to a slight discrepancy with the actual tip data
#'   and estimated ancestral character values.
#'
#'   True ML trait estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   Alternatively to a unique `contMap`, the user can provide a full set of continuous stochastic maps
#'   as `contMaps`. Each map then represents an independent evolutionary history of ancestral trait evolution,
#'   conditioned to observed trait data and model fit (i.e., stochastic maps). The 'most likely' trait values
#'   will be extracted as the mean of values observed across the stochastic maps.
#'   This quantity closely approximates the ML estimates as typically provided with a `contMap`.
#'
#'   For categorical trait and biogeographic data:
#'
#'   `densityMaps` are the default input. Most likely states/ranges are directly extracted from
#'   the posterior probabilities displayed in the `densityMaps`.
#'   The states/ranges with the highest probability is assigned to each tip and cut branches at `focal_time`.
#'
#'   True ML states/ranges will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   Alternatively to `densityMaps`, the user can provide a full set of stochastic maps as `simmaps`.
#'   Each map then represents an independent evolutionary history of ancestral trait evolution,
#'   conditioned to observed trait data and model fit (i.e., stochastic maps).
#'   The 'most likely' states / ranges will be extracted as the most frequent state/range observed across the stochastic maps.
#'   This quantity is fully equivalent to what is derived from the `densityMaps`.
#'
#'   ----- Update the `contMap(s)`/`densityMaps`/`simmaps` -----
#'
#'   To obtain updated `contMap(s)`/`densityMaps`/`simmaps` alongside the trait data, set `update_Map = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The mapping in `contMap(s)`/`densityMaps` (`$tree$maps` and `$tree$mapped.edge`) or `simmaps` (`$maps` and `$mapped.edge`)
#'   is updated accordingly by removing mapping associated with the cut off branches.
#'
#'   To extract all trait data across multiple stochastic maps, and not just the most likely trait value/state/range,
#'   see this associated function: [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named numerical vector with ML trait values found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous", "categorical", or "biogeographic". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_Map = TRUE`, the output also contains updated mapped phylogenies as: `$contMap` and/or `$contMaps` or `$densityMaps` and/or `$simmaps`.
#'
#'   For continuous trait data:
#'
#'   * `$contMap` An object of class `"contMap"` that contains the updated `contMap` with  branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$contMap$tree` element.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#'   * `$contMaps` List of `"contMap"` objects as described above.
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
#'   * `$simmaps` A list of objects of classes `"phylo"` and `"simmap"`, that contains the updated stochastic maps `simmap`
#'     for each simulation of trait/range evolution.
#'     The same useful sub-elements as describe above are also added to each `simmap` object.
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMap_for_focal_time()] [deepSTRAPP::cut_densityMaps_for_focal_time()]
#'
#' Equivalent function to extract all data for multiple stochastic maps
#' [deepSTRAPP::extract_all_trait_values_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Continuous trait ----- #
#'
#' ## Prepare data
#'
#' # Load eel data from the R package phytools
#' # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract body size
#' eel_data <- setNames(eel.data$Max_TL_cm,
#'                      rownames(eel.data))
#'
#' \donttest{ # (May take several minutes to run)
#' ## Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' eel_ACE <- phytools::fastAnc(tree = eel.tree, x = eel_data)
#'
#' ## Interpolate ML estimates of trait values based on a Brownian Motion model
#' # to interpolate values along branches and obtain a "contMap" object
#' eel_contMap <- phytools::contMap(eel.tree, x = eel_data,
#'                                  res = 100, # Number of time steps
#'                                  plot = FALSE)
#'
#' # Set focal time to 50 Mya
#' focal_time <- 50
#'
#' ## Extract trait data and update contMap for the given focal_time
#'
#' # Extract from the contMap (values are not exact ML estimates)
#' eel_cont_50 <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = eel_contMap,
#'    trait_data_type = "continuous",
#'    focal_time = focal_time,
#'    update_Map = TRUE)
#' # Extract from tip data and ML estimates of ancestral characters (values are true ML estimates)
#' eel_cont_50 <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = eel_contMap,
#'    ace = eel_ACE, tip_data = eel_data,
#'    trait_data_type = "continuous",
#'    focal_time = focal_time,
#'    update_Map = TRUE)
#'
#' ## Visualize outputs
#'
#' # Print trait data
#' eel_cont_50$trait_data
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(eel_contMap, fsize = c(0.5, 1))
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(eel_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' plot(eel_cont_50$contMap)
#' ape::nodelabels(text = eel_cont_50$contMap$tree$initial_nodes_ID) }
#'
#'
#' # ----- Example 2: Categorical trait ----- #
#'
#' \donttest{ # (May take several minutes to run)
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
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_cat_3lvl_data_10My <- extract_most_likely_trait_values_for_focal_time(
#'    densityMaps = eel_cat_3lvl_data$densityMaps,
#'    trait_data_type = "categorical",
#'    focal_time = focal_time,
#'    update_Map = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_3lvl_data_10My, 1)
#' eel_cat_3lvl_data_10My$trait_data
#'
#' ## Plot density maps as overlay of all state posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' plot_densityMaps_overlay(densityMaps = eel_cat_3lvl_data$densityMaps)
#' abline(v = max(phytools::nodeHeights(eel_cat_3lvl_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(eel_cat_3lvl_data_10My$densityMaps) }
#'
#'
#'# ----- Example 3: Biogeographic ranges ----- #
#'
#' \donttest{ # (May take several minutes to run)
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
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_biogeo_data_10My <- extract_most_likely_trait_values_for_focal_time(
#'    densityMaps = eel_biogeo_data$densityMaps,
#'    # ace = eel_biogeo_data$ace,
#'    trait_data_type = "biogeographic",
#'    focal_time = focal_time,
#'    update_Map = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data
#'
#' ## Plot density maps as overlay of all range posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' plot_densityMaps_overlay(densityMaps = eel_biogeo_data$densityMaps)
#' abline(v = max(phytools::nodeHeights(eel_biogeo_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(eel_biogeo_data_10My$densityMaps) }
#'

extract_most_likely_trait_values_for_focal_time <- function (contMap = NULL,
                                                             contMaps = NULL,
                                                             densityMaps = NULL,
                                                             simmaps = NULL,
                                                             ace = NULL,
                                                             tip_data = NULL,
                                                             trait_data_type,
                                                             focal_time,
                                                             update_Map = FALSE, # Change for update_Map in the wrapper. Adjust to update_contMap / update_densityMaps in the sub-functions
                                                             keep_tip_labels = TRUE)
{
  ### Check input validity

  ## trait_data_type
  # trait_data_type must be "continuous", categorical" or "biogeographic"
  if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
  {
    stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
  }

  ## Check that what is provided contMap OR densityMaps match the trait_data_type
  if ( (trait_data_type %in% c("categorical", "biogeographic")) & (!is.null(contMap) | !is.null(contMaps)) )
  {
    stop(paste0("You provided a 'contMap' or 'contMaps' but selected '",trait_data_type,"' as 'trait_data_type'. contMaps are used to map continuous traits.\n",
                "If you wish to extract trait values for a continuous trait, provide 'trait_data_type = continuous'.\n",
                "If you wish to extract trait states/ranges for ",trait_data_type," data, provide 'densityMaps' or 'simmaps' as input instead of 'contMap(s)'"))
  }
  if ( (trait_data_type == "continuous") & (!is.null(densityMaps) | !is.null(simmaps)) )
  {
    stop(paste0("You provided 'densityMaps' or 'simmaps' but selected 'trait_data_type = continuous'. 'densityMaps/simmaps' are used to map categorical or biogeographic data.\n",
                "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                "If you wish to extract trait values for a continuous trait, provide a 'contMap' or 'contMaps' as input instead of 'densityMaps/simmaps'.\n",
                "\tYou can also provide 'ace' and tip_data' in complement to a 'contMap' to use accurate ML estimates of trait values at nodes and tips."))
  }

  # Must provide a contMap or contMaps for continuous traits
  if ( (trait_data_type == "continuous") & (is.null(contMap) & is.null(contMaps)) )
  {
    stop(paste0("You must provide a 'contMap' or 'contMaps' for continuous traits.\n",
                "See ?deepSTRAPP::prepare_trait_data() and ?phytools::contMap() to learn how to generate those objects."))
  }

  # Must provide densityMaps or simmaps for categorical and biogeographic traits
  if ( (trait_data_type %in% c("categorical", "biogeographic")) & (is.null(densityMaps) & is.null(simmaps)) )
  {
    stop(paste0("You must provide 'densityMaps' or 'simmaps' for categorical/biogeographic traits.\n",
                "See ?deepSTRAPP::prepare_trait_data(), ?phytools::densityMap(), and ?phytools::make.simmaps() to learn how to generate those objects."))
  }

  ## Compute the appropriate sub-function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous = { # Case for continuous data

           if (!is.null(contMap))
           {
             # Input = contMap. ML trait values are interpolated along branches.
             trait_data_extract <- extract_most_likely_trait_values_from_contMap_for_focal_time(
               contMap = contMap,
               ace = ace,
               tip_data = tip_data,
               focal_time = focal_time,
               update_contMap = update_Map,
               keep_tip_labels = keep_tip_labels)

             # Special case when both contMap and contMaps are provided and update_Map = TRUE
             # Need to update contMaps too and add them to the contMap output
             if (!is.null(contMaps) & update_Map)
             {
               ## Cut contMap$tree at focal time and update trait mapping in contMap$tree$maps and contMap$tree$mapped.edge for all contMaps in the list
               updated_contMaps <- cut_contMaps_for_focal_time(contMaps = contMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
               trait_data_extract$contMaps <- updated_contMaps
             }

           } else {
             # Input = contMaps. Mean trait values from simulations are aggregated along branches.
             trait_data_extract <- extract_most_likely_trait_values_from_contMaps_for_focal_time(
               contMaps = contMaps,
               tip_data = tip_data,
               focal_time = focal_time,
               update_contMaps = update_Map,
               keep_tip_labels = keep_tip_labels)
           }
         },
         categorical = { # Case for categorical data

           if (!is.null(densityMaps))
           {
             # Input = densityMaps. Based on posterior densities from stochastic simulations of trait states.
             trait_data_extract <- extract_most_likely_states_from_densityMaps_for_focal_time(
               densityMaps = densityMaps,
               ace = ace,
               tip_data = tip_data,
               focal_time = focal_time,
               update_densityMaps = update_Map,
               keep_tip_labels = keep_tip_labels)

             # Special case when both densityMaps and simmaps are provided and update_Map = TRUE
             # Need to update simmaps too and add them to the densityMaps output
             if (!is.null(simmaps) & update_Map)
             {
               ## Cut simmaps at focal time and update trait mapping in simmaps$maps and simmaps$mapped.edge for all simmaps in the list
               updated_simmaps <- cut_simmaps_for_focal_time(simmaps = simmaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
               trait_data_extract$simmaps <- updated_simmaps
             }
           } else {
             # Input = simmaps. Based on stochastic simulations of trait states.
             trait_data_extract <- extract_most_likely_states_from_simmaps_for_focal_time(
               simmaps = simmaps,
               ace = ace,
               tip_data = tip_data,
               focal_time = focal_time,
               update_simmaps = update_Map,
               keep_tip_labels = keep_tip_labels)
           }
         },
         biogeographic = { # Case for biogeographic data

           if (!is.null(densityMaps))
           {
             # Input = densityMaps. Based on posterior densities from stochastic simulations of geographic ranges.
             trait_data_extract <- extract_most_likely_ranges_from_densityMaps_for_focal_time(
               densityMaps = densityMaps,
               ace = ace,
               tip_data = tip_data,
               focal_time = focal_time,
               update_densityMaps = update_Map,
               keep_tip_labels = keep_tip_labels)

             # Special case when both densityMaps and simmaps are provided and update_Map = TRUE
             # Need to update simmaps too and add them to the densityMaps output
             if (!is.null(simmaps) & update_Map)
             {
               ## Cut simmaps at focal time and update trait mapping in simmaps$maps and simmaps$mapped.edge for all simmaps in the list
               updated_simmaps <- cut_simmaps_for_focal_time(simmaps = simmaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
               trait_data_extract$simmaps <- updated_simmaps
             }
           } else {
             # Input = simmaps. Based on stochastic simulations of geographic ranges.
             trait_data_extract <- extract_most_likely_ranges_from_simmaps_for_focal_time(
               simmaps = simmaps,
               ace = ace,
               tip_data = tip_data,
               focal_time = focal_time,
               update_simmaps = update_Map,
               keep_tip_labels = keep_tip_labels)
           }
         }
  )

  ## Export the output
  return(invisible(trait_data_extract))
}


### Sub-function for continuous trait data ####

#' @title Extract continuous trait data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely trait values found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogeny (`contMap`) such as
#'   branches overlapping the `focal_time` are shorten to the `focal_time`, and
#'   the continuous trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param contMap Object of class `"contMap"`, typically generated with [deepSTRAPP::prepare_trait_data()]
#'   or [phytools::contMap()], that contains a phylogenetic tree and associated continuous trait mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param ace Named numeric vector (Optional). Ancestral Character Estimates (ACE) of the internal nodes,
#'   typically generated with [phytools::fastAnc()], [phytools::anc.ML()], or [ape::ace()].
#'   Names are nodes_ID of the internal nodes. Values are ACE of the trait.
#'   Needed to provide accurate estimates of trait values.
#' @param tip_data Named numeric vector (Optional). Tip values of the trait.
#'   Names are nodes_ID of the internal nodes.
#'   Needed to provide accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_contMap Logical. Specify whether the mapped phylogeny (`contMap`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated contMap. Default is `TRUE`.
#'   Used only if `update_contMap = TRUE`.
#'
#' @importFrom phytools nodeHeights plot.contMap
#' @importFrom ape nodelabels
#'
#' @details The mapped phylogeny (`contMap`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   If providing only the `contMap` trait values at tips and internal nodes will be extracted from
#'   the mapping of the `contMap` leading to a slight dependency with the actual tip data
#'   and estimated ancestral character values.
#'
#'   True ML estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `contMap` -----
#'
#'   To obtain an updated `contMap` alongside the trait data, set `update_contMap = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The continuous trait mapping `contMap` (`$tree$maps` and `$tree$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named numerical vector with ML trait values found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_contMap = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$contMap`.
#'   * `$contMap` An object of class that contains the updated `contMap` with  branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$contMap$tree` element.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMap_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Only extent taxa (Ultrametric tree) ----- #
#'
#' ## Prepare data
#'
#' # Load eel data from the R package phytools
#' # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract body size
#' eel_data <- setNames(eel.data$Max_TL_cm,
#'                      rownames(eel.data))
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' eel_ACE <- phytools::fastAnc(tree = eel.tree, x = eel_data)
#'
#' # Interpolate ML values along branches based on a Brownian Motion model
#' # and obtain a "contMap" object
#' eel_contMap <- phytools::contMap(eel.tree, x = eel_data,
#'                                  res = 100, # Number of time steps
#'                                  plot = FALSE)
#'
#' # Set focal time to 50 Mya
#' focal_time <- 50
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update contMap for the given focal_time
#'
#' # Extract from the contMap (values are not exact ML estimates)
#' eel_test <- extract_most_likely_trait_values_from_contMap_for_focal_time(
#'    contMap = eel_contMap,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#' # Extract from tip data and ML estimates of ancestral characters (values are true ML estimates)
#' eel_test <- extract_most_likely_trait_values_from_contMap_for_focal_time(
#'    contMap = eel_contMap,
#'    ace = eel_ACE, tip_data = eel_data,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#'
#' ## Visualize outputs
#'
#' # Print trait data
#' eel_test$trait_data
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(eel_contMap, fsize = c(0.5, 1))
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(eel_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' plot(eel_test$contMap)
#' ape::nodelabels(text = eel_test$contMap$tree$initial_nodes_ID) }
#'
#' # ----- Example 2: Include fossils (Non-ultrametric tree) ----- #
#'
#' ## Test with non-ultrametric trees like mammals in motmot
#'
#' ## Prepare data
#'
#' # Load mammals phylogeny and data from the R package motmot included within deepSTRAPP
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#' data("mammals", package = "deepSTRAPP")
#'
#' mammals_tree <- mammals$mammal.phy
#' mammals_data <- setNames(object = mammals$mammal.mass$mean,
#'                          nm = row.names(mammals$mammal.mass))[mammals_tree$tip.label]
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' mammals_ACE <- phytools::fastAnc(tree = mammals_tree, x = mammals_data)
#'
#' # Run a Stochastic Mapping based on a Brownian Motion model
#' # to interpolate values along branches and obtain a "contMap" object
#' mammals_contMap <- phytools::contMap(mammals_tree, x = mammals_data,
#'                                      res = 100, # Number of time steps
#'                                      plot = FALSE)
#'
#' # Set focal time to 80 Mya
#' focal_time <- 80
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update contMap for the given focal_time
#'
#' # Extract from the contMap (values are not exact ML estimates)
#' mammals_test <- extract_most_likely_trait_values_from_contMap_for_focal_time(
#'    contMap = mammals_contMap,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#' # Extract from tip data and ML estimates of ancestral characters (values are true ML)
#' mammals_test <- extract_most_likely_trait_values_from_contMap_for_focal_time(
#'    contMap = mammals_contMap,
#'    ace = mammals_ACE, tip_data = mammals_data,
#'    focal_time = focal_time,
#'    update_contMap = TRUE)
#'
#' ## Visualize outputs
#'
#' # Print trait data
#' mammals_test$trait_data
#'
#' # Plot node labels on initial stochastic map with cut-off
#' phytools::plot.contMap(mammals_contMap, fsize = c(0.5, 1))
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(mammals_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' phytools::plot.contMap(mammals_test$contMap, fsize = c(0.8, 1))
#' ape::nodelabels(text = mammals_test$contMap$tree$initial_nodes_ID) }
#'


extract_most_likely_trait_values_from_contMap_for_focal_time <- function (
    contMap,
    ace = NULL,
    tip_data = NULL,
    focal_time,
    update_contMap = FALSE,
    keep_tip_labels = TRUE)
{

  ### Check input validity

  {
    ## contMap
    # Must provide a contMap for continuous traits
    if (is.null(contMap))
    {
      stop(paste0("You must provide a 'contMap' for continuous traits.\n",
                  "See ?deepSTRAPP::prepare_trait_data() and ?phytools::contMap() to learn how to generate those objects."))
    }
    # contMap must be a "contMap" class object
    if (!("contMap" %in% class(contMap)))
    {
      stop("'contMap' must have the 'contMap' class. See ?phytools::contMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects.")
    }
    # contMap$tree must have a $maps element
    if (is.null(contMap$tree$maps))
    {
      stop(paste0("'contMap' must have a $tree$maps element that provides the mapping of the evolution of the continuous trait on the phylogeny.\n",
                  "See ?phytools::contMap() and ?deepSTRAPP::prepare_trait_data() to learn how to generate those objects."))
    }

    ## ace
    if (!is.null(ace))
    {
      # ace must be a named numerical vector
      if (!is.numeric(ace))
      {
        stop(paste0("For continuous traits, 'ace' must be a named numerical vector that provides trait values for internal nodes.\n",
                    "The object you provided is not a numerical vector."))
      }
      # ace should have many values as there are internal nodes in the contMap$tree
      if (length(ace) != contMap$tree$Nnode)
      {
        stop(paste0("'ace' should have as many values as there are internal nodes in the contMap$tree.\n",
                    "Number of values in 'ace' = ",length(ace),"; number of internal nodes in the contMap$tree = ",contMap$tree$Nnode,"."))
      }
      # names(ace) = internal node IDs
      if (!all(as.numeric(names(ace)) %in% (length(contMap$tree$tip.label) + 1):(length(contMap$tree$tip.label) + contMap$tree$Nnode)))
      {
        stop(paste0("'names(ace)' should match numerical ID of internal nodes in the contMap$tree."))
      }
      if (!all(as.numeric(names(ace)) == (length(contMap$tree$tip.label) + 1):(length(contMap$tree$tip.label) + contMap$tree$Nnode)))
      {
        warning(paste0("Values in 'ace' are not ordered in increasing numerical ID of internal nodes.\n",
                       "They were reordered to follow the numerical ID of internal nodes."))
      }
    }

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named numerical vector
      if (!is.numeric(tip_data))
      {
        stop(paste0("For continuous traits, 'tip_data' must be a named numerical vector that provides trait values for tips.\n",
                    "The object you provided is not a numerical vector."))
      }
      # tip_data should have many values as there are tips in the contMap$tree
      if (length(tip_data) != length(contMap$tree$tip.label))
      {
        stop(paste0("'tip_data' should have as many values as there are tips in the contMap$tree.\n",
                    "Number of values in 'tip_data' = ",length(tip_data),"; number of tips in the contMap$tree = ",length(contMap$tree$tip.label),"."))
      }
      # names(tip_data) = contMap$tree$tip.label
      if (!all(names(tip_data) %in% contMap$tree$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the contMap$tree$tip.label."))
      }
      if (!all(names(tip_data) == contMap$tree$tip.label))
      {
        warning(paste0("Values in 'tip_data' are not ordered as tip labels in the contMap$tree.\n",
                       "They were reordered to follow tip labels."))
      }
    }

    ## focal_time

    # Extract root age
    root_age <- max(phytools::nodeHeights(contMap$tree)[,2])

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

  ## Warn against not providing ace and tip_data
  if (is.null(ace))
  {
    cat(paste0("WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using values interpolated in the contMap instead.\n"))
  }
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using values interpolated in the contMap instead.\n"))
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
  all_edges_df <- identify_edges_at_focal_time(phylo = contMap$tree, focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- contMap$tree$edge[which.min(contMap$tree$edge[, 1] %in% contMap$tree$edge[, 2]), 1]

  # If no edge present, send warning
  if (sum(all_edges_df$edge_present) == 0)
  {
    warning(paste0("No branch is present at focal time = ", focal_time, ". Return a NULL object.\n"))

    # Return a NULL object for trait_data
    trait_data <- NULL

    if (!update_contMap)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "continuous"))
    } else {
      # Return a NULL object for contMap
      updated_contMap <- NULL
      return(list(trait_data = trait_data, focal_time = focal_time, data_type = "continuous", contMap = updated_contMap))
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

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ACE_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
    }


    ## Update contMap if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_contMap)
    {
      ## Cut contMap$tree at focal time and update trait mapping in contMap$tree$maps and contMap$tree$mapped.edge
      updated_contMap <- cut_contMap_for_focal_time(contMap = contMap, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_contMap)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "continuous"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "continuous", contMap = updated_contMap))
    }
  }
}


#' @title Extract continuous trait data mapped across stochastic maps at a given time in the past
#'
#' @description Extracts the most likely trait values found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogeny (`contMap`) such as
#'   branches overlapping the `focal_time` are shorten to the `focal_time`, and
#'   the continuous trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param contMaps For continuous trait data. List of objects of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   each containing a phylogenetic tree and associated continuous trait mapping that
#'   represents an independent evolutionary history of ancestral trait evolution,
#'   conditioned to observed trait data and model fit (i.e., stochastic maps).
#' @param tip_data Named numeric vector (Optional). Tip values of the trait.
#'   Names are nodes_ID of the internal nodes.
#'   Needed to provide accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_contMaps Logical. Specify whether the continuous stochastic maps (`contMaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated contMap. Default is `TRUE`.
#'   Used only if `update_contMaps = TRUE`.
#'
#' @importFrom phytools nodeHeights plot.contMap
#' @importFrom ape nodelabels
#'
#' @details The mapped phylogeny (`contMap`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   If providing only the `contMap` trait values at tips and internal nodes will be extracted from
#'   the mapping of the `contMap` leading to a slight dependency with the actual tip data
#'   and estimated ancestral character values.
#'
#'   True ML estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `contMap` -----
#'
#'   To obtain an updated `contMap` alongside the trait data, set `update_contMap = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The continuous trait mapping `contMap` (`$tree$maps` and `$tree$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named numerical vector with ML trait values found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_contMap = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$contMap`.
#'   * `$contMap` An object of class that contains the updated `contMap` with  branches and mapping that are younger than the `focal_time` cut off.
#'      The function also adds multiple useful sub-elements to the `$contMap$tree` element.
#'     + `$root_age` Integer. Stores the age of the root of the tree.
#'     + `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` to the `initial_node_ID`. Each row is a node.
#'     + `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'     + `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` to the `initial_edge_ID`. Each row is an edge/branch.
#'     + `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMap_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
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
#'  ## Extract the most likely (mean) values across contMaps for focal time = 10 Mya
#'  extract_contMaps_trait_data_10My <- extract_most_likely_trait_values_from_contMaps_for_focal_time(
#'     contMaps = eel_contMaps, focal_time = 10, update_contMaps = T)
#'
#'  head(extract_contMaps_trait_data_10My$trait_data)
#'
#'  ## Compare to ML estimates from the contMap
#'
#'  # Extract the ML values across contMap for focal time = 10 Mya
#'  extract_contMap_trait_data_10My <- extract_most_likely_trait_values_from_contMap_for_focal_time(
#'     contMap = eel_cont_data$contMap, focal_time = 10, update_contMap = T)
#'
#'  head(extract_contMap_trait_data_10My$trait_data)
#'
#'  ## Plot updated contMaps
#'
#'  updated_contMaps_10My <- extract_contMaps_trait_data_10My$contMaps
#'  plot_contMap(contMap = updated_contMaps_10My[[1]])
#'  plot_contMap(contMap = updated_contMaps_10My[[10]])
#'  plot_contMap(contMap = updated_contMaps_10My[[100]])
#'  }
#' }
#'

extract_most_likely_trait_values_from_contMaps_for_focal_time <- function (
    contMaps,
    tip_data = NULL,
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
                  "See ?deepSTRAPP::prepare_trait_data(), ?contsimmap::make.contsimmap(), and ?deepSTRAPP::convert_contsimmap_to_contMaps() to learn how to generate those objects."))
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

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named numerical vector
      if (!is.numeric(tip_data))
      {
        stop(paste0("For continuous traits, 'tip_data' must be a named numerical vector that provides trait values for tips.\n",
                    "The object you provided is not a numerical vector."))
      }
      # tip_data should have many values as there are tips in the contMaps[[1]]$tree
      if (length(tip_data) != length(contMaps[[1]]$tree$tip.label))
      {
        stop(paste0("'tip_data' should have as many values as there are tips in the contMaps[[1]]$tree.\n",
                    "Number of values in 'tip_data' = ",length(tip_data),"; number of tips in the contMaps[[1]]$tree = ",length(contMaps[[1]]$tree$tip.label),"."))
      }
      # names(tip_data) = contMaps[[1]]$tree$tip.label
      if (!all(names(tip_data) %in% contMaps[[1]]$tree$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the contMaps[[1]]$tree$tip.label."))
      }
      if (!all(names(tip_data) == contMaps[[1]]$tree$tip.label))
      {
        warning(paste0("Values in 'tip_data' are not ordered as tip labels in the contMaps[[1]]$tree.\n",
                       "They were reordered to follow tip labels."))
      }
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

  ## Warn against using contMaps to extract average data
  cat(paste0("WARNING: You are using 'contMaps' to extract the most likely trait values recorded at focal time = ",focal_time,".\n",
             "This will effectively returns the average trait values observed across all continuous stochastic maps.\n",
             "While this quantity approximates the ML estimates, it is recommended to provide directly the 'contMap' as input if you aim to extract ML ancestral trait estimates.\n"))

  ## Warn against not providing ace and tip_data
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using values interpolated in the contMaps instead.\n"))
  }

  ## Aggregate all contMaps into a mean contMap
  mean_contMap <- aggregate_contMaps(contMaps = contMaps,
                                     fun = "mean", # To select the aggregating function. Options are mean and median
                                     display_plot = FALSE, # Whether to plot the resulting aggregated contMap
                                     verbose = FALSE)

  ## Identify edges present at focal time
  all_edges_df <- identify_edges_at_focal_time(phylo = mean_contMap$tree, focal_time = focal_time, tolerance = 10^-5)

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- mean_contMap$tree$edge[which.min(mean_contMap$tree$edge[, 1] %in% mean_contMap$tree$edge[, 2]), 1]

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

    # Create vector to convert scaled values from 0 to 1000 in the mean_contMap to initial values
    trans <- 0:1000/1000 * (mean_contMap$lims[2] - mean_contMap$lims[1]) + mean_contMap$lims[1]
    names(trans) <- 0:1000

    # Loop per edge
    for (i in 1:nrow(present_edges_df))
    {
      # i <- 1

      ## Extract scaled ACE values at nodes from the mean_contMap

      # Extract edge ID
      edge_ID_i <- as.numeric(present_edges_df$edge_ID[i])

      # Extract associated edge mapping
      edge_map_i <- mean_contMap$tree$maps[[edge_ID_i]]

      # Extract rootward node scaled ACE values as the first mapped values on the edge
      # Discrepancy with actual rootward node ACE values as this is the expected value for the mean age of the first segment of the edge...
      present_edges_df$rootward_node_scaled_ACE[i] <- as.numeric(names(edge_map_i)[1])

      # Extract tipward node scaled ACE values as the first mapped values on the edge
      # Discrepancy with actual tipward node ACE values as this is the expected value for the mean age of the last segment of the edge...
      present_edges_df$tipward_node_scaled_ACE[i] <- as.numeric(names(edge_map_i)[length(edge_map_i)])

      ## Convert the scaled ACE values at nodes into initial values

      present_edges_df$rootward_node_ACE[i] <- trans[as.character(present_edges_df$rootward_node_scaled_ACE[i])]
      present_edges_df$tipward_node_ACE[i] <- trans[as.character(present_edges_df$tipward_node_scaled_ACE[i])]

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

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ACE_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
    }

    ## Update contMaps if needed
    # Not needed for STRAPP test. Useful only for visualization.
    if (update_contMaps)
    {
      ## Cut contMap$tree at focal time and update trait mapping in contMap$tree$maps and contMap$tree$mapped.edge
      updated_contMaps <- cut_contMaps_for_focal_time(contMaps = contMaps, focal_time = focal_time, keep_tip_labels = keep_tip_labels)
    }

    ## Export outputs
    if (!update_contMaps)
    {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "continuous"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "continuous", contMaps = updated_contMaps))
    }
  }
}


### Sub-function for categorical trait data ####

#' @title Extract categorical trait data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely states found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, the mapped trait data are summarized as posterior probabilities of each state in `densityMaps`.
#'   Optionally, the function can update the`densityMaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements in each `densityMap`.
#'
#' @param densityMaps List of objects of class `"densityMap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param ace (Optional) Numerical matrix that record the posterior probabilities of ancestral states at internal nodes,
#'   obtained with [deepSTRAPP::prepare_trait_data()] as output in the `$ace` slot.
#'   Rows are internal nodes_ID. Columns are states. Values are posterior probabilities of each state per node.
#'   Needed to provide accurate estimates of ancestral states.
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
#' @importFrom phytools nodeHeights plot.densityMap
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The mapped phylogeny (`densityMaps`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Most likely states are extracted from the posterior probabilities displayed in the `densityMaps`.
#'   The state with the highest probability is assigned to each tip and cut branches at `focal_time`.
#'
#'   True ML estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
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
#'   The categorical trait mapping in `densityMap` (`$tree$maps` and `$tree$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named character string vector with ML states found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
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
#' Associated main function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
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
#' eel_cat_3lvl_data_10My <- extract_most_likely_states_from_densityMaps_for_focal_time(
#'    densityMaps = eel_cat_3lvl_data$densityMaps,
#'    # ace = eel_cat_3lvl_data$ace,
#'    focal_time = focal_time,
#'    update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_3lvl_data_10My, 1)
#' eel_cat_3lvl_data_10My$trait_data
#'
#' ## Plot density maps as overlay of all state posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' plot_densityMaps_overlay(densityMaps = eel_cat_3lvl_data$densityMaps, fsize = 0.5)
#' abline(v = max(phytools::nodeHeights(eel_cat_3lvl_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(eel_cat_3lvl_data_10My$densityMaps, fsize = 0.7) }
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
#'                                        plot_map = FALSE)
#'
#' # Set focal time
#' focal_time <- 80
#'
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' mammals_cat_data_80My <- extract_most_likely_states_from_densityMaps_for_focal_time(
#'     densityMaps = mammals_cat_data$densityMaps,
#'     focal_time = focal_time,
#'     update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(mammals_cat_data_80My, 1)
#' mammals_cat_data_80My$trait_data
#'
#' ## Plot density maps as overlay of all state posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' plot_densityMaps_overlay(densityMaps = mammals_cat_data$densityMaps, fsize = 0.7)
#' abline(v = max(phytools::nodeHeights(mammals_cat_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(mammals_cat_data_80My$densityMaps, fsize = 0.8) }
#'

extract_most_likely_states_from_densityMaps_for_focal_time <- function (
    densityMaps,
    ace = NULL,
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

    ## ace
    if (!is.null(ace))
    {
      # ace must be a numerical matrix
      if (!is.matrix(ace))
      {
        stop(paste0("For categorical traits, 'ace' must be a numerical matrix that provides posterior probability of each state per internal nodes.\n",
                    "The object you provided is not a matrix."))
      }
      # ace should have as many rows as there are internal nodes in the densityMaps$tree
      if (nrow(ace) != densityMaps[[1]]$tree$Nnode)
      {
        stop(paste0("'ace' should have as many rows as there are internal nodes in the densityMaps[[i]]$tree.\n",
                    "Number of rows in 'ace' = ",nrow(ace),"; number of internal nodes in the densityMaps[[1]]$tree = ",densityMaps[[1]]$tree$Nnode,"."))
      }
      internal_nodes_ID <- (length(densityMaps[[1]]$tree$tip.label) + 1):(length(densityMaps[[1]]$tree$tip.label) + densityMaps[[1]]$tree$Nnode)
      # row.names(ace) = internal node IDs
      if (!all(as.numeric(row.names(ace)) %in% internal_nodes_ID))
      {
        stop(paste0("'row.names(ace)' should match numerical ID of internal nodes in the densityMaps[[i]]$tree."))
      }
      if (!all(as.numeric(row.names(ace)) == internal_nodes_ID))
      {
        warning(paste0("Rows in 'ace' are not ordered in increasing numerical ID of internal nodes.\n",
                       "They were reordered to follow the numerical ID of internal nodes."))
      }
      # ace should have as many columns as there are densityMaps associated to each state
      if (ncol(ace) != length(densityMaps))
      {
        stop(paste0("'ace' should have as many columns as there are states = objects in the densityMaps.\n",
                    "Number of columns in 'ace' = ",ncol(ace),"; number of states = objects in the 'densityMaps' = ",length(densityMaps),"."))
      }
      # ace columns should match states
      if (!all(colnames(ace) %in% states_list))
      {
        stop(paste0("'ace' column names should match the states in the densityMaps.\n",
                    "Column names in 'ace' = ",paste(colnames(ace), collapse = ", "),".\n",
                    "States in 'densityMaps' = ",paste(states_list, collapse = ", "),"."))
      }
      # ace columns should match ordered states
      if (!identical(colnames(ace), states_list))
      {
        warning(paste0("'ace' columns should match the order of states in the densityMaps.\n",
                       "They were reordered to follow the order of states in the densityMaps."))
      }
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

  ## Warn against not providing ace and tip_data
  if (is.null(ace))
  {
    cat(paste0("WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using most likely states extracted from the densityMaps instead.\n"))
  }
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

  ## Extract node states if provided with 'ace'
  if (!is.null(ace))
  {
    internal_nodes_ID <- (length(densityMaps[[1]]$tree$tip.label) + 1):(length(densityMaps[[1]]$tree$tip.label) + densityMaps[[1]]$tree$Nnode)
    # Reorder ace rows according to their internal node index
    ace <- ace[match(x = row.names(ace), table = internal_nodes_ID), ]

    # Reorder ace columns as states in densityMaps
    ace <- ace[ , states_list]

    # Extract most likely state
    ace_state_ID <- apply(X = ace, MARGIN = 1, FUN = which.max)
    ace_states <- states_list[ace_state_ID]
    names(ace_states) <- names(ace_state_ID)

    # Use them only if a focal_time match exactly a node age
    node_data_is_provided <- T
  } else {
    node_data_is_provided <- F
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

    # Initiate field for ACE = state with the highest posterior probability at focal time
    present_edges_df$ML_state_at_focal_time <- NA

    # Loop per edge
    for (i in 1:nrow(present_edges_df))
    {
      # i <- 1

      ## Extract posterior probabilities at focal time from densityMaps

      # Extract edge ID
      edge_ID_i <- as.numeric(present_edges_df$edge_ID[i])

      # Extract associated edge mappings across states
      edge_maps_i <- lapply(X = densityMaps, FUN = function (x) { x$tree$maps[[edge_ID_i]] } )

      # Compute rootward ages of segments
      segment_rootward_ages_i <- rev(cumsum(rev(edge_maps_i[[1]])) + present_edges_df$tipward_node_age[i])
      # Identify segment matching the given focal time
      if (all(!(segment_rootward_ages_i < focal_time)))
      {
        # Case where all rootward ages are lower than focal_time, then focal_segment is the last one
        focal_segment_ID <- length(segment_rootward_ages_i)
      } else {
        # Otherwise focal_segment is the last to have a rootward age > to focal_time
        # focal_segment_ID <- which.max(segment_rootward_ages_i < focal_time) - 1
        focal_segment_ID <- which.min(segment_rootward_ages_i >= focal_time) - 1
      }

      # Extract posterior probability for focal segments
      edge_PP_i <- as.numeric(unlist(lapply(X = edge_maps_i, FUN = function (x) { names(x)[focal_segment_ID] } )))
      # Extract ML states as the state with the highest posterior probabilities
      ML_state_i <- states_list[which.max(edge_PP_i)]

      # Export ML state in present_edges_df
      present_edges_df$ML_state_at_focal_time[i] <- ML_state_i
    }

    ## Match states from ace and tip_data if needed to correct for possible discrepancy from the densityMaps

    if (tip_data_is_provided | node_data_is_provided)
    {
      # Build df for tips/nodes to adjust
      if (tip_data_is_provided)
      {
        tip_data_df <- as.data.frame(tip_data)
        tip_data_df$node_label <- row.names(tip_data_df)
        names(tip_data_df) <- c("state", "node_label")
      }
      if (node_data_is_provided)
      {
        ace_states_df <- as.data.frame(ace_states)
        ace_states_df$node_label <- row.names(ace_states_df)
        names(ace_states_df) <- c("state", "node_label")
      }
      if (tip_data_is_provided)
      {
        if (node_data_is_provided)
        {
          # Case with both tip and node data
          accurate_states_df <- rbind(tip_data_df[, c("node_label", "state")], ace_states_df[, c("node_label", "state")])
        } else {
          # Case with only tip data
          accurate_states_df <- tip_data_df[, c("node_label", "state")]
        }
      } else {
        # Case with only node data
        accurate_states_df <- ace_states_df[, c("node_label", "state")]
      }
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
        # Extract only matched node/tips
        accurate_states_df_to_patch <- accurate_states_df[(abs(accurate_states_df$tipward_node_age - focal_time) < 1e-05), ]
        # Replace ML_state_at_focal_time with provided tip/node data
        present_edges_df$ML_state_at_focal_time[match(x = accurate_states_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_states_df_to_patch$state
      }
    }

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ML_state_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
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
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "categorical"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "categorical", densityMaps = updated_densityMaps))
    }
  }
}


#' @title Extract categorical trait data mapped across stochastic maps at a given time in the past
#'
#' @description Extracts the most likely states found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, the mapped trait data are provided as stochastic maps simulating trait evolution in `simmaps`.
#'   Optionally, the function can update the`simmaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$maps` and `$mapped.edge` elements in each `simmap`.
#'
#' @param simmaps List of objects of classes `"phylo"` and `"simmap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains multiple evolutionary histories mapping state evolution along branches.
#'   Each object (i.e., `simmap`) corresponds to an independent evolutionary histories simulated conditioned to observed data and model fit.
#' @param ace (Optional) Numerical matrix that record the posterior probabilities of ancestral states at internal nodes,
#'   obtained with [deepSTRAPP::prepare_trait_data()] as output in the `$ace` slot.
#'   Rows are internal nodes_ID. Columns are states. Values are posterior probabilities of each state per node.
#'   Needed to provide fully accurate estimates of ancestral states.
#' @param tip_data (Optional) Named character string vector of tip states.
#'   Names are nodes_ID of the internal nodes. Needed to provide fully accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_simmaps Logical. Specify whether the mapped phylogeny (`simmaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated simmaps. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @importFrom phytools nodeHeights
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The stochastic maps (`simmaps`) are cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Most likely states are extracted as the most frequent states displayed in the `simmaps`.
#'
#'   True ML estimates of state probabilities will be used if `tip_data` and/or `ace` are provided as optional inputs.
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
#'   The categorical trait mapping in `simmaps` (`$maps` and `$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named character string vector with most frequent states found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
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
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_simmaps_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
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
#' ## Extract trait data and update simmaps for the given focal_time
#'
#' # Extract from the simmaps
#' eel_cat_3lvl_data_10My <- extract_most_likely_states_from_simmaps_for_focal_time(
#'    simmaps = eel_cat_3lvl_data$simmaps,
#'    focal_time = focal_time,
#'    ace = eel_cat_3lvl_data$ace,
#'    update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_3lvl_data_10My, 1)
#' eel_cat_3lvl_data_10My$trait_data
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
#' ## Produce simmaps using stochastic character mapping based on an equal-rates (ER) Mk model
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
#' ## Extract trait data and update simmaps for the given focal_time
#'
#' # Extract from the simmaps
#' mammals_cat_data_80My <- extract_most_likely_states_from_simmaps_for_focal_time(
#'     simmaps = mammals_cat_data$simmaps,
#'     tip_data = mammals_data,
#'     ace = mammals_cat_data$ace,
#'     focal_time = focal_time,
#'     update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(mammals_cat_data_80My, 1)
#' mammals_cat_data_80My$trait_data
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

extract_most_likely_states_from_simmaps_for_focal_time <- function (
    simmaps,
    ace = NULL,
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

    # Extract state list
    all_maps <- unlist(lapply(X = simmaps, FUN = function (x) { lapply(X = x$maps, FUN = function (x) { names(x) } ) } ))
    states_list <- unique(all_maps)
    states_list <- states_list[order(states_list)]

    ## ace
    if (!is.null(ace))
    {
      # ace must be a numerical matrix
      if (!is.matrix(ace))
      {
        stop(paste0("For categorical traits, 'ace' must be a numerical matrix that provides posterior probability of each state per internal nodes.\n",
                    "The object you provided is not a matrix."))
      }
      # ace should have as many rows as there are internal nodes in the simmaps
      if (nrow(ace) != simmaps[[1]]$Nnode)
      {
        stop(paste0("'ace' should have as many rows as there are internal nodes in the simmaps.\n",
                    "Number of rows in 'ace' = ",nrow(ace),"; number of internal nodes in the simmaps[[1]] = ",simmaps[[1]]$Nnode,"."))
      }
      internal_nodes_ID <- (length(simmaps[[1]]$tip.label) + 1):(length(simmaps[[1]]$tip.label) + simmaps[[1]]$Nnode)
      # row.names(ace) = internal node IDs
      if (!all(as.numeric(row.names(ace)) %in% internal_nodes_ID))
      {
        stop(paste0("'row.names(ace)' should match numerical ID of internal nodes in the simmaps."))
      }
      if (!all(as.numeric(row.names(ace)) == internal_nodes_ID))
      {
        warning(paste0("Rows in 'ace' are not ordered in increasing numerical ID of internal nodes.\n",
                       "They were reordered to follow the numerical ID of internal nodes."))
      }
      # ace should have as many columns as there are simmaps associated to each state
      if (ncol(ace) != length(states_list))
      {
        stop(paste0("'ace' should have as many columns as there are in the simmaps.\n",
                    "Number of columns in 'ace' = ",ncol(ace),"; number of states in the 'simmaps' = ",length(states_list),"."))
      }
      # ace columns should match states
      if (!all(colnames(ace) %in% states_list))
      {
        stop(paste0("'ace' column names should match the states in the simmaps.\n",
                    "Column names in 'ace' = ",paste(colnames(ace), collapse = ", "),".\n",
                    "States in 'simmaps' = ",paste(states_list, collapse = ", "),"."))
      }
      # ace columns should match ordered states
      if (!identical(colnames(ace), states_list))
      {
        warning(paste0("'ace' columns should match the order of states in the simmaps.\n",
                       "They were reordered to follow the order of states in the simmaps."))
      }
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
      # tip_data should have many states as there are tips in the simmaps
      if (length(tip_data) != length(simmaps[[1]]$tip.label))
      {
        stop(paste0("'tip_data' should have as many states as there are tips in the simmaps.\n",
                    "Number of states in 'tip_data' = ",length(tip_data),"; number of tips in the simmaps = ",length(simmaps[[1]]$tip.label),"."))
      }
      # names(tip_data) = simmaps[[i]]$tip.label
      if (!all(names(tip_data) %in% simmaps[[1]]$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the simmaps[[i]]$tip.label."))
      }
      if (!all(names(tip_data) == simmaps[[1]]$tip.label))
      {
        warning(paste0("States in 'tip_data' are not ordered as tip labels in the simmaps.\n",
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

  ## Warn against not providing ace and tip_data
  if (is.null(ace))
  {
    cat(paste0("WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using most frequent states extracted from the simmaps instead.\n"))
  }
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

  ## Extract node states if provided with 'ace'
  if (!is.null(ace))
  {
    internal_nodes_ID <- (length(simmaps[[1]]$tip.label) + 1):(length(simmaps[[1]]$tip.label) + simmaps[[1]]$Nnode)
    # Reorder ace rows according to their internal node index
    ace <- ace[match(x = row.names(ace), table = internal_nodes_ID), ]

    # Reorder ace columns as states in simmaps
    ace <- ace[ , states_list]

    # Extract most likely state
    ace_state_ID <- apply(X = ace, MARGIN = 1, FUN = which.max)
    ace_states <- states_list[ace_state_ID]
    names(ace_states) <- names(ace_state_ID)

    # Use them only if a focal_time match exactly a node age
    node_data_is_provided <- T
  } else {
    node_data_is_provided <- F
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

    # Compute node distances to focal time
    present_edges_df$rootward_node_dist <- abs(present_edges_df$rootward_node_age - focal_time)
    present_edges_df$tipward_node_dist <- abs(present_edges_df$tipward_node_age - focal_time)

    # Initiate field for ACE = state with the highest frequency at focal time
    present_edges_df$ML_state_at_focal_time <- NA

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

      # Identify the most frequent state across all stochastic maps for edge n°k
      ML_state_k <- names(table(states_edge_k)[order(table(states_edge_k), decreasing = T)])[1]
      # Record the most frequent state across all stochastic maps for edge n°k
      present_edges_df$ML_state_at_focal_time[k] <- ML_state_k
    }

    ## Match states from ace and tip_data if needed to correct for possible discrepancy from the simmaps

    if (tip_data_is_provided | node_data_is_provided)
    {
      # Build df for tips/nodes to adjust
      if (tip_data_is_provided)
      {
        tip_data_df <- as.data.frame(tip_data)
        tip_data_df$node_label <- row.names(tip_data_df)
        names(tip_data_df) <- c("state", "node_label")
      }
      if (node_data_is_provided)
      {
        ace_states_df <- as.data.frame(ace_states)
        ace_states_df$node_label <- row.names(ace_states_df)
        names(ace_states_df) <- c("state", "node_label")
      }
      if (tip_data_is_provided)
      {
        if (node_data_is_provided)
        {
          # Case with both tip and node data
          accurate_states_df <- rbind(tip_data_df[, c("node_label", "state")], ace_states_df[, c("node_label", "state")])
        } else {
          # Case with only tip data
          accurate_states_df <- tip_data_df[, c("node_label", "state")]
        }
      } else {
        # Case with only node data
        accurate_states_df <- ace_states_df[, c("node_label", "state")]
      }
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
        # Extract only matched node/tips
        accurate_states_df_to_patch <- accurate_states_df[(abs(accurate_states_df$tipward_node_age - focal_time) < 1e-05), ]
        # Replace ML_state_at_focal_time with provided tip/node data
        present_edges_df$ML_state_at_focal_time[match(x = accurate_states_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_states_df_to_patch$state
      }
    }

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ML_state_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
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
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "categorical"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "categorical", simmaps = updated_simmaps))
    }
  }
}


### Sub-function for biogeographic range data ####

#' @title Extract biogeographic range data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely ranges found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, the mapped range data are summarized as posterior probabilities of each range in `densityMaps`.
#'   Optionally, the function can update the`densityMaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the range mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements in each `densityMap`..
#'
#' @param densityMaps List of objects of class `"densityMap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a range. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param ace (Optional) Numerical matrix that record the posterior probabilities of ancestral ranges at internal nodes,
#'   obtained with [deepSTRAPP::prepare_trait_data()] as output in the `$ace` slot.
#'   Rows are internal nodes_ID. Columns are ranges. Values are posterior probabilities of each range per node.
#'   Needed to provide accurate estimates of ancestral ranges.
#' @param tip_data (Optional) Named character string vector of tip ranges.
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
#' @importFrom phytools nodeHeights plot.densityMap
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The mapped phylogeny (`densityMaps`) is cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Most likely ranges are extracted from the posterior probabilities displayed in the `densityMaps`.
#'   The range with the highest probability is assigned to each tip and cut branches at `focal_time`.
#'
#'   True ML estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
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
#'   The ancestral range mapping in `densityMap` (`$tree$maps` and `$tree$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named character string vector with ML ranges found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "biogeographic". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_densityMaps = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$densityMaps`.
#'
#'   * `$densityMaps` A list of objects of class `"densityMap"` that contains the updated `densityMap` of each range/range,
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
#' Associated main function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
#' @examples
#'
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
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_biogeo_data_10My <- extract_most_likely_ranges_from_densityMaps_for_focal_time(
#'    densityMaps = eel_biogeo_data$densityMaps,
#'    # ace = eel_biogeo_data$ace,
#'    focal_time = focal_time,
#'    update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data
#'
#' ## Plot density maps as overlay of all range posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' plot_densityMaps_overlay(densityMaps = eel_biogeo_data$densityMaps, fsize = 0.7)
#' abline(v = max(phytools::nodeHeights(eel_biogeo_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(eel_biogeo_data_10My$densityMaps, fsize = 0.7) }
#'
#' # ----- Example 2: Using all ranges ----- #
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update densityMaps_all_ranges for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_biogeo_data_10My <- extract_most_likely_ranges_from_densityMaps_for_focal_time(
#'   densityMaps = eel_biogeo_data$densityMaps_all_ranges,
#'   # ace = eel_biogeo_data$ace_all_ranges,
#'   focal_time = focal_time,
#'   update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data
#'
#' ## Plot density maps as overlay of all range posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' root_age <- max(phytools::nodeHeights(eel_biogeo_data$densityMaps_all_ranges[[1]]$tree)[,2])
#' plot_densityMaps_overlay(densityMaps = eel_biogeo_data$densityMaps_all_ranges, fsize = 0.7)
#' abline(v =  root_age - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(eel_biogeo_data_10My$densityMaps, fsize = 0.7) }
#'

extract_most_likely_ranges_from_densityMaps_for_focal_time <- function (
    densityMaps,
    ace = NULL,
    tip_data = NULL,
    focal_time,
    update_densityMaps = FALSE,
    keep_tip_labels = TRUE)
{
  ### Check input validity

  {
    ## densityMaps
    # Must provide densityMaps for biogeographic data
    if (is.null(densityMaps))
    {
      stop(paste0("You must provide 'densityMaps' for biogeographic data).\n",
                  "See ?deepSTRAPP::prepare_trait_data(), ?deepSTRAPP:::BSM_to_phytools_simmap(), and ?phytools::densityMap() to learn how to generate those objects."))
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
    range_list <- names(densityMaps)
    range_list <- str_remove(string = range_list, pattern = "Density_map_")

    ## ace
    if (!is.null(ace))
    {
      # ace must be a numerical matrix
      if (!is.matrix(ace))
      {
        stop(paste0("For biogeographic data, 'ace' must be a numerical matrix that provides posterior probability of each range per internal nodes.\n",
                    "The object you provided is not a matrix."))
      }
      # ace should have as many rows as there are internal nodes in the densityMaps$tree
      if (nrow(ace) != densityMaps[[1]]$tree$Nnode)
      {
        stop(paste0("'ace' should have as many rows as there are internal nodes in the densityMaps[[i]]$tree.\n",
                    "Number of rows in 'ace' = ",nrow(ace),"; number of internal nodes in the densityMaps[[1]]$tree = ",densityMaps[[1]]$tree$Nnode,"."))
      }
      internal_nodes_ID <- (length(densityMaps[[1]]$tree$tip.label) + 1):(length(densityMaps[[1]]$tree$tip.label) + densityMaps[[1]]$tree$Nnode)
      # row.names(ace) = internal node IDs
      if (!all(as.numeric(row.names(ace)) %in% internal_nodes_ID))
      {
        stop(paste0("'row.names(ace)' should match numerical ID of internal nodes in the densityMaps[[i]]$tree."))
      }
      if (!all(as.numeric(row.names(ace)) == internal_nodes_ID))
      {
        warning(paste0("Rows in 'ace' are not ordered in increasing numerical ID of internal nodes.\n",
                       "They were reordered to follow the numerical ID of internal nodes."))
      }
      # ace should have as many columns as there are densityMaps associated to each range
      if (ncol(ace) != length(densityMaps))
      {
        stop(paste0("'ace' should have as many columns as there are ranges = objects in the densityMaps.\n",
                    "Number of columns in 'ace' = ",ncol(ace),"; number of ranges = objects in the 'densityMaps' = ",length(densityMaps),"."))
      }
      # ace columns should match ranges
      if (!all(colnames(ace) %in% range_list))
      {
        stop(paste0("'ace' column names should match the ranges in the densityMaps.\n",
                    "Column names in 'ace' = ",paste(colnames(ace), collapse = ", "),".\n",
                    "Ranges in 'densityMaps' = ",paste(range_list, collapse = ", "),"."))
      }
      # ace columns should match ordered ranges
      if (!identical(colnames(ace), range_list))
      {
        warning(paste0("'ace' columns should match the order of ranges in the densityMaps.\n",
                       "They were reordered to follow the order of ranges in the densityMaps."))
      }
    }

    ## tip_data
    if (!is.null(tip_data))
    {
      # tip_data must be a named character string vector
      if (!is.character(tip_data))
      {
        stop(paste0("For biogeographic data, 'tip_data' must be a character string vector that provides ranges for tips.\n",
                    "The object you provided is not a character string vector."))
      }
      # tip_data should have many ranges as there are tips in the densityMaps[[i]]$tree
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

  ## Warn against not providing ace and tip_data
  if (is.null(ace))
  {
    cat(paste0("WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using most likely ranges extracted from the densityMaps instead.\n"))
  }
  if (is.null(tip_data))
  {
    cat(paste0("WARNING: No tip data have been provided. Using ranges extracted from the densityMaps instead.\n"))
  }

  ## Split multi-area ranges in tip_data if densityMaps only have unique areas
  if (all(nchar(range_list) == 1))
  {
    # Update tip_data to host only unique areas
    if (!is.null(tip_data))
    {
      # Split ranges at tips in unique areas
      unique_areas_in_tip_data <- strsplit(x = tip_data, split = "")
      # Select randomly an area
      tip_data_split <- unlist(lapply(X = unique_areas_in_tip_data, FUN = function (x) { y <- sample(x = 1:length(x), size = 1) ; z <- x[y] ; return(z) }))
      # Replace tip_data
      tip_data <- tip_data_split
    }
    # No need to update 'ace' as the validity check ensure ranges are the same as in the densityMaps
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

  ## Extract node ranges if provided with 'ace'
  if (!is.null(ace))
  {
    internal_nodes_ID <- (length(densityMaps[[1]]$tree$tip.label) + 1):(length(densityMaps[[1]]$tree$tip.label) + densityMaps[[1]]$tree$Nnode)
    # Reorder ace rows according to their internal node index
    ace <- ace[match(x = row.names(ace), table = internal_nodes_ID), ]

    # Reorder ace columns as ranges in densityMaps
    ace <- ace[ , range_list]

    # Extract most likely range
    ace_range_ID <- apply(X = ace, MARGIN = 1, FUN = which.max)
    ace_ranges <- range_list[ace_range_ID]
    names(ace_ranges) <- names(ace_range_ID)

    # Use them only if a focal_time match exactly a node age
    node_data_is_provided <- T
  } else {
    node_data_is_provided <- F
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

    # Initiate field for ACE = range with the highest posterior probability at focal time
    present_edges_df$ML_range_at_focal_time <- NA

    # Loop per edge
    for (i in 1:nrow(present_edges_df))
    {
      # i <- 1

      ## Extract posterior probabilities at focal time from densityMaps

      # Extract edge ID
      edge_ID_i <- as.numeric(present_edges_df$edge_ID[i])

      # Extract associated edge mappings across ranges
      edge_maps_i <- lapply(X = densityMaps, FUN = function (x) { x$tree$maps[[edge_ID_i]] } )

      # Compute rootward ages of segments
      segment_rootward_ages_i <- rev(cumsum(rev(edge_maps_i[[1]])) + present_edges_df$tipward_node_age[i])
      # Identify segment matching the given focal time
      if (all(!(segment_rootward_ages_i < focal_time)))
      {
        # Case where all rootward ages are lower than focal_time, then focal_segment is the last one
        focal_segment_ID <- length(segment_rootward_ages_i)
      } else {
        # Otherwise focal_segment is the last to have a rootward age > to focal_time
        # focal_segment_ID <- which.max(segment_rootward_ages_i < focal_time) - 1
        focal_segment_ID <- which.min(segment_rootward_ages_i >= focal_time) - 1
      }

      # Extract posterior probability for focal segments
      edge_PP_i <- as.numeric(unlist(lapply(X = edge_maps_i, FUN = function (x) { names(x)[focal_segment_ID] } )))
      # Extract ML ranges as the range with the highest posterior probabilities
      ML_range_i <- range_list[which.max(edge_PP_i)]

      # Export ML range in present_edges_df
      present_edges_df$ML_range_at_focal_time[i] <- ML_range_i
    }

    ## Match ranges from ace and tip_data if needed to correct for possible discrepancy from the densityMaps

    if (tip_data_is_provided | node_data_is_provided)
    {
      # Build df for tips/nodes to adjust
      if (tip_data_is_provided)
      {
        tip_data_df <- as.data.frame(tip_data)
        tip_data_df$node_label <- row.names(tip_data_df)
        names(tip_data_df) <- c("range", "node_label")
      }
      if (node_data_is_provided)
      {
        ace_ranges_df <- as.data.frame(ace_ranges)
        ace_ranges_df$node_label <- row.names(ace_ranges_df)
        names(ace_ranges_df) <- c("range", "node_label")
      }
      if (tip_data_is_provided)
      {
        if (node_data_is_provided)
        {
          # Case with both tip and node data
          accurate_ranges_df <- rbind(tip_data_df[, c("node_label", "range")], ace_ranges_df[, c("node_label", "range")])
        } else {
          # Case with only tip data
          accurate_ranges_df <- tip_data_df[, c("node_label", "range")]
        }
      } else {
        # Case with only node data
        accurate_ranges_df <- ace_ranges_df[, c("node_label", "range")]
      }
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
        # Extract only matched node/tips
        accurate_ranges_df_to_patch <- accurate_ranges_df[(abs(accurate_ranges_df$tipward_node_age - focal_time) < 1e-05), ]
        # Replace ML_range_at_focal_time with provided tip/node data
        present_edges_df$ML_range_at_focal_time[match(x = accurate_ranges_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_ranges_df_to_patch$range
      }
    }

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ML_range_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
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
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "biogeographic"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "biogeographic", densityMaps = updated_densityMaps))
    }
  }
}

#' @title Extract range data mapped across biogeographic stochastic maps at a given time in the past
#'
#' @description Extracts the most likely ranges found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   As input, the mapped trait data are provided as stochastic maps simulating trait evolution in `simmaps`.
#'   Optionally, the function can update the`simmaps`
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$maps` and `$mapped.edge` elements in each `simmap`.
#'
#' @param simmaps List of objects of classes `"phylo"` and `"simmap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains multiple biogeographic histories mapping range evolution along branches.
#'   Each object (i.e., `simmap`) corresponds to an independent evolutionary histories simulated conditioned to observed data and model fit.
#' @param ace (Optional) Numerical matrix that record the posterior probabilities of ancestral ranges at internal nodes,
#'   obtained with [deepSTRAPP::prepare_trait_data()] as output in the `$ace` slot.
#'   Rows are internal nodes_ID. Columns are ranges. Values are posterior probabilities of each range per node.
#'   Needed to provide fully accurate estimates of ancestral ranges.
#' @param tip_data (Optional) Named character string vector of tip ranges.
#'   Names are nodes_ID of the internal nodes. Needed to provide fully accurate tip values.
#' @param focal_time Integer. The time, in terms of time distance from the present,
#'   at which the tree and mapping must be cut. It must be smaller than the root age of the phylogeny.
#' @param update_simmaps Logical. Specify whether the mapped phylogeny (`simmaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated simmaps. Default is `TRUE`.
#'   Used only if `update_Map = TRUE`.
#'
#' @importFrom phytools nodeHeights
#' @importFrom ape nodelabels
#' @importFrom dplyr left_join join_by
#'
#' @details The stochastic maps (`simmaps`) are cut at a specific time in the past
#'   (i.e. the `focal_time`) and the current trait values of the overlapping edges/branches are extracted.
#'
#'   ----- Extract `trait_data` -----
#'
#'   Most likely ranges are extracted as the most frequent ranges displayed in the `simmaps`.
#'
#'   True ML estimates of range probabilities will be used if `tip_data` and/or `ace` are provided as optional inputs.
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
#'   The categorical trait mapping in `simmaps` (`$maps` and `$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named character string vector with most frequent ranges found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
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
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_simmaps_for_focal_time()]
#'
#' Associated main function: [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
#' @examples
#' ## Load categorical trait data mapped on a phylogeny
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_biogeo_data, 1)
#' length(eel_biogeo_data$simmaps) # 100 stochastic maps: one per simulated evolutionary history
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' \donttest{ # (May take several minutes to run)
#' ## Extract trait data and update simmaps for the given focal_time
#'
#' # Extract from the simmaps
#' eel_biogeo_data_10My <- extract_most_likely_ranges_from_simmaps_for_focal_time(
#'    simmaps = eel_biogeo_data$simmaps,
#'    focal_time = focal_time,
#'    ace = eel_biogeo_data$ace_all_ranges,
#'    update_simmaps = TRUE)
#'
#' ## Print trait data
#' str(eel_biogeo_data_10My, 1)
#' eel_biogeo_data_10My$trait_data
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

extract_most_likely_ranges_from_simmaps_for_focal_time <- function (
    simmaps,
    ace = NULL,
    tip_data = NULL,
    focal_time,
    update_simmaps = FALSE,
    keep_tip_labels = TRUE)
{
  ### Check input validity

  {
    ## simmaps
    # Must provide simmaps for biogeographic data
    if (is.null(simmaps))
    {
      stop(paste0("You must provide 'simmaps' for biogeographic data).\n",
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
      stop(paste0("'simmaps' objects must have a $maps element that provides the mapping of the biogeographic range evolution on the phylogeny.\n",
                  "?deepSTRAPP::prepare_trait_data() and ?phytools::make.simmap() to learn how to generate those objects."))
    }

    # Extract range list
    all_maps <- unlist(lapply(X = simmaps, FUN = function (x) { lapply(X = x$maps, FUN = function (x) { names(x) } ) } ))
    ranges_list <- unique(all_maps)
    ranges_list <- ranges_list[order(ranges_list)]

    ## ace
    if (!is.null(ace))
    {
      # ace must be a numerical matrix
      if (!is.matrix(ace))
      {
        stop(paste0("For biogeographic data, 'ace' must be a numerical matrix that provides posterior probability of each range per internal nodes.\n",
                    "The object you provided is not a matrix."))
      }
      # ace should have as many rows as there are internal nodes in the simmaps
      if (nrow(ace) != simmaps[[1]]$Nnode)
      {
        stop(paste0("'ace' should have as many rows as there are internal nodes in the simmaps.\n",
                    "Number of rows in 'ace' = ",nrow(ace),"; number of internal nodes in the simmaps[[1]] = ",simmaps[[1]]$Nnode,"."))
      }
      internal_nodes_ID <- (length(simmaps[[1]]$tip.label) + 1):(length(simmaps[[1]]$tip.label) + simmaps[[1]]$Nnode)
      # row.names(ace) = internal node IDs
      if (!all(as.numeric(row.names(ace)) %in% internal_nodes_ID))
      {
        stop(paste0("'row.names(ace)' should match numerical ID of internal nodes in the simmaps."))
      }
      if (!all(as.numeric(row.names(ace)) == internal_nodes_ID))
      {
        warning(paste0("Rows in 'ace' are not ordered in increasing numerical ID of internal nodes.\n",
                       "They were reordered to follow the numerical ID of internal nodes."))
      }
      # ace should have as many columns as there are simmaps associated to each range
      if (ncol(ace) != length(ranges_list))
      {
        stop(paste0("'ace' should have as many columns as there are in the simmaps.\n",
                    "Number of columns in 'ace' = ",ncol(ace),"; number of ranges in the 'simmaps' = ",length(ranges_list),"."))
      }
      # ace columns should match ranges
      if (!all(colnames(ace) %in% ranges_list))
      {
        stop(paste0("'ace' column names should match the ranges in the simmaps.\n",
                    "Column names in 'ace' = ",paste(colnames(ace), collapse = ", "),".\n",
                    "ranges in 'simmaps' = ",paste(ranges_list, collapse = ", "),"."))
      }
      # ace columns should match ordered ranges
      if (!identical(colnames(ace), ranges_list))
      {
        warning(paste0("'ace' columns should match the order of ranges in the simmaps.\n",
                       "They were reordered to follow the order of ranges in the simmaps."))
      }
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
          stop(paste0("For biogeographic data, 'tip_data' must be a character string vector that provides ranges for tips.\n",
                      "The object you provided is not a character string vector."))
        }
      }
      # tip_data should have many ranges as there are tips in the simmaps
      if (length(tip_data) != length(simmaps[[1]]$tip.label))
      {
        stop(paste0("'tip_data' should have as many ranges as there are tips in the simmaps.\n",
                    "Number of ranges in 'tip_data' = ",length(tip_data),"; number of tips in the simmaps = ",length(simmaps[[1]]$tip.label),"."))
      }
      # names(tip_data) = simmaps[[i]]$tip.label
      if (!all(names(tip_data) %in% simmaps[[1]]$tip.label))
      {
        stop(paste0("'names(tip_data)' should match tip labels in the simmaps[[i]]$tip.label."))
      }
      if (!all(names(tip_data) == simmaps[[1]]$tip.label))
      {
        warning(paste0("Ranges in 'tip_data' are not ordered as tip labels in the simmaps.\n",
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

  ## Warn against not providing ace and tip_data
  if (is.null(ace))
  {
    cat(paste0("WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using most frequent ranges extracted from the simmaps instead.\n"))
  }
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

  ## Extract node ranges if provided with 'ace'
  if (!is.null(ace))
  {
    internal_nodes_ID <- (length(simmaps[[1]]$tip.label) + 1):(length(simmaps[[1]]$tip.label) + simmaps[[1]]$Nnode)
    # Reorder ace rows according to their internal node index
    ace <- ace[match(x = row.names(ace), table = internal_nodes_ID), ]

    # Reorder ace columns as ranges in simmaps
    ace <- ace[ , ranges_list]

    # Extract most likely range
    ace_range_ID <- apply(X = ace, MARGIN = 1, FUN = which.max)
    ace_ranges <- ranges_list[ace_range_ID]
    names(ace_ranges) <- names(ace_range_ID)

    # Use them only if a focal_time match exactly a node age
    node_data_is_provided <- T
  } else {
    node_data_is_provided <- F
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

    # Compute node distances to focal time
    present_edges_df$rootward_node_dist <- abs(present_edges_df$rootward_node_age - focal_time)
    present_edges_df$tipward_node_dist <- abs(present_edges_df$tipward_node_age - focal_time)

    # Initiate field for ACE = range with the highest frequency at focal time
    present_edges_df$ML_range_at_focal_time <- NA

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

      # Identify the most frequent range across all stochastic maps for edge n°k
      ML_range_k <- names(table(ranges_edge_k)[order(table(ranges_edge_k), decreasing = T)])[1]
      # Record the most frequent range across all stochastic maps for edge n°k
      present_edges_df$ML_range_at_focal_time[k] <- ML_range_k
    }

    ## Match ranges from ace and tip_data if needed to correct for possible discrepancy from the simmaps

    if (tip_data_is_provided | node_data_is_provided)
    {
      # Build df for tips/nodes to adjust
      if (tip_data_is_provided)
      {
        tip_data_df <- as.data.frame(tip_data)
        tip_data_df$node_label <- row.names(tip_data_df)
        names(tip_data_df) <- c("range", "node_label")
      }
      if (node_data_is_provided)
      {
        ace_ranges_df <- as.data.frame(ace_ranges)
        ace_ranges_df$node_label <- row.names(ace_ranges_df)
        names(ace_ranges_df) <- c("range", "node_label")
      }
      if (tip_data_is_provided)
      {
        if (node_data_is_provided)
        {
          # Case with both tip and node data
          accurate_ranges_df <- rbind(tip_data_df[, c("node_label", "range")], ace_ranges_df[, c("node_label", "range")])
        } else {
          # Case with only tip data
          accurate_ranges_df <- tip_data_df[, c("node_label", "range")]
        }
      } else {
        # Case with only node data
        accurate_ranges_df <- ace_ranges_df[, c("node_label", "range")]
      }
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
        # Extract only matched node/tips
        accurate_ranges_df_to_patch <- accurate_ranges_df[(abs(accurate_ranges_df$tipward_node_age - focal_time) < 1e-05), ]
        # Replace ML_range_at_focal_time with provided tip/node data
        present_edges_df$ML_range_at_focal_time[match(x = accurate_ranges_df_to_patch$node_label, table = present_edges_df$tip.label)] <- accurate_ranges_df_to_patch$range
      }
    }

    ## Format "trait_data" output = named vector of most likely values at focal time
    trait_data <- present_edges_df$ML_range_at_focal_time
    # names(trait_data) <- present_edges_df$edge_ID
    if (keep_tip_labels) # Names = tip.labels of tipward nodes
    {
      names(trait_data) <- present_edges_df$tip.label
    } else { # Names = tipward nodes ID
      names(trait_data) <- present_edges_df$tipward_node_ID
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
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "biogeographic"))

    } else {
      return(list(trait_data = trait_data, focal_time = focal_time, trait_data_type = "biogeographic", simmaps = updated_simmaps))
    }
  }
}


### Possible update: Make it work with non-dichotomous trees!!!

## Make unit tests for ultrametric (eel.tree / eel_contMap) and non-ultrametric trees (mammals$mammals.phy / mammals_contMap)

## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0


