## Functions to extract trait data from a mapped phylogeny at a given focal time
# One master function to select the proper pipeline according to data type
# Three sub-functions mapping trait evolution according to data type

### Master function to select the proper sub-function according to data type ####

#' @title Extracts trait data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely trait values found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogeny (`contMap` or `densityMaps`)
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
#'
#' @param contMap For continuous trait data. Object of class `"contMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()] or [phytools::contMap()],
#'   that contains a phylogenetic tree and associated continuous trait mapping.
#'   The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
#' @param densityMaps For categorical trait or biogeographic data. List of objects of class `"densityMap"`,
#'   typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state/range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state/range. The phylogenetic tree must be rooted and fully resolved/dichotomous,
#'   but it does not need to be ultrametric (it can includes fossils).
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
#' @param update_map Logical. Specify whether the mapped phylogeny (`contMap` or `densityMaps`)
#'   provided as input should be updated for visualization and returned among the outputs. Default is `FALSE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#' @param keep_tip_labels Logical. Specify whether terminal branches with a single descendant tip
#'   must retained their initial `tip.label` on the updated contMap. Default is `TRUE`.
#'   Used only if `update_map = TRUE`.
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
#'   If providing only the `contMap` trait values at tips and internal nodes will be extracted from
#'   the mapping of the `contMap` leading to a slight discrepancy with the actual tip data
#'   and estimated ancestral character values.
#'
#'   True ML trait estimates will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   For categorical trait and biogeographic data:
#'
#'   Most likely states/ranges are extracted from the posterior probabilities displayed in the `densityMaps`.
#'   The states/ranges with the highest probability is assigned to each tip and cut branches at `focal_time`.
#'
#'   True ML states/ranges will be used if `tip_data` and/or `ace` are provided as optional inputs.
#'   In practice the discrepancy is negligible.
#'
#'   ----- Update the `contMap`/`densityMaps` -----
#'
#'   To obtain an updated `contMap`/`densityMaps` alongside the trait data, set `update_map = TRUE`.
#'   The update consists in cutting off branches and mapping that are younger than the `focal_time`.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = TRUE`,
#'       the leaf left is labeled with the tip.label of the unique descendant tip.
#'   * When a branch with a single descendant tip is cut and `keep_tip_labels = FALSE`,
#'     the leaf left is labeled with the node ID of the unique descendant tip.
#'   * In all cases, when a branch with multiple descendant tips (i.e., a clade) is cut,
#'     the leaf left is labeled with the node ID of the MRCA of the cut-off clade.
#'
#'   The mapping in `contMap`/`densityMaps` (`$tree$maps` and `$tree$mapped.edge`) is updated accordingly by removing mapping associated with the cut off branches.
#'
#'   A specific sub-function (that can be used independently) is called according to the type of trait data:
#'   * For continuous traits: [deepSTRAPP::extract_most_likely_trait_values_from_contMap_for_focal_time()]
#'   * For categorical traits: [deepSTRAPP::extract_most_likely_states_from_densityMaps_for_focal_time()]
#'   * For biogeographic ranges: [deepSTRAPP::extract_most_likely_ranges_from_densityMaps_for_focal_time()]
#'
#' @return By default, the function returns a list with three elements.
#'
#'   * `$trait_data` A named numerical vector with ML trait values found along branches overlapping the `focal_time`. Names are the tip.label/tipward node ID.
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted.
#'   * `$trait_data_type` Character string. Define the type of trait data as "continuous", "categorical", or "biogeographic". Used in downstream analyses to select appropriate statistical processing.
#'
#'   If `update_map = TRUE`, the output is a list with four elements: `$trait_data`, `$focal_time`, `$trait_data_type`, and `$contMap` or `$densityMaps`.
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
#' @author Maël Doré
#'
#' @seealso [deepSTRAPP::cut_phylo_for_focal_time()] [deepSTRAPP::cut_contMap_for_focal_time()] [deepSTRAPP::cut_densityMaps_for_focal_time()]
#'
#' Associated sub-functions per type of trait data:
#'
#' [deepSTRAPP::extract_most_likely_trait_values_from_contMap_for_focal_time()]
#' [deepSTRAPP::extract_most_likely_states_from_densityMaps_for_focal_time()]
#' [deepSTRAPP::extract_most_likely_ranges_from_densityMaps_for_focal_time()]
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
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' eel_ACE <- phytools::fastAnc(tree = eel.tree, x = eel_data)
#'
#' # Run a Stochastic Mapping based on a Brownian Motion model
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
#' eel_test <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = eel_contMap,
#'    focal_time = focal_time,
#'    update_map = TRUE)
#' # Extract from tip data and ML estimates of ancestral characters (values are true ML estimates)
#' eel_test <- extract_most_likely_trait_values_for_focal_time(
#'    contMap = eel_contMap,
#'    ace = eel_ACE, tip_data = eel_data,
#'    focal_time = focal_time,
#'    update_map = TRUE)
#'
#' ## Visualize outputs
#'
#' # Print trait data
#' eel_test$trait_data
#'
#' # Plot node labels on initial stochastic map with cut-off
#' plot(eel_contMap)
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(eel_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' plot(eel_test$contMap)
#' ape::nodelabels(text = eel_test$contMap$tree$initial_nodes_ID)
#'

extract_most_likely_trait_values_for_focal_time <- function (contMap = NULL,
                                                             densityMaps = NULL,
                                                             ace = NULL,
                                                             tip_data = NULL,
                                                             trait_data_type = "continuous", # Remove entry from the wrapper. Remove argument from the sub-functions
                                                             focal_time,
                                                             update_map = FALSE, # Change for update_map in the wrapper. Adjust to update_contMap / update_densityMaps in the sub-functions
                                                             keep_tip_labels = TRUE)
{
  ### Check input validity

  ## contMap OR densityMaps
  if ((!is.null(contMap) & !is.null(densityMaps)))
  {
    stop(paste0("You must provide a 'contMap' (for continuous traits) OR a 'densityMaps' (for categorical and biogeographic traits) according to the type of trait data.\n",
                "'ace' and tip_data' can also be provided in complement to use accurate ML estimates of trait values at nodes and tips."))
  }

  ## trait_data_type
  # trait_data_type must be "continuous", categorical" or "biogeographic"
  if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
  {
    stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
  }

  ## Check that what is provided contMap OR densityMaps match the trait_data_type
  if ((!is.null(contMap) & trait_data_type != "continuous"))
  {
    stop(paste0("You provided a 'contMap' but selected '",trait_data_type,"' as 'trait_data_type'. contMaps are used to map continuous traits.\n",
                "If you wish to extract trait values for a continuous trait, provide 'trait_data_type = continuous'.\n",
                "If you wish to extract trait states/ranges for ",trait_data_type," data, provide 'densityMaps' as input instead of 'contMap'"))
  }
  if ((!is.null(densityMaps) & !(trait_data_type %in% c("categorical", "biogeographic"))))
  {
    stop(paste0("You provided 'densityMaps' but selected '",trait_data_type,"' as 'trait_data_type'. densityMaps are used to map categorical or biogeographic data.\n",
                "If you wish to extract trait states/ranges for categorical or biogeographic data, provide 'trait_data_type = categorical' or 'trait_data_type = biogeographic' accordingly.\n",
                "If you wish to extract trait values for a continuous trait, provide 'contMap' as input instead of 'densityMaps'.\n",
                "\tYou can also provide 'ace' and tip_data' in complement to a 'contMap' to use accurate ML estimates of trait values at nodes and tips."))
  }

  ## Compute the appropriate sub-function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous = { # Case for continuous data
           # Input = contMap. Trait values are interpolated along branches.
           trait_data_extract <- extract_most_likely_trait_values_from_contMap_for_focal_time(
             contMap = contMap,
             ace = ace,
             tip_data = tip_data,
             focal_time = focal_time,
             update_contMap = update_map,
             keep_tip_labels = keep_tip_labels
           )
         },
         categorical = { # Case for categorical data
           # Input = densityMaps. Based on stochastic simulations of trait states.
           trait_data_extract <- extract_most_likely_states_from_densityMaps_for_focal_time(
             densityMaps = densityMaps,
             ace = ace,
             tip_data = tip_data,
             focal_time = focal_time,
             update_densityMaps = update_map, # Change for update_map in the wrapper. Adjust to update_contMap / update_densityMaps in the sub-functions
             keep_tip_labels = keep_tip_labels
           )
         },
         biogeographic = { # Case for biogeographic data
           # Input = densityMaps. Based on stochastic simulations of geographic ranges.
           trait_data_extract <- extract_most_likely_ranges_from_densityMaps_for_focal_time(
             densityMaps = densityMaps,
             ace = ace,
             tip_data = tip_data,
             focal_time = focal_time,
             update_densityMaps = update_map, # Change for update_map in the wrapper. Adjust to update_contMap / update_densityMaps in the sub-functions
             keep_tip_labels = keep_tip_labels
           )
         }
  )

  ## Export the output
  return(invisible(trait_data_extract))
}


### Sub-function for continuous trait data ####

#' @title Extracts continuous trait data mapped on a phylogeny at a given time in the past
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
#' @export
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
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_most_likely_states_from_densityMaps_for_focal_time()]
#' [deepSTRAPP::extract_most_likely_ranges_from_densityMaps_for_focal_time()]
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
#' # Run a Stochastic Mapping based on a Brownian Motion model
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
#' plot(eel_contMap)
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(eel_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' plot(eel_test$contMap)
#' ape::nodelabels(text = eel_test$contMap$tree$initial_nodes_ID)
#'
#' # ----- Example 2: Include fossils (Non-ultrametric tree) ----- #
#'
#' ## Test with non-ultrametric trees like mammals in motmot
#'
#' ## Prepare data
#'
#' # Load mammals phylogeny and data from the R package motmot
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#'
#' library(motmot)
#' data("mammals")
#' force(mammals)
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
#' phytools::plot.contMap(mammals_contMap)
#' ape::nodelabels()
#' abline(v = max(phytools::nodeHeights(mammals_contMap$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated contMap with initial node labels
#' phytools::plot.contMap(mammals_test$contMap)
#' ape::nodelabels(text = mammals_test$contMap$tree$initial_nodes_ID)
#'

## Once datasets are included in my package, remove motmot from dependencies

## Currently, For continuous traits
# Input = contMap

## Make a different function for each type of data, then make a wrapper function for all types
# extract_most_likely_trait_values_for_focal_time() = wrapper
   # Detect type of data based on $trait_data_type, but check it match the input (contMap or densityMaps) and the type of data in ace and tip_data
# extract_most_likely_trait_values_from_contMap_for_focal_time() = For continuous traits
# extract_most_likely_state_values_from_densityMaps_for_focal_time() = For categorical traits
# extract_most_likely_range_values_from_densityMaps_for_focal_time() = For biogeographic traits

## See if better to extract most likely range/state from densityMaps instead of simmaps
# Check my scripts in preparing trait data for the Macroevolution Course

## Update the reference to extract_most_likely_trait_values_for_focal_time() to link to the wrapper function!
# Ex: in the doc of compute_STRAPP_test_for_focal_time.R

# Note that the use of ace and tip_data is useful only for contMap! densityMaps should provide accurate ace and tip_data
# Adjust doc, to explain better what is mandatory (contMap or densityMaps), and what is optional (ace and tip_data, for continuous traits)

### Possible update: Make it work with non-dichotomous trees!!!


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
                  "See ?BAMMtools::prepare_trait_data() and ?phytools::contMap() to learn how to generate those objects."))
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

  # Edge, rootward_node, tipward_node, length (once cut)

  # Get node ages per edge (no root edge)
  all_edges_df <- phytools::nodeHeights(contMap$tree)
  root_age <- max(phytools::nodeHeights(contMap$tree)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 5)) # # May be an issue for trees with very short time span
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- contMap$tree$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  # Assign tip.labels. Use tipward_node_ID for internal edges
  all_edges_df$tip.label <- contMap$tree$tip.label[match(x = all_edges_df$tipward_node_ID, 1:length(contMap$tree$tip.label))]
  all_edges_df$tip.label[is.na(all_edges_df$tip.label)] <- all_edges_df$tipward_node_ID[is.na(all_edges_df$tip.label)]

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
    present_edges_df <- all_edges_df[all_edges_df$time_test, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

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


### Sub-function for categorical trait data ####

#' @title Extracts categorical trait data mapped on a phylogeny at a given time in the past
#'
#' @description Extracts the most likely states found along branches
#'   at a specific time in the past (i.e. the `focal_time`).
#'   Optionally, the function can update the mapped phylogeny (`densityMaps`)
#'   such as branches overlapping the `focal_time` are shorten to the `focal_time`,
#'   and the trait mapping for the cut off branches are removed
#'   by updating the `$tree$maps` and `$tree$mapped.edge` elements.
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
#'   Used only if `update_map = TRUE`.
#'
#' @export
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
#'   The states with the highest probability is assigned to each tip and cut branches at `focal_time`.
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
#'   * `$densityMaps` A list of objects of class `"densityMap"` that contains the updated `densityMap` of each state/range,
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
#' Sub-functions for other types of trait data:
#'
#' [deepSTRAPP::extract_most_likely_trait_values_from_contMap_for_focal_time()]
#' [deepSTRAPP::extract_most_likely_ranges_from_densityMaps_for_focal_time()]
#'
#' @examples
#' # ----- Example 1: Only extent taxa (Ultrametric tree) ----- #
#'
#' ## Load categorical trait data mapped on a phylogeny
#' data(eel_cat_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_cat_data, 1)
#' eel_cat_data$densityMaps # Three density maps: one per state
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' ## Extract trait data and update densityMaps for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_cat_data_10My <- extract_most_likely_states_from_densityMaps_for_focal_time(
#'    densityMaps = eel_cat_data$densityMaps,
#'    focal_time = focal_time,
#'    update_densityMaps = TRUE)
#'
#' ## Print trait data
#' str(eel_cat_data_10My, 1)
#' eel_cat_data_10My$trait_data
#'
#' ## Plot density maps as overlay of all state posterior probabilities
#'
#' # Plot initial density maps with ACE pies
#' plot_densityMaps_overlay(densityMaps = eel_cat_data$densityMaps)
#' abline(v = max(phytools::nodeHeights(eel_cat_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(eel_cat_data_10My$densityMaps)
#'
#'
#' # ----- Example 2: Include fossils (Non-ultrametric tree) ----- #
#' ## Test with non-ultrametric trees like mammals in motmot
#'
#' ## Prepare data
#'
#' # Load mammals phylogeny and data from the R package motmot
#' # Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
#' library(motmot)
#'
#' data("mammals")
#' force(mammals)
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
#' ## Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
#' mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
#'                                        trait_data_type = "categorical",
#'                                        evolutionary_models = "ER",
#'                                        nb_simulations = 100)
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
#' plot_densityMaps_overlay(densityMaps = mammals_cat_data$densityMaps)
#' abline(v = max(phytools::nodeHeights(mammals_cat_data$densityMaps[[1]]$tree)[,2]) - focal_time,
#'        col = "red", lty = 2, lwd = 2)
#'
#' # Plot updated densityMaps with ACE pies
#' plot_densityMaps_overlay(mammals_cat_data_80My$densityMaps)
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
                  "See ?BAMMtools::prepare_trait_data(), ?phytools::make.simmap(), and ?phytools::densityMap() to learn how to generate those objects."))
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
    state_list <- names(densityMaps)
    state_list <- str_remove(string = state_list, pattern = "Density_map_")

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
                    "Number of rows in 'ace' = ",nrow(ace),"; number of internal nodes in the densityMaps[[1]]$tree = ",densityMaps[[1]]$Nnode,"."))
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
      if (!all(colnames(ace) %in% state_list))
      {
        stop(paste0("'ace' column names should match the states in the densityMaps.\n",
                    "Column names in 'ace' = ",paste(colnames(ace), collapse = ", "),".\n",
                    "States in 'densityMaps' = ",paste(state_list, collapse = ", "),"."))
      }
      # ace columns should match ordered states
      if (!identical(colnames(ace), state_list))
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
        stop(paste0("For categorical traits, 'tip_data' must be a character string vector that provides states for tips.\n",
                    "The object you provided is not a character string vector."))
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
  if (!is.null(ace))
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
    ace <- ace[ , state_list]

    # Extract most likely state
    ace_state_ID <- apply(X = ace, MARGIN = 1, FUN = which.max)
    ace_states <- state_list[ace_state_ID]
    names(ace_states) <- names(ace_state_ID)

    # Use them only if a focal_time match exactly a node age
    node_data_is_provided <- T
  } else {
    node_data_is_provided <- F
  }

  ## Identify edges present at focal time

  # Edge, rootward_node, tipward_node, length (once cut)

  # Get node ages per edge (no root edge)
  all_edges_df <- phytools::nodeHeights(densityMaps[[1]]$tree)
  root_age <- max(phytools::nodeHeights(densityMaps[[1]]$tree)[,2])
  all_edges_df <- as.data.frame(round(root_age - all_edges_df, 5)) # # May be an issue for trees with very short time span
  names(all_edges_df) <- c("rootward_node_age", "tipward_node_age")
  all_edges_df$edge_ID <- row.names(all_edges_df)

  # Get nodes ID per edge
  all_edges_ID_df <- densityMaps[[1]]$tree$edge
  colnames(all_edges_ID_df) <- c("rootward_node_ID", "tipward_node_ID")
  all_edges_df <- cbind(all_edges_df, all_edges_ID_df)
  all_edges_df <- all_edges_df[, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age")]

  # Assign tip.labels. Use tipward_node_ID for internal edges
  all_edges_df$tip.label <- densityMaps[[1]]$tree$tip.label[match(x = all_edges_df$tipward_node_ID, 1:length(densityMaps[[1]]$tree$tip.label))]
  all_edges_df$tip.label[is.na(all_edges_df$tip.label)] <- all_edges_df$tipward_node_ID[is.na(all_edges_df$tip.label)]

  # # Detect root node ID as the only rootward node that is not also the tipward node of any edge
  # root_node_ID <- densityMaps[[1]]$tree$edge[which.min(densityMaps[[1]]$tree$edge[, 1] %in% densityMaps[[1]]$tree$edge[, 2]), 1]

  # Identify edges present at the focal time
  all_edges_df$rootward_test <- all_edges_df$rootward_node_age > focal_time
  all_edges_df$tipward_test <- all_edges_df$tipward_node_age <= focal_time
  all_edges_df$time_test <- all_edges_df$rootward_test & all_edges_df$tipward_test

  # If no edge present, send warning
  if (sum(all_edges_df$time_test) == 0)
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
    present_edges_df <- all_edges_df[all_edges_df$time_test, c("edge_ID", "rootward_node_ID", "tipward_node_ID", "rootward_node_age", "tipward_node_age", "tip.label")]

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
      ML_state_i <- state_list[which.max(edge_PP_i)]

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


## Make unit tests for ultrametric (eel.tree / eel_contMap) and non-ultrametric trees (mammals$mammals.phy / mammals_contMap)

## Make unit tests for edge cases: focal_time > root_age; focal_time = root_age; focal_time = 0


