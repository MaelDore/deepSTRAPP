#### R scripts used to provide documentation for datasets ####

### 1/ Ponerinae_tree ####

#' @title Dataset providing the extensive time-calibrated phylogeny of extant ponerine ants
#'
#' @description A `phylo` object describing the time-calibrated phylogeny of the 1534 extant ponerine ants (Ponerinae subfamily).
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_tree)
#' @format A `phylo` object with 4 elements.
#'
#' @details A time-calibrated phylogeny as a `phylo` object with 4 elements.
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_tree
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  \url{https://doi_of_Paper_to_provide.html}
#'
"Ponerinae_tree"


### 2/ Ponerinae_BAMM_object ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant ponerine ants
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant ponerine ants (Ponerinae subfamily) modeled with BAMM.
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_BAMM_object)
#' @format A list with 24 elements.
#'
#' @details An object of class `"bammdata"` containing information of diversification dynamics
#'   of extant ponerine ants (Ponerinae subfamily) modeled with BAMM.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'
#'   BAMM internal elements used for tree exploration:
#'   * `$begin` Vector of numerical. Absolute time since root of edge/branch start (rootward).
#'   * `$end` Vector of numerical.  Absolute time since root of edge/branch end (tipward).
#'   * `$downseq` Vector of integers. Order of node visits when using a pre-order tree traversal.
#'   * `$lastvisit` ID of the last node visited when starting from the node in the corresponding position in `$downseq`.
#'
#'   BAMM elements summarizing diversification data:
#'   * `$numberEvents` Vector of integer. Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts.
#'   * `$eventData` List of data.frames. One per posterior sample. Records shift events and macroevolutionary regimes parameters. 1st line = Background root regime.
#'   * `$eventVectors` List of integer vectors. One per posterior sample. Record regime ID per branches.
#'   * `$tipStates` List of integer vectors. One per posterior sample. Record regime ID per tips.
#'   * `$tipLambda` List of numerical vectors. One per posterior sample. Record speciation rates per tips.
#'   * `$tipMu` List of numerical vectors. One per posterior sample. Record extinction rates per tips.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of numerical. Mean tip speciation rates across all posterior configurations of tips.
#'   * `$meanTipMu` Vector of numerical. Mean tip extinction rates across all posterior configurations of tips.
#'   * `$type` Character string. Set the type of data modeled with BAMM. Here, type = "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_BAMM_object
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  \url{https://doi_of_Paper_to_provide.html}
#' @seealso BAMM software website: \url{http://bamm-project.org/}
#'
"Ponerinae_BAMM_object"


### 3/ Ponerinae_BAMM_object_10My ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant ponerine ants updated to 10 My
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant ponerine ants (Ponerinae subfamily) modeled with BAMM and cut-off for a `focal_time` of 10 My.
#'  This object was obtained with [deepSTRAPP::update_rates_and_regimes_for_focal_time()].
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_BAMM_object_10My)
#' @format A list with 32 elements.
#'
#' @details An object of class `"bammdata"` containing information of diversification dynamics
#'   of extant ponerine ants (Ponerinae subfamily) modeled with BAMM.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips. If an initial extant clade was cut-off, the tip.label is the tipward edge ID of the cut branch.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'
#'   BAMM internal elements used for tree exploration:
#'   * `$begin` Vector of numerical. Absolute time since root of edge/branch start (rootward).
#'   * `$end` Vector of numerical.  Absolute time since root of edge/branch end (tipward).
#'   * `$downseq` Vector of integers. Order of node visits when using a pre-order tree traversal.
#'   * `$lastvisit` ID of the last node visited when starting from the node in the corresponding position in `$downseq`.
#'
#'   BAMM elements summarizing diversification data updated for a `focal_time` of 10 My:
#'   * `$numberEvents` Vector of integer. Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts.
#'   * `$eventData` List of data.frames. One per posterior sample. Records shift events and macroevolutionary regimes parameters. 1st line = Background root regime.
#'   * `$eventVectors` List of integer vectors. One per posterior sample. Record regime ID per branches.
#'   * `$tipStates` List of integer vectors. One per posterior sample. Record regime ID per tips.
#'   * `$tipLambda` List of numerical vectors. One per posterior sample. Record speciation rates per tips.
#'   * `$tipMu` List of numerical vectors. One per posterior sample. Record extinction rates per tips.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of numerical. Mean tip speciation rates across all posterior configurations of tips.
#'   * `$meanTipMu` Vector of numerical. Mean tip extinction rates across all posterior configurations of tips.
#'   * `$type` Character string. Set the type of data modeled with BAMM. Here, type = "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'      whose origin is older that `focal_time`.
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history, updated for the `focal_time`.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history, updated for the `focal_time`.
#'
#'   New elements added to provide updated information:
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` in the cut tree to the `initial_node_ID` in the extant tree. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` in the cut tree to the `initial_edge_ID` in the extant tree. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'   * `$dtrates` List of three elements.
#'     + 1/ `$dtrates$tau` Numerical. Resolution factor describing the fraction of each segment length used in [deepSTRAPP::plot_BAMM_rates()]
#'       compared to the full depth of the initial tree (i.e., the root_age)
#'     + 2/ `$dtrates$rates` List of two numerical vectors. Speciation and extinction rates along segments used by [deepSTRAPP::plot_BAMM_rates()].
#'     + 3/ `$dtrates$tmat` Matrix of numerical. Start and end times of segments in term of distance to the root.
#'   * `$initial_colorbreaks` List of three vectors of numerical. Rate values of the percentiles delimiting the bins for mapping rates to colors with [BAMMtools::plot.bammdata()].
#'     Each element provides values for different type of rates (`$speciation`, `$extinction`, `$net_diversification`).
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the rates/regimes were extracted and the tree was eventually cut. Here: 10 My.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_BAMM_object_10My
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  \url{https://doi_of_Paper_to_provide.html}
#' @seealso [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#'
#'  BAMM software website: \url{http://bamm-project.org/}
#'
"Ponerinae_BAMM_object_10My"

### 4/ Ponerinae_BAMM_object_25My ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant ponerine ants updated to 25 My
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant ponerine ants (Ponerinae subfamily) modeled with BAMM and cut-off for a `focal_time` of 25 My.
#'  This object was obtained with [deepSTRAPP::update_rates_and_regimes_for_focal_time()].
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_BAMM_object_25My)
#' @format A list with 32 elements.
#'
#' @details An object of class `"bammdata"` containing information of diversification dynamics
#'   of extant ponerine ants (Ponerinae subfamily) modeled with BAMM.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips. If an initial extant clade was cut-off, the tip.label is the tipward edge ID of the cut branch.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'
#'   BAMM internal elements used for tree exploration:
#'   * `$begin` Vector of numerical. Absolute time since root of edge/branch start (rootward).
#'   * `$end` Vector of numerical.  Absolute time since root of edge/branch end (tipward).
#'   * `$downseq` Vector of integers. Order of node visits when using a pre-order tree traversal.
#'   * `$lastvisit` ID of the last node visited when starting from the node in the corresponding position in `$downseq`.
#'
#'   BAMM elements summarizing diversification data updated for a `focal_time` of 25 My:
#'   * `$numberEvents` Vector of integer. Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts.
#'   * `$eventData` List of data.frames. One per posterior sample. Records shift events and macroevolutionary regimes parameters. 1st line = Background root regime.
#'   * `$eventVectors` List of integer vectors. One per posterior sample. Record regime ID per branches.
#'   * `$tipStates` List of integer vectors. One per posterior sample. Record regime ID per tips.
#'   * `$tipLambda` List of numerical vectors. One per posterior sample. Record speciation rates per tips.
#'   * `$tipMu` List of numerical vectors. One per posterior sample. Record extinction rates per tips.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of numerical. Mean tip speciation rates across all posterior configurations of tips.
#'   * `$meanTipMu` Vector of numerical. Mean tip extinction rates across all posterior configurations of tips.
#'   * `$type` Character string. Set the type of data modeled with BAMM. Here, type = "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'      whose origin is older that `focal_time`.
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history, updated for the `focal_time`.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history, updated for the `focal_time`.
#'
#'   New elements added to provide updated information:
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` in the cut tree to the `initial_node_ID` in the extant tree. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` in the cut tree to the `initial_edge_ID` in the extant tree. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'   * `$dtrates` List of three elements.
#'     + 1/ `$dtrates$tau` Numerical. Resolution factor describing the fraction of each segment length used in [deepSTRAPP::plot_BAMM_rates()]
#'       compared to the full depth of the initial tree (i.e., the root_age)
#'     + 2/ `$dtrates$rates` List of two numerical vectors. Speciation and extinction rates along segments used by [deepSTRAPP::plot_BAMM_rates()].
#'     + 3/ `$dtrates$tmat` Matrix of numerical. Start and end times of segments in term of distance to the root.
#'   * `$initial_colorbreaks` List of three vectors of numerical. Rate values of the percentiles delimiting the bins for mapping rates to colors with [BAMMtools::plot.bammdata()].
#'     Each element provides values for different type of rates (`$speciation`, `$extinction`, `$net_diversification`).
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the rates/regimes were extracted and the tree was eventually cut. Here: 25 My.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_BAMM_object_25My
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  \url{https://doi_of_Paper_to_provide.html}
#' @seealso [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#'
#'  BAMM software website: \url{http://bamm-project.org/}
#'
"Ponerinae_BAMM_object_25My"


### 5/ Ponerinae_trait_data ####

#' @title Dataset providing head width trait data for extant ponerine ants
#'
#' @description A data.frame of head width measurements covering the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'
#'  Source: TBA
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_trait_data)
#' @format A data.frame with 1534 rows and 5 columns.
#'
#' @details A data.frame of head width measurements covering the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'   * `$Taxa` Character string. Names of ponerinae ant taxa.
#'   * `$HW` Numeric. Head width measurements in mm. Includes missing taxa with `NA`.
#'   * `$ln_HW` Numeric. Head width measurements in log-scale. Includes missing taxa with `NA`.
#'   * `$sim_ln_HW` Numeric. Head width measurements in log-scale. Missing data are replaced with randomly simulated data generated from a Brownian Motion model.
#'   * `$sim_HW` Numeric. Head width measurements in mm. Missing data are replaced with randomly simulated data generated from a Brownian Motion model.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_trait_data
#'
#' @references TBA
#'  \url{https://doi_of_Paper_to_provide.html}
#'
"Ponerinae_trait_data"


### 6/ Ponerinae_trait_data_10My ####

#' @title Dataset providing head width trait data for extant ponerine ants updated to 10 My
#'
#' @description A list containing head width trait data of extant ponerine ants (Ponerinae subfamily)
#'  mapped on the phylogeny, modeled with [phytools::contMap] and cut-off for a `focal_time` of 10 My.
#'  This updated object was obtained with [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()].
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_trait_data_10My)
#' @format A list with 4 elements.
#'
#' @details A list containing head width trait data of extant ponerine ants (Ponerinae subfamily)
#'  mapped on the phylogeny, modeled with [phytools::contMap()] and cut-off for a `focal_time` of 10 My.
#'   * `$trait_data` A named vector of numerical. Contains the head width data found on the ponerine phylogeny at time = 10 My.
#'   * `$focal_time` Numeric. Time in the past at which the trait were extracted. Here, 10 My.
#'   * `$data_type` Character string. The type of trait data. Here, "continuous".
#'   * `$contMap` An object of class that contains the updated `contMap` with branches and mapping that are younger than the 10 My cut off.
#'      It can be plotted with [phytools::plot.contMap()] to visualize the mapping of head width (log scale) evolution on
#'      the phylogeny of ponerine ants cut at time = 10 My.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_trait_data_10My
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  \url{https://doi_of_Paper_to_provide.html}
#' @seealso [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'
"Ponerinae_trait_data_10My"


### 7/ Ponerinae_binary_range_table ####

#' @title Dataset providing biogeographic range data for extant ponerine ants
#'
#' @description A data.frame of range location for the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_binary_range_table)
#' @format A data.frame with 1534 rows and 10 columns.
#'
#' @details A data.frame of range locations covering the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'   * `$Taxa` Character string. Names of ponerinae ant taxa.
#'   * `$Afrotropics` Logical. Whether the range of the taxa extends to Afrotropics.
#'   * `$Australasia` Logical. Whether the range of the taxa extends to Australasia.
#'   * `$Indomalaya` Logical. Whether the range of the taxa extends to Indomalaya.
#'   * `$Nearctic` Logical. Whether the range of the taxa extends to Nearctic.
#'   * `$Neotropics` Logical. Whether the range of the taxa extends to Neotropics.
#'   * `$Eastern_Palearctic` Logical. Whether the range of the taxa extends to Eastern Palearctic.
#'   * `$Western_Palearctic` Logical. Whether the range of the taxa extends to Western Palearctic.
#'   * `$Old_World` Logical. Whether the range of the taxa extends to the Old World: encompassing any bioregion among
#'     Afrotropics, Australasia, Indomalaya, Eastern Palearctic, or Western Palearctic.
#'   * `$New_World` Logical. Whether the range of the taxa extends to the New World: encompassing any bioregion among
#'     Nearctic, or Neotropics.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_binary_range_table
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  \url{https://doi_of_Paper_to_provide.html}
#'
"Ponerinae_binary_range_table"


### 8/ Temporary test sets as output from run_deepSTRAPP_over_time to use for plotting p_values histogram ####

#' @title Temporary test set as output from run_deepSTRAPP_over_time to use for plotting p_values histogram
#'
#' @description Temporary test set as output from run_deepSTRAPP_over_time to use for plotting p_values histogram.
#'
#' @usage data(STRAPP_tests_over_time_temp_example)
#' @format A list with X elements.
#'
#' @docType data
#' @keywords datasets
#' @name STRAPP_tests_over_time_temp_example
#'
"STRAPP_tests_over_time_temp_example"

#' @title Temporary test set as output from run_deepSTRAPP_over_time to use for plotting rates_through_time
#'
#' @description Temporary test set as output from run_deepSTRAPP_over_time to use for plotting rates_through_time.
#'
#' @usage data(STRAPP_tests_over_time_temp_example_2)
#' @format A list with X elements.
#'
#' @docType data
#' @keywords datasets
#' @name STRAPP_tests_over_time_temp_example_2
#'
"STRAPP_tests_over_time_temp_example_2"

### 9/ Temporary test sets as output from run_deepSTRAPP_over_time on continuous trait to use for plotting stuff ####

#' @title Temporary test set as output from run_deepSTRAPP_over_time on continuous trait to use for plotting rates_through_time
#'
#' @description Temporary test set as output from run_deepSTRAPP_over_time to use for plotting rates_through_time.
#'
#' @usage data(Ponerinae_deepSTRAPP_cont_0_40)
#' @format A list with X elements.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_deepSTRAPP_cont_0_40
#'
"Ponerinae_deepSTRAPP_cont_0_40"

### 10/ Temporary test sets as output from run_deepSTRAPP_over_time on categorical trait to use for plotting stuff ####

#' @title Temporary test set as output from run_deepSTRAPP_over_time on categorical trait to use for plotting rates_through_time
#'
#' @description Temporary test set as output from run_deepSTRAPP_over_time to use for plotting rates_through_time.
#'
#' @usage data(Ponerinae_deepSTRAPP_cat_0_40)
#' @format A list with X elements.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_deepSTRAPP_cat_0_40
#'
"Ponerinae_deepSTRAPP_cat_0_40"

### 11/ Temporary test sets as output from run_deepSTRAPP_over_time on biogeographic ranges to use for plotting stuff ####

#' @title Temporary test set as output from run_deepSTRAPP_over_time on categorical trait to use for plotting rates_through_time
#'
#' @description Temporary test set as output from run_deepSTRAPP_over_time to use for plotting rates_through_time.
#'
#' @usage data(Ponerinae_deepSTRAPP_biogeo_0_40)
#' @format A list with X elements.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_deepSTRAPP_biogeo_0_40
#'
"Ponerinae_deepSTRAPP_biogeo_0_40"

### 12/ BAMM output for whale phylogeny ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant whales
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant whales (Cetacea order) modeled with BAMM.
#'
#'  Source: Steeman, M. E., M. B. Hebsgaard, R. E. Fordyce, S. Y. W. Ho, D. L. Rabosky, R. Nielsen, C. Rahbek, H. Glenner, M. V. Sorensen, and E. Willerslev (2009)
#'  Radiation of extant cetaceans driven by restructuring of the oceans. Systematic Biology, 58, 573-585.
#'
#' @usage data(whale_BAMM_object)
#' @format A list with 24 elements.
#'
#' @details An object of class `"bammdata"` containing information of diversification dynamics
#'   of extant ponerine ants (Ponerinae subfamily) modeled with BAMM.
#'
#'   Phylogeny-related elements used to plot a phylogeny with [ape::plot.phylo()]:
#'   * `$edge` Matrix of integers. Defines the tree topology by providing rootward and tipward node ID of each edge.
#'   * `$Nnode` Integer. Number of internal nodes.
#'   * `$tip.label` Vector of character strings. Labels of all tips.
#'   * `$edge.length` Vector of numerical. Length of edges/branches.
#'
#'   BAMM internal elements used for tree exploration:
#'   * `$begin` Vector of numerical. Absolute time since root of edge/branch start (rootward).
#'   * `$end` Vector of numerical.  Absolute time since root of edge/branch end (tipward).
#'   * `$downseq` Vector of integers. Order of node visits when using a pre-order tree traversal.
#'   * `$lastvisit` ID of the last node visited when starting from the node in the corresponding position in `$downseq`.
#'
#'   BAMM elements summarizing diversification data:
#'   * `$numberEvents` Vector of integer. Number of events/macroevolutionary regimes (k+1) recorded in each posterior configuration. k = number of shifts.
#'   * `$eventData` List of data.frames. One per posterior sample. Records shift events and macroevolutionary regimes parameters. 1st line = Background root regime.
#'   * `$eventVectors` List of integer vectors. One per posterior sample. Record regime ID per branches.
#'   * `$tipStates` List of integer vectors. One per posterior sample. Record regime ID per tips.
#'   * `$tipLambda` List of numerical vectors. One per posterior sample. Record speciation rates per tips.
#'   * `$tipMu` List of numerical vectors. One per posterior sample. Record extinction rates per tips.
#'   * `$eventBranchSegs` List of matrix of numerical. One per posterior sample. Record regime ID per segments of branches.
#'   * `$meanTipLambda` Vector of numerical. Mean tip speciation rates across all posterior configurations of tips.
#'   * `$meanTipMu` Vector of numerical. Mean tip extinction rates across all posterior configurations of tips.
#'   * `$type` Character string. Set the type of data modeled with BAMM. Here, type = "diversification".
#'
#'   Additional elements providing key information for downstream analyses:
#'   * `$expectedNumberOfShifts` Integer. The expected number of regime shifts used to set the prior in BAMM.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_indices` Vector of integers. The indices of the Maximum A Posteriori probability (MAP) configurations among the posterior samples.
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history.
#'   * `$MSC_indices` Vector of integers. The indices of the Maximum Shift Credibility (MSC) configurations among the posterior samples.
#'   * `$MSC_BAMM_object` List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum Shift Credibility (MSC) configurations. All BAMM elements summarizing diversification data holds a single entry describing
#'      this mean diversification history.
#'
#' @docType data
#' @keywords datasets
#' @name whale_BAMM_object
#'
#' @references Steeman, M. E., M. B. Hebsgaard, R. E. Fordyce, S. Y. W. Ho, D. L. Rabosky, R. Nielsen, C. Rahbek, H. Glenner, M. V. Sorensen, and E. Willerslev (2009)
#'  Radiation of extant cetaceans driven by restructuring of the oceans. Systematic Biology, 58, 573-585.
#' @seealso BAMM software website: \url{http://bamm-project.org/}
#'
"whale_BAMM_object"


### 13/ Template file for BAMM diversification analyses ####

#' @title Template file for BAMM diversification analyses
#'
#' @description Template file for BAMM diversification analyses provided as character strings.
#'
#'  Source: bamm-2.5.0
#'  References: http://bamm-project.org/; https://github.com/macroevolution/bamm
#'
#' @usage data(BAMM_template_diversification)
#' @format A vector of 260 character strings.
#'
#' @details A vector of 260 character strings that can be displayed with `print(BAMM_template_diversification)`.
#'   This is the template used to generate the 'config_file.txt' controlling settings used for a BAMM run.
#'   It provides detailed explanations of the `additional_BAMM_settings` that can be used in [deepSTRAPP::prepare_diversification_data()]
#'   to customize the BAMM run.
#'   It is called internally by [deepSTRAPP::prepare_diversification_data()] to produce
#'   the custom 'config_file' used in the subsequent BAMM run.
#'
#' @docType data
#' @keywords datasets
#' @name BAMM_template_diversification
#'
#' @references Authors: Daniel Rabosky (BAMM) & Pascal Title ([BAMMtools]). Modified by Maël Doré for deepSTRAPP.
#' @seealso BAMM software website: \url{http://bamm-project.org/}
#'
"BAMM_template_diversification"


### 14/ Categorical trait evolution data for eel using 3-level factor ####

#' @title Data summarizing the evolution of feeding habits in eels using a 3-level factor as categorical trait
#'
#' @description A list containing feeding habits data of eels mapped on the phylogeny,
#'  modeled with [geiger::fitDiscrete]. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  Initial data was altered arbitrarily to create three categories, adding a "kiss" feeding habit to the initial
#'  "bite" and "suction" data. This is NOT real biological data. Please refer to the initial article for real data.
#'
#'  Original data source: Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \href{https://doi.org/10.1038/ncomms6505}{https://doi.org/10.1038/ncomms6505}
#'
#' @usage data(eel_cat_data)
#' @format A list with 5 elements.
#'
#' @details A list of five objects containing information on the evolution of feeding habits in eels.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state/range along branches. The list contains one `"densityMap` per state/range found in the `tip_data`.
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "categorical".
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral states/ranges (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   * `$best_model_fit` List that provides the output of the best fitting model (Here: ER model).
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name eel_cat_data
#'
#' @references Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \href{https://doi.org/10.1038/ncomms6505}{https://doi.org/10.1038/ncomms6505}
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"eel_cat_data"

### 15/ Categorical trait evolution data for Ponerinae ants using 3-level factor ####

#' @title Data summarizing the evolution of head width in Ponerinae ants using a 3-level factor as categorical trait
#'
#' @description A list containing head width data of Ponerinae ants mapped on the phylogeny,
#'  modeled with [geiger::fitDiscrete]. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  Initial data was altered arbitrarily to create three categories from a continuous trait (head width).
#'
#' @usage data(Ponerinae_cat_data)
#' @format A list with 5 elements.
#'
#' @details A list of five objects containing information on the evolution of feeding habits in eels.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state/range along branches. The list contains one `"densityMap` per state/range found in the `tip_data`.
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "categorical".
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral states/ranges (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   * `$best_model_fit` List that provides the output of the best fitting model (Here: ER model).
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_cat_data
#'
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"Ponerinae_cat_data"

### 16/ Biogeographic range evolution data for eel ####

#' @title Data summarizing the evolution of geographic ranges in eels
#'
#' @description A list containing (fake) geographic ranges data of eels mapped on the phylogeny,
#'  modeled with R package `BioGeoBEARS`. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  Initial data based on feeding habits was altered to be transformed into range "A" and "B", and then adding arbitrarily multi-area "AB" ranges.
#'  This is NOT real biogeographic data. Please refer to the initial article for real data.
#'
#'  Original data source: Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \href{https://doi.org/10.1038/ncomms6505}{https://doi.org/10.1038/ncomms6505}
#'
#' @usage data(eel_biogeo_data)
#' @format A list with 9 elements.
#'
#' @details A list of 9 elements containing information on the evolution of geographic ranges in eels.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given range along branches. The list contains only a `"densityMap` per unique areas because `split_multi_area_ranges` was set to TRUE.
#'   * `$densityMaps_all_ranges` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given range along branches. The list contains one `"densityMap` per range found along branches during the simulated biogeographic histories.
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "biogeographic".
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral ranges estimated at internal nodes.
#'     Only unique areas are considered among the ranges. Multi-area ranges have been split among unique ranges.
#'     Rows are internal nodes. Columns are ranges. Values are posterior probabilities of each range per node.
#'   * `$ace_all_ranges` Numerical matrix that record the posterior probabilities of ancestral ranges estimated at internal nodes.
#'     All ranges observed along branches during the simulated biogeographic histories are present.
#'     Rows are internal nodes. Columns are ranges. Values are posterior probabilities of each range per node.
#'   * `$BSM_output` List of two lists that contains summary information of cladogenetic (`$RES_caldo_events_tables`) and anagenetic (`$RES_ana_events_tables`) events
#'     recording across the 1000 simulations of biogeographic histories performed during Biogeographic Stochastic Mapping (BSM).
#'     Each element of the list is a data.frame recording events occurring during one simulation.
#'   * `$simmaps` List of 1000 objects of class `"simmap"`.
#'     Each simmap object is a phylogeny with one simulated biogeographic history (i.e., transitions in geographic ranges) mapped along branches.
#'   * `$best_model_fit` List that provides the output of the best fitting model (Here: DEC+J model).
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name eel_biogeo_data
#'
#' @references Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \href{https://doi.org/10.1038/ncomms6505}{https://doi.org/10.1038/ncomms6505}
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"eel_biogeo_data"


### 17/ Biogeographic range evolution for Ponerinae ants between Old World (O) and New World (N) ####

#' @title Data summarizing the evolution of geographic ranges in Ponerinae ants
#'
#' @description A list containing geographic ranges data of Ponerinae mapped on the phylogeny,
#'  modeled with R package `BioGeoBEARS`. Ranges are labeled between "Old World" (O) and New World (N).
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
#'  [https://doi_of_Paper_to_provide.html]
#'
#' @usage data(Ponerinae_biogeo_data)
#' @format A list with 6 elements.
#'
#' @details A list of five objects containing information on the evolution of feeding habits in eels.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given range along branches. The list contains only a `"densityMap` per unique areas because `split_multi_area_ranges` was set to TRUE.
#'     In this case, unique areas are "N" (= "New World") and "O" (= "Old World)
#'   * `$densityMaps_all_ranges` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given range along branches. The list contains one `"densityMap` per range found along branches during the simulated biogeographic histories.
#'     Here those ranges are "N" (= "New World"), "O" (= "Old World), and "NO" for multi-area ranges encompassing both regions.
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "biogeographic".
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral ranges estimated at internal nodes.
#'     Only unique areas are considered among the ranges. Multi-area ranges have been split among unique ranges.
#'     Rows are internal nodes. Columns are ranges. Values are posterior probabilities of each range per node.
#'   * `$ace_all_ranges` Numerical matrix that record the posterior probabilities of ancestral ranges estimated at internal nodes.
#'     All ranges observed along branches during the simulated biogeographic histories are present.
#'     Rows are internal nodes. Columns are ranges. Values are posterior probabilities of each range per node.
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_biogeo_data
#'
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"Ponerinae_biogeo_data"

