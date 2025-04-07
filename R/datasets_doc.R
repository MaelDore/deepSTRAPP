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
#' @format A list with 21 elements.
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
#'   * `$MSC_index` Integer. The index of the Maximum Shift Credibility configuration among the posterior samples.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history.
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
#' @format A list with 29 elements.
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
#'   * `$MSC_index` Integer. The index of the Maximum Shift Credibility configuration among the posterior samples.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history.
#'
#'   New elements added to provide updated information:
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` in the cut tree to the `initial_node_ID` in the extant tree. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` in the cut tree to the `initial_edge_ID` in the extant tree. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'   * `$dtrates` List of three elements.
#'     + 1/ `$dtrates$tau` Numerical. Resolution factor describing the fraction of each segment length used in [BAMMtools::plot.bammdata()]
#'       compared to the full depth of the initial tree (i.e., the root_age)
#'     + 2/ `$dtrates$rates` List of two numerical vectors. Speciation and extinction rates along segments used by [BAMMtools::plot.bammdata()].
#'     + 3/ `$dtrates$tmat` Matrix of numerical. Start and end times of segments in term of distance to the root.
#'   * `$initial_colorbreaks` Vector of numerical. Diversification rate values of the percentiles delimiting the bins for mapping rates to colors with [BAMMtools::plot.bammdata()].
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
#' @format A list with 29 elements.
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
#'   * `$MSC_index` Integer. The index of the Maximum Shift Credibility configuration among the posterior samples.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history.
#'
#'   New elements added to provide updated information:
#'   * `$root_age` Integer. Stores the age of the root of the tree.
#'   * `$nodes_ID_df` Data.frame with two columns. Provides the conversion from the `new_node_ID` in the cut tree to the `initial_node_ID` in the extant tree. Each row is a node.
#'   * `$initial_nodes_ID` Vector of character strings. Provides the initial ID of internal nodes. Used to plot internal node IDs as labels with [ape::nodelabels()].
#'   * `$edges_ID_df` Data.frame with two columns. Provides the conversion from the `new_edge_ID` in the cut tree to the `initial_edge_ID` in the extant tree. Each row is an edge/branch.
#'   * `$initial_edges_ID` Vector of character strings. Provides the initial ID of edges/branches. Used to plot edge/branch IDs as labels with [ape::edgelabels()].
#'   * `$dtrates` List of three elements.
#'     + 1/ `$dtrates$tau` Numerical. Resolution factor describing the fraction of each segment length used in [BAMMtools::plot.bammdata()]
#'       compared to the full depth of the initial tree (i.e., the root_age)
#'     + 2/ `$dtrates$rates` List of two numerical vectors. Speciation and extinction rates along segments used by [BAMMtools::plot.bammdata()].
#'     + 3/ `$dtrates$tmat` Matrix of numerical. Start and end times of segments in term of distance to the root.
#'   * `$initial_colorbreaks` Vector of numerical. Diversification rate values of the percentiles delimiting the bins for mapping rates to colors with [BAMMtools::plot.bammdata()].
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
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
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
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B., (2025),
#'  Timing is everything: Evolution of ponerine ants highlights how dispersal history shapes modern biodiversity, Nature Communications.
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

### 6/ Temporary test set as output from run_deepSTRAPP_over_time to use for plotting p_values histogram

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

### 7/ Temporary BAMM output for whale phylogeny ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant whales
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant whales (Cetacea order) modeled with BAMM.
#'
#'  Source: Steeman, M. E., M. B. Hebsgaard, R. E. Fordyce, S. Y. W. Ho, D. L. Rabosky, R. Nielsen, C. Rahbek, H. Glenner, M. V. Sorensen, and E. Willerslev (2009)
#'  Radiation of extant cetaceans driven by restructuring of the oceans. Systematic Biology, 58, 573-585.
#'
#' @usage data(whale_BAMM_object)
#' @format A list with 21 elements.
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
#'   * `$MSC_index` Integer. The index of the Maximum Shift Credibility configuration among the posterior samples.
#'   * `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating information from the Phylogeny-related elements above,
#'      except `$MSP_tree$edge.length` is recording the Marginal Shift Probability of each branch (i.e., the probability of a regime shift to occur along each branch)
#'   * `$MAP_BAMM_object`. List of 18 elements of class `"bammdata" recording the mean rates and regime shift locations found across
#'      the Maximum A Posteriori probability (MAP) configuration. All BAMM elements summarizing diversification data holds a single entry describing this
#'      the mean diversification history.
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


### 8/ Template file for BAMM diversification analyses ####

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


