#### R scripts used to provide documentation for datasets ####

##### Phylogenies #####

### 1/ Ponerinae_tree ####

#' @title Dataset providing the extensive time-calibrated phylogeny of extant ponerine ants
#'
#' @description A `phylo` object describing the time-calibrated phylogeny of the 1534 extant ponerine ants (Ponerinae subfamily).
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
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
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
"Ponerinae_tree"

### 2/ Ponerinae_tree_old_calib ####

#' @title Dataset providing the extensive time-calibrated phylogeny of extant ponerine ants using an old calibration for illustrative purposes
#'
#' @description A `phylo` object describing the time-calibrated phylogeny of the 1534 extant ponerine ants (Ponerinae subfamily).
#'  This is NOT a properly time-calibrated phylogeny. It uses an ill-designed old calibration for illustrative purposes.
#'  For a proper time-calibrated phylogeny of ponerine ants, see the `Ponerinae_tree` object in deepSTRAPP.
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
#' @usage data(Ponerinae_tree_old_calib)
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
#' @name Ponerinae_tree_old_calib
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
"Ponerinae_tree_old_calib"



##### Trait datasets #####

### 3/ Mammals data from motmot R package ####

#' @title Phylogeny and body mass data for extant and extinct mammals families/genera from Slater, 2013
#'
#' @description A list containing two elements:
#'  * `$mammal.mass` A data.frame with mean and standard error of mammal body masses in ln(mass in g).
#'  * `$mammal.phy` A time-calibrated phylogeny of extinct and extant mammal families/genera.
#'
#'  Source: Slater, G. J. (2013). Phylogenetic evidence for a shift in the mode of mammalian body size evolution at the Cretaceous-Palaeogene boundary.
#'  Methods in Ecology and Evolution, 4(8), 734-744. \doi{10.1111/2041-210X.12084}
#'
#' @usage data(mammals)
#' @format A list with 2 elements.
#'  * `$mammal.mass` A data.frame of 213 observations and two columns. Values are numerical.
#'  * `$mammal.phy` A time-calibrated phylogeny of 211 tips.
#'
#' @details Initial dataset from Slater, G. J. (2013). The R object was imported from R package `motmot` to include in deepSTRAPP.
#'
#' @docType data
#' @keywords datasets
#' @name mammals
#'
"mammals"


### 4/ Ponerinae_trait_tip_data ####

#' @title Dataset providing fake trait data for extant ponerine ants for illustrative purposes
#'
#' @description A data.frame of fake trait data covering the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'  This is NOT real biological/ecological data. They were designed for illustrative purposes only.
#'
#' @usage data(Ponerinae_trait_tip_data)
#' @format A data.frame with 1534 rows and 4 columns.
#'
#' @details A data.frame of fake trait data covering the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'   * `$Taxa` Character string. Names of ponerinae ant taxa.
#'   * `$fake_cont_tip_data` Numeric. Fake continuous trait data.
#'   * `$fake_cat_2lvl_tip_data` Vector of character strings. Fake categorical size data with two levels: "large" and "small".
#'   * `$fake_cat_3lvl_tip_data` Vector of character strings. Fake categorical habitat data with three levels: "arboreal", "subterranean", and "terricolous".
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_trait_tip_data
#'
"Ponerinae_trait_tip_data"


### 5/ Ponerinae_binary_range_table ####

#' @title Dataset providing biogeographic range data for extant ponerine ants
#'
#' @description A data.frame of range location for the 1534 extant ponerine ant taxa (Ponerinae subfamily).
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
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
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
"Ponerinae_binary_range_table"


##### Trait evolution data #####

### 6/ Ponerinae_trait_cont_tip_data_10My ####

#' @title Data summarizing the evolution of a fake continuous trait in Ponerinae ants extracted for 10 Mya
#'
#' @description A list containing estimated values for a fake continuous trait mapped on the Ponerinae ant phylogeny,
#'  modeled with [geiger::fitContinuous]. This object was obtained with [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()]
#'  applied on trait evolution data obtained with [deepSTRAPP::prepare_trait_data()].
#'  The phylogeny used is NOT a properly time-calibrated phylogeny. It uses an ill-designed old calibration for illustrative purposes.
#'
#' @usage data(Ponerinae_trait_cont_tip_data_10My)
#' @format A list with 4 elements.
#'
#' @details  A list of four elements containing information on the evolution of a continuous trait in Ponerinae ants extracted for 10 Mya.
#'
#'   * `$trait_data` Named vector of numerical values. Names are the taxa or internal tipward node ID associated with the values.
#'     Values are the continuous trait data estimated along branches for 10 Mya.
#'   * `$focal_time` Numerical. Time in the past at which the trait data where extracted.
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "continuous".
#'   * `$contMap` A phylogenetic tree and associated mapping of estimated trait values. It was updated such as the tips correspond to
#'     lineages found 10 Mya (i.e., at the focal time in the past).
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_trait_cont_tip_data_10My
#'
"Ponerinae_trait_cont_tip_data_10My"


### 7/ Categorical trait evolution data for eel using 3-level factor ####

#' @title Data summarizing the evolution of feeding habits in eels using a 3-level factor as categorical trait
#'
#' @description A list containing feeding habits data of eels mapped on the phylogeny,
#'  modeled with [geiger::fitDiscrete]. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  Initial data was altered arbitrarily to create three categories, adding a "kiss" feeding habit to the initial
#'  "bite" and "suction" data. This is NOT real biological data. Please refer to the initial article for real data.
#'
#'  Original data source: Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \doi{10.1038/ncomms6505}
#'
#' @usage data(eel_cat_3lvl_data)
#' @format A list with 6 elements.
#'
#' @details A list of five objects containing information on the evolution of feeding habits in eels.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state/range along branches. The list contains one `"densityMap` per state/range found in the `tip_data`.
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "categorical".
#'   * `$simmaps` List of 100 stochastic mapping simulations for trait evolution. Each element is a `"simmap"` object (see [phytools::make.simmap])
#'     representing a possible evolutionary history that fits states observed on tips, infered ACE at internal nodes, and transition rates as estimated from the best fit model.
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral states/ranges (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   * `$best_model_fit` List that provides the output of the best fitting model (Here: ER model).
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name eel_cat_3lvl_data
#'
#' @references Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \doi{10.1038/ncomms6505}
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"eel_cat_3lvl_data"


### 8/ Categorical trait evolution data for Ponerinae ants using 2-level factor ####

#' @title Data summarizing the evolution of fake size data in Ponerinae ants using a 2-level factor as categorical trait
#'
#' @description A list containing fake size data of Ponerinae ants mapped on the phylogeny,
#'  modeled with [geiger::fitDiscrete]. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  This is NOT real biological/ecological data. They were designed for illustrative purposes only.
#'  The phylogeny used is also NOT a properly time-calibrated phylogeny. It uses an ill-designed old calibration for illustrative purposes.
#'
#' @usage data(Ponerinae_cat_2lvl_data_old_calib)
#' @format A list with 5 elements.
#'
#' @details A list of five objects containing information on the evolution of fake size data in Ponerinae ants.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of two objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state along branches. The list contains one `"densityMap` per state found in the `tip_data` (i.e., "large" and "small").
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "categorical".
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral states (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states (i.e., "large" and "small"). Values are posterior probabilities of each state per node.
#'   * `$best_model_fit` List that provides the output of the best fitting model (Here: ARD model).
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_cat_2lvl_data_old_calib
#'
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"Ponerinae_cat_2lvl_data_old_calib"


### 9/ Categorical trait evolution data for Ponerinae ants using 3-level factor ####

#' @title Data summarizing the evolution of fake habitat data in Ponerinae ants using a 3-level factor as categorical trait
#'
#' @description A list containing fake habitat data of Ponerinae ants mapped on the phylogeny,
#'  modeled with [geiger::fitDiscrete]. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  This is NOT real biological/ecological data. They were designed for illustrative purposes only.
#'  The phylogeny used is also NOT a properly time-calibrated phylogeny. It uses an ill-designed old calibration for illustrative purposes.
#'
#' @usage data(Ponerinae_cat_3lvl_data_old_calib)
#' @format A list with 5 elements.
#'
#' @details A list of five objects containing information on the evolution of fake habitat data in Ponerinae ants.
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'
#'   * `$densityMaps` List of three objects of class `"densityMap` that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state along branches. The list contains one `"densityMap` per state found in the `tip_data` (i.e., "arboreal", "subterranean", and "terricolous").
#'   * `$trait_data_type` Character string. Record the type of trait data. Here: "categorical".
#'   * `$ace` Numerical matrix that record the posterior probabilities of ancestral states (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states (i.e., "large" and "small"). Values are posterior probabilities of each state per node.
#'   * `$best_model_fit` List that provides the output of the best fitting model (Here: ARD model).
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_cat_3lvl_data_old_calib
#'
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"Ponerinae_cat_3lvl_data_old_calib"


### 10/ Biogeographic range evolution data for eel ####

#' @title Data summarizing the evolution of geographic ranges in eels
#'
#' @description A list containing (fake) geographic ranges data of eels mapped on the phylogeny,
#'  modeled with R package `BioGeoBEARS`. This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  Initial data based on feeding habits was altered to be transformed into range "A" and "B", and then adding arbitrarily multi-area "AB" ranges.
#'  This is NOT real biogeographic data. Please refer to the initial article for real data.
#'
#'  Original data source: Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. Nature Communications, 5, 5505.
#'  \doi{10.1038/ncomms6505}
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
#'  \doi{10.1038/ncomms6505}
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
"eel_biogeo_data"


### 11/ Biogeographic range evolution for Ponerinae ants between Old World (O) and New World (N) ####

#' @title Data summarizing the evolution of geographic ranges in Ponerinae ants using an old ill-calibrated phylogeny for illustrative purposes
#'
#' @description A list containing geographic ranges data of Ponerinae mapped on the old ill-calibrated phylogeny,
#'  modeled with R package `BioGeoBEARS`. Ranges are labeled between "Old World" (O) and New World (N).
#'  This object was obtained with [deepSTRAPP::prepare_trait_data()].
#'  The phylogeny used is also NOT a properly time-calibrated phylogeny. It uses an ill-designed old calibration for illustrative purposes.
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
#' @usage data(Ponerinae_biogeo_data_old_calib)
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
#'     Only unique areas (i.e., "N" and "O") are considered among the ranges. Multi-area ranges (i.e., "NO") have been split among unique ranges.
#'     Rows are internal nodes. Columns are ranges. Values are posterior probabilities of each range per node.
#'   * `$ace_all_ranges` Numerical matrix that record the posterior probabilities of ancestral ranges estimated at internal nodes.
#'     All ranges observed along branches during the simulated biogeographic histories are present (i.e., "N", "O", and "NO").
#'     Rows are internal nodes. Columns are ranges. Values are posterior probabilities of each range per node.
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_biogeo_data_old_calib
#'
#' @seealso [deepSTRAPP::prepare_trait_data()]
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
"Ponerinae_biogeo_data_old_calib"


##### BAMM data #####

### 12/ Template file for BAMM diversification analyses ####

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
#' @references Authors: Daniel Rabosky (BAMM) & Pascal Title (`{BAMMtools}`). Modified by Maël Doré for deepSTRAPP.
#' @seealso BAMM software website: \url{http://bamm-project.org/}
#'
"BAMM_template_diversification"


### 13/ BAMM output for whale phylogeny ####

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


### 14/ Ponerinae_BAMM_object ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant ponerine ants
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant ponerine ants (Ponerinae subfamily) modeled with BAMM.
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
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
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#' @seealso BAMM software website: \url{http://bamm-project.org/}
#'
"Ponerinae_BAMM_object"


### 15/ Ponerinae_BAMM_object_old_calib ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant ponerine ants based on an old time-calibration for illustrative purposes
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant ponerine ants (Ponerinae subfamily) modeled with BAMM based on an old time-calibration for illustrative purposes.
#'  For a BAMM based on the proper time-calibrated phylogeny of ponerine ants, see the `Ponerinae_BAMM_object` in deepSTRAPP.
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#'
#' @usage data(Ponerinae_BAMM_object_old_calib)
#' @format A list with 24 elements.
#'
#' @details An object of class `"bammdata"` containing information of diversification dynamics
#'   of extant ponerine ants (Ponerinae subfamily) modeled with BAMM based on an old time-calibration for illustrative purposes.
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
#' @name Ponerinae_BAMM_object_old_calib
#'
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#' @seealso BAMM software website: \url{http://bamm-project.org/}
#'
"Ponerinae_BAMM_object_old_calib"


### 16/ Ponerinae_BAMM_object_10My ####

#' @title Dataset summarizing 1000 posterior samples of BAMM for extant ponerine ants updated to 10 My
#'
#' @description An object of class `"bammdata"` containing information of diversification dynamics
#'  of extant ponerine ants (Ponerinae subfamily) modeled with BAMM and cut-off for a `focal_time` of 10 My.
#'  This object was obtained with [deepSTRAPP::update_rates_and_regimes_for_focal_time()].
#'
#'  Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
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
#' @references Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
#'  Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications, 16, 8297.
#'  \doi{10.1038/s41467-025-63709-3}
#' @seealso [deepSTRAPP::update_rates_and_regimes_for_focal_time()]
#'
#'  BAMM software website: \url{http://bamm-project.org/}
#'
"Ponerinae_BAMM_object_10My"


##### deepSTRAPP data #####

# ### 17/ Output from run_deepSTRAPP_over_time on continuous trait ####
#
# #' @title deepSTRAPP output for fake continuous trait data for old time-calibrated ponerine ant phylogeny over 0 to 40 Mya
# #'
# #' @description deepSTRAPP output for fake continuous trait data for old time-calibrated ponerine ant phylogeny over 0 to 40 Mya.
# #'   This object is the typical result of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
# #'   It summaries STRAPP test results and can optionally include trait and diversification rate data frames, STRAPP test outputs,
# #'   updated trait data objects, and updated BAMM objects, that are required fro downstream plots and analyses.
# #'
# #' @usage data(Ponerinae_deepSTRAPP_cont_old_calib_0_40)
# #' @format A list with 10 elements.
# #'
# #' @details deepSTRAPP output summarizing results of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
# #'
# #'   * `$pvalues_summary_df` Data.frame with three columns providing test stat `$estimate` and `$p_value` obtained for each time step (i.e., `$focal_time`),
# #'     that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
# #'   * `$time_steps` Numerical vector. Time steps at which the STRAPP tests were carried out in the same order as the objects returned in the output lists.
# #'   * `$trait_data_type` Character string. Specify the type of trait data. Here: "continuous".
# #'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. Here: "continuous"".
# #'   * `$rate_type` Character string. The type of diversification rates used in the tests. Here: "net_diversification".
# #'
# #'   Optional melted data.frames:
# #'   * `$trait_data_df_over_time` Data.frame with three columns providing `$trait_value` associated with each `$tip_ID` found along each time step (i.e., `$focal_time`).
# #'   * `$diversification_data_df_over_time` Data.frame with six columns providing diversification regimes (`$regime_ID`) and `$rates` sorted by `$rate_type` along tips (`$tip_ID`)
# #'     found across all posterior samples (`$BAMM_sample_ID`) over each time step (i.e., `$focal_time`).
# #'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
# #'     the evolution diversification rates across trait values over time.
# #'
# #'   Optional objects generated for each time step (i.e., `focal_time`) and ordered as in `$time_steps`:
# #'   * `$STRAPP_results_over_time` List of objects summarizing the results of the STRAPP tests
# #'     See [compute_STRAPP_test_for_focal_time()] for a detailed description of the elements in each object.
# #'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
# #'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
# #'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
# #'   * `$updated_trait_data_with_Map_over_time` List of objects containing trait data and updated `contMap`.
# #'     Updated `contMap` can be plotted with [deepSTRAPP::plot_contMap()] to display a phylogeny mapped with trait values with branches cut at each `focal_time`.
# #'   * `$updated_BAMM_objects_over_time` List of objects containing rates and regimes ID mapped on phylogeny.
# #'     Updated `BAMM_object` can be plotted with [deepSTRAPP::plot_BAMM_rates()] to display a phylogeny mapped with
# #'     diversification rates with branches cut at each `focal_time`.
# #'
# #' @docType data
# #' @keywords datasets
# #' @name Ponerinae_deepSTRAPP_cont_old_calib_0_40
# #'
# "Ponerinae_deepSTRAPP_cont_old_calib_0_40"

### 18/ Output from run_deepSTRAPP_over_time on categorical trait with 2-levels ####

#' @title deepSTRAPP output for fake categorical size data with 2-levels mapped on old time-calibrated ponerine ant phylogeny over 0 to 40 Mya
#'
#' @description deepSTRAPP output for fake categorical size data with 2-levels mapped on old time-calibrated ponerine ant phylogeny over 0 to 40 Mya.
#'   This object is the typical result of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
#'   It summaries STRAPP test results and can optionally include trait and diversification rate data frames, STRAPP test outputs,
#'   updated trait data objects, and updated BAMM objects, that are required fro downstream plots and analyses.
#'
#' @usage data(Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40)
#' @format A list with 10 elements.
#'
#' @details deepSTRAPP output summarizing results of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'   * `$pvalues_summary_df` Data.frame with three columns providing test stat `$estimate` and `$p_value` obtained for each time step (i.e., `$focal_time`),
#'     that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
#'   * `$time_steps` Numerical vector. Time steps at which the STRAPP tests were carried out in the same order as the objects returned in the output lists.
#'   * `$trait_data_type` Character string. Specify the type of trait data. Here: "categorical".
#'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. Here: "binary".
#'   * `$rate_type` Character string. The type of diversification rates used in the tests. Here: "net_diversification".
#'
#'   Optional summary df for multinominal data, if `posthoc_pairwise_tests = TRUE`:
#'   * `$pvalues_summary_df_for_posthoc_pairwise_tests` Data.frame with four or five columns providing test stat `$estimate`, `$p_value`, and `$p_value_adjusted`
#'     (if `p.adjust_method` used is not "none") for each `$pair` of states involved in post hoc Dunn's tests obtained for each time step (i.e., `$focal_time`).
#'      This data.frame can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the post hoc test results across time.
#'
#'   Optional melted data.frames:
#'   * `$trait_data_df_over_time` Data.frame with three columns providing `$trait_value` associated with each `$tip_ID` found along each time step (i.e., `$focal_time`).
#'     Set `extract_trait_data_melted_df = TRUE` to include it in the output.
#'   * `$diversification_data_df_over_time` Data.frame with six columns providing diversification regimes (`$regime_ID`) and `$rates` sorted by `$rate_type` along tips (`$tip_ID`)
#'     found across all posterior samples (`$BAMM_sample_ID`) over each time step (i.e., `$focal_time`).
#'     Set `extract_diversification_data_melted_df = TRUE` to include it in the output.
#'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
#'     the evolution diversification rates across trait values over time.
#'
#'   Optional objects generated for each time step (i.e., `focal_time`) and ordered as in `$time_steps`:
#'   * `$STRAPP_results_over_time` List of objects summarizing the results of the STRAPP tests
#'     See [compute_STRAPP_test_for_focal_time()] for a detailed description of the elements in each object.
#'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
#'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
#'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
#'   * `$updated_trait_data_with_Map_over_time` List of objects containing trait data and updated `contMap`/`densityMaps`.
#'     Updated `densityMaps` can be plotted with [deepSTRAPP::plot_densityMaps_overlay()] to display a phylogeny mapped with inferred states with branches cut at each `focal_time`.
#'   * `$updated_BAMM_objects_over_time` List of objects containing rates and regimes ID mapped on phylogeny.
#'     Updated `BAMM_object` can be plotted with [deepSTRAPP::plot_BAMM_rates()] to display a phylogeny mapped with
#'     diversification rates with branches cut at each `focal_time`.
#'
#' @docType data
#' @keywords datasets
#' @name Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40
#'
"Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40"

# ### 19/ Output from run_deepSTRAPP_over_time on categorical trait with 3-levels ####
#
# #' @title deepSTRAPP output for fake categorical habitat data with 3-levels mapped on old time-calibrated ponerine ant phylogeny over 0 to 40 Mya
# #'
# #' @description deepSTRAPP output for fake categorical habitat data with 3-levels mapped on old time-calibrated ponerine ant phylogeny over 0 to 40 Mya.
# #'   This object is the typical result of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
# #'   It summaries STRAPP test results and can optionally include trait and diversification rate data frames, STRAPP test outputs,
# #'   updated trait data objects, and updated BAMM objects, that are required fro downstream plots and analyses.
# #'
# #' @usage data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40)
# #' @format A list with 11 elements.
# #'
# #' @details deepSTRAPP output summarizing results of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
# #'
# #'   * `$pvalues_summary_df` Data.frame with three columns providing test stat `$estimate` and `$p_value` obtained for each time step (i.e., `$focal_time`),
# #'     that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
# #'   * `$time_steps` Numerical vector. Time steps at which the STRAPP tests were carried out in the same order as the objects returned in the output lists.
# #'   * `$trait_data_type` Character string. Specify the type of trait data. Here: "categorical".
# #'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. Here: "multinominal".
# #'   * `$rate_type` Character string. The type of diversification rates used in the tests. Here: "net_diversification".
# #'
# #'   Optional summary df for multinominal data, if `posthoc_pairwise_tests = TRUE`:
# #'   * `$pvalues_summary_df_for_posthoc_pairwise_tests` Data.frame with five columns providing test stat `$estimate`, `$p_value`, and `$p_value_adjusted`
# #'      for each `$pair` of states involved in post hoc Dunn's tests obtained for each time step (i.e., `$focal_time`).
# #'      This data.frame can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the post hoc test results across time.
# #'
# #'   Optional melted data.frames:
# #'   * `$trait_data_df_over_time` Data.frame with three columns providing `$trait_value` associated with each `$tip_ID` found along each time step (i.e., `$focal_time`).
# #'   * `$diversification_data_df_over_time` Data.frame with six columns providing diversification regimes (`$regime_ID`) and `$rates` sorted by `$rate_type` along tips (`$tip_ID`)
# #'     found across all posterior samples (`$BAMM_sample_ID`) over each time step (i.e., `$focal_time`).
# #'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
# #'     the evolution diversification rates across trait values over time.
# #'
# #'   Optional objects generated for each time step (i.e., `focal_time`) and ordered as in `$time_steps`:
# #'   * `$STRAPP_results_over_time` List of objects summarizing the results of the STRAPP tests
# #'     See [compute_STRAPP_test_for_focal_time()] for a detailed description of the elements in each object.
# #'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
# #'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
# #'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
# #'   * `$updated_trait_data_with_Map_over_time` List of objects containing trait data and updated `contMap`/`densityMaps`.
# #'     Updated `densityMaps` can be plotted with [deepSTRAPP::plot_densityMaps_overlay()] to display a phylogeny mapped with inferred states with branches cut at each `focal_time`.
# #'   * `$updated_BAMM_objects_over_time` List of objects containing rates and regimes ID mapped on phylogeny.
# #'     Updated `BAMM_object` can be plotted with [deepSTRAPP::plot_BAMM_rates()] to display a phylogeny mapped with
# #'     diversification rates with branches cut at each `focal_time`.
# #'
# #' @docType data
# #' @keywords datasets
# #' @name Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40
# #'
# "Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40"

# ### 20/ Output from run_deepSTRAPP_over_time on biogeographic ranges ####
#
# #' @title deepSTRAPP output for biogeographic ranges mapped on old time-calibrated ponerine ant phylogeny over 0 to 40 Mya
# #'
# #' @description deepSTRAPP output for biogeographic ranges mapped on old time-calibrated ponerine ant phylogeny over 0 to 40 Mya.
# #'   This object is the typical result of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
# #'   It summaries STRAPP test results and can optionally include trait and diversification rate data frames, STRAPP test outputs,
# #'   updated trait data objects, and updated BAMM objects, that are required fro downstream plots and analyses.
# #'
# #' @usage data(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40)
# #' @format A list with 10 elements.
# #'
# #' @details deepSTRAPP output summarizing results of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()].
# #'
# #'   * `$pvalues_summary_df` Data.frame with three columns providing test stat `$estimate` and `$p_value` obtained for each time step (i.e., `$focal_time`),
# #'     that can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the test results across time.
# #'   * `$time_steps` Numerical vector. Time steps at which the STRAPP tests were carried out in the same order as the objects returned in the output lists.
# #'   * `$trait_data_type` Character string. Specify the type of trait data. Here: "biogeographic".
# #'   * `$trait_data_type_for_stats` Character string. The type of trait data used to select statistical method. Here: "binary".
# #'   * `$rate_type` Character string. The type of diversification rates used in the tests. Here: "net_diversification".
# #'
# #'   Optional summary df for multinominal data, if `posthoc_pairwise_tests = TRUE`:
# #'   * `$pvalues_summary_df_for_posthoc_pairwise_tests` Data.frame with four or five columns providing test stat `$estimate`, `$p_value`, and `$p_value_adjusted`
# #'     (if `p.adjust_method` used is not "none") for each `$pair` of states involved in post hoc Dunn's tests obtained for each time step (i.e., `$focal_time`).
# #'      This data.frame can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot showing the evolution of the post hoc test results across time.
# #'
# #'   Optional melted data.frames:
# #'   * `$trait_data_df_over_time` Data.frame with three columns providing `$trait_value` associated with each `$tip_ID` found along each time step (i.e., `$focal_time`).
# #'   * `$diversification_data_df_over_time` Data.frame with six columns providing diversification regimes (`$regime_ID`) and `$rates` sorted by `$rate_type` along tips (`$tip_ID`)
# #'     found across all posterior samples (`$BAMM_sample_ID`) over each time step (i.e., `$focal_time`).
# #'   * Those data.frames can be passed down to [deepSTRAPP::plot_rates_through_time()] to generate a plot showing
# #'     the evolution diversification rates across trait values over time.
# #'
# #'   Optional objects generated for each time step (i.e., `focal_time`) and ordered as in `$time_steps`:
# #'   * `$STRAPP_results_over_time` List of objects summarizing the results of the STRAPP tests
# #'     See [compute_STRAPP_test_for_focal_time()] for a detailed description of the elements in each object.
# #'     Combined with `return_perm_data = TRUE`, it allows to plot the histograms of the null distributions
# #'     used to assess significance of the tests with [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()].
# #'     (for a single `focal_time`) and [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] (for multiple `time_steps`).
# #'   * `$updated_trait_data_with_Map_over_time` List of objects containing trait data and updated `contMap`/`densityMaps`.
# #'     Updated `densityMaps` can be plotted with [deepSTRAPP::plot_densityMaps_overlay()] to display a phylogeny mapped with inferred states with branches cut at each `focal_time`.
# #'   * `$updated_BAMM_objects_over_time` List of objects containing rates and regimes ID mapped on phylogeny.
# #'     Updated `BAMM_object` can be plotted with [deepSTRAPP::plot_BAMM_rates()] to display a phylogeny mapped with
# #'     diversification rates with branches cut at each `focal_time`.
# #'
# #' @docType data
# #' @keywords datasets
# #' @name Ponerinae_deepSTRAPP_biogeo_old_calib_0_40
# #'
# "Ponerinae_deepSTRAPP_biogeo_old_calib_0_40"




