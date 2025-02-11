#### R scripts used to provide documentation for datasets ####

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
#' @format A list with 18 elements.
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
