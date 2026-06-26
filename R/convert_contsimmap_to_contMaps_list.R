
## Function to convert a contsimmap object into a list of contMaps

#' @title Convert a contsimmap object into a list of contMaps
#'
#' @description Convert a contsimmap object generated with [contsimmap::make.contsimmap()]
#'   into a list of `contMaps` as produced by the [phytools::contMap()] function.
#'
#' @param contsimmap List of class `"contsimmap"`, typically generated with [contsimmap::make.contsimmap()],
#'   that summarizes a set of continuous stochastic maps. Each map represents a simulation of the evolution of a continuous trait
#'   compatible with the associated evolutionary model fit on observed trait data.
#' @param verbose Logical. Whether to display conversion progress every 10 maps.
#'
#' @export
#' @importFrom grDevices rainbow
#'
#' @details A contsimmap object is typically produced by the [contsimmap::make.contsimmap()] function.
#'    It stores results of a continuous stochastic mapping simulation that represent
#'    independent evolutionary history of a continuous trait mapped along a time-calibrated phylogeny,
#'    and compatible with the associated evolutionary model fit on observed trait data.
#'    The method is described in Martin & Weber, 2026.
#'
#'    A typical `"contsimmap"` object store simulated trait values across maps in a 3D array,
#'    recording values for all time-points used to represent the continuous evolution.
#'    While `"contsimmap"` objects can record evolution of multivariate traits,
#'    deepSTRAPP currently works only with simulation of univariate trait evolution.
#'    See the [contsimmap::make.contsimmap()] help section for a detailed description of the structure of the object.
#'
#'    A `"contMap"` object typically produced by the [phytools::contMap()] function also records the evolution of a continuous
#'    trait alongside branches of a time-calibrated phylogeny. Trait values are stored in `$maps` as list of named vectors recording
#'    trait values (as names scaled between 0 and 1000) and length of the associated branch segment between two time-points (as values).
#'    Therefore, a unique `"contsimmap"` object is converted in a list of `"contMaps"`, with each item being a `"contMap"` that represents
#'    a unique evolutionary history simulation conditioned to the observed trait data and model fit.
#'
#' @return Returns a list of `"contMap"` objects.
#'    Each `"contMap"` represents a unique evolutionary history simulation conditioned to the observed trait data and model fit
#'    (i.e., a continuous stochastic map).
#'
#'    A `"contMap"` typically contains:
#'    * `$tree` A list of classes `"simmap"` and `"phylo"`. The mapped phylogeny including:
#'      * `$maps` A list of named numerical vectors. Provides the mapping of trait values along each edge (scaled between 0 and 1000).
#'      * `$mapped.edge` A numeric matrix. Provides the evolutionary time spent across trait values (columns) along the edges (rows).
#'    * `$cols` A named vector mapping values to colors used to display trait evolution along the branches.
#'    * `$lims` A numeric vector. Minimum and maximum trait values recorded during the simulation.
#'
#' @author Maël Doré
#'
#' @seealso [contsimmap::make.contsimmap()] [phytools::contMap()]
#'  [deepSTRAPP::prepare_trait_data()]
#'
#' @references For continuous stochastic mapping: Martin, B. S., & Weber, M. G. (2026). Stochastic character mapping of continuous traits on phylogenies.
#'  Systematic Biology, syag031. \doi{10.1093/sysbio/syag031}.
#'
#' @examples
#' #' if (deepSTRAPP::is_dev_version())
#' {
#'  ## The R package 'contsimmap' is needed for this example to work.
#'  # Please install it manually from: https://github.com/bstaggmartin/contsimmap.
#'
#'  # ----- Prepare data ----- #
#'
#'  # Load eel phylogeny and tip data from the R package phytools
#'  # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#'  data("eel.tree", package = "phytools")
#'  data("eel.data", package = "phytools")
#'
#'  # Extract body size
#'  eel_tip_data <- stats::setNames(eel.data$Max_TL_cm,
#'                                  rownames(eel.data))
#'
#'  # ----- Run continuous stochastic mapping ----- #
#'
#'  \donttest{ # (May take several seconds to run)
#'  eel_contsimmap <- contsimmap::make.contsimmap(
#'     tree = eel.tree,
#'     trait.data = eel_tip_data ,
#'     nsim = 100, res = 100,
#'     verbose = TRUE)
#'
#'  dim(eel_contsimmap) # 121 edges, 1 trait, 100 simulations
#'
#'  # ----- Convert to a list of contMaps ----- #
#'
#'  eel_contMaps_list <- convert_contsimmap_to_contMaps_list(
#'     contsimmap = eel_contsimmap, verbose = TRUE)
#'
#'  # ----- Explore contMaps ----- #
#'
#'  # Plot contMap from simulation n°1
#'  plot_contMap(eel_contMaps_list[[1]])
#'  # Plot contMap from simulation n°50
#'  plot_contMap(eel_contMaps_list[[50]])
#'  }
#' }
#'

convert_contsimmap_to_contMaps_list <- function (contsimmap, verbose = FALSE)
{
  ## Control for contsimmap install
  if (!requireNamespace("contsimmap", quietly = TRUE))
  {
    stop("Package 'contsimmap' is required for continuous stochastic mapping.
       Please install it manually from: https://github.com/bstaggmartin/contsimmap;
       or from the alterative deepSTRAPP repository (https://maeldore.github.io/drat).
       For instructions, see the Dependencies section on deepSTRAPP homepage (https://github.com/MaelDore/deepSTRAPP).")
  }

  ### Check input validity
  {
    ## contsimmap
    # Check that object is of class contsimmap
    if (!inherits(contsimmap, "contsimmap"))
    {
      stop(paste0("'contsimmap' must be of class 'contsimmap'.\n",
                  "Check 'contsimmap::make.contsimmap()' for information on how to produce this kind of object."))
    }

    # Check that attr(contsimmap, "maps") is available
    if (is.null(attr(contsimmap, "maps")))
    {
      stop(paste0("'contsimmap' must have a 'maps' attribute.\n",
                  "Check 'contsimmap::make.contsimmap()' for information on how to produce this kind of object."))
    }

    # Check that the contsimmap contains data for a unique trait
    if (dim(contsimmap)[2] != 1)
    {
      stop(paste0("'contsimmap' must contain continuous stochastic mapping for a single univariate trait.\n",
                  "Check 'contsimmap::make.contsimmap()' for information on how to produce this kind of object."))
    }
  }

  # Initiate final list of contMaps
  contMaps_list <- list()

  # Extract phylo
  tree <- attr(contsimmap, "tree")[[1]]
  # str(tree, 1)

  # # Need to rebuild $maps and $mapped.edge
  # tree$mapped.edge
  # tree$maps

  # Extract nb of simulations
  nb_sim <- dim(contsimmap)[3]
  nb_edges <- dim(contsimmap)[1] - 1

  # Build color default contMap scheme
  colors <- grDevices::rainbow(1001, start = 0, end = 0.7)
  names(colors) <- 0:1000

  ## Loop per simulations
  for (k in 1:nb_sim)
  {
    # k <- 2

    # Extract length of time intervals on edges
    maps_k <- attr(contsimmap, "maps")[ , 1, "dts"]
    names(maps_k) <- NULL

    # Add zero length at the beginning of each edge to match with rootward values recorded in the contsimmap
    maps_k <- lapply(X = maps_k, FUN = function (x) { c(0, x) } )

    # Extract trait values to use as names
    raw_maps_k <- contsimmap[ , 1, k][-(nb_edges + 1)]
    # str(raw_maps_k, 1)
    for (i in 1:nb_edges)
    {
      # i <- 985
      # names(maps_k[[i]]) <- raw_maps_k[[i]]
      names(maps_k[[i]]) <- raw_maps_k[[i]]$values
    }

    # Detect limit values across the simulation
    min_value_k <- min(as.numeric(names(unlist(maps_k))))
    max_value_k <- max(as.numeric(names(unlist(maps_k))))

    # Adjust values recorded as names along a 0 to 1000 scale (as in contMap from phytools)
    for (i in 1:nb_edges)
    {
      # i <- 985
      names(maps_k[[i]]) <- scale_0_1000(x = as.numeric(names(maps_k[[i]])),
                                         min_val = min_value_k, max_val = max_value_k)
      names(maps_k[[i]]) <- round(x = as.numeric(names(maps_k[[i]])), digits = 0)
    }

    # Build $mapped.edge based on $maps and $edge data
    mapped.edge_k <- makeMappedEdge(edge = tree$edge, maps = maps_k)
    mapped.edge_k <- mapped.edge_k[, order(as.numeric(colnames(mapped.edge_k)))]

    # Update $mapped.edge and $maps for this simulation k
    tree_k <- tree
    class(tree_k) <- c("simmap", setdiff(class(tree_k), "simmap"))
    tree_k$mapped.edge <- mapped.edge_k
    tree_k$maps <- maps_k
    attr(tree_k, "map.order") <- "right-to-left"

    # Build full contMap object for simulation k
    contMap_k <- list(tree = tree_k, cols = colors, lims = c(min_value_k, max_value_k))
    class(contMap_k) <- "contMap"

    # Plot contMap for simulation k
    # plot(contMap_k)

    # Build final list of contMaps
    contMaps_list <- append(x = contMaps_list, values = list(contMap_k))
    names(contMaps_list) <- paste0("contMap_", 1:length(contMaps_list))

    ## Print progress
    if (verbose & (k %% 10 == 0))
    {
      cat(paste0(Sys.time(), " - Continuous stochastic map n\u00B0",k,"/",nb_sim," converted in contMap\n"))
    }
  }

  ## Export final list of contMaps
  return(contMaps_list)
}


## Helper function to standardized trait value between 0 to 1000 (as in a phytools::contMap)
scale_0_1000 <- function (x, min_val, max_val)
{
  (x - min_val) / (max_val - min_val) * 1000
}

## Helper function to reverse standardization from 0 to 1000 scaling
unscale_0_1000 <- function (x_scaled, min_val, max_val)
{
  min_val + (x_scaled / 1000) * (max_val - min_val)
}


## Function to aggregate trait values across a list of contMap

#' @title Aggregate a list of contMaps into a unique mean/median contMap
#'
#' @description Aggregate a list of contMapsinto a unique mean/median contMap
#'   as produced by the [phytools::contMap()] function.
#'
#' @param contMaps_list List of objects of class `"contMap"` that represent independent simulations of the evolution of a continuous trait
#'   (i.e., continuous stochastic maps).
#' @param fun Character string. Select the aggregating function. Available options are `"mean"` and `"median"`. Default = `"mean"`.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   showing the evolution of the continuous trait. From lowest values to highest values.
#'   Default (`color_scale = NULL`) is using the color palette recorded in the `contMap$cols` item. If none was provided, the `rainbow()` palette is used.
#' @param display_plot Logical. Whether to plot the resulting aggregated contMap. Default = `TRUE`.
#' @param ... List of named arguments. Additional arguments to be passed down to `deepSTRAPP::plot_contMap()`.
#' @param verbose Logical. Whether to display progress every 100 edges.
#'
#' @export
#' @importFrom phytools setMap
#'
#' @details The function is primarily designed to average trait values across multiple continuous stochastic maps
#'    produced from the same simulation, typically with [contsimmap::make.contsimmap()] and then converted to a list of contMaps
#'    (phytools format) with the [convert_contsimmap_to_contMaps_list] function.
#'
#'    The result is a unique `contMap` which should be consistent with the `contMap` produced by [phytools::contMap()]
#'    when interpolating maximum likelihood estimates of ancestral trait values between nodes.
#'
#' @return Returns a unique `"contMap"` object representing the average trait evolutionary history
#'    recorded across all continuous stochastic maps provided as input in `contMaps_list`.
#'
#'    The resulting aggregated `"contMap"` is a list with three items:
#'    * `$tree` A list of classes `"simmap"` and `"phylo"`. The mapped phylogeny including:
#'      * `$maps` A list of named numerical vectors. Provides the mapping of aggregated trait values along each edge (scaled between 0 and 1000).
#'      * `$mapped.edge` A numeric matrix. Provides the evolutionary time spent across trait values (columns) along the edges (rows).
#'    * `$cols` A named vector mapping values to colors used to display trait evolution along the branches.
#'    * `$lims` A numeric vector. Minimum and maximum aggregated trait values recorded.
#'
#' @author Maël Doré
#'
#' @seealso [contsimmap::make.contsimmap()] [phytools::contMap()]
#'  [deepSTRAPP::plot_contMap()]
#'
#' @references For continuous stochastic mapping: Martin, B. S., & Weber, M. G. (2026). Stochastic character mapping of continuous traits on phylogenies.
#'  Systematic Biology, syag031. \doi{10.1093/sysbio/syag031}.
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  ## The R package 'contsimmap' is needed for this example to work.
#'  # Please install it manually from: https://github.com/bstaggmartin/contsimmap.
#'
#'  # ----- Prepare data ----- #
#'
#'  # Load eel phylogeny and tip data from the R package phytools
#'  # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
#'
#'  data("eel.tree", package = "phytools")
#'  data("eel.data", package = "phytools")
#'
#'  # Extract body size
#'  eel_tip_data <- stats::setNames(eel.data$Max_TL_cm,
#'                                  rownames(eel.data))
#'
#'  # ----- Run continuous stochastic mapping ----- #
#'
#'  \donttest{ # (May take several seconds to run)
#'  eel_contsimmap <- contsimmap::make.contsimmap(
#'     tree = eel.tree,
#'     trait.data = eel_tip_data ,
#'     nsim = 100, res = 100,
#'     verbose = TRUE)
#'
#'  dim(eel_contsimmap) # 121 edges, 1 trait, 100 simulations
#'
#'  # ----- Convert to a list of contMaps ----- #
#'
#'  eel_contMaps_list <- convert_contsimmap_to_contMaps_list(
#'     contsimmap = eel_contsimmap, verbose = TRUE)
#'
#'  # ----- Explore contMaps ----- #
#'
#'  # Plot contMap from simulation n°1
#'  plot_contMap(eel_contMaps_list[[1]])
#'  # Plot contMap from simulation n°50
#'  plot_contMap(eel_contMaps_list[[50]])
#'
#'  # ----- Aggregate all simulations ----- #
#'
#'  eel_aggregated_contMap <- aggregate_contMaps_list(
#'     contMaps_list = eel_contMaps_list,
#'     verbose = TRUE, display_plot = FALSE)
#'
#'  plot_contMap(eel_aggregated_contMap)
#'
#'  # ----- Compare to interpolated contMap from phytools ----- #
#'
#'  # Produce interpolated contMap with phytools
#'  eel_contMap <- phytools::contMap(tree = eel.tree, x = eel_tip_data,
#'                                   res = 100, # Number of time steps
#'                                   plot = FALSE)
#'
#'  # Plot the interpolated contMap mapping Maximum Likelihood estimates
#'  # of ancestral trait values interpolated between nodes
#'  plot_contMap(eel_contMap)
#'  # Plot the aggregated contMap representing average trait values
#'  # from true continuous stochastic mapping simulations
#'  plot_contMap(eel_aggregated_contMap)
#'  # Both should be visually similar
#'  }
#' }
#'

aggregate_contMaps_list <- function (contMaps_list,
                                     fun = "mean", # To select the aggregating function. Options are mean and median
                                     color_scale = NULL,
                                     display_plot = TRUE, # Whether to plot the resulting aggregated contMap
                                     ..., # Arguments to pass down to plot_contMap() if ploting is requested
                                     verbose = TRUE) # Whether to display progress every 1000 edges
{
  ### Check input validity
  {
    ## contMaps_list
    # Check that object is a list of contMap objects
    if (!all(lapply(X = contMaps_list, FUN = class) == "contMap"))
    {
      stop(paste0("'contMaps_list' must be a list of objects of class 'contMap'.\n",
                  "Check 'deepSTRAPP::convert_contsimmap_to_contMaps_list()' to learn how to produce this kind of object from 'contsimmap::make.contsimmap()' outputs."))
    }

    # All contMap objects must be built on the same phylogeny
    all_trees <- lapply(X = contMaps_list, FUN = function (x) { x$tree[c("edge", "edge.length", "Nnode", "tip.label")] })
    if (!all(unlist(lapply(X = all_trees, FUN = identical, all_trees[[1]]))))
    {
      stop(paste0("'contMaps_list' must be built on the same phylogeny.\n",
                  "Check 'deepSTRAPP::convert_contsimmap_to_contMaps_list()' to learn how to produce this kind of object from 'contsimmap::make.contsimmap()' outputs."))
    }

    # All contMap objects must be built with the same time points
    all_maps <- lapply(X = contMaps_list, FUN = function (x) { unname(unlist(x$tree$maps)) })
    if (!all(unlist(lapply(X = all_maps, FUN = identical, all_maps[[1]]))))
    {
      stop(paste0("'contMaps_list' must be built with the same time points.\n",
                  "Check 'deepSTRAPP::convert_contsimmap_to_contMaps_list()' to learn how to produce this kind of object from 'contsimmap::make.contsimmap()' outputs."))
    }

    ## fun
    # fun must be either mean or median
    if (!(fun %in% c("mean", "median")))
    {
      stop(paste0("'fun' must be either 'mean' or 'median'."))
    }
  }

  ## Extract number of edges
  nb_edges <- length(contMaps_list[[1]]$tree$maps)

  ## Loop across edges to record trait values
  aggregated_map <- list()
  for (k in 1:nb_edges)
  {
    # k <- 1

    # Initiate df to record all maps for edge k
    edge_maps_df_k <- data.frame()

    ## Loop across maps
    for (i in seq_along(contMaps_list))
    {
      # i <- 1

      # Extract trait values
      edge_k_map_i <- as.numeric(names(contMaps_list[[i]]$tree$maps[[k]]))
      # Convert trait values back to raw scale
      edge_k_map_i <- unscale_0_1000(x_scaled = edge_k_map_i,
                                     min_val = contMaps_list[[i]]$lims[1],
                                     max_val = contMaps_list[[i]]$lims[2])
      # Store values
      edge_maps_df_k <- rbind(edge_maps_df_k, edge_k_map_i)
    }

    ## Aggregate trait values using selected function
    if (fun == "mean")
    {
      aggregated_map[[k]] <- apply(X = edge_maps_df_k, MARGIN = 2, FUN = mean, na.rm = TRUE)
    } else {
      aggregated_map[[k]] <- apply(X = edge_maps_df_k, MARGIN = 2, FUN = median, na.rm = TRUE)
    }

    ## Print progress
    if (verbose & (k %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Trait values aggregated across egde n\u00B0",k,"/",nb_edges,"\n"))
    }
  }

  ## Extract new limits from aggregated trait values
  min_value <- min(unlist(aggregated_map))
  max_value <- max(unlist(aggregated_map))

  ## Rescale trait values to new limits
  aggregated_map_recaled <- lapply(X = aggregated_map, FUN = scale_0_1000,
                                   min_val = min_value, max_val = max_value)

  ## Remap averaged values across initial map
  new_map <- contMaps_list[[1]]$tree$maps
  for (k in 1:nb_edges)
  {
    # k <- 1

    names(new_map[[k]]) <- round(aggregated_map_recaled[[k]], 0)
  }

  ## Extract phylo
  tree <- contMaps_list[[1]]$tree
  # str(tree, 1)

  ## Replace the $maps
  tree$maps <- new_map

  ## Build $mapped.edge based on $maps and $edge data
  mapped.edge <- makeMappedEdge(edge = tree$edge, maps = tree$maps)
  mapped.edge <- mapped.edge[, order(as.numeric(colnames(mapped.edge)))]
  tree$mapped.edge <- mapped.edge

  # Build the final aggregated contMap object
  aggregated_contMap <- list(tree = tree, cols = contMaps_list[[1]]$cols, lims = c(min_value, max_value))
  class(aggregated_contMap) <- "contMap"

  ## Update color palette if requested
  if (!is.null(color_scale))
  {
    aggregated_contMap <- phytools::setMap(x = aggregated_contMap, colors = color_scale)
  }

  ## Plot contMap if requested
  if (display_plot)
  {
    plot_contMap(contMap = aggregated_contMap, color_scale = color_scale, display_plot = TRUE, ...)
  }

  ## Export contMap
  return(aggregated_contMap)
}



