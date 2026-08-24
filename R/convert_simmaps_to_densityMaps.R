
#' @title Convert simmaps into densityMaps suitable for deepSTRAPP
#'
#' @description Convert `simmaps` representing independent simulations of evolutionary history
#'   into a `densityMaps` object suitable to use as input for a deepSTRAPP run.
#'
#'   `simmaps` are mapping discrete character/geographic evolutionary history as
#'   transitions in character states/geographic ranges across nodes and branches of a phylogeny.
#'
#'   `densityMaps` are recording the posterior probability/frequencies of being in a given state/range along branches.
#'   In deepSTRAPP format, `densityMaps` is a list of objects of class `"densityMap"` where each object corresponds to
#'   the frequency of presence/absence of a state/range.
#'
#' @param simmaps List of objects of class `"simmap"`, typically generated with [deepSTRAPP::prepare_trait_data()] or [phytools::make.simmap()],
#'   that represent discrete character/geographic evolutionary history
#'   (i.e., transitions in character states/geographic ranges) mapped along branches.
#' @param colors_per_levels Named character string. To set the colors to use to map each state/range posterior probabilities. Names = states/ranges; values = colors.
#'   If `NULL` (default), the `rainbow()` color scale will be used.
#' @param tol Positive numerical. To set the tolerance used to match node ages and time steps (i.e., consider them equal). Default = 1e-5.
#' @param verbose Logical. Whether to display progress every 100 edges. Default = `TRUE`.
#'
#' @export
#' @importFrom grDevices rainbow colorRampPalette col2rgb rgb rgb2hsv hsv
#' @importFrom methods hasArg
#' @importFrom stats setNames
#'
#' @details The function is a wrapper of the original [phytools::densityMap()] by Liam Revell.
#'   Although, it does not produce a single `densityMap` for binary states,
#'   but rather can handle `simmaps` mapping any number of states/ranges, and produce a `densityMaps` list
#'   that summarizes the frequency of presence/absence of each state/range in subsequent `densityMap` objects.
#'
#'   Both `simmaps` and `densityMaps` objects can be used as inputs for a deepSTRAPP run
#'   with [deepSTRAPP::run_deepSTRAPP_for_focal_time()] or [deepSTRAPP::run_deepSTRAPP_over_time()].
#'
#'  `simmaps` retain the identity of each simulated history,
#'   allowing deepSTRAPP to keep track of which simulation generated each set of trait values (i.e., states/ranges).
#'   However, retaining all simulated histories can require substantial RAM, particularly for large phylogenies or many simulations.
#'
#'   Alternatively, `densityMaps` summarize the frequency of states/ranges across simulations.
#'   They require substantially less RAM and can be used to visualize the overall uncertainty
#'   in trait evolution with [deepSTRAPP::plot_densityMaps_overlay()].
#'
#'   The trade-off is that `densityMaps` discard the identity of individual simulations and
#'   therefore cannot be used to track which simulated history generated a given set of trait values.
#'
#' @return The function returns a `densityMaps` as a list objects of class `"densityMap"`,
#'   where each object is mapping the frequency of presence/absence of a state/range.
#'
#'   The number of objects depends on the number of states/ranges observed across the `simmaps` provided for conversion.
#'   Each `densityMap` is named as "Density_map_X" with X being each of the state/range recorded across the `simmaps`.
#'
#'   Each `densityMap` is a list of three elements:
#'   * `$tree` List of classes `"simmap"` and `"phylo"` that contains the phylogeny in [ape] format
#'      and the mapping of states/ranges frequencies along branches in `$tree$maps`.
#'   * `$cols` Named character strings. Colors mapped to the 0 to 1000 scale used to record frequencies in `$tree$maps`.
#'   * `$states` Character string with two values. First entry is the absence of the state/range recorded as "Not X".
#'     Second entry is the presence of the state/range recorded as "X", X being the state/range name.
#'
#' @author Maël Doré
#'
#' @seealso [phytools::densityMap()] [deepSTRAPP::plot_densityMaps_overlay()] [deepSTRAPP::prepare_trait_data()]
#'
#' @examples
#' ## Load data
#'
#' # Load trait df
#' data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#' # Load phylogeny
#' data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
#'
#' # Extract categorical data with 3-levels
#' Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#'                                          nm = Ponerinae_trait_tip_data$Taxa)
#' table(Ponerinae_cat_3lvl_tip_data)
#'
#' # Select color scheme for states
#' colors_per_states <- c("forestgreen", "sienna", "goldenrod")
#' names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
#'
#' \donttest{ # (May take several minutes to run)
#' ## Produce densityMaps using stochastic character mapping based on an ER Mk model
#' Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     phylo = Ponerinae_tree_old_calib,
#'     trait_data_type = "categorical",
#'     colors_per_levels = colors_per_states,
#'     evolutionary_models = "ER", # Use ER model
#'     run_stochastic_maps = TRUE,
#'     nb_simulations = 100, # Reduce number of simulations to save time
#'     seed = 1234, # Set seed for reproducibility
#'     return_simmaps = TRUE, # Return simmaps in the output
#'     plot_map = FALSE)
#'
#' # Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
#' # but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
#' # using [deepSTRAPP::convert_simmaps_to_densityMaps()].
#'
#' ## Convert simmaps to densityMaps as input format for deepSTRAPP
#' Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
#'    simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#'
#' # Plot densityMaps one by one
#' plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")
#' plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")
#' plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")
#'
#' # Plot overlay of all densityMaps
#' plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)
#' }
#'


convert_simmaps_to_densityMaps <- function (simmaps,
                                            colors_per_levels = NULL,
                                            tol = 1e-5,
                                            verbose = TRUE)
{
  # Get list of states
  states_list <- unique(names(unlist(lapply(simmaps, FUN = function (x) { x$maps }))))
  states_list <- states_list[order(states_list)]

  ### Check input validity
  {
    ## colors_per_levels
    if (!is.null(colors_per_levels))
    {
      # Check that the color scale match the states
      if (!all(states_list %in% names(colors_per_levels)))
      {
        missing_states <- states_list[!(states_list %in% names(colors_per_levels))]
        stop(paste0("Not all states are found in 'colors_per_levels'.\n",
                    "Missing: ", paste(missing_states, collapse = ", "), "."))
      }
      # Check whether all colors are valid
      if (!all(is_color(colors_per_levels)))
      {
        invalid_colors <- colors_per_levels[!is_color(colors_per_levels)]
        stop(paste0("Some color names in 'colors_per_levels' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }

  ## If not provided, define a color per state/range
  if (is.null(colors_per_levels))
  {
    colors_per_levels <- grDevices::rainbow(n = length(states_list), s = 0.8) # With HSV color space
    # colors_per_levels <- colorspace::rainbow_hcl(n = length(states_list), c = 90, l = 65) # With HLC color space
    names(colors_per_levels) <- states_list
  } else {
    # If provided, ensure it is properly ordered
    colors_per_levels <- colors_per_levels[match(x = states_list, table = names(colors_per_levels))]
  }

  # Initiate list of densityMaps
  densityMaps_all_states <- list()

  ## Loop per state
  for (i in seq_along(states_list))
  {
    # i <- 1

    # Extract state
    state_i <- states_list[i]
    # Define other states by contrast
    other_states <- states_list[states_list != state_i]

    ## Binarize states in simmaps
    simmaps_binary_i <- simmaps
    simmaps_binary_i <- lapply(X = simmaps_binary_i, FUN = phytools::mergeMappedStates, old.states = other_states, new.state = "0")
    simmaps_binary_i <- lapply(X = simmaps_binary_i, FUN = phytools::mergeMappedStates, old.states = state_i, new.state = "1")

    class(simmaps_binary_i) <- c("list", "multiSimmap", "multiPhylo")

    # Check that remaining states are all binary
    # unique(unlist(lapply(X = simmaps_binary_i, FUN = function (x) { lapply(X = x$maps, FUN = names) })))
    # simmaps_binary_i[[1]]$maps

    ## Estimate the posterior probabilities of states along all branches (from the set of simulated maps)

    densityMap_state_i <- densityMap_custom(trees = simmaps_binary_i,
                                            tol = 1e-5, verbose = verbose,
                                            col_scale = NULL,
                                            plot = FALSE)

    ## Update color gradient

    # Set color gradient from grey to focal color
    focal_color <- colors_per_levels[i]
    col_fn <- grDevices::colorRampPalette(colors = c("grey90", focal_color))
    col_scale <- col_fn(n = 1001)

    # Update color gradient
    densityMap_state_i <- phytools::setMap(densityMap_state_i, c("grey90", focal_color))

    ## Update state names
    densityMap_state_i$states <- c(paste0("Not ", state_i), state_i)

    # plot(densityMap_state_i)

    ## Store in final object with all density maps
    densityMaps_all_states[[i]] <- densityMap_state_i

    ## Print progress for each state
    if (verbose)
    {
      cat(paste0(Sys.time(), " - Posterior probabilities computed for State = ",state_i," - n\u00B0", i, "/", length(states_list),"\n"))
    }
  }
  names(densityMaps_all_states) <- paste0("Density_map_", states_list)

  ## Export output
  return(densityMaps_all_states)
}



### Helper function to produce densityMap from simmaps ####
# Original function written by Liam Revell, 2012
# Input = simmaps

#' @title Plot posterior density of stochastic mapping on a tree
#'
#' @description Visualize posterior probability density from stochastic mapping using a color gradient on the tree.
#'   Original function written by Liam Revell, 2012 in the [phytools] package: [phytools::densityMap()].
#'
#' @inheritParams phytools::densityMap
#' @param tol Positive numerical. To set the tolerance used to match node ages and time steps (i.e., consider them equal). Default = 1e-5.
#' @param verbose Logical. To display or progress every 100 edges. Default = `TRUE`.
#' @param col_scale Character string vector. To set the color scale manually. Need to provide 1001 colors for the scale.
#'   If `NULL` (the default), the `rainbow()` color scale will be used.
#'
#' @return The function plots a tree with mapped trait probability densities and returns an object of class `densityMap` invisibly.
#'   A `densityMap` is a list with three elements.
#'     * `$tree` List of at least 8 elements. Includes the phylogeny, the trait evolution model data from the simmaps, and the newly mapped trait posterior densities.
#'       * `$maps` List of N elements, one per edge. Each list comprises a named numerical vector that represent changes in posterior probability density of the focal state along segments of equal time.
#'         Named are posterior probabilities scaled from 0 to 1000. Values are length of the segments. Segments are ordered from root to tips.
#'       * `$mapped.edge` Matrix of edge per posterior probability summarizing the overall length of each edge attributed to a specific posterior probabiliy value.
#'       * `$Q` Numerical square matrix summarizing instantaneous transition rates between states as estimated from the evolutionary model.
#'       Rows = initial states. Cols = final states.
#'       * `$logL` Numerical. Log-likelihood of the data as optimized when estimated model parameters.
#'     * `$col` Named character string vector. Color scale used to map posterior probabilities. Names are the posterior probabilities scaled from 0 to 1000. Values are the colors.
#'     * `$states` Character string. The name of the states.
#'
#' @details Wrapped function of [phytools::densityMap()].
#'   Additions to the initial function:
#'   * Can modify manually the tolerance to handle issue with mismatch between node ages and time steps used.
#'   * Can print progress across egdes
#'   * Can provide a manual color scale to replace the default rainbow scale. The color scale must have 1001 colors.
#'
#' @author Maël Doré. Initial function by Liam Revell, 2012 in the [phytools] package.
#'
#' @seealso [phytools::densityMap()]
#'
#' @noRd
#'

densityMap_custom <- function (trees, res = 100, fsize = NULL, ftype = NULL, lwd = 3,
                               tol = 1e-5, verbose = T, col_scale = NULL,
                               check = FALSE, legend = NULL, outline = FALSE, type = "phylogram",
                               direction = "rightwards", plot = TRUE, ...)
{
  # Set graphical parameters
  if (methods::hasArg(mar))
    mar <- list(...)$mar
  else mar <- rep(0.3, 4)
  if (methods::hasArg(offset))
    offset <- list(...)$offset
  else offset <- NULL
  if (methods::hasArg(states))
    states <- list(...)$states
  else states <- NULL
  if (methods::hasArg(hold))
    hold <- list(...)$hold
  else hold <- TRUE
  if (length(lwd) == 1)
    lwd <- rep(lwd, 2)
  else if (length(lwd) > 2)
    lwd <- lwd[1:2]

  # Adjust tolerance
  tol <- tol

  # Check validity of the class of the input object
  if (!inherits(trees, "multiPhylo") && inherits(trees, "phylo"))
  {
    stop("trees not \"multiPhylo\" object; just use plotSimmap.")
  }
  if (!inherits(trees, "multiPhylo"))
  {
    stop("trees should be an object of class \"multiPhylo\".")
  }

  # Extract root age
  h <- sapply(unclass(trees), function(x) max(phytools::nodeHeights(x)))

  # Define time steps
  steps <- 0:res/res * max(h)

  # Rescale trees to ensure they all have the same root age
  trees <- phytools::rescaleSimmap(trees, totalDepth = max(h))

  # Check that phylogeny topology and branch length are equal
  if (check)
  {
    X <- matrix(FALSE, length(trees), length(trees))
    for (i in 1:length(trees)) X[i, ] <- sapply(trees, ape::all.equal.phylo,
                                                current = trees[[i]])
    if (!all(X))
      stop("some of the trees don't match in topology or relative branch lengths")
  }

  # Extract first tree as reference
  tree <- trees[[1]]

  # Remove class
  trees <- unclass(trees)

  # Extract all states from the first tree (dangerous if some states are not present in this tree but in other!)
  if (is.null(states))
  {
    ss <- sort(unique(c(phytools::getStates(tree, "nodes"), phytools::getStates(tree, "tips"))))
  }  else {
    ss <- states
  }

  # If states are not binary, rename the first two states as "0" and "1" (dangerous as if there are more states, will lead to errors)
  if (!all(ss == c("0", "1")))
  {
    c1 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    c2 <- paste(sample(c(letters, LETTERS), 6), collapse = "")
    trees <- lapply(trees, phytools::mergeMappedStates, ss[1], c1)
    trees <- lapply(trees, phytools::mergeMappedStates, ss[2], c2)
    trees <- lapply(trees, phytools::mergeMappedStates, c1, "0")
    trees <- lapply(trees, phytools::mergeMappedStates, c2, "1")
  }

  # Extract all node ages per edge
  H <- phytools::nodeHeights(tree)
  # message("sorry - this might take a while; please be patient")

  # Reinitiate the map of the reference tree with NULL data
  tree$maps <- vector(mode = "list", length = nrow(tree$edge))

  # Loop per edge/item in tree$maps
  for (i in 1:nrow(tree$edge))
  {
    # i <- 5
    # i <- 983

    # YY = Matrix of ages to use as time step along edge i
    # Include start (0) and end (edge length)
    # One raw = one interval
    # Columns = Start and End relative age
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]),
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) - H[i, 1]

    # Initiate vector of final time step values for edge i
    ZZ <- rep(0, nrow(YY))

    # Loop per trees/simmaps
    for (j in 1:length(trees))
    {
      # j <- 1

      # XX = a matrix of the states detected for edge i on simmap j
      # One row = one state
      # Colmuns = start and end relative age
      XX <- matrix(data = 0,
                   nrow = length(trees[[j]]$maps[[i]]),
                   ncol = 2,
                   dimnames = list(names(trees[[j]]$maps[[i]]),
                                   c("start", "end")))
      # Fill the first raw with information on start and end relative age of the first state
      XX[1, 2] <- trees[[j]]$maps[[i]][1]

      # Case with multiple states: fill information for other states
      if (length(trees[[j]]$maps[[i]]) > 1)
      {
        for (k in 2:length(trees[[j]]$maps[[i]]))
        {
          XX[k, 1] <- XX[k - 1, 2]
          XX[k, 2] <- XX[k, 1] + trees[[j]]$maps[[i]][k]
        }
      }

      # Loop per time interval wanted for the density mapping
      for (k in 1:nrow(YY))
      {
        # k <- 1

        # Detect which state start before the k time step
        lower <- which(XX[, 1] <= YY[k, 1])
        lower <- lower[length(lower)] # Take the last one as the last state recorded before the beginning of the time step

        # Detect which state end after the k time step
        upper <- which(XX[, 2] >= (YY[k, 2] - tol))[1]


        AA <- 0
        names(lower) <- names(upper) <- NULL
        if (!all(XX == 0))
        {
          # Case for internal edge (end time > 0)

          # Loop per states on the time interval
          for (l in lower:upper)
          {
            # Compute weighted mean of the time interval
            AA <- AA + (min(XX[l, 2], YY[k, 2]) - max(XX[l, 1], YY[k, 1]))/(YY[k, 2] - YY[k, 1]) * as.numeric(rownames(XX)[l])
          }
        } else {
          # Case for tips (or null branches) (start and end time = 0)
          AA <- as.numeric(rownames(XX)[1]) # Use the tip state
        }
        # Increment the final value of the edge i
        ZZ[k] <- ZZ[k] + AA/length(trees)
      }
    }
    # Record time steps length
    tree$maps[[i]] <- YY[, 2] - YY[, 1]

    # Record time steps continuous value as names of the maps of edge i
    names(tree$maps[[i]]) <- round(ZZ * 1000) # Convert proportion in a scale from 0 to 1000

    # Print progress every 100 edges
    if ((i %% 100 == 0) & verbose)
    {
      cat(paste0(Sys.time(), " - Posterior probability computed for edge n\u00B0", i, "/", nrow(tree$edge),"\n"))
    }
  }

  # Create color scale
  if (is.null(col_scale))
  {
    # cols <- grDevices::rainbow(n =1001, start = 0.7, end = 0) # Initial HSV color scale from phytools::contMap
    cols <- grDevices::rainbow(n = 1001, s = 0.8, start = 0.7, end = 0) # HSV color space with less saturation
    # cols <- colorspace::rainbow_hcl(n = 1001, c = 90, l = 65, start = 240) # With HCL color space
    names(cols) <- 0:1000
  } else {
    cols <- col_scale
  }

  # Recreate map using continuous values
  tree$mapped.edge <- makeMappedEdge(tree$edge, tree$maps)
  tree$mapped.edge <- tree$mapped.edge[, order(as.numeric(colnames(tree$mapped.edge)))]

  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  attr(tree, "map.order") <- "right-to-left"
  x <- list(tree = tree, cols = cols, states = ss)
  class(x) <- "densityMap"

  # Plot
  if (plot)
    plot(x, fsize = fsize, ftype = ftype, lwd = lwd,
         legend = legend, outline = outline, type = type,
         mar = mar, direction = direction, offset = offset,
         hold = hold)

  # Return final simmap
  invisible(x)
}

