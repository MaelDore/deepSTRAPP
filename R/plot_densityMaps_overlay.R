

#' @title Plot posterior probabilities of states/ranges on phylogeny from densityMaps
#'
#' @description Plot on a time-calibrated phylogeny the evolution of a categorical trait/biogeographic ranges
#'   summarized from `densityMaps` typically generated with [deepSTRAPP::prepare_trait_data()].
#'   Each branch is colored according to the posterior probability of being in a given state/range.
#'   Color for each state/range are overlaid using transparency to produce a single plot for all states/ranges.
#'
#' @param densityMaps List of objects of class `"densityMap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state/range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state/range. If no color is provided for multi-area ranges, they will be interpolated.
#' @param colors_per_levels Named character string. To set the colors to use to map each state/range posterior probabilities. Names = states/ranges; values = colors.
#'   If `NULL` (default), the color scale provided `densityMaps` will be used.
#' @param add_ACE_pies Logical. Whether to add pies of posterior probabilities of states/ranges at internal nodes on the mapped phylogeny. Default = `TRUE`.
#' @param cex_pies Numerical. To adjust the size of the ACE pies. Default = `0.5`.
#' @param ace Numerical matrix. To provide the posterior probabilities of ancestral states/ranges (characters) estimates (ACE) at internal nodes
#'   used to plot the ACE pies. Rows are internal nodes. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'   Typically generated with [deepSTRAPP::prepare_trait_data()] in the `$ace` slot.
#'   If `NULL` (default), the ACE are extracted from the `densityMaps` with a possible slight discrepancy with the actual tip states
#'   and estimated posterior probabilities of ancestral states.
#' @param ... Additional arguments to pass down to [phytools::plotSimmap()] to control plotting.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#'
#' @export
#' @importFrom graphics par
#' @importFrom phytools add.simmap.legend
#' @importFrom ape nodelabels
#' @importFrom grDevices pdf dev.off
#'
#' @return If `display_plot = TRUE`, the function plots a time-calibrated phylogeny displaying the evolution of a categorical trait/biogeographic ranges.
#' If `PDF_file_path` is provided, the function exports the plot into a PDF file.
#'
#' @author Maël Doré
#' @author Original functions by Liam Revell in R package `{phytools}`. Contact: \email{liam.revell@umb.edu}
#'
#' @seealso [phytools::plot.densityMap()] [phytools::plotSimmap()]
#'
#' @examples
#'
#' # Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Transform feeding mode data into a 3-level factor
#' eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
#' eel_data <- as.character(eel_data)
#' eel_data[c(1, 5, 6, 7, 10, 11, 15, 16, 17, 24, 25, 28, 30, 51, 52, 53, 55, 58, 60)] <- "kiss"
#' eel_data <- stats::setNames(eel_data, rownames(eel.data))
#' table(eel_data)
#'
#' # Manually define a Q_matrix for rate classes of state transition to use in the 'matrix' model
#' # Does not allow transitions from state 1 ("bite") to state 2 ("kiss") or state 3 ("suction")
#' # Does not allow transitions from state 3 ("suction") to state 1 ("bite")
#' # Set symmetrical rates between state 2 ("kiss") and state 3 ("suction")
#' Q_matrix = rbind(c(NA, 0, 0), c(1, NA, 2), c(0, 2, NA))
#'
#' # Set colors per states
#' colors_per_levels <- c("limegreen", "orange", "dodgerblue")
#' names(colors_per_levels) <- c("bite", "kiss", "suction")
#'
#' \donttest{ # (May take several minutes to run)
#' # Run evolutionary models to prepare trait data
#' eel_cat_3lvl_data <- prepare_trait_data(tip_data = eel_data, phylo = eel.tree,
#'     trait_data_type = "categorical",
#'     colors_per_levels = colors_per_levels,
#'     evolutionary_models = c("ER", "SYM", "ARD", "meristic", "matrix"),
#'     Q_matrix = Q_matrix,
#'     nb_simulations = 1000,
#'     plot_map = TRUE,
#'     plot_overlay = TRUE,
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE) }
#'
#' # Load directly output
#' data(eel_cat_3lvl_data, package = "deepSTRAPP")
#'
#' # Plot densityMaps one by one
#' plot(eel_cat_3lvl_data$densityMaps[[1]]) # densityMap for state n°1 ("bite")
#' plot(eel_cat_3lvl_data$densityMaps[[2]]) # densityMap for state n°1 ("kiss")
#' plot(eel_cat_3lvl_data$densityMaps[[3]]) # densityMap for state n°1 ("suction")
#'
#' # Plot overlay of all densityMaps
#' plot_densityMaps_overlay(densityMaps = eel_cat_3lvl_data$densityMaps)
#'


plot_densityMaps_overlay <- function (
    densityMaps,
    colors_per_levels = NULL,
    add_ACE_pies = TRUE,
    cex_pies = 0.5,
    ace = NULL,
    ..., # To allow to pass down arguments in the plotSimmap() function
    display_plot = TRUE,
    PDF_file_path = NULL)
{
  # Get list of states
  states_list <- unname(unlist(lapply(densityMaps, FUN = function (x) { x$states[2] })))

  # Get number of tips
  nb_tips <- length(densityMaps[[1]]$tree$tip.label)
  # Get number of nodes
  nb_nodes <- nb_tips + densityMaps[[1]]$tree$Nnode

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

    ## ace
    if (!is.null(ace))
    {
      # Check that ace have proper rows (internal nodes)
      if (!(identical(row.names(ace), as.character((nb_tips+1):nb_nodes))))
      {
        stop(paste0("Row.names in 'ace' do not match internal node IDs."))
      }
      # Check that ace have proper columns (states)
      if (!(identical(colnames(ace), states_list)))
      {
        missing_states <- states_list[!(states_list %in% colnames(ace))]
        stop(paste0("Not all states are found in 'ace'.\n",
                    "Missing: ", paste(missing_states, collapse = ", "), "."))
      }
    }

    ## display_plot OR PDF_file_path
    # Check that at least one type of output is requested
    if (!display_plot & is.null(PDF_file_path))
    {
      stop(paste0("You must request at least one option between displaying the plot (`display_plot` = TRUE), or producing a PDF (fill the `PDF_file_path` argument)."))
    }

    ## PDF_file_path
    # If provided, PDF_file_path must end with ".pdf"
    if (!is.null(PDF_file_path))
    {
      if (length(grep(pattern = "\\.pdf$", x = PDF_file_path)) != 1)
      {
        stop("'PDF_file_path' must end with '.pdf'")
      }
    }
  }

  ## Filter list of additional arguments to ensure default values are used if not provided
  add_args <- list(...)
  # Extract additional args for phytools::plotSimmap()
  args_names_for_plotSimmap <- c("fsize", "ftype", "lwd", "pts", "node.numbers", "mar", "offset", "direction",
                                  "type", "setEnv", "part", "xlim", "ylim", "nodes", "tips", "maxY", "hold",
                                  "split.vertical", "lend", "asp", "outline", "underscore", "arc_height")
  add_args_for_plotSimmap <- add_args[names(add_args) %in% args_names_for_plotSimmap]
  # Set default value if not provided
  if ("fsize" %in% names(add_args_for_plotSimmap)) { fsize <- add_args_for_plotSimmap$fsize } else { fsize <- 0.7 }
  if ("ftype" %in% names(add_args_for_plotSimmap)) { ftype <- add_args_for_plotSimmap$ftype } else { ftype <- "reg" }
  if ("lwd" %in% names(add_args_for_plotSimmap)) { lwd <- add_args_for_plotSimmap$lwd } else { lwd <- 2 }
  if ("mar" %in% names(add_args_for_plotSimmap)) { mar <- add_args_for_plotSimmap$mar } else { mar <- graphics::par()$mar }
  if ("tips" %in% names(add_args_for_plotSimmap)) { tips <- add_args_for_plotSimmap$tips } else { tips <- stats::setNames(object = 1:nb_tips, nm = densityMaps[[1]]$tree$tip.label) }
  add_args_for_plotSimmap <- add_args_for_plotSimmap[!(names(add_args_for_plotSimmap) %in% c("fsize", "ftype", "lwd", "mar", "tips"))]

  ## Retrieve colors_per_levels if not provided
  if (is.null(colors_per_levels))
  {
    colors_per_levels <- unname(unlist(lapply(densityMaps, FUN = function (x) { x$cols[length(x$cols)] })))
    names(colors_per_levels) <- states_list
  }

  ### Display plots
  if (display_plot)
  {
    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    ## Loop per state
    for (i in seq_along(states_list))
    {
      # i <- 1

      # Extract densityMap
      densityMap_state_i <- densityMaps[[i]]

      # Set color gradient from transparent to focal color
      focal_color <- colors_per_levels[i]
      focal_color_rgb <- grDevices::col2rgb(focal_color, alpha = TRUE)
      focal_color_hexa0 <- grDevices::rgb(red = focal_color_rgb[1,1], green = focal_color_rgb[2,1], blue = focal_color_rgb[3,1], alpha = 0, maxColorValue = 255)
      focal_color_hexa1 <- grDevices::rgb(red = focal_color_rgb[1,1], green = focal_color_rgb[2,1], blue = focal_color_rgb[3,1], alpha = 255, maxColorValue = 255)
      col_fn <- grDevices::colorRampPalette(colors = c(focal_color_hexa0, focal_color_hexa1), alpha = TRUE)
      col_scale <- col_fn(n = 1001)

      # Update color gradient
      densityMap_state_i <- phytools::setMap(densityMap_state_i, c(focal_color_hexa0, focal_color_hexa1), alpha = TRUE)

      # plot(densityMap_state_i, legend = FALSE)

      if (i == 1)
      {
        add_plot <- FALSE
      } else {
        add_plot <- TRUE
      }

      # # Plot each densityMap as a Simmap with transparent colors
      # plot(x = densityMap_state_i$tree, colors = densityMap_state_i$cols,
      #      fsize = fsize, ftype = ftype, lwd = lwd, add = add_plot,
      #      mar = graphics::par()$mar, tips = tips,
      #      plot = TRUE, ...)

      do.call(what = plot, # phytools::plotSimmap
              args = c(list(x = densityMap_state_i$tree, colors = densityMap_state_i$cols,
                            fsize = fsize, ftype = ftype, lwd = lwd, add = add_plot,
                            mar = graphics::par()$mar, tips = tips,
                            plot = TRUE),
                       add_args_for_plotSimmap))

    }

    # Add node pies of ACE if requested
    if (add_ACE_pies)
    {
      # Compute ACE if not provided
      if (is.null(ace))
      {
        # Get root ID
        # root_ID <- nb_tips + 1
        root_ID <- (1:nb_nodes)[!(1:nb_nodes %in% densityMaps[[1]]$tree$edge[,2])]

        # Initiate matrix of state posterior probability per nodes
        PP_per_nodes <- matrix(data = NA, nrow = nb_nodes, ncol = length(densityMaps))

        ## Loop per states
        for (i in 1:length(densityMaps))
        {
          # i <- 1

          # Extract densityMap of state i
          maps_i <- densityMaps[[i]]$tree$maps

          # Get last posterior frequency per edges
          last_freq_edges_i <- unlist(lapply(X = maps_i, FUN = function (x) { names(x)[length(x)] }))

          # Match edges with tipward nodes
          freqs_per_nodes <- data.frame(node = 1:nb_nodes, state = NA)
          freqs_per_nodes$freq <- last_freq_edges_i[match(x = freqs_per_nodes$node, table = densityMaps[[1]]$tree$edge[,2])]

          # Extract state for the root
          root_edges_ID <- which(densityMaps[[1]]$tree$edge[,1] == root_ID)
          root_state <- names(maps_i[root_edges_ID][[1]])[1] # Take the initial state of the first descending edge (should be equal among both descending edges)
          freqs_per_nodes$freq[root_ID] <- root_state

          # Store freqs in summary matrix
          PP_per_nodes[, i] <- as.numeric(freqs_per_nodes$freq) / 1000
        }

        ## Format in posterior probabilities for each state
        row.names(PP_per_nodes) <- 1:nb_nodes
        colnames(PP_per_nodes) <- states_list

        # Rearrange with internal nodes followed by tips
        ace_matrix <- PP_per_nodes
        ace_matrix[1:(nb_nodes-nb_tips), ] <- PP_per_nodes[(nb_tips+1):nb_nodes, ]
        ace_matrix[(nb_nodes-nb_tips+1):nb_nodes, ] <- PP_per_nodes[1:nb_tips, ]
        row.names(ace_matrix) <- c((nb_nodes-nb_tips+2):nb_nodes, densityMaps[[1]]$tree$tip.label)

        ## Display ACE posterior probabilities
        # print(PP_per_nodes)


      } else {
        # Add tip in ACE matrix
        tip_matrix <- matrix(data = NA, nrow = nb_tips, ncol = ncol(ace))
        ace_matrix <- rbind(ace, tip_matrix)
        row.names(ace_matrix) <- c((nb_nodes-nb_tips+2):nb_nodes, densityMaps[[1]]$tree$tip.label)
      }

      # Add ACE pies
      ape::nodelabels(pie = ace_matrix, piecol = colors_per_levels, cex = cex_pies)
    }

    # Add legend
    phytools::add.simmap.legend(colors = colors_per_levels, x = par()$usr[1] + 0.05 * (par()$usr[2] - par()$usr[1]), y = par()$usr[3] - 0.01 * (par()$usr[4] - par()$usr[3]),
                                vertical = FALSE,
                                prompt = FALSE)

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Save PDF
  if (!is.null(PDF_file_path))
  {
    # Adjust width/height according to the nb of tips
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10

    ## Open PDF
    grDevices::pdf(file = file.path(PDF_file_path),
                   width = width, height = height)

    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    ## Loop per state
    for (i in seq_along(states_list))
    {
      # i <- 1

      # Extract densityMap
      densityMap_state_i <- densityMaps[[i]]

      # Set color gradient from transparent to focal color
      focal_color <- colors_per_levels[i]
      focal_color_rgb <- grDevices::col2rgb(focal_color, alpha = TRUE)
      focal_color_hexa0 <- grDevices::rgb(red = focal_color_rgb[1,1], green = focal_color_rgb[2,1], blue = focal_color_rgb[3,1], alpha = 0, maxColorValue = 255)
      focal_color_hexa1 <- grDevices::rgb(red = focal_color_rgb[1,1], green = focal_color_rgb[2,1], blue = focal_color_rgb[3,1], alpha = 255, maxColorValue = 255)
      col_fn <- grDevices::colorRampPalette(colors = c(focal_color_hexa0, focal_color_hexa1), alpha = TRUE)
      col_scale <- col_fn(n = 1001)

      # Update color gradient
      densityMap_state_i <- phytools::setMap(densityMap_state_i, c(focal_color_hexa0, focal_color_hexa1), alpha = TRUE)

      # plot(densityMap_state_i, legend = FALSE)

      if (i == 1)
      {
        add_plot <- FALSE
      } else {
        add_plot <- TRUE
      }

      # # Plot each densityMap as a Simmap with transparent colors
      # plot(x = densityMap_state_i$tree, colors = densityMap_state_i$cols,
      #      fsize = fsize, ftype = ftype, lwd = lwd, add = add_plot,
      #      mar = graphics::par()$mar, tips = tips,
      #      plot = TRUE, ...)

      do.call(what = plot, # phytools::plotSimmap
              args = c(list(x = densityMap_state_i$tree, colors = densityMap_state_i$cols,
                            fsize = fsize, ftype = ftype, lwd = lwd, add = add_plot,
                            mar = graphics::par()$mar, tips = tips,
                            plot = TRUE),
                       add_args_for_plotSimmap))

    }

    # Add node pies of ACE if requested
    if (add_ACE_pies)
    {
      # Compute ACE if not provided
      if (is.null(ace))
      {
        # Get root ID
        # root_ID <- nb_tips + 1
        root_ID <- (1:nb_nodes)[!(1:nb_nodes %in% densityMaps[[1]]$tree$edge[,2])]

        # Initiate matrix of state posterior probability per nodes
        PP_per_nodes <- matrix(data = NA, nrow = nb_nodes, ncol = length(densityMaps))

        ## Loop per states
        for (i in 1:length(densityMaps))
        {
          # i <- 1

          # Extract densityMap of state i
          maps_i <- densityMaps[[i]]$tree$maps

          # Get last posterior frequency per edges
          last_freq_edges_i <- unlist(lapply(X = maps_i, FUN = function (x) { names(x)[length(x)] }))

          # Match edges with tipward nodes
          freqs_per_nodes <- data.frame(node = 1:nb_nodes, state = NA)
          freqs_per_nodes$freq <- last_freq_edges_i[match(x = freqs_per_nodes$node, table = densityMaps[[1]]$tree$edge[,2])]

          # Extract state for the root
          root_edges_ID <- which(densityMaps[[1]]$tree$edge[,1] == root_ID)
          root_state <- names(maps_i[root_edges_ID][[1]])[1] # Take the initial state of the first descending edge (should be equal among both descending edges)
          freqs_per_nodes$freq[root_ID] <- root_state

          # Store freqs in summary matrix
          PP_per_nodes[, i] <- as.numeric(freqs_per_nodes$freq) / 1000
        }

        ## Format in posterior probabilities for each state
        row.names(PP_per_nodes) <- 1:nb_nodes
        colnames(PP_per_nodes) <- states_list

        # Rearrange with internal nodes followed by tips
        ace_matrix <- PP_per_nodes
        ace_matrix[1:(nb_nodes-nb_tips), ] <- PP_per_nodes[(nb_tips+1):nb_nodes, ]
        ace_matrix[(nb_nodes-nb_tips+1):nb_nodes, ] <- PP_per_nodes[1:nb_tips, ]
        row.names(ace_matrix) <- c((nb_nodes-nb_tips+2):nb_nodes, densityMaps[[1]]$tree$tip.label)

        ## Display ACE posterior probabilities
        # print(PP_per_nodes)


      } else {
        # Add tip in ACE matrix
        tip_matrix <- matrix(data = NA, nrow = nb_tips, ncol = ncol(ace))
        ace_matrix <- rbind(ace, tip_matrix)
        row.names(ace_matrix) <- c((nb_nodes-nb_tips+2):nb_nodes, densityMaps[[1]]$tree$tip.label)
      }

      # Add ACE pies
      ape::nodelabels(pie = ace_matrix, piecol = colors_per_levels, cex = cex_pies)
    }

    # Add legend
    phytools::add.simmap.legend(colors = colors_per_levels, x = par()$usr[1] + 0.05 * (par()$usr[2] - par()$usr[1]), y = par()$usr[3] - 0.01 * (par()$usr[4] - par()$usr[3]),
                                vertical = FALSE,
                                prompt = FALSE)

    # Reset $xpd to initial values
    par(xpd = xpd_init)

    ## Close PDF
    invisible(grDevices::dev.off())
  }
}


