
#' @title Plot continuous trait evolution on the tree
#'
#' @description Plot on a time-calibrated phylogeny the evolution of a continuous trait as
#'   summarized in a `contMap` object typically generated with [deepSTRAPP::prepare_trait_data()].
#'
#'   This function is a wrapper of original functions from the R package [phytools]:
#'   * Step 1: Use [phytools::setMap()] to update the color scale if requested.
#'   * Step 2: Use [phytools::plot.contMap()] to plot the mapped phylogeny.
#'
#' @param contMap List of class `"contMap"`, typically generated with [deepSTRAPP::prepare_trait_data()],
#'   that contains a phylogenetic tree and associated posterior probability of being in a given state/range along branches.
#'   Each object (i.e., `densityMap`) corresponds to a state/range. If no color is provided for multi-area ranges, they will be interpolated.
#' @param color_scale Vector of character string. List of colors to use to build the color scale with [grDevices::colorRampPalette()]
#'   showing the evolution of a continuous trait. From lowest values to highest values.
#' @param ... Additional arguments to pass down to [phytools::plot.contMap()] to control plotting.
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#'
#' @export
#' @importFrom phytools setMap
#' @importFrom grDevices pdf dev.off
#'
#' @details This function is a wrapper of original functions [phytools::setMap()] and [phytools::plot.contMap()]. Additions are listed below:
#'   * The color scale can be controlled directly with the argument `color_scale`.
#'   * The plot can be exported in PDF using `PDF_file_path` to define the output file.
#'
#' @return If `display_plot = TRUE`, the function plots the time-calibrated phylogeny displaying the evolution of a continuous trait.
#' If `PDF_file_path` is provided, the function exports the plot into a PDF file.
#'
#' An object of class `"contMap"` with an (optionally) updated color scale (`$cols`) is returned invisibly.
#'
#' @author Maël Doré
#' @author Original functions by Liam Revell in R package [phytools]. Contact: \email{liam.revell@umb.edu}
#'
#' @seealso [phytools::plot.contMap()] [deepSTRAPP::plot_densityMaps_overlay()]
#'
#' @examples
#' # Load phylogeny
#' data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#' # Load trait df
#' data(Ponerinae_tree, package = "deepSTRAPP")
#'
#' ## Prepare trait data
#'
#' # Extract continuous trait data as a named vector
#' Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#'                                     nm = Ponerinae_trait_tip_data$Taxa)
#'
#' # Get Ancestral Character Estimates based on a Brownian Motion model
#' # To obtain values at internal nodes
#' Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree, x = Ponerinae_cont_tip_data)
#'
#' # Infer Ancestral Character Estimates based on a Brownian Motion model
#' # and run a Stochastic Mapping to interpolate values along branches and obtain a "contMap" object
#' Ponerinae_contMap <- phytools::contMap(Ponerinae_tree, x = Ponerinae_cont_tip_data,
#'                                        res = 100, # Number of time steps
#'                                        plot = FALSE)
#'
#' # Plot contMap with the phytools method
#' plot(x = Ponerinae_contMap)
#'
#' # Plot contMap with an updated color scale
#' plot_contMap(contMap = Ponerinae_contMap,
#'              color_scale = c("darkgreen", "limegreen", "orange", "red"))
#'              # PDF_file_path = "Ponerinae_contMap.pdf")
#'

plot_contMap <- function (contMap,
                          color_scale = NULL,
                          ..., # To pass down to BAMMtools::plot.bammdata(), BAMMtools::addBAMMshifts(), and par()
                          display_plot = TRUE,
                          PDF_file_path = NULL)
{
  ### Check input validity
  {
    ## color_scale
    # Check whether all colors are valid
    if (!is.null(color_scale))
    {
      if (!all(is_color(color_scale)))
      {
        invalid_colors <- color_scale[!is_color(color_scale)]
        stop(paste0("Some color names in 'color_scale' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
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

  ## Update color scale if requested
  if (!is.null(color_scale))
  {
    # Update color palette in contMap
    updated_contMap <- phytools::setMap(x = contMap, colors = color_scale)
    # print(contMap$cols)
  } else {
    updated_contMap <- contMap
  }

  ## Plot the contMap
  phytools::plot.contMap(x = contMap,
                         plot = display_plot,
                         ...)

  ## Export PDF if requested
  if (!is.null(PDF_file_path))
  {
    # Adjust width and height according to phylo
    nb_tips <- length(contMap$tree$tip.label)
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10

    # Open PDF
    grDevices::pdf(file = file.path(PDF_file_path),
                   width = width, height = height)

    ## Plot the contMap
    phytools::plot.contMap(x = contMap,
                           plot = TRUE,
                           ...)
    # Close PDF
    grDevices::dev.off()
  }

  # Return output of phytools::plot.contMap()
  return(invisible(updated_contMap))
}
