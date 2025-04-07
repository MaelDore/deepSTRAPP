##### plot_BAMM_rates() #####

# Explain that it needs the expectedNumberOfShifts to be recorded in the BAMM_object, as in the output of deepSTRAPP::prepare_diversification_data()

# See if it may be more efficient to compute the MAP and MSC in advance in the deepSTRAPP::prepare_diversification_data() and store it in the output?
  # Add in doc of deepSTRAPP::prepare_diversification_data() and update_rates_and_regimes_for_focal_time()
  # Add it to whale/Ponerinae/Ponerine_10My
# See if adding regime shifts. From what?

# Add in the seeAlso of prepare_diversification_data(), after Step 5 in the details section, and in the examples

# Input = BAMM_object with expectedNumberOfShifts, typically from deepSTRAPP::prepare_diversification_data()

# Argument = configuration_type = c("MAP", "MSC")
# Odd-ratio threshold (See if it affects shifts displayed by addBAMMshifts => Does it show only the 'core'-shifts?)
# Argument = method in plot.bammdata: "polar" or "phylogram"
# With regime shifts or not

## Explain the function is a wrapper around BAMMtools::plot.bammdata, BAMMtools::addBAMMshifts and BAMMtools::marginalShiftProbsTree that let
# plot rates and regimes also from BAMM_object updated for a given focal-time in the past (with cut phylogeny)
# Default bg color changed for a more neutral "grey"

# # Import grDevices::gray
#
# ?plot.bammdata
#
# data("Ponerinae_BAMM_object")
#
# # Explore output
# str(Ponerinae_BAMM_object, 1)
#
# # Plot mean net diversification rates on the phylogeny
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object,
#                          labels = FALSE, legend = TRUE)
# BAMMtools::addBAMMshifts(Ponerinae_BAMM_object,
#                          index = Ponerinae_BAMM_object$MSC_index,
#                          cex = 2)
# # Replace by plot_BAMM_rates()
#
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object, lwd = 1.25)
# BAMMtools::addBAMMshifts(MAP, cex=2)
#
# # Need only $eventData, and $tip.label (only used by getShiftNodesFromIndex to extract the root ID based on length($tip.label))
#
# ## Check how the $eventData is generated. Is it simply extracting a random posterior sample or is it somewhat a mean of all samples fiting the configuration?
# # If it is simply a random posterior, I can just record the MAP ID and make it works like MSC.
# # If it is a mean, I need to record the whole MAP bammdata object (so it can be plotted, and used to extract regime shifts)
#
# # Only found in 18 posterior samples across 1000!
# str(MAP, max.level = 1)
#
#
#
# ## NEED to update MAP bamm object in update_rates_and_regimes_for_focal_time too !
# # And to pass down MSC_index and eventually MAP$EventData
#
# ### Make a version that works with MAP too
#
# plot_BAMM_rates(Ponerinae_BAMM_object, show_regime_shifts = TRUE,
#                 configuration_type = "MSC",
#                 regimes_size = 2, bg = "coral")
#
# plot_BAMM_rates(Ponerinae_BAMM_object, show_regime_shifts = TRUE,
#                 configuration_type = "MAP",
#                 regimes_size = 2, bg = "black")

### Replace BAMMtools::plot.bammdata() in all docs!

### Currently, the MSC option is selecting the first sample among all MSC configs.
# It will be better if I used an average approach just like in MAP
# Make my own averaging function based on [BAMMtools::getBestShiftConfiguration()])
# that subsampled all the posterior samples detected with BAMMtools::maximumShiftCredibility() => $bestconfigs[[1]] and not just $sampleindex


plot_BAMM_rates <- function (BAMM_object,
                             rate_type = "net_diversification",
                             configuration_type = "MSC", # MSC, MAP, or index
                             sample_index = 1,
                             show_regime_shifts = TRUE,
                             adjust_size_to_prob = TRUE, # To adjust size of points representing regime shifts to their marginal posterior probabilities
                             method = "phylogram",
                             regimes_fill = "grey", # Replace 'bg' argument in BAMMtools::addBAMMshifts()
                             regimes_size = 1, # Replace 'cex' argument in BAMMtools::addBAMMshifts()
                             regimes_pch = 21, # Replace 'pch' argument in BAMMtools::addBAMMshifts()
                             regimes_border_col = "black", # Replace 'col' argument in BAMMtools::addBAMMshifts()
                             regimes_border_width = 1, # Replace 'lwd' argument in BAMMtools::addBAMMshifts()
                             ...) # To pass down to BAMMtools::plot.bammdata(), BAMMtools::addBAMMshifts(), and par()
{
  ### Check input validity
  {
    ## BAMM_object is a bammdata object

    ## rate_type must be either 'net_diversification', 'speciation', 'extinction'

    ## Need $MSC_index
    # Ask to run prepare_diversification_data, or get it from BAMMtools::maximumShiftCredibility()

    ## Need $MSP_tree if using 'adjust_size_to_prob = TRUE'

    ## configuration_type is either 'MSC' or 'MAP' or 'index'

    ## sample index must be an integer from 1 to length of EventData

    ## method is either 'phylogram' or 'polar'
  }

  ## Convert 'rate_type' into 'spex'
  if (rate_type == "net_diversification") { spex <- "netdiv" }
  if (rate_type == "speciation") { spex <- "s" }
  if (rate_type == "extinction") { spex <- "e" }

  ## Filter list of additional arguments to avoid warnings from par()
  add_args <- list(...)
  args_names_for_plot <- c("tau", "xlim", "ylim", "vtheta", "rbf", "show", "labels", "legend",
                           "spex", "pal", "mask", "mask.color", "colorbreaks", "logcolor",
                           "breaksmethod", "color.interval", "JenksSubset",
                           "par.reset", "direction")
  args_names_for_addBAMMshifts <- c("shiftnodes", "par.reset")

  add_args_for_plot <- add_args[names(add_args) %in% args_names_for_plot]
  add_args_for_par <- add_args[!(names(add_args) %in% c(args_names_for_plot, args_names_for_addBAMMshifts))]

  ## Retrieve named arguments for plot.bammdata
  if ("tau" %in% names(add_args_for_plot)) { tau <- add_args_for_plot$tau } else { tau <- 0.01 }
  if ("xlim" %in% names(add_args_for_plot)) { xlim <- add_args_for_plot$xlim } else { xlim <- NULL }
  if ("ylim" %in% names(add_args_for_plot)) { ylim <- add_args_for_plot$ylim } else { ylim <- NULL }
  if ("vtheta" %in% names(add_args_for_plot)) { vtheta <- add_args_for_plot$vtheta } else { vtheta <- 5 }
  if ("rbf" %in% names(add_args_for_plot)) { rbf <- add_args_for_plot$rbf } else { rbf <- 0.001 }
  if ("show" %in% names(add_args_for_plot)) { show <- add_args_for_plot$show } else { show <- TRUE }
  if ("labels" %in% names(add_args_for_plot)) { labels <- add_args_for_plot$labels } else { labels <- FALSE }
  if ("legend" %in% names(add_args_for_plot)) { legend <- add_args_for_plot$legend } else { legend <- FALSE }
  # if ("spex" %in% names(add_args_for_plot)) { spex <- add_args_for_plot$spex } else { spex <- "s" }
  if ("lwd" %in% names(add_args_for_plot)) { lwd <- add_args_for_plot$lwd } else { lwd <- 1 }
  if ("cex" %in% names(add_args_for_plot)) { cex <- add_args_for_plot$cex } else { cex <- 1 }
  if ("pal" %in% names(add_args_for_plot)) { pal <- add_args_for_plot$pal } else { pal <- "RdYlBu" }
  if ("mask" %in% names(add_args_for_plot)) { mask <- add_args_for_plot$mask } else { mask <- integer(0) }
  if ("mask.color" %in% names(add_args_for_plot)) { mask.color <- add_args_for_plot$mask.color } else { mask.color <- grDevices::gray(0.5) }
  if ("colorbreaks" %in% names(add_args_for_plot)) { colorbreaks <- add_args_for_plot$colorbreaks } else { colorbreaks <- NULL }
  if ("logcolor" %in% names(add_args_for_plot)) { logcolor <- add_args_for_plot$logcolor } else { logcolor <- FALSE }
  if ("breaksmethod" %in% names(add_args_for_plot)) { breaksmethod <- add_args_for_plot$breaksmethod } else { breaksmethod <- "linear" }
  if ("color.interval" %in% names(add_args_for_plot)) { color.interval <- add_args_for_plot$color.interval } else { color.interval <- NULL }
  if ("JenksSubset" %in% names(add_args_for_plot)) { JenksSubset <- add_args_for_plot$JenksSubset } else { JenksSubset <- 20000 }
  if ("par.reset" %in% names(add_args_for_plot)) { par.reset <- add_args_for_plot$par.reset } else { par.reset <- TRUE } # Set to TRUE to avoid affecting next plots
  if ("direction" %in% names(add_args_for_plot)) { direction <- add_args_for_plot$direction } else { direction <- "rightwards" }

  # ## Plot rates
  # output <- BAMMtools::plot.bammdata(x = BAMM_object,
  #                                    method = method,
  #                                    ...)

  ## Plot rates while separating names arguments from additional arguments for par() in the ellipsis (...)
  output <- do.call(what = BAMMtools::plot.bammdata,
                    args = c(list(x = BAMM_object, method = method, tau = tau, xlim = xlim, ylim = ylim, vtheta = vtheta,
                                  rbf = rbf, show = show, labels = labels, legend = legend, spex = spex, lwd = lwd,
                                  cex = cex, pal = pal, mask = mask, mask.color = mask.color, colorbreaks = colorbreaks,
                                  logcolor = logcolor, breaksmethod = breaksmethod, color.interval = color.interval,
                                  JenksSubset = JenksSubset, par.reset = par.reset, direction = direction),
                             add_args_for_par))

  ## Plot regimes if requested
  if (show_regime_shifts)
  {
    ## Filter list of additional arguments to avoid conflicts and warnings from par()
    add_args <- list(...)
    add_args_for_addBAMMshifts <- add_args[names(add_args) %in% args_names_for_addBAMMshifts]
    if ("shiftnodes" %in% names(add_args_for_addBAMMshifts)) { shiftnodes <- add_args_for_addBAMMshifts$shiftnodes } else { shiftnodes <- NULL }
    if ("par.reset" %in% names(add_args_for_addBAMMshifts)) { par.reset <- add_args_for_addBAMMshifts$par.reset } else { par.reset <- TRUE }
    add_args_for_par <- add_args_for_par[!(names(add_args_for_par) %in% c("bg", "cex", "pch", "col", "lwd"))]

    # Provide Marginal Shift Probabilities to adjust size if requested
    if (adjust_size_to_prob)
    {
      msp <- BAMM_object$MSP_tree
    } else {
      msp <- NULL
    }

    # Case for 'index' => Plotting the configuraiton for a given posterior sample
    if (configuration_type == "index")
    {
      # ## Add regime shifts on the plot
      # BAMMtools::addBAMMshifts(ephy = BAMM_object,
      #                          index = sample_index,
      #                          method = method,
      #                          msp = msp,
      #                          bg = regimes_fill,
      #                          cex = regimes_size,
      #                          pch = regimes_pch,
      #                          col = regimes_border_col,
      #                          lwd = regimes_border_width
      #                          ...)

      ## Add regime shifts on the plot while separating names arguments from additional arguments for points() and par() in the ellipsis (...)
      do.call(what = addBAMMshifts_custom,
              args = c(list(ephy = BAMM_object,
                            index = sample_index,
                            method = method,
                            msp = msp,
                            regimes_fill = regimes_fill,
                            regimes_size = regimes_size,
                            regimes_pch = regimes_pch,
                            regimes_border_col = regimes_border_col,
                            regimes_border_width = regimes_border_width,
                            shiftnodes = shiftnodes,
                            par.reset = par.reset),
                       add_args_for_par))

    }

    # Case for Maximum Shift Credibility (MSC) configuration
    if (configuration_type == "MSC")
    {
      # ## Add regime shifts on the plot
      # BAMMtools::addBAMMshifts(ephy = BAMM_object,
      #                          index = BAMM_object$MSC_index,
      #                          method = method,
      #                          msp = msp,
      #                          bg = regimes_fill,
      #                          cex = regimes_size,
      #                          pch = regimes_pch,
      #                          col = regimes_border_col,
      #                          lwd = regimes_border_width
      #                          ...)

      ## Add regime shifts on the plot while separating names arguments from additional arguments for points() and par() in the ellipsis (...)
      do.call(what = addBAMMshifts_custom,
              args = c(list(ephy = BAMM_object,
                            index = BAMM_object$MSC_index,
                            method = method,
                            msp = msp,
                            regimes_fill = regimes_fill,
                            regimes_size = regimes_size,
                            regimes_pch = regimes_pch,
                            regimes_border_col = regimes_border_col,
                            regimes_border_width = regimes_border_width,
                            shiftnodes = shiftnodes,
                            par.reset = par.reset),
                       add_args_for_par))

    }

    # Case for Maximum A Posteriori probability (MAP) configuration
    if (configuration_type == "MAP")
    {
      ## Use the BAMM_object$MAP_BAMM_object to get locations of shifts
      do.call(what = addBAMMshifts_custom,
              args = c(list(ephy = BAMM_object$MAP_BAMM_object,
                            index = 1,
                            method = method,
                            msp = msp,
                            regimes_fill = regimes_fill,
                            regimes_size = regimes_size,
                            regimes_pch = regimes_pch,
                            regimes_border_col = regimes_border_col,
                            regimes_border_width = regimes_border_width,
                            shiftnodes = shiftnodes,
                            par.reset = par.reset),
                       add_args_for_par))
    }
  }

  # Return output of BAMMtools::plot.bammdata()
  return(invisible(output))
}



### Modified version of BAMMtools::addBAMMshifts

## Handle adjustment of regime shift point size controlled by both cex and msp
## Fix issue with conflicting parameter names between the main function and par()

# Source: BAMMtools::addBAMMshifts()
# Author: Mike Grundler

#' @importFrom graphics par points
#' @importFrom ape as.phylo branching.times .PlotPhyloEnv
#' @importFrom BAMMtools getShiftNodesFromIndex

addBAMMshifts_custom <- function (ephy, index = 1, method = "phylogram", regimes_size = 1, regimes_pch = 21,
                                  regimes_border_col = 1, regimes_border_width = 1, regimes_fill = 2,
                                  msp = NULL, shiftnodes = NULL, par.reset = TRUE, ...)
{
  if (!inherits(ephy, "bammdata"))
    stop("Object ephy must be of class bammdata")
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  if (par.reset) {
    op <- graphics::par(no.readonly = TRUE)
    graphics::par(lastPP$pp)
  }
  if (length(ephy$eventData) == 1) {
    index <- 1
  }
  if (is.null(shiftnodes))
    shiftnodes <- BAMMtools::getShiftNodesFromIndex(ephy, index)
  isShift <- ephy$eventData[[index]]$node %in% shiftnodes
  times <- ephy$eventData[[index]]$time[isShift]
  if (!is.null(msp)) {
    # Key change: regimes_size (cex) is update, not replaced by the msp size
    regimes_size <- 0.5 + regimes_size * 2.5 * msp$edge.length[msp$edge[, 2] %in%
                                        shiftnodes]
  }
  if (method == "phylogram") {
    XX <- times
    YY <- lastPP$yy[shiftnodes]
  }
  else if (method == "polar") {
    rb <- lastPP$rb
    XX <- (rb + times/max(ape::branching.times(ape::as.phylo(ephy)))) *
      cos(lastPP$theta[shiftnodes])
    YY <- (rb + times/max(ape::branching.times(ape::as.phylo(ephy)))) *
      sin(lastPP$theta[shiftnodes])
  }
  graphics::points(XX, YY, pch = regimes_pch, cex = regimes_size, col = regimes_border_col,
                   bg = regimes_fill, lwd = regimes_border_width, ...)
  if (par.reset) {
    graphics::par(op)
  }
}
