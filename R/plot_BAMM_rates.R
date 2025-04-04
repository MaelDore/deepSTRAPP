##### plot_BAMM_rates() #####

# Explain that it needs the expectedNumberOfShifts to be recorded in the BAMM_object, as in the output of deepSTRAPP::prepare_diversification_data()

# See if it may be more efficient to compute the MPP and MSC in advance in the deepSTRAPP::prepare_diversification_data() and store it in the output?
  # Add in doc of deepSTRAPP::prepare_diversification_data() and update_rates_and_regimes_for_focal_time()
  # Add it to whale/Ponerinae/Ponerine_10My
# See if adding regime shifts. From what?

# Add in the seeAlso of prepare_diversification_data(), after Step 5 in the details section, and in the examples

# Input = BAMM_object with expectedNumberOfShifts, typically from deepSTRAPP::prepare_diversification_data()

# Argument = configuration_type = c("MPP", "MSC")
# Odd-ratio threshold (See if it affects shifts displayed by addBAMMshifts => Does it show only the 'core'-shifts?)
# Argument = method in plot.bammdata: "polar" or "phylogram"
# With regime shifts or not

## Explain the function is a wrapper around BAMMtools::plot.bammdata, BAMMtools::addBAMMshifts and BAMMtools::marginalShiftProbsTree that let
# plot rates and regimes also from BAMM_object updated for a given focal-time in the past (with cut phylogeny)
# Default bg color changed for a more neutral "grey"

# Import grDevices::gray

# ?plot.bammdata
#
# data("Ponerinae_BAMM_object")
#
# # Explore output
# str(Ponerinae_BAMM_object, 1)
#
# # Plot mean net diversification rates on the phylogeny
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object,
#                          labels = TRUE, legend = TRUE)
# BAMMtools::addBAMMshifts(Ponerinae_BAMM_object,
#                          index = Ponerinae_BAMM_object$MSC_index,
#                          cex = 2)
# # Replace by plot_BAMM_rates()
#
# ## Try to cut Ponerinae deep enough so some events should disappear
# MSC_detection <- BAMMtools::maximumShiftCredibility(Ponerinae_BAMM_object)
# MSC_detection$sampleindex
#
# branch_marg_posterior_probs <- BAMMtools::marginalShiftProbsTree(Ponerinae_BAMM_object)
#
# ## NEED to update MSP tree in update_rates_and_regimes_for_focal_time too !
# # And to pass down MSC_index and eventually MAP$EventData
#
# ### Make a version that works with MPP too
#
# plot_BAMM_rates(Ponerinae_BAMM_object,
#                 regimes_size = 2, bg = "black")
#
# MSP_tree <- BAMMtools::marginalShiftProbsTree(Ponerinae_BAMM_object)
# Ponerinae_BAMM_object$MSP_tree <- MSP_tree


plot_BAMM_rates <- function (BAMM_object,
                             configuration_type = "MSC", # MSC or MPP
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

    ## Need $MSC_index
    # Ask to run prepare_diversification_data, or get it from BAMMtools::maximumShiftCredibility()

    ## Need $MSP_tree if using 'adjust_size_to_prob = TRUE'

    ## configuration_type is either 'MSC' or 'MPP'

    ## method is either 'phylogram' or 'polar'
  }

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
  if ("spex" %in% names(add_args_for_plot)) { spex <- add_args_for_plot$spex } else { spex <- "s" }
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

    # Provide Marginal Shift Probabilities to adjust size if requested
    if (adjust_size_to_prob)
    {
      msp <- BAMM_object$MSP_tree
    } else {
      msp <- NULL
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
                            bg = regimes_fill,
                            cex = regimes_size,
                            pch = regimes_pch,
                            col = regimes_border_col,
                            lwd = regimes_border_width,
                            shiftnodes = shiftnodes,
                            par.reset = par.reset),
                       add_args_for_par))

    } else {
      # Case for Maximum Posterior Probability (MPP) configuration
    }
  }

  # Return output of BAMMtools::plot.bammdata()
  return(invisible(output))
}

# ##### 7/ Identify the most likely configuration ####
#
# ### 7.1/ Maximum A posteriori Probability (MPP) configuration ###
#
# ## Get the most frequent configuration = maximum a posteriori probability (MPP) shift configuration = the single configuration of shift location showing up the most in the posterior sample
#
# # Idea: find the posterior configuration that is the most frequent in the posterior samples
#
# MPP <- BAMMtools::getBestShiftConfiguration(Ponerinae_BAMM_object,
#    expectedNumberOfShifts = Ponerinae_BAMM_object$expectedNumberOfShifts,
#    threshold = 0) # Odd-ratio threshold used to select core-shifts based on odd-ratio of prior/posterior probabilities
#
#
# ## Why does the threshold influence the position along branchs of the shifts!? (And not just the selection of 'core-shifts' ?)
# # Use zero by default as you want to show all shifts selected for the MPP.
#
# x <-  BAMMtools::credibleShiftSet(Ponerinae_BAMM_object,
#                                   expectedNumberOfShifts = Ponerinae_BAMM_object$expectedNumberOfShifts,
#                                   threshold = 5, set.limit = 0.95)
# x$indices[[1]]
#
# str(MPP, 1)
#
# plot(Ponerinae_BAMM_object)
# plot(MPP)
#
# BAMMtools::plot.bammdata(MPP, lwd = 1.25)
# BAMMtools::addBAMMshifts(MPP, cex=2)
#
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object, lwd = 1.25)
# BAMMtools::addBAMMshifts(MPP, cex=2)
#
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object, lwd = 1.25)
# BAMMtools::addBAMMshifts(Ponerinae_BAMM_object, index = x$indices[[1]][1], cex=2)
#
# for (i in seq_along(x$indices[[1]]))
# {
#   BAMMtools::plot.bammdata(Ponerinae_BAMM_object, lwd = 1.25, main = paste0("MPP n°", i))
#   BAMMtools::addBAMMshifts(Ponerinae_BAMM_object, index = x$indices[[1]][i], cex=2)
#
# }
#
# # Only useful thing is the eventData table. Which can be filtered to remove shifts beyond focal-time!
# # or I could simply use MSC only...
#
# cset <- BAMMtools::credibleShiftSet(whale_BAMM_object,
#                          expectedNumberOfShifts = whale_BAMM_object$expectedNumberOfShifts, # Initial parameter used in the config file to set the exponential prior of the lambda parameter for number of shifts
#                          threshold = 5) # Significance threshold of odd-ratios between prior and posterior probabilities of rate shift on branches
#
# # Only found in 18 posterior samples across 1000!
# str(cset$indices, max.level = 1)
#
#
# ### Think about how to handle that on cut phylo!
#  # Extract the MPP and include it in the main object (during prepare_diversification_data), to be cut too?
#  # More efficent:
#     # Identify its ID from the BAMMtools::getBestShiftConfiguration
#     # And then plot! Adding shifts with BAMMtools::addBAMMshifts using the proper index !!!!
#
# # Do the same for MSC !
#
# # Identify the MSC configs
# MSC_detection <- BAMMtools::maximumShiftCredibility(Ponerinae_BAMM_object)
# MSC_detection$sampleindex
#
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object, lwd = 1.25)
# BAMMtools::addBAMMshifts(Ponerinae_BAMM_object, index = MSC_detection$sampleindex, cex = 2)
#
# BAMMtools::plot.bammdata(Ponerinae_BAMM_object_10My, lwd = 1.25)
# BAMMtools::addBAMMshifts(Ponerinae_BAMM_object_10My, index = MSC_detection$sampleindex, cex = 2)
#
# # Too make it work with MPP => Need to save the $eventData of the MPP output in BAMM_object.
# # When running update regime_rates, also update the $eventData of the MPP
# # When plotting, create a temporary copy and replace the $eventData of the all posterior samples with the $eventData of the MPP,
# # so the BAMMtools::addBAMMshifts function will use it to plot regimes.
#
# # It means, the MPP and MSC ID must be recorded on the non-updated BAMM_object !!!
# # So maybe, I need to include this in deepSTRAPP::prepare_diversification_data() directly...
#   # Then need to update doc for deepSTRAPP::prepare_diversification_data() and update_rates_and_regimes
#   # And update the whales, Ponerinae, and Ponerinae_10My object so they have it recorded!
#
#
# ### 7.2/ Maximum Shift Credibility configuration (MSC) ####
#
# ## Get the most likely configuration = maximum shift credibility configuration (MSC)
#
# # Idea: find the posterior configuration that has the highest probability according to the marginal probability of the branches and number of shifts
#
# # Useful for big phylogenies when there are no clear MPP because of the many possible configurations
#
# # Identify the MSC configs
# MSC_detection <- maximumShiftCredibility(BAMM_posterior_samples_data)
# # Extract the MSC config
# MSC_tree <- subsetEventData(BAMM_posterior_samples_data, index = MSC_detection$sampleindex)
# plot.bammdata(MSC_tree, lwd = 1.25)
# addBAMMshifts(MSC_tree, par.reset = FALSE, cex = 2)
#
#
# ## Plot branch rates on the phylogeny: polar
#
# # pdf(file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/Ponerinae_diversity_dynamics_polar.pdf", height = 7, width = 6)
# pdf(file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/Ponerinae_diversity_dynamics_polar.pdf", height = 7, width = 6)
#
# par(mfrow = c(1,1), xpd = TRUE)
#
# # tau = length of segments used to discretize the continuous rates
# # tau = 0.001 => each segment is 1/1000 of the tree height (root age)
#
# # Quantile method to discretize rates in categories of equal frequencies
# branch_div_rates_plot <- plot.bammdata(BAMM_posterior_samples_data,
#                                        method = "polar",
#                                        tau = 0.001,
#                                        # breaksmethod = 'quantile',
#                                        breaksmethod = 'jenks',
#                                        lwd = 2)
#
# # Need to identify a configuration with the credible shifts!
# addBAMMshifts(BAMM_posterior_samples_data, index = MSC_detection$sampleindex,
#               method = "polar",
#               col = "black",
#               # bg = MSC_shifts_df$trend_col[-1], # Color shift by trend fo the new regime
#               # bg = MSC_shifts_df$rate_shift_type_col[-1], # Color shift by discrete rate shift on the edge
#               bg = "grey20",
#               par.reset = FALSE, cex = 2)
# title(main = 'Ponerinae diversity dynamics', cex = 2, line = 0)
# # text(x = 0.1, y = 0.9, labels = "Net Div. Rates", font = 2)
# addBAMMlegend(branch_div_rates_plot, location = "topleft", longFrac = 0.2, font = 2)
#
# dev.off()



### Modified version of BAMMtools::addBAMMshifts

## Handle adjustment of regime shift point size controlled by both cex and msp

# Source: BAMMtools::addBAMMshifts()
# Author: Mike Grundler

#' @importFrom graphics par points
#' @importFrom ape as.phylo branching.times .PlotPhyloEnv
#' @importFrom BAMMtools getShiftNodesFromIndex

addBAMMshifts_custom <- function (ephy, index = 1, method = "phylogram", cex = 1, pch = 21,
          col = 1, bg = 2, msp = NULL, shiftnodes = NULL, par.reset = TRUE,
          ...)
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
    # Key change: cex is update, not replaced by the msp size
    cex <- 0.5 + cex * 2.5 * msp$edge.length[msp$edge[, 2] %in%
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
  graphics::points(XX, YY, pch = pch, cex = cex, col = col, bg = bg,
         ...)
  if (par.reset) {
    graphics::par(op)
  }
}
