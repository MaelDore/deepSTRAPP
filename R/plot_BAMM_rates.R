
#' @title Plot diversification rates and regime shifts from BAMM on phylogeny
#'
#' @description Plot on a time-calibrated phylogeny the evolution of diversification rates and
#'   the location of regime shifts estimated from a BAMM (Bayesian Analysis of Macroevolutionary Mixtures).
#'   Each branch is colored accroding to the estimated rates of speciation, extinction, or net diversification
#'   stored in an object of class `bammdata`. Rates can vary along time, thus colors evolved along individual branches.
#'
#'   This function is a wrapper of original functions from the R package `{BAMMtools}`:
#'
#'   * Step 1: Use [BAMMtools::plot.bammdata()] to map rates on the phylogeny.
#'   * Step 2: Add the location of regime shifts with [BAMMtools::addBAMMshifts()] (if `add_regime_shifts = TRUE`).
#'
#' @param BAMM_object Object of class `"bammdata"`, typically generated with [deepSTRAPP::prepare_diversification_data()],
#'   that contains a phylogenetic tree and associated diversification rate mapping across selected posterior samples.
#'   It works also for `BAMM_object` updated for a specific `focal_time` using [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'   or the deepSTRAPP workflow with [deepSTRAPP::run_deepSTRAPP_for_focal_time] and [deepSTRAPP::run_deepSTRAPP_over_time()].
#' @param rate_type A character string specifying the type of diversification rates to plot.
#'   Must be one of 'speciation', 'extinction' or 'net_diversification' (default).
#' @param method A character string indicating the method for plotting the phylogenetic tree.
#'   * `method = "phylogram"` (default) plots the phylogenetic tree using rectangular coordinates.
#'   * `method = "polar"` plots the phylogenetic tree using polar coordinates.
#' @param add_regime_shifts Logical. Whether to add the location of regime shifts on the phylogeny (Step 2). Default is `TRUE`.
#' @param configuration_type A character string specifying how to select the location of regime shifts across posterior samples.
#'   * `configuration_type = "MAP"`: Use the average locations recorded in posterior samples with the Maximum A Posteriori probability (MAP) configuration.
#'     This regime shift configuration is the most frequent configuration among the posterior samples (See [BAMMtools::getBestShiftConfiguration()]).
#'     This is the default option.
#'   * `configuration_type = "MSC"`: Use the average locations recorded in posterior samples with the Maximum Shift Credibility (MSC) configuration.
#'     This regime shift configuration has the highest product of marginal probabilities across branches (See [BAMMtools::maximumShiftCredibility()]).
#'   * `configuration_type = "index"`: Use the configuration of a unique posterior sample those index is provided in `sample_index`.
#' @param sample_index Integer. Index of the posterior samples to use to plot the location of regime shifts.
#'   Used only if `configuration_type = index`. Default = `1`.
#' @param adjust_size_to_prob Logical. Whether to scale the size of the symbols showing the location of regime shifts according to
#'   the marginal shift probability of the shift happening on each location/branch. This will only works if there is an `$MSP_tree` element
#'   summarizing the marginal shift probabilities across branches in the `BAMM_object`. Default is `TRUE`.
#' @param regimes_fill Character string. Set the color of the background of the symbols showing the location of regime shifts.
#'   Equivalent to the `bg` argument in [BAMMtools::addBAMMshifts()]. Default is `"grey"`.
#' @param regimes_size Numerical. Set the size of the symbols showing the location of regime shifts.
#'   Equivalent to the `cex` argument in [BAMMtools::addBAMMshifts()]. Default is `1`.
#' @param regimes_pch Integer. Set the shape of the symbols showing the location of regime shifts.
#'   Equivalent to the `pch` argument in [BAMMtools::addBAMMshifts()]. Default is `21`.
#' @param regimes_border_col Character string. Set the color of the border of the symbols showing the location of regime shifts.
#'   Equivalent to the `col` argument in [BAMMtools::addBAMMshifts()]. Default is `"black"`.
#' @param regimes_border_width Numerical. Set the width of the border of the symbols showing the location of regime shifts.
#'   Equivalent to the `lwd` argument in [BAMMtools::addBAMMshifts()]. Default is `1`.
#' @param ... Additional graphical arguments to pass down to [BAMMtools::plot.bammdata()], [BAMMtools::addBAMMshifts()], and [par()].
#' @param display_plot Logical. Whether to display the plot generated in the R console. Default is `TRUE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#'
#' @export
#' @importFrom grDevices pdf dev.off gray
#' @importFrom BAMMtools plot.bammdata addBAMMshifts
#'
#' @details The main input `BAMM_object` is the typical output of [deepSTRAPP::prepare_diversification_data()].
#'   It provides information on rates and regimes shifts across the posterior samples of a BAMM.
#'
#'   `$MAP_BAMM_object` and `$MSC_BAMM_object` elements are required to plot regime shift locations following the
#'   "MAP" or "MSC" `configuration_type` respectively.
#'   A `$MSP_tree` element is required to scale the size of the symbols showing the location of regime shifts according marginal shift probabilities.
#'   (If `adjust_size_to_prob = TRUE`).
#'
#'   The default option to display regime shift is to use the average locations from the posterior samples with the Maximum A Posteriori probability (MAP) configuration.
#'   However, sometimes, multiple configurations have similarly high frequency in the posterior samples (See [BAMMtools::credibleShiftSet()].
#'   An alternative is to use the average locations from posterior samples with the Maximum Shift Credibility (MSC) configuration instead.
#'   This regime shift configuration has the highest product of marginal probabilities across branches where a shift is estimated.
#'   It may differ from the MAP configuration. (See [BAMMtools::maximumShiftCredibility()]).
#'
#' @return The function returns (invisibly) a list with three three elements similarly to [BAMMtools::plot.bammdata()].
#'  * `$coords`: A matrix of plot coordinates. Rows correspond to branches. Columns 1-2 are starting (x,y) coordinates of each branch and columns 3-4 are ending (x,y) coordinates of each branch. If method = "polar" a fifth column gives the angle(in radians) of each branch.
#'  * `$colorbreaks`: A vector of percentiles used to group macroevolutionary rates into color bins.
#'  * `$colordens`: A matrix of the kernel density estimates (column 2) of evolutionary rates (column 1) and the color (column 3) corresponding to each rate value.
#'
#' @author Maël Doré
#' @author Original functions by Mike Grundler & Pascal Title in R package `{BAMMtools}`.
#'
#' @seealso Initial functions in BAMMtools: [BAMMtools::plot.bammdata()] [BAMMtools::addBAMMshifts()]
#'
#' Associated functions in deepSTRAPP: [deepSTRAPP::prepare_diversification_data()] [deepSTRAPP::update_rates_and_regimes_for_focal_time()] [deepSTRAPP::run_deepSTRAPP_for_focal_time()] [deepSTRAPP::run_deepSTRAPP_over_time()]
#'
#' @examples
#' # Load BAMM output
#' data(whale_BAMM_object, package = "deepSTRAPP")
#'
#' ## Plot overall mean rates with MAP configuration for regime shifts
#' # (rates are averaged across all posterior samples)
#' plot_BAMM_rates(whale_BAMM_object, # Use the full object to plot mean overall rates
#'                 add_regime_shifts = TRUE,
#'                 configuration_type = "MAP",
#'                 regimes_size = 3, bg = "black")
#' ## Plot overall mean rates with MSC configuration for regime shifts
#' # (rates are averaged across all posterior samples)
#' plot_BAMM_rates(whale_BAMM_object, add_regime_shifts = TRUE,
#'                 configuration_type = "MSC",
#'                 regimes_size = 3, bg = "black")
#'
#' ## Plot mean MAP rates with regime shifts
#' # (rates are averaged only across MAP samples)
#' plot_BAMM_rates(whale_BAMM_object$MAP_BAMM_object, # Use the MAP object to plot mean MAP rates
#'                 add_regime_shifts = TRUE,
#'                 configuration_type = "index",
#'                 # Set to index to use the regime shift location from MAP configuration
#'                 regimes_size = 3, bg = "black")
#' ## Plot mean MSC rates with regime shifts
#' # (rates averaged only across MSC samples)
#' plot_BAMM_rates(whale_BAMM_object$MSC_BAMM_object, # Use the MSC object to plot mean MSC rates
#'                 add_regime_shifts = TRUE,
#'                 configuration_type = "index",
#'                 # Set to index to use the regime shift data from MSC configuration
#'                 regimes_size = 3, bg = "black")
#'


plot_BAMM_rates <- function (BAMM_object,
                             rate_type = "net_diversification",
                             method = "phylogram",
                             add_regime_shifts = TRUE,
                             configuration_type = "MAP", # MAP, MSC, or index
                             sample_index = 1,
                             adjust_size_to_prob = TRUE, # To adjust size of points representing regime shifts to their marginal posterior probabilities
                             regimes_fill = "grey", # Replace 'bg' argument in BAMMtools::addBAMMshifts()
                             regimes_size = 1, # Replace 'cex' argument in BAMMtools::addBAMMshifts()
                             regimes_pch = 21, # Replace 'pch' argument in BAMMtools::addBAMMshifts()
                             regimes_border_col = "black", # Replace 'col' argument in BAMMtools::addBAMMshifts()
                             regimes_border_width = 1, # Replace 'lwd' argument in BAMMtools::addBAMMshifts()
                             ..., # To pass down to BAMMtools::plot.bammdata(), BAMMtools::addBAMMshifts(), and par()
                             display_plot = TRUE,
                             PDF_file_path = NULL)
{
  ### Check input validity
  {
    ## BAMM_object
    # BAMM_object must be of class "bammdata"
    if (!("bammdata" %in% class(BAMM_object)))
    {
      stop("'BAMM_object' must have the 'bammdata' class. See ?deepSTRAPP::prepare_diversification_data() to learn how to generate those objects.")
    }

    ## rate_type must be either "speciation", "extinction" or "net_diversification"
    if (!(rate_type %in% c("speciation", "extinction", "net_diversification")))
    {
      stop("'rate_type' can only be 'speciation', 'extinction', or 'net_diversification'.")
    }

    ## method is either 'phylogram' or 'polar'
    if (!(method %in% c("phylogram", "polar")))
    {
      stop("'method' can only be 'phylogram', 'polar'.")
    }

    ## configuration_type is either 'MAP', 'MSC' or 'index'
    if (!(configuration_type %in% c("MAP", "MSC", "index")))
    {
      stop("'configuration_type' should be either 'MAP', 'MSC', or 'index'.")
    }
    # Need $MAP_BAMM_object if configuration_type = 'MAP'
    if ((configuration_type == "MAP") & !("MAP_BAMM_object" %in% names(BAMM_object)))
    {
      stop(paste0("'BAMM_object' must have a '$MAP_BAMM_object' to be able to plot regime shifts from samples with the MAP configuration.\n",
                  "See ?deepSTRAPP::prepare_diversification_data() to learn how to generate those objects."))
    }
    # Need $MSC_BAMM_object if configuration_type = 'MSC'
    if ((configuration_type == "MSC") & !("MSC_BAMM_object" %in% names(BAMM_object)))
    {
      stop(paste0("'BAMM_object' must have a '$MSC_BAMM_object' to be able to plot regime shifts from samples with the MSC configuration.\n",
                  "See ?deepSTRAPP::prepare_diversification_data() to learn how to generate those objects."))
    }

    ## sample_index must be an integer from 1 to the length of EventData
    # sample_index must be a positive integer.
    if ((sample_index != abs(sample_index)) | (sample_index != round(sample_index)))
    {
      stop(paste0("'sample_index' must be a positive integer defining the ID of the posterior sample used to plot regime shifts."))
    }
    if (!(sample_index %in% 1:length(BAMM_object$eventData)))
    {
      stop(paste0("'sample_index' must be compatible with the posterior samples available in 'BAMM_object'.\n",
                  "The 'BAMM_object' provided has diversification histories from ",length(BAMM_object$eventData)," posterior samples.\n",
                  "Current 'sample_index' is: ",sample_index,"."))
    }

    ## Need $MSP_tree if using 'adjust_size_to_prob = TRUE'
    # BAMM_object must be of class "bammdata"
    if (!("bammdata" %in% class(BAMM_object)))
    {
      stop("'BAMM_object' must have the 'bammdata' class. See ?deepSTRAPP::prepare_diversification_data() to learn how to generate those objects.")
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

  ## Save initial par() and reassign them on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ## Convert 'rate_type' into 'spex'
  if (rate_type == "net_diversification") { spex <- "netdiv" }
  if (rate_type == "speciation") { spex <- "s" }
  if (rate_type == "extinction") { spex <- "e" }

  ## Filter list of additional arguments to avoid warnings from par()
  add_args <- list(...)
  args_names_for_plot <- c("tau", "xlim", "ylim", "vtheta", "rbf", "show", "labels", "legend",
                           "spex", "lwd", "cex", "pal", "mask", "mask.color", "colorbreaks", "logcolor",
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
  # if ("show" %in% names(add_args_for_plot)) { show <- add_args_for_plot$show } else { show <- TRUE }
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

  ## Display plot if requested
  if (display_plot)
  {
    # ## Plot rates
    # output <- BAMMtools::plot.bammdata(x = BAMM_object,
    #                                    method = method,
    #                                    ...)

    ## Plot rates while separating names arguments from additional arguments for par() in the ellipsis (...)
    output <- do.call(what = BAMMtools::plot.bammdata,
                      args = c(list(x = BAMM_object, method = method, tau = tau, xlim = xlim, ylim = ylim, vtheta = vtheta,
                                    rbf = rbf, show = display_plot, labels = labels, legend = legend, spex = spex, lwd = lwd,
                                    cex = cex, pal = pal, mask = mask, mask.color = mask.color, colorbreaks = colorbreaks,
                                    logcolor = logcolor, breaksmethod = breaksmethod, color.interval = color.interval,
                                    JenksSubset = JenksSubset, par.reset = par.reset, direction = direction),
                               add_args_for_par))

    ## Plot regimes if requested
    if (add_regime_shifts)
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

      # Case for Maximum Shift Credibility (MSC) configuration
      if (configuration_type == "MSC")
      {
        # ## Add regime shifts on the plot
        # BAMMtools::addBAMMshifts(ephy = BAMM_object$MSC_BAMM_object,
        #                          index = 1,
        #                          method = method,
        #                          msp = msp,
        #                          bg = regimes_fill,
        #                          cex = regimes_size,
        #                          pch = regimes_pch,
        #                          col = regimes_border_col,
        #                          lwd = regimes_border_width
        #                          ...)

        ## Use the BAMM_object$MSC_BAMM_object to get locations of shifts
        do.call(what = addBAMMshifts_custom,
                args = c(list(ephy = BAMM_object$MSC_BAMM_object,
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
    }
  }

  ## Save plot if requested
  if (!is.null(PDF_file_path))
  {
    # Adjust width and height according to phylo
    nb_tips <- length(BAMM_object$tip.label)
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10

    ## Force the plot to be displayed in the exported graphics device
    display_plot <- TRUE

    ## Open PDF
    grDevices::pdf(file = file.path(PDF_file_path),
                   width = width, height = height)

    # ## Plot rates
    # output <- BAMMtools::plot.bammdata(x = BAMM_object,
    #                                    method = method,
    #                                    ...)

    ## Plot rates while separating names arguments from additional arguments for par() in the ellipsis (...)
    output <- do.call(what = BAMMtools::plot.bammdata,
                      args = c(list(x = BAMM_object, method = method, tau = tau, xlim = xlim, ylim = ylim, vtheta = vtheta,
                                    rbf = rbf, show = display_plot, labels = labels, legend = legend, spex = spex, lwd = lwd,
                                    cex = cex, pal = pal, mask = mask, mask.color = mask.color, colorbreaks = colorbreaks,
                                    logcolor = logcolor, breaksmethod = breaksmethod, color.interval = color.interval,
                                    JenksSubset = JenksSubset, par.reset = par.reset, direction = direction),
                               add_args_for_par))

    ## Plot regimes if requested
    if (add_regime_shifts)
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

      # Case for Maximum Shift Credibility (MSC) configuration
      if (configuration_type == "MSC")
      {
        # ## Add regime shifts on the plot
        # BAMMtools::addBAMMshifts(ephy = BAMM_object$MSC_BAMM_object,
        #                          index = 1,
        #                          method = method,
        #                          msp = msp,
        #                          bg = regimes_fill,
        #                          cex = regimes_size,
        #                          pch = regimes_pch,
        #                          col = regimes_border_col,
        #                          lwd = regimes_border_width
        #                          ...)

        ## Use the BAMM_object$MSC_BAMM_object to get locations of shifts
        do.call(what = addBAMMshifts_custom,
                args = c(list(ephy = BAMM_object$MSC_BAMM_object,
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
    }

    ## Close PDF
    grDevices::dev.off()
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
