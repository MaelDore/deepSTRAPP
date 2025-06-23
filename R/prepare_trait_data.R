## Functions to prepare trait data by mapping their evolution on the phylogeny
# One master function to select the proper pipeline according to data type
# Three sub-functions mapping trait evolution according to data type

#' @title Map trait evolution on a time-calibrated phylogeny
#'
#' @description Map trait evolution on a time-calibrated phylogeny in several steps:
#'
#'   * Step 1: Fit evolutionary models to trait data using Maximum Likelihood.
#'   * Step 2: Select the best fitting model comparing AICc.
#'   * Step 3: Infer ancestral characters estimates (ACE) at nodes.
#'   * Step 4: Run stochastic mapping simulations to generate evolutionary histories
#'     compatible with the best model and inferred ACE. (Only for categorical and biogeographic data)
#'   * Step 5: Infer ancestral states along branches.
#'     - For continuous traits: use interpolation to produce a `contMap`.
#'     - For categorical and biogeographic data: compute posterior frequencies of each state/range
#'       to produce a `densityMap` for each state/range.
#'
#' @param tip_data Named numerical or character string vector of trait values/states/ranges at tips.
#'   Names should be ordered as the tip labels in the phylogeny foudn in `phylo$tip.label`.
#' @param trait_data_type Character string. Type of trait data. Either: "continuous", "categorical" or "biogeographic".
#' @param phylo Time-calibrated phylogeny. Object of class `"phylo"` as defined in [ape].
#'   Tip labels (`phylo$tip.label`) should match names in `tip_data`.
#' @param evolutionary_models (Vector of) character string(s). To provide the set of evolutionary models to fit on the data.
#'   * Models available for continuous data are detailed in [geiger::fitContinuous()].
#'   * Models available for categorical data are detailed in [geiger::fitDiscrete()].
#'   * Models for biogeographic data are fit with `BioGeoBEARS`.
#'   * See list in "Details" section.
#' @param Q_matrix Custom Q-matrix for categorical data representing transition classes between states.
#'   Elements that are zero signify rates that are fixed to zero (i.e., impossible transition).
#' @param ... Additional arguments to be passed down to the functions used to fit models (See `evolutionary_models`) and produce simmaps (See [phytools::make.simmap()]).
#' @param res Integer. Define the number of time steps used to interpolate/estimate trait value/state/range in `contMap`/`densityMaps`.
#' @param nb_simulations Integer. Define the number of simulations generated for stochastic mapping. Default = 1000. Only for "categorical" and "biogeographic" data.
#' @param colors_per_states Named character string. To set the colors to use to map each state posterior probabilities. Names = states; values = colors.
#'  If `NULL` (default), the `rainbow()` color scale will be used. Only for categorical and biogeographic data.
#' @param plot_map Logical. Whether to plot or not the phylogeny with mapped trait evolution.
#' @param plot_overlay Logical. If `TRUE` (default), plot a unique `densityMap` with overlapping states/ranges using transparency.
#'    If `FALSE`, plot a `densityMap` per state/range. Only for "categorical" and "biogeographic" data.
#' @param add_ACE_pies Logical. Whether to add pies of posterior probabilities of states/ranges at internal nodes on the mapped phylogeny. Default = `TRUE`.
#'    Only for categorical and biogeographic data.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with '.pdf'.
#' @param return_ace Logical. Whether the named vector of ancestral characters estimates (ACE) at internal nodes should be returned in the output. Default = `TRUE`.
#' @param return_simmaps Logical. Whether the evolutionary histories simulated during stochastic mapping (i.e., `simmaps`) should be returned in the output.
#'  Default = `TRUE`. Only for "categorical" and "biogeographic" data.
#' @param return_best_model_fit Logical. Whether to include the output of the best fitting model in the function output. Default = `FALSE`.
#' @param return_model_selection_df Logical. Whether to include the data.frame summarizing model comparisons used to select the best fitting model should be returned in the output. Default = `FALSE`.
#' @param verbose Logical. Should progression be displayed? A message will be printed for every steps in the process. Default is `FALSE`.
#'
#' @export
#' @importFrom stats logLik
#' @importFrom geiger fitContinuous
#' @importFrom ape all.equal.phylo
#' @importFrom phytools rescale fastAnc contMap densityMap rescaleSimmap as.Qmatrix setMap
#' @importFrom grDevices pdf dev.off rainbow colorRampPalette col2rgb rgb
#' @importFrom methods hasArg
#' @importFrom stats setNames
#'
#' @details Map trait evolution on a time-calibrated phylogeny in several steps:
#'
#'  Step 1: Models are fit using Maximum Likelihood approach:
#'    * For "continuous" data models are fit with [geiger::fitContinuous()]: "BM", "OU", "EB", "rate_trend", "lambda", "kappa", "delta".
#'    * For "categorical" data models are fit with [geiger::fitDiscrete()]: "ER", "SYM", "ARD".
#'    * For "biogeographic" data models are fit with R package `BioGeoBEARS`: "BAYAREALIKE", "DIVALIKE", "DEC", "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J".
#'
#'  Step 2: Best model is identified among the list of `evolutionary_models` by comparing the corrected AIC (AICc)
#'    and selecting the  model with lowest AICc.
#'
#'  Step 3: For continuous traits: Ancestral characters estimates (ACE) are inferred with [phytools::fastAnc] on a tree
#'    with modified branch lengths scaled to reflect the evolutionary rates estimated from the best model using [phytools::rescale()].
#'
#'  Step 4: Stochastic Mapping.
#'
#'    For categorical and biogeographic data, stochastic mapping simulations are performed to generate evolutionary histories
#'    compatible with the best model and inferred ACE. Node states/ranges are drawn from the scaled marginal likelihoods of ACE,
#'    and states/ranges shifts along branches are simulated according to the transition matrix Q estimated from the best fitting model.
#'
#'  Step 5: Infer ancestral states along branches.
#'    * For continuous traits: ancestral trait values along branches are interpolated with [phytools::contMap()].
#'      This provides quick estimates of trait value at any point in time, but it does not provide accurate ML estimates in
#'      case of models that are time or trait-value dependent (such as "EB" or "OU") as the interpolation used to built the contMap is assuming
#'      a constant rate along each branch. However, ancestral trait values at nodes remain accurate
#'    * For categorical and biogeographic data: compute posterior frequencies of each state/range among the simulated evolutionary histories (`simmaps`)
#'      to produce a `densityMap` for each state/range that reflects the changes along branches in probability of harboring a given state/range.
#'
#'  # Note on macroevolutionary models of trait evolution
#'
#'  This function provides an easy solution to map trait evolution on a time-calibrated phylogeny
#'  and obtain the `contMap`/`densityMaps` objects needed to run the deepSTRAPP workflow ([run_deepSTRAPP_for_focal_time], [run_deepSTRAPP_over_time]).
#'  However, it does not explore the most complex options for trait evolution. You may need to explore more complex models to capture the dynamics of trait evolution.
#'  such as trait-dependent multi-rate models ([phytools::brownie.lite()], [OUwie::OUwie]), Bayesian reversible jump MCMC implementations allowing a thorough exploration
#'  of location and number of regime shifts (Ex: BayesTraits, BAMM), or RRphylo for a penalized phylogenetic ridge regression approach that allows regime shifts across all branches.
#'
#' @return The function returns a list with at least two elements.
#'
#'   * `$contMap` (For "continuous" data) Object of class `"contMap"`, typically generated with [phytools::contMap()],
#'     that contains a phylogenetic tree and associated continuous trait mapping.
#'   * `$densityMaps` (For "categorical" and "biogeographic" data) List of objects of class `"densityMap`,
#'     typically generated with [phytools::densityMap()], that contains a phylogenetic tree and associated mapping of probability
#'     to harbor a given state/range along branches. The list contains one `"densityMap` per state/range found in the `tip_data`.
#'   * `$trait_data_type` Character string. Record the type of trait data. Either: "continuous", "categorical" or "biogeographic".
#'
#'   If `return_ace = TRUE`,
#'   * `$ace` For continuous traits: Named vector that record the ancestral characters estimates (ACE) at internal nodes.
#'     For categorical and biogeographic data: Matrix that record the posterior probabilities of ancestral states/ranges (characters) estimates (ACE) at internal nodes.
#'     Rows are internal nodes. Columns are states/ranges. Values are posterior probabilities of each state per node.
#'
#'   If `return_best_model_fit = TRUE`,
#'   * `$best_model_fit` List that provides the output of the best fitting model.
#'
#'   If `model_selection_df = TRUE`,
#'   * `$model_selection_df` Data.frame that summarizes model comparisons used to select the best fitting model.
#'
#' @author Maël Doré
#'
#' @seealso [geiger::fitContinuous()] [geiger::fitDiscrete()] [phytools::contMap()] [phytools::densityMap()]
#'
#' @references For macroevolutionary models in geiger: Pennell, M. W., Eastman, J. M., Slater, G. J., Brown, J. W., Uyeda, J. C., FitzJohn, R. G., ... & Harmon, L. J. (2014).
#'  geiger v2. 0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics, 30(15), 2216-2218..
#'  \url{https://doi.org/10.1093/bioinformatics/btu181}.
#'
#'  For BioGeoBEARS: Matzke, Nicholas J. (2018). BioGeoBEARS: BioGeography with Bayesian (and likelihood) Evolutionary Analysis with R Scripts.
#'    version 1.1.1, published on GitHub on November 6, 2018. DOI: \url{http://dx.doi.org/10.5281/zenodo.1478250}. Website: \url{http://phylo.wikidot.com/}.
#'
#' @examples
#' # ----- Example 1: Continuous data ----- #
#'
#' # Load phylogeny and tip data
#' library(phytools)
#' data(eel.tree)
#' data(eel.data)
#'
#' # Extract body size
#' eel_data <- stats::setNames(eel.data$Max_TL_cm,
#'                             rownames(eel.data))
#'
#' # Map trait evolution on the phylogeny, selecting among four models ("BM", "OU", "lambda", "kappa")
#' mapped_cont_traits <- prepare_trait_data(
#'    tip_data = eel_data,
#'    trait_data_type = "continuous",
#'    phylo = eel.tree,
#'    evolutionary_models = c("BM", "OU", "lambda", "kappa"),
#'    # Example of an additional argument ('control') that can be provided to geiger::fitContinuous()
#'    control = list(niter = 200),
#'    plot_map = FALSE,
#'    return_best_model_fit = TRUE,
#'    return_model_selection_df = TRUE,
#'    verbose = TRUE)
#'
#' # Explore output
#' plot(mapped_cont_traits$contMap) # contMap with interpolated trait values
#' mapped_cont_traits$model_selection_df # Summary of model selection
#' # Parameter estimates and optimization summary of the best model
#' # (Here, the best model is Pagel's lambda)
#' mapped_cont_traits$best_model_fit$opt
#' mapped_cont_traits$ace # Ancestral character estimates at internal nodes
#'
#' # ----- Example 2: Categorical data ----- #
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
#' colors_per_states <- c("limegreen", "orange", "dodgerblue")
#' names(colors_per_states) <- c("bite", "kiss", "suction")
#'
#' mapped_cat_traits <- prepare_trait_data(tip_data = eel_data, phylo = eel.tree,
#'                                         trait_data_type = "categorical",
#'                                         colors_per_states = colors_per_states,
#'                                         evolutionary_models = c("ER", "ARD", "matrix"),
#'                                         Q_matrix = Q_matrix,
#'                                         nb_simulations = 10, # Set to 10 to save time.
#'                                         # But recommended value = 1000.
#'                                         plot_map = TRUE,
#'                                         plot_overlay = TRUE,
#'                                         return_best_model_fit = TRUE,
#'                                         return_model_selection_df = TRUE)
#'
#' # Explore output
#' plot(mapped_cat_traits$densityMaps[[1]]) # densityMap for state n°1 ("bite")
#' mapped_cat_traits$model_selection_df # Summary of model selection
#' # Parameter estimates and optimization summary of the best model
#' # (Here, the best model is ER)
#' print(mapped_cat_traits$best_model_fit)$ # Summary of the best evolutionary model
#' mapped_cat_traits$ace # Posterior probabilities of each state (= ACE) at internal nodes
#'
#'
#' # ----- Example 3: Biogeographic data ----- #
#'
#' # TBA
#'


## Make a different function for each type of data, then make a wrapper function for all types
# prepare_trait_data() = wrapper
# Detect type of data based on $trait_data_type, but check it match the input (tip_data)
# prepare_trait_data_for_continuous_data() = For continuous traits
# prepare_trait_data_for_categorical_data() = For categorical traits
# prepare_trait_data_for_biogeographic_data() = For biogeographic traits

# Add an option to get ACE for categorical and biogeographic traits too (as scaled marginal likelihoods)

# Default for categorical = "ARD". "All" = ER, SYM, ARD.
# May need to use more complex models to capture the dynamics of trait evolution. Custom constrains in the Q matrix, Multi-regimes models, ... (cite packages)

# Default for Biogeographic = "DEC+J". "All" = DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J.
# May need to use more complex models to capture the dynamics of trait evolution. +W, +X, DECX, ... (cite packages/refs)

# Check reference to contMap in the doc => add DensityMaps


### Master function to prepare data and select the proper test function according to data type ####

prepare_trait_data <- function (
    tip_data,
    trait_data_type,
    phylo,
    evolutionary_models = NULL, # Default = "BM" for continuous data; "ARD" for categorical; "DEC+J" for biogeographic
    Q_matrix = NULL, # Custom Q-matrix for categorical data
    ..., # To allow to pass down arguments in the functions used to fit the models
    res = 100, # Number of time steps used to interpolate trait value in the contMap
    nb_simulations = 1000, # Only for categorical and biogeographic data
    colors_per_states = NULL, # Only for categorical and biogeographic data. To set the colors to use to map each state posterior probabilities
    plot_map = TRUE,
    plot_overlay = TRUE, # Only for categorical and biogeographic data
    add_ACE_pies = TRUE, # Only for categorical and biogeographic data
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_simmaps = TRUE, # Only for categorical and biogeographic data
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE)
{
  ### Check input validity
  {
    ## tip_data
    # tip_data must have the same names as the tip.label in the phylogeny
    if (!all((names(tip_data) %in% phylo$tip.label)))
    {
      stop(paste0("Names in 'tip_data' must match with '$tip.label' in 'phylo'."))
    }
    # Check if tip_data need to be reorder to match phylo$tip.label. If so, do it, but send a warning
    if (!all(names(tip_data) == phylo$tip.label))
    {
      tip_data <- tip_data[match(phylo$tip.label, table = names(tip_data))]
      warning(paste0("Entries in 'tip_data' were reordered to match 'phylo$tip.label."))
    }

    ## trait_data_type
    # trait_data_type must be "continuous", categorical" or "biogeographic"
    if (!(trait_data_type %in% c("continuous", "categorical", "biogeographic")))
    {
      stop("'trait_data_type' can only be 'continuous', 'categorical', or 'biogeographic'.")
    }
    # Check that trait_data_type is coherent with what is found in tip_data
    if ((trait_data_type == "continuous") & !is.numeric(tip_data))
    {
      stop("For 'trait_data_type = continuous', 'tip_data' must be a named numeric vector.")
    }
    if ((trait_data_type == "categorical") & !is.character(tip_data))
    {
      stop("For 'trait_data_type = categorical', 'tip_data' must be a named vector of character strings providing tip states.")
    }
    if ((trait_data_type == "biogeographic") & !is.character(tip_data))
    {
      stop("For 'trait_data_type = biogeographic', 'tip_data' must be a named vector of character strings providing tip ranges.")
    }

    ## phylo
    # phylo must be a "phylo" class object
    if (!("phylo" %in% class(phylo)))
    {
      stop("'phylo' must have the 'phylo' class. See ?ape::read.tree() and ?ape::read.nexus() to learn how to import phylogenies in R.")
    }
    # phylo must be rooted
    if (!(ape::is.rooted(phylo)))
    {
      stop(paste0("'phylo' must be a rooted phylogeny."))
    }
    # phylo must be fully resolved/dichotomous
    if (!(ape::is.binary(phylo)))
    {
      stop(paste0("'phylo' must be a fully resolved/dichotomous/binary phylogeny."))
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

    ## Other checks are carried in dedicated sub-functions
  }

  # ## Catch additional arguments
  # add_args <- list(...)

  ## Compute the appropriate internal function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous =   { # Case for continuous data
           # Output = contMap. Trait values are interpolated along branches.
           trait_data_output <- prepare_trait_data_for_continuous_data(
             tip_data = tip_data,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "BM" for continuous data
             ..., # Additional arguments for geiger::fitContinuous()
             res = res,
             plot_map = plot_map,
             PDF_file_path = PDF_file_path,
             return_ace = return_ace,
             return_best_model_fit = return_best_model_fit,
             return_model_selection_df = return_model_selection_df,
             verbose = verbose
           )
         },
         categorical =  { # Case for categorical data
           # Output = densityMap (+ simmaps). Include simulations of trait states.
           trait_data_output <- prepare_trait_data_for_categorical_data(
             tip_data = tip_data,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "ARD" for categorical data
             Q_matrix = Q_matrix, # Custom Q-matrix for categorical data
             ..., # Additional arguments for geiger::fitDiscrete()
             res = res,
             nb_simulations = nb_simulations, # Only for categorical and biogeographic data
             colors_per_states = colors_per_states, # Only for categorical and biogeographic data
             plot_map = plot_map,
             plot_overlay = plot_overlay, # Only for categorical and biogeographic data
             add_ACE_pies = add_ACE_pies, # Only for categorical and biogeographic data
             PDF_file_path = PDF_file_path,
             return_ace = return_ace,
             return_simmaps = return_simmaps, # Only for categorical and biogeographic data
             return_best_model_fit = return_best_model_fit,
             return_model_selection_df = return_model_selection_df,
             verbose = verbose
           )
         },
         biogeographic = { # Case for biogeographic data
           # Output = densityMap (+ simmaps). Include simulations of geographic ranges.
           trait_data_output <- prepare_trait_data_for_biogeographic_data(
             tip_data = tip_data,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "DEC+J" for biogeographic data
             ..., # Additional arguments for BioGeoBEARS functions
             nb_simulations = nb_simulations, # Only for categorical and biogeographic data
             res = res,
             plot_map = plot_map,
             plot_overlay = plot_overlay, # Only for categorical and biogeographic data
             PDF_file_path = PDF_file_path,
             return_ace = return_ace,
             return_simmaps = return_simmaps, # Only for categorical and biogeographic data
             return_best_model_fit = return_best_model_fit,
             return_model_selection_df = return_model_selection_df,
             verbose = verbose
           )
         }
  )

  ## Export the output
  return(invisible(trait_data_output))
}


### Sub-function to handle continuous data ####

prepare_trait_data_for_continuous_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "BM", # Default = "BM" for continuous data
    ..., # Additional arguments for geiger::fitContinuous()
    res = 100,
    plot_map = TRUE,
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE
    )
{
  ### Check input validity
  {
    ## evolutionary_models
    if (!is.null(evolutionary_models))
    {
      if (!all(evolutionary_models %in% c("BM", "OU", "EB", "rate_trend", "lambda", "kappa", "delta")))
      {
        stop(paste0("For 'trait_data_type = continuous', 'evolutionary_models' must be selected among: 'BM', 'OU', 'EB', 'rate_trend', 'lambda', 'kappa', 'delta'.\n",
                    "See details in ?geiger::fitContinuous()."))
      }
    }

    ## res
    # res must be a positive integer.
    if ((res != abs(res)) | (res != round(res)))
    {
      stop(paste0("'res' must be a positive integer defining the number of time steps used to interpolate trait value in the contMap."))
    }

  }

  # Set default model (BM) if absent
  if (is.null(evolutionary_models)) { evolutionary_models <- "BM" }

  ##	Run all models and store their results in a list
  # ?geiger::fitContinuous

  nb_models <- length(evolutionary_models)
  if (verbose) { cat(paste0(Sys.time(), " - Fit ",nb_models," evolutionary model(s): ", paste(evolutionary_models, collapse = ", "), ".\n\n")) }

  # Initiate list to store models outputs
  models_fits <- list()

  # Fit BM model
  if (any(evolutionary_models == "BM"))
  {
    BM_ID <- which(evolutionary_models == "BM")
    BM_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "BM", ...)

    if (verbose) { print(BM_fit) ; cat("\n") }
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[BM_ID]] <- BM_fit
  }

  # Fit OU model
  if (any(evolutionary_models == "OU"))
  {
    OU_ID <- which(evolutionary_models == "OU")
    OU_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "OU", ...)

    if (verbose) { print(OU_fit) ; cat("\n") }
    # sigsq = the drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space
      # It is not the evolutionary rate since such rate quantifying the ratio of movement in trait space per unit of time results of the combine action of the drift parameter and the strength of selection (alpha))
    # alpha = parameter describing the strength of selection pushing trait evolution toward the optimum
    # z0 = ancestral root state

    # Store output
    models_fits[[OU_ID]] <- OU_fit
  }

  # Fit EB model
  if (any(evolutionary_models == "EB"))
  {
    EB_ID <- which(evolutionary_models == "EB")
    EB_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "EB", ...)

    if (verbose) { print(EB_fit) ; cat("\n") }
    # sigsq = the initial evolutionary rate at t = 0
    # a = time-dependency parameter describing the shape of the exponential decrease of sigsq though time
      # If a < 0 => the model is an EB with rates decaying in time
      # If a > 0 => the model is a late burst with rates increasing exponentially in time
    # z0 = ancestral root state

    # Store output
    models_fits[[EB_ID]] <- EB_fit
  }

  # Fit rate_trend model
  if (any(evolutionary_models == "rate_trend"))
  {
    rate_trend_ID <- which(evolutionary_models == "rate_trend")
    rate_trend_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "rate_trend", ...)

    if (verbose) { print(rate_trend_fit) ; cat("\n") }
    # slope = coefficient of the linear trend in rates through time
      # If slope < 0 => rates decrease linearly in time
      # If slope > 0 => rates increase linearly in time
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[rate_trend_ID]] <- rate_trend_fit
  }

  # Fit Pagel's lambda model
  if (any(evolutionary_models == "lambda"))
  {
    lambda_ID <- which(evolutionary_models == "lambda")
    lambda_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "lambda", ...)

    if (verbose) { print(lambda_fit) ; cat("\n") }
    # lambda = Internal branch length multiplier = modulates strength of the phylogenetic signal
      # Close to 0 => star-like phylogeny = no phylogenetic structure/signal
      # Close to 1 => BM = "perfect" phylogenetic signal
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[lambda_ID]] <- lambda_fit
  }

  # Fit Pagel's kappa model
  if (any(evolutionary_models == "kappa"))
  {
    kappa_ID <- which(evolutionary_models == "kappa")
    kappa_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "kappa", ...)

    if (verbose) { print(kappa_fit) ; cat("\n") }
    # kappa = Branch length exponent = modulates effect of cladogenesis on trait evolution (punctuational/speciational) model of trait evolution by raising all branch length to power kappa
      # Close to 0 => all branches have equal length, thus only the number of speciation events affect trait evolution.
      # Close to 1 => No effect of cladogenesis = BM.
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[kappa_ID]] <- kappa_fit
  }

  # Fit Pagel's delta model
  if (any(evolutionary_models == "delta"))
  {
    delta_ID <- which(evolutionary_models == "delta")
    delta_fit <- geiger::fitContinuous(phy = phylo, dat = tip_data, model = "delta", ...)

    if (verbose) { print(delta_fit) ; cat("\n") }
    # delta = Node height exponent = time-dependent model that adjusts the relative contributions of early versus late evolution in the tree.
      # delta < 1 => relative node heights are reduced, particularly for shallower (closer to the tips) branches. Equivalent to an Early Burst model when most trait evolution happens in early times.
      # delta = 1 => No time-dependency = BM.
      # delta > 1 => relative node heights are increased, particularly for shallower (closer to the tips) branches. Equivalent to a Late Burst model when most trait evolution happens in recent times.
    # sigsq = drift parameter quantifying the ability of trait to evolve quickly and drift in the trait space = trait evolutionary rate
    # z0 = ancestral root state

    # Store output
    models_fits[[delta_ID]] <- delta_fit
  }

  # Name model fits according to the associated models
  names(models_fits) <- evolutionary_models

  ## Compare model fits and select best model

  # Use it to compare model fits with AICc and Akaike's weights
  models_comparison <- select_best_trait_model_from_geiger(list_model_fits = models_fits)
  best_model_name <- models_comparison$best_model_name
  best_model_fit <- models_comparison$best_model_fit

  # Display result of model comparison
  if (verbose)
  {
    cat(paste0(Sys.time(), " - Compare model fits.\n\n"))
    print(models_comparison$models_comparison_df)
  }

  ### Get ACE (useful in all cases to build the contMap)

  if (verbose) { cat(paste0("\n", Sys.time(), " - Infer Ancestral Character Estimates from the best fitting model: ",best_model_name,".\n")) }

  ## Rescale the tree such as the relative branch length account for trait rates in the best fitting model
  # ?geiger::rescale.phylo

  # Adjust model name
  if (best_model_name == "rate_trend") { best_model_name <- "trend" }
  # Create function to rescale phylo
  rescaled_phylo_fn <- phytools::rescale(x = phylo, model = best_model_name)
  # Extract parameters from best model fit
  fun_args <- formals(fun = rescaled_phylo_fn)
  best_model_args <- best_model_fit$opt[names(fun_args)]
  # Rescale phylogeny using the estimated parameters from the best model
  rescaled_phylo <- do.call(what = rescaled_phylo_fn, args = best_model_args)

  ## Run ACE inference with BM on transformed tree
  # ?phytools::fastAnc

  ACE_output <- phytools::fastAnc(tree = rescaled_phylo,
                                  x = tip_data,
                                  vars = FALSE, # Compute the variance of ancestral state estimates
                                  CI = FALSE) # Compute the 95% CI of ancestral state estimates
  # ACE_output

  ## Create contMap using ACE computed with the best fitting model

  if (verbose) { cat(paste0("\n", Sys.time(), " - Create contMap by interpolating values along branches.\n\n")) }

  contMap <- phytools::contMap(tree = phylo,
                               method = "user",
                               x = tip_data,
                               anc.states = ACE_output,
                               res = res, # Number of time steps
                               plot = plot_map) # Plot only if requested

  ## Export contMap in PDF if requested
  if (!is.null(PDF_file_path))
  {
    nb_tips <- length(phylo$tip.label)
    grDevices::pdf(file = file.path(PDF_file_path),
                   width = nb_tips/60*8, height = nb_tips/60*10)

    plot(contMap)

    grDevices::dev.off()

  }

  ## Build output
  output <- list(contMap = contMap,
                 trait_data_type = "continuous")
  # Include ACE if requested
  if(return_ace) { output$ace <- ACE_output }
  # Include output of best model if requested
  if(return_best_model_fit) { output$best_model_fit <- best_model_fit }
  # Include df for model comparison if requested
  if(return_model_selection_df) { output$model_selection_df <- models_comparison$models_comparison_df }

  ## Return output
  return(invisible(output))

}


### Sub-function to handle categorical data ####

prepare_trait_data_for_categorical_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "ARD", # Default = "ARD" for categorical
    Q_matrix = NULL, # Custom Q-matrix for categorical data
    ..., # Additional arguments for geiger::fitDiscrete()
    res = 100,
    nb_simulations = 1000, # Only for categorical and biogeographic data
    colors_per_states = NULL,
    plot_map = TRUE,
    plot_overlay = TRUE, # Only for categorical and biogeographic data
    add_ACE_pies = TRUE, # Only for categorical and biogeographic data
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_simmaps = TRUE, # Only for categorical and biogeographic data
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE)
{
  ### Check input validity
  {
    ## evolutionary_models
    if (!is.null(evolutionary_models))
    {
      if (!all(evolutionary_models %in% c("ER", "SYM", "ARD", "meristic", "matrix")))
      {
        stop(paste0("For 'trait_data_type = categorical', 'evolutionary_models' must be selected among: 'ER', 'SYM', 'ARD', 'meristic', 'matrix'.\n",
                    "See details in ?geiger::fitDiscrete()."))
      }
      ## Q_matrix
      if ("matrix" %in% evolutionary_models)
      {
        if (is.null(Q_matrix))
        {
         stop(paste0("For 'evolutionary_models = matrix', you must provide a 'Q_matrix' defining rate transition classes.\n",
              "See details in ?geiger::fitDiscrete()."))
        }
      }
    }

    ## res
    # res must be a positive integer.
    if ((res != abs(res)) | (res != round(res)))
    {
      stop(paste0("'res' must be a positive integer defining the number of time steps used to interpolate trait value in the densityMaps."))
    }

    ## nb_simulations
    # Check that nb_simulations is a positive integer.
    if ((nb_simulations != abs(nb_simulations)) | (nb_simulations != round(nb_simulations)))
    {
      stop(paste0("'nb_simulations' must be a positive integer defining the number of simulations generated for stochastic mapping."))
    }
    #  Send warnings if higher than 10000 and lower than 100
    if ((nb_simulations >= 10000))
    {
      cat(paste0("WARNING: 'nb_simulations' is set to ",nb_simulations,". High number of simulations may be time-consuming and only improve marginally the robustness of the tests.\n"))
    }
    if ((nb_simulations <= 100))
    {
      cat(paste0("WARNING: 'nb_simulations' is set to ",nb_simulations,". Low number of simulations may provide biased estimates of states/ranges and affect test outputs.\n"))
    }

    ## colors_per_states
    if (!is.null(colors_per_states))
    {
      # Check that the color scale match the states
      states_list <- levels(as.factor(tip_data))
      if (!all(states_list %in% names(colors_per_states)))
      {
        missing_states <- states_list[!(states_list %in% names(colors_per_states))]
        stop(paste0("Not all states are found in 'colors_per_states'.\n",
                    "Missing: ", paste(missing_states, collapse = ", "), "."))
      }
      # Check whether all colors are valid
      if (!all(is_color(colors_per_states)))
      {
        invalid_colors <- colors_per_states[!is_color(colors_per_states)]
        stop(paste0("Some color names in 'colors_per_states' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
  }

  # Set default model (ARD) if absent
  if (is.null(evolutionary_models)) { evolutionary_models <- "ARD" }

  # Get number of tips
  nb_tips <- length(phylo$tip.label)
  # Get number of nodes
  nb_nodes <- nb_tips + phylo$Nnode

  ##	Run all models and store their results in a list
  # ?geiger::fitDiscrete

  nb_models <- length(evolutionary_models)
  if (verbose) { cat(paste0(Sys.time(), " - Fit ",nb_models," evolutionary model(s): ", paste(evolutionary_models, collapse = ", "), ".\n\n")) }

  # Initiate list to store models outputs
  models_fits <- list()

  # Fit ARD model
  if (any(evolutionary_models == "ARD"))
  {
    ARD_ID <- which(evolutionary_models == "ARD")
    ARD_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "ARD", ...)

    if (verbose) { cat("------ ARD model ------ \n\n") ; print(ARD_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(ARD_fit$opt)
    # q12 = transition rate from state 1 to state 2
    # q21 = transition rate from state 2 to state 1
    # q13 = transition rate from state 1 to state 3
    # All transition rates can differ

    # Store output
    models_fits[[ARD_ID]] <- ARD_fit
  }

  # Fit ER model
  if (any(evolutionary_models == "ER"))
  {
    ER_ID <- which(evolutionary_models == "ER")
    ER_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "ER", ...)

    if (verbose) { cat("------ ER model ------ \n\n") ; print(ER_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(ER_fit$opt)
    # q12 = q21 = q13 = q31 = q23 = q32
    # Transition rates between states are all equals

    # Store output
    models_fits[[ER_ID]] <- ER_fit
  }

  # Fit SYM model
  if (any(evolutionary_models == "SYM"))
  {
    SYM_ID <- which(evolutionary_models == "SYM")
    SYM_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "SYM", ...)

    if (verbose) { cat("------ SYM model ------ \n\n") ; print(SYM_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(SYM_fit$opt)
    # q12 = q21
    # Transition rates from state 1 to state 2 are equal in both direction
    # But transition between different states can differ: q12 ≠ q13

    # Store output
    models_fits[[SYM_ID]] <- SYM_fit
  }

  # Fit 'meristic' model = step-wise transitions => 1 <-> 2 <-> 3
  if (any(evolutionary_models == "meristic"))
  {
    meristic_ID <- which(evolutionary_models == "meristic")
    meristic_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = "meristic", ...)
    # Can define if the Q-matrix is symmetrical or not with symmetric = TRUE/FALSE. Default is TRUE.

    if (verbose) { cat("------ Meristic model ------ \n\n") ; print(meristic_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(meristic_fit$opt)
    # q13 = q31 = 0
    # No transition possible between non-sequential states

    # Store output
    models_fits[[meristic_ID]] <- meristic_fit
  }

  # Fit 'matrix' model = User-defined Q-matrix defining state transition classes
  if (any(evolutionary_models == "matrix"))
  {
    matrix_ID <- which(evolutionary_models == "matrix")
    matrix_fit <- geiger::fitDiscrete(phy = phylo, dat = tip_data, model = Q_matrix, ...)
    # Can define if the Q-matrix is symmetrical or not with symmetric = TRUE/FALSE. Default is TRUE.

    if (verbose) { cat("------ User-defined matrix model ------ \n\n") ; print(matrix_fit) ; cat("\n") }
    # fitted Q matrix = transition parameters defining instantaneous rates of transitions between states in nb of events / time

    # print(matrix_fit$opt)

    # Store output
    models_fits[[matrix_ID]] <- matrix_fit
  }

  # Name model fits according to the associated models
  names(models_fits) <- evolutionary_models

  ## Compare model fits and select best model

  # Compare model fits with AICc and Akaike's weights
  models_comparison <- select_best_trait_model_from_geiger(list_model_fits = models_fits)

  # Extract best model
  best_model_name <- models_comparison$best_model_name
  best_model_fit <- models_comparison$best_model_fit

  # Display result of model comparison
  if (verbose)
  {
    cat(paste0(Sys.time(), " - Compare model fits.\n\n"))
    print(models_comparison$models_comparison_df)
  }

  ### Run stochastic mapping = obtain simmaps

  # Display steps
  if (verbose)
  {
    cat(paste0("\n", Sys.time(), " - Run simulations for stochastic mapping.\n\n"))
  }

  # Extract Q-matrix from best model
  best_Q_matrix <- phytools::as.Qmatrix(best_model_fit)

  simmaps <- phytools::make.simmap(tree = phylo, x = tip_data, model = best_model_name,
                                   nsim = nb_simulations, # Run at least 100 simulations because you want to observe the mean trend
                                   Q = best_Q_matrix,
                                   ...) # Can provide additional arguments for root states

  ### Get ACE at nodes based on posterior sampling of the simmaps (if return_ace)

  # Display steps
  if (verbose)
  {
    cat(paste0("\n", Sys.time(), " - Extract ACE as posterior sampling from stochastic mapping.\n"))
  }

  # Use phytools summary function
  simmaps_summary_obj <- phytools::describe.simmap(tree = simmaps, plot = FALSE)
  ace_matrix <- simmaps_summary_obj$ace

  ### Create densityMaps from simmaps

  if (verbose) { cat(paste0("\n", Sys.time(), " - Create densityMaps by summarizing simulations of evolutionary history (simmaps).\n\n")) }

  # Works only for binary traits.
  # Need to binarize every state as 0/1 with 0 = Not the focal state; 1 = Focal state.

  ## If not provided, define a color per state
  if (is.null(colors_per_states))
  {
    colors_per_states <- grDevices::rainbow(n = length(states_list))
    names(colors_per_states) <- states_list
  } else {
    # If provided, ensure it is properly ordered
    colors_per_states <- colors_per_states[match(x = states_list, table = names(colors_per_states))]
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
    focal_color <- colors_per_states[i]
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

  ## Plot densityMaps (if plot_map)
  if (plot_map)
  {
    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    if (!plot_overlay)
    {
      ## Plot one densityMap per state
      for (i in seq_along(densityMaps_all_states))
      {
        plot(densityMaps_all_states[[i]])
      }
    } else {

      ## Plot the overlay of densityMaps with alpha
      plot_densityMaps_overlay(densityMaps = densityMaps_all_states,
                               add_ACE_pies = TRUE,
                               ace = ace_matrix[1:(nrow(ace_matrix)-nb_tips), ])
    }

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Export densityMap(s) in PDF if requested
  if (!is.null(PDF_file_path))
  {
    # Allows plotting outside of figure range
    xpd_init <- par()$xpd
    par(xpd = TRUE)

    if (!plot_overlay)
    {
      ## Plot one densityMap per state
      grDevices::pdf(file = file.path(PDF_file_path),
                     width = nb_tips/60*8, height = nb_tips/60*10)

      for (i in seq_along(densityMaps_all_states))
      {
        plot(densityMaps_all_states[[i]])
      }

      grDevices::dev.off()

    } else {

      ## Plot the overlay of densityMaps with alpha
      grDevices::pdf(file = file.path(PDF_file_path),
                     width = nb_tips/60*8, height = nb_tips/60*10)

      plot_densityMaps_overlay(densityMaps = densityMaps_all_states,
                               add_ACE_pies = TRUE,
                               ace = ace_matrix[1:(nrow(ace_matrix)-nb_tips), ])

      grDevices::dev.off()
    }

    # Reset $xpd to initial values
    par(xpd = xpd_init)
  }

  ## Build output
  output <- list(densityMaps = densityMaps_all_states,
                 trait_data_type = "categorical")
  # Include simmaps if requested
  if(return_simmaps) { output$simmaps <- simmaps } # Filter to include only internal nodes (to be consistent with continuous trait)
  # Include ACE if requested
  if(return_ace) { output$ace <- ace_matrix[1:(nrow(ace_matrix)-nb_tips), ] } # Filter to include only internal nodes (to be consistent with continuous trait)
  # Include output of best model if requested
  if(return_best_model_fit) { output$best_model_fit <- best_model_fit }
  # Include df for model comparison if requested
  if(return_model_selection_df) { output$model_selection_df <- models_comparison$models_comparison_df }

  ## Return output
  return(invisible(output))
}


### Sub-function to handle biogeographic data ####

prepare_trait_data_for_biogeographic_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "DEC+J", # Default = "DEC+J" for biogeographic
    ..., # Additional arguments for BioGeoBEARS functions
    res = 100,
    nb_simulations = 1000, # Only for categorical and biogeographic data
    plot_map = TRUE,
    plot_overlay = TRUE, # Only for categorical and biogeographic data
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_simmaps = TRUE, # Only for categorical and biogeographic data
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE)
{
  ### Check input validity
  {
    ## evolutionary_models
    if (!is.null(evolutionary_models))
    {
      if (!all(evolutionary_models %in% c("BAYAREALIKE", "DIVALIKE", "DEC", "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J")))
      {
        stop(paste0("For 'trait_data_type = biogeographic', 'evolutionary_models' must be selected among: 'BAYAREALIKE', 'DIVALIKE', 'DEC', 'BAYAREALIKE+J', 'DIVALIKE+J', 'DEC+J'.\n",
                    "See details in documentation from R package `BioGeoBEARS`."))
      }
    }

    ## res
    # res must be a positive integer.
    if ((res != abs(res)) | (res != round(res)))
    {
      stop(paste0("'res' must be a positive integer defining the number of time steps used to interpolate trait value in the densityMaps."))
    }

    ## nb_simulations
    # Check that nb_simulations is a positive integer.
    if ((nb_simulations != abs(nb_simulations)) | (nb_simulations != round(nb_simulations)))
    {
      stop(paste0("'nb_simulations' must be a positive integer defining the number of simulations generated for stochastic mapping."))
    }
    #  Send warnings if higher than 10000 and lower than 100
    if ((nb_simulations >= 10000))
    {
      warning(paste0("'nb_simulations' is set to ",nb_simulations,". High number of simulations may be time-conusming and only improve marginally the robustness of the tests."))
    }
    if ((nb_simulations <= 100))
    {
      warning(paste0("'nb_simulations' is set to ",nb_simulations,". Low number of simulations may provide bias estimates of states/ranges and affect test outputs."))
    }
  }

  # If is.null(evolutionary_models), use "DEC+J" as default
  # If NULL, send a message to inform of the model used as default

  #	Run all models and store their results in a list

  #	select_best_trait_model() # Use it to compare model fits with AICc applying on the list of model results
  # Generate a df for model comparison

  # Get ACE (if return_ace)

  # Run stochastic mapping = obtain BSM outputs from BioGeoBEARS

  # Convert BSM outputs to simmaps
  # convert_BSM_to_simmaps()

  # Create densityMaps from simmaps

  # Plot densityMaps (if plot_map)
  # If plot_overlay, plot the overlay of densityMaps with alpha
  # If !plot_overlay,  plot one densityMap per state

  # Export densityMaps in PDF (if !is.null(PDF_file_path))

  # Build output
  # Include contMap by default
  # Include trait_data_type <- "biogeographic" by default
  # Include ACE if return_ace
  # Include simmaps if return_simmaps
  # Include output of best model if best_model_fit
  # Include df for model comparison if return_model_selection_df

  # Return output
}

### Helper function to compare model fits from geiger::fitContinuous or geiger::fitDiscrete with AICc and Akaike's weights ####
# Input = list of models fit with geiger::fitContinuous or geiger::fitDiscrete

#' @title Compare model fits with AICc and Akaike's weights
#'
#' @description Compare models fit with [geiger::fitContinuous()] or [geiger::fitDiscrete()] using AICc and Akaike's weights.
#'   Generate a data.frame summarizing information. Identify the best model and extract its results.
#'
#' @param list_model_fits Named list with the results of a model fit with [geiger::fitContinuous()] or [geiger::fitDiscrete()] in each element.
#'
#' @return The function returns a list with three elements.
#' * `$model_comparison_df` Data.frame summarizing information to compare model fits. It includes the model name (`$model`),
#'   the log-likelihood (`$logLik`), the number of free-parameters (`$k`), the corrected AIC (`$AICc`),
#'   the Akaike weights (`$Akaike_weights`), and their rank based on AICc (`$rank`).
#' * `$best_model_name` Character string. Name of the best model.
#' * `$best_model_fit` List containing the output of [geiger::fitContinuous()] or [geiger::fitDiscrete()] for the model with the best fit.
#'
#' @author Maël Doré
#'
#' @seealso [geiger::fitContinuous()] [geiger::fitDiscrete()]
#'
#' @noRd
#'

select_best_trait_model_from_geiger <- function (list_model_fits)
{
  # Homemade function to extract the number of free-parameters from geiger model outputs
  extract_k <- function (x)
  {
    k <- x$opt$k
    return(k)
  }

  # Homemade function to extract AICc from geiger model outputs
  extract_AICc <- function (x)
  {
    AICc <- x$opt$aicc
    return(AICc)
  }

  # Extract ln-likelihood, number of parameters, and AICc from models' list
  models_comparison_df <- data.frame(model = names(list_model_fits),
                                     logL = sapply(X = list_model_fits, FUN = stats::logLik),
                                     k = sapply(X = list_model_fits, FUN = extract_k),
                                     AICc = sapply(X = list_model_fits, FUN = extract_AICc))

  # Compute Akaike's weights based on AICc
  models_comparison_df$Akaike_weights <- round(phytools::aic.w(models_comparison_df$AICc), 3) * 100

  # Compute ranks based on AICc (The lowest = the best)
  models_comparison_df$rank <- rank(models_comparison_df$AICc, na.last = "keep", ties.method = "first")

  # Extract name of best model
  best_model_ID <- which(models_comparison_df$rank == 1)
  best_model_name <- names(list_model_fits)[best_model_ID]

  # Extract output of best model
  best_model_fit <- list_model_fits[[best_model_ID]]

  # Export output without printing best model fit
  return(invisible(list(models_comparison_df = models_comparison_df,
                        best_model_name = best_model_name,
                        best_model_fit = best_model_fit)))
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
#' @param tol Positive numerical. To set the tolerance used to match node ages and time steps (i.e., onsider them equal). Default = 1e-5.
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
    cols <- grDevices::rainbow(1001, start = 0.7, end = 0)
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

### Utility function from phytools used to produce $mapped.edge in densityMap ####
# Original function written by Liam Revell, 2012

makeMappedEdge <- function (edge, maps)
{
  st <- sort(unique(unlist(sapply(maps, function(x) names(x)))))
  mapped.edge <- matrix(0, nrow(edge), length(st))
  rownames(mapped.edge) <- apply(edge, 1, function(x) paste(x, collapse = ","))
  colnames(mapped.edge) <- st

  for (i in 1:length(maps))
  {
    for (j in 1:length(maps[[i]]))
    {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }

  return(mapped.edge)
}

### Utility function to check that an object is a valid color according to grDevices ####

is_color <- function(x, null_ok = FALSE)
{
  if (is.null(x) && null_ok)
  {
    return(TRUE)
  }
  vapply(x, function(i)
  {
    tryCatch({
      is.matrix(grDevices::col2rgb(i))
    }, error = function(e) FALSE)
  }, TRUE)
}



# ### Utility functions for tree transformation ####
#
# Used to copy the S3 method rescale.phylo() from geiger
#
# # Wrapper function for tree transformation from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# rescale.phylo <- function(phy,
#                           model = c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"),
#                           ...)
# {
#
#   model = match.arg(model, c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda",
#                              "kappa", "delta", "white", "depth"))
#
#   if (!("phylo" %in% class(phy))) stop("supply 'phy' as a 'phylo' object")
#
#   FUN = switch(model,
#                BM = .bm.phylo(phy),
#                OU = .ou.phylo(phy),
#                EB = .eb.phylo(phy),
#                nrate = .nrate.phylo(phy),
#                lrate = .lrate.phylo(phy),
#                trend = .trend.phylo(phy),
#                lambda = .lambda.phylo(phy),
#                kappa = .kappa.phylo(phy),
#                delta = .delta.phylo(phy),
#                white = .white.phylo(phy),
#                depth = .depth.phylo(phy)
#   )
#   class(FUN) = c("transformer", "function")
#   dots = list(...)
#   if(!missing(...))
#   {
#     if(!all(names(dots) %in% argn(FUN)))
#       stop(paste("The following parameters are expected:\n\t", paste(argn(FUN),
#                                                                      collapse = "\n\t", sep = ""), sep = ""))
#     return(FUN(...))
#   } else {
#     return(FUN)
#   }
# }
#
# # Tree transformation for BM from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .bm.phylo = function(phy)
# {
#   el = phy$edge.length
#   z = function(sigsq)
#   {
#     phy$edge.length = el * sigsq
#     phy
#   }
#   attr(z, "argn") = "sigsq"
#   return(z)
# }
#
# # Tree transformation for OU from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .ou.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t1 = Tmax - ht$end[phy$edge[mm,1]]
#   ht$t2 = ht$start - ht$end + ht$t1
#   z = function(alpha, sigsq=1)
#   {
#     if(alpha<0) stop("'alpha' must be positive valued")
#     bl = (1/(2*alpha)) * exp(-2*alpha * (Tmax-ht$t2)) * (1 - exp(-2 * alpha * ht$t2)) - (1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t1)) * (1 - exp(-2 * alpha * ht$t1))
#     phy$edge.length = bl[phy$edge[,2]]
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("alpha", "sigsq")
#   return(z)
# }
#
# # Tree transformation for EB from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .eb.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t1 = Tmax - ht$end[phy$edge[mm,1]]
#   ht$t2 = ht$start - ht$end + ht$t1
#
#   z = function(a, sigsq = 1)
#   {
#     if(a == 0) return(phy)
#     bl = (exp(a*ht$t2) - exp(a*ht$t1)) / (a)
#     phy$edge.length = bl[phy$edge[,2]]
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("a", "sigsq")
#   return(z)
# }
#
# # Tree transformation for multi-regime BM with rates split by time from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .nrate.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t = Tmax - ht$end
#   ht$e = ht$start - ht$end
#   ht$a = ht$t - ht$e
#   ht$rS = ht$a/Tmax
#   ht$rE = ht$t/Tmax
#
#   dd = phy$edge[,2]
#
#   relscale.brlen = function(start, end, len, dat)
#   {
#     ss = start <= dat[,"time"]
#     strt = min(which(ss))
#
#     ee = dat[,"time"] < end
#     etrt = max(which(ee)) + 1
#
#     bl = numeric()
#     fragment = numeric()
#     marker = start
#     for(i in strt:etrt)
#     {
#       fragment = c(fragment, (nm <- (min(c(end, dat[i, "time"])))) - marker)
#       bl = c(bl,dat[i, "rate"])
#       marker = nm
#     }
#     fragment = fragment/(sum(fragment))
#     sclbrlen = numeric()
#     for (i in 1:length(bl)) sclbrlen = c(sclbrlen, len * fragment[i] * bl[i])
#     sc = structure(as.numeric(sclbrlen), names = strt:etrt)
#     return(sc)
#   }
#
#
#   z = function(time, rate, sigsq=1)
#   {
#     if(any(time > 1) | any(time < 0)) stop("supply 'time' as a vector of relative time:\n\tvalues should be in the range 0 (root) to 1 (present)")
#     if(any(rate < 0)) stop("'rate' must consist of positive values")
#     if(length(time) != length(rate)) stop("'time' and 'rate' must be of equal length")
#     phy$edge.length = phy$edge.length * sigsq
#     ordx = order(time)
#     time = time[ordx]
#     rate = rate[ordx]
#     dat = cbind(time=c(0,time, 1), rate = (c(1, 1, rate)))
#     rs = sapply(dd, function(x) as.numeric(sum(relscale.brlen(ht$rS[x], ht$rE[x], ht$e[x], dat))))
#     phy$edge.length = rs
#     phy
#   }
#   attr(z, "argn") = c("time", "rate", "sigsq")
#   return(z)
# }
#
# # Tree transformation for multi-regime BM with rates split by clades from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .lrate.phylo = function(phy)
# {
#   N = Ntip(phy)
#   n = Nnode(phy)
#   cache = .cache.tree(phy)
#   cache$phy = phy
#   vv = rep(1, N+n-1)
#
#   z = function(node, rate, sigsq = 1)
#   {
#     shifts = c(sort(node[node > N]), node[node <= N])
#     mm = match(shifts, node)
#     rates = rate[mm]
#
#     phy$edge.length = phy$edge.length * sigsq
#     for(i in 1:length(shifts)) vv = .assigndescendants(vv, shifts[i], rates[i], exclude = shifts, cache = cache)
#     phy$edge.length = vv * phy$edge.length
#     phy
#   }
#
#   attr(z, "argn") = c("node", "rate", "sigsq")
#   return(z)
# }
#
# # Tree transformation for "rate_trend" model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .trend.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$head = Tmax - ht$end[phy$edge[mm,1]] # age
#   ht$tail = ht$head + (ht$start - ht$end)
#
#   z = function(slope, sigsq = 1)
#   {
#     # begin (age): head
#     # end: tail
#     ht$br = 1 + ht$head * slope
#     ht$er = 1 + ht$tail * slope
#     scl = sapply(1:nrow(ht),
#                  function(idx)
#                  {
#                    if(idx == N+1) return(NA)
#                    if(ht$br[idx] > 0 & ht$er[idx] > 0)
#                    {
#                      return((ht$br[idx]+ht$er[idx])/2)
#                    } else if (ht$br[idx] < 0 & ht$er[idx] < 0) {
#                      return(0)
#                    } else {
#                      si = -1/slope
#                      return(ht$br[idx] * (si-ht$head[idx]) / (2 * (ht$tail[idx] - ht$head[idx])))
#                    }
#                  })
#
#     phy$edge.length = phy$edge.length * scl[phy$edge[,2]]
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("slope", "sigsq")
#   return(z)
# }
#
# # Tree transformation for Pagel's lambda model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .lambda.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:N, phy$edge[,2])
#   ht$e = ht$start - ht$end
#   paths = .paths.phylo(phy)
#
#   z = function(lambda, sigsq = 1)
#   {
#     if(lambda < 0) stop("'lambda' must be positive valued")
#
#     bl = phy$edge.length * lambda
#     bl[mm] = bl[mm] + (paths - (paths * lambda))
#     phy$edge.length = bl
#     if(any(phy$edge.length < 0))
#     {
#       warning("negative branch lengths encountered:\n\tlambda may be too large")
#     }
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("lambda", "sigsq")
#   return(z)
# }
#
# # Tree transformation for Pagel's kappa model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .kappa.phylo = function(phy)
# {
#   z = function(kappa, sigsq = 1)
#   {
#     if(kappa < 0) stop("'kappa' must be positive valued")
#
#     phy$edge.length = phy$edge.length^kappa
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("kappa", "sigsq")
#   return(z)
# }
#
# # Tree transformation for Pagel's delta model from R package geiger
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .delta.phylo = function(phy)
# {
#   ht = heights.phylo(phy)
#   N = Ntip(phy)
#   Tmax = ht$start[N+1]
#   mm = match(1:nrow(ht), phy$edge[,2])
#   ht$t = Tmax - ht$end
#   ht$e = ht$start - ht$end
#   ht$a = ht$t - ht$e
#   if(sum(ht$a < -0.1) > 0) stop("Calculation error; contact developers (LJ Harmon and JM Eastman).")
#   ht$a[ht$a < 0] <- 0
#
#   z = function(delta, sigsq = 1, rescale=TRUE)
#   {
#     if(delta < 0) stop("'delta' must be positive valued")
#     bl = (ht$a + ht$e)^delta - ht$a^delta
#     phy$edge.length = bl[phy$edge[,2]]
#
#     if(rescale)
#     {
#       scl = Tmax^delta
#       phy$edge.length = (phy$edge.length/scl) * Tmax
#     }
#     phy$edge.length = phy$edge.length * sigsq
#     phy
#   }
#   attr(z, "argn") = c("delta", "sigsq")
#   return(z)
# }
#
# # Tree transformation for a non-phylogenetic "white" model from R package geiger
# # that assumes there is no covariance structure among species
# # Internal branches are set to a length of zero
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .white.phylo = function(phy)
# {
#   N = Ntip(phy)
#   phy$edge.length[] = 0
#   phy$edge.length[phy$edge[,2] <= N] = 1
#   el = phy$edge.length
#   z = function(sigsq = 1)
#   {
#     phy$edge.length = el * sigsq
#     phy
#   }
#   attr(z, "argn") = "sigsq"
#   return(z)
# }
#
# # Tree transformation "depth" model from R package geiger
# # Adjusts the total depth of the tree
# # which in return increase/decrease rates of evolution
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .depth.phylo = function(phy)
# {
#   orig = max(heights.phylo(phy))
#   z = function(depth)
#   {
#     phy$edge.length <- (phy$edge.length/orig) * depth
#     if(!is.null(phy$root.edge)) phy$root.edge = (phy$root.edge/orig) * depth
#     phy
#   }
#   attr(z, "argn") = "depth"
#   z
# }
#
# # Internal function from R package geiger
# # Compute path length from root to tip
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# .paths.phylo = function(phy, ...)
# {
#
#   ## much from ape:::vcv.phylo()
#   phy <- reorder(phy, "postorder")
#
#   FUN = function(vcv = FALSE)
#   {
#     n <- length(phy$tip.label)
#     pp <- .cache.descendants(phy)$tips
#     e1 <- phy$edge[, 1]
#     e2 <- phy$edge[, 2]
#     EL <- phy$edge.length
#     xx <- numeric(n + phy$Nnode)
#     if(vcv) vmat = matrix(0, n, n)
#     for (i in length(e1):1)
#     {
#       var.cur.node <- xx[e1[i]]
#       xx[e2[i]] <- var.cur.node + EL[i]
#       if(vcv)
#       {
#         j <- i - 1L
#         while (e1[j] == e1[i] && j > 0)
#         {
#           left = pp[[e2[j]]]
#           right = pp[[e2[i]]]
#           vmat[left, right] <- vmat[right, left] <- var.cur.node
#           j <- j - 1L
#         }
#       }
#     }
#     if(vcv)
#     {
#       diags <- 1 + 0:(n - 1) * (n + 1)
#       vmat[diags] <- xx[1:n]
#       colnames(vmat) <- rownames(vmat) <- phy$tip.label
#       return(vmat)
#     } else {
#       return(xx[1:n])
#     }
#   }
#
#   FUN(...)
# }
#
# ## Internal function from R package geiger
# # Associated C scripts in /src/utilities.cpp
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
#
# # Declare the use of compiled C code in the package
# # Add an import in NAMESPACE
# # Need to add things listed in Rcpp::Rcpp.package.skeleton() to get the necessary changes for the package to include Rcpp
# # Need to add CallEntries in init.c such as R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
# #' @useDynLib deepSTRAPP
# #' @importFrom Rcpp evalCpp
#
# .cache.descendants = function(phy)
# {
#   # Fetches all tips subtended by each internal node
#
#   N = as.integer(Ntip(phy))
#   n = as.integer(Nnode(phy))
#
#   phy = reorder(phy, "postorder")
#
#   zz = list(N=N,
#             MAXNODE=N+n,
#             ANC=as.integer(phy$edge[,1]),
#             DES=as.integer(phy$edge[,2])
#   )
#
#   res = .Call("cache_descendants", phy = zz, PACKAGE = "geiger")
#   return(res)
# }
#
#
# # Internal function from R package geiger
# # Compute heights of nodes in phylo
# # Source: geiger/R/utilities-phylo.R
# # Authors: LJ Harmon and JM Eastman
# heights.phylo = function(x)
# {
#   phy = x
#   phy <- reorder(phy, "postorder")
#   n <- length(phy$tip.label)
#   n.node <- phy$Nnode
#   xx <- numeric(n + n.node) # ending times
#   for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
#   root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
#   labs = c(phy$tip.label, phy$node.label)
#   depth = max(xx)
#   tt = depth - xx # time to 'present day' of branch starts
#   idx = 1:length(tt)
#   #dd = phy$edge.length[idx]
#   mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
#   dd = c(phy$edge.length, root)[mm] # reordered bls
#   ss = tt + dd
#   res = cbind(ss, tt)
#   rownames(res) = idx
#   colnames(res) = c("start", "end")
#   res = data.frame(res)
#   res
# }

### Dummy function to force the use of Rcpp ####
#' @useDynLib deepSTRAPP
#' @importFrom Rcpp evalCpp
dummy_Rcpp <- function ()
{
  Rcpp::evalCpp(code = "__cplusplus" )
}
