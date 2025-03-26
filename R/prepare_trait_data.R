## Functions to prepare trait data by mapping their evolution on the phylogeny
# One master function to select the proper pipeline according to data type
# Three sub-functions mapping trait evolution according to data type

#' @title Map trait evolution on a time-calibrated phylogeny
#'


## Currently, For continuous traits
# Input = contMap
# Output = contMap + plot of contMap + parameters of best model of continuous traits (BM, EB/ACDC, OU1, Pagel's lambda ?)

## Make a different function for each type of data, then make a wrapper function for all types
# prepare_trait_data() = wrapper
# Detect type of data based on $trait_data_type, but check it match the input (tip_data)
# prepare_trait_data_for_continuous_data() = For continuous traits
# prepare_trait_data_for_categorical_data() = For categorical traits
# prepare_trait_data_for_biogeographic_data() = For biogeographic traits

# Add an option to get ACE for categorical and biogeographic traits to (as scaled marginal likelihoods)

## Explain that this only compare and fit 'classic' simple models using geiger::fitContinuous in a ML approach
# Default for continuous = "BM". "All" = BM, EB/ACDC, OU, Pagel's lambda, kappa, delta, and rate_trend as available in geiger::fitContinuous
# May need to use more complex models to capture the dynamics of trait evolution. BayesTraits. BAMM. Multi-regimes models (BM and OU)
# Trait-dependent multi-rate models (use a simmap to assign branches to a specific regime)
# ?phytools::brownie.lite()
# ?OUwie::OUwie with model = "BMS" (Multi-BM) and "OUM", "OUMV", "OUA", "OUMVA" (Multi-OU with sigma, theta and alpha varying or not across regimes)
# BayesTraits and BAMM for a thorough exploration of location and number of regime shifts independent from any other trait/extrinsic variable
# RRphylo for a penalized phylogenetic ridge regression approach that allows regime shifts across all branches

## Note that interpolation of ancestral trait values along branches is not accurate for time-dependent models,
# as the interpolation used to built the contMap is assuming a constant rate along each branch.
# However, ancestral trait values at nodes are accurate.

# Default for categorical = "ARD". "All" = ER, SYM, ARD.
# May need to use more complex models to capture the dynamics of trait evolution. Custom constrains in the Q matrix, Multi-regimes models, ... (cite packages)

# Default for Biogeographic = "DEC+J". "All" = DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J.
# May need to use more complex models to capture the dynamics of trait evolution. +W, +X, DECX, ... (cite packages/refs)

## plot_overlay for categorical and biogeographic data
# Need plot_map = TRUE. If plot_overlay = TRUE: Plot the overlay of Density maps with alpha. If plot_overlay = FALSE, plot one densityMap per state/range

### Master function to prepare data and select the proper test function according to data type ####

prepare_trait_data <- function (
    tip_data,
    trait_data_type,
    phylo,
    evolutionary_models = NULL, # Default = "BM" for continuous data; "ARD" for categorical; "DEC+J" for biogeographic
    res = 100, # Number of time steps used to interpolate trait value in the contMap
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

  # tip_data must have names %in% phylo$tip.label
  # Check if tip_data need to be reorder to match phylo$tip.label. If so, do it, but send a warning

  # Check that trait_data_type is one of "continuous", "categorical", "biogeographic"
  # Check that trait_data_type is coherent with what is found in tip_data

  # Check that phylo is of class "phylo"
  # Check that the phylo is fully resolved (is it needed for those models?)
  # Check that the phylo is ultrametric (is it needed for those models?)

  # Check that models are compatible with trait_data_type. If not, send error showing which one is wrong what are the valid options for the given trait_data_type
  # No error if is.null() as the default will be used

  # Only for categorical and biogeographic data
  # Check that nb_simulations is a positive integer. Send warnings if higher than 10000 and lower than 100
    # Higher than 10000 => May be time-consuming
    # Lower than 100 => May be bias

  # If provided, PDF_file_path must end with ".pdf"

  ## Compute the appropriate internal function depending on the type of data

  switch(EXPR = trait_data_type,
         continuous =   { # Case for continuous data
           # Output = contMap. Trait values are interpolated along branches.
           trait_data_output <- prepare_trait_data_for_continuous_data(
             tip_data = tip_data,
             trait_data_type = trait_data_type,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "BM" for continuous data
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
             trait_data_type = trait_data_type,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "ARD" for categorical data
             res = res,
             nb_simulations = nb_simulations, # Only for categorical and biogeographic data
             plot_map = plot_map,
             plot_overlay = plot_overlay, # Only for categorical and biogeographic data
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
             trait_data_type = trait_data_type,
             phylo = phylo,
             evolutionary_models = evolutionary_models, # Default = "DEC+J" for biogeographic data
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


library(phytools)
data(eel.tree)
data(eel.data)

# Extract body size
eel_data <- setNames(eel.data$Max_TL_cm,
                     rownames(eel.data))

tip_data <- eel_data
phylo <- eel.tree
evolutionary_models <- c("EB", "rate_trend", "BM")

test <- prepare_trait_data_for_continuous_data(tip_data = eel_data,
                                               phylo = eel.tree,
                                               evolutionary_models = c("EB", "BM", "lambda"),
                                               plot_map = TRUE,
                                               return_best_model_fit = TRUE,
                                               return_model_selection_df = TRUE,
                                               verbose = FALSE)
str(test, 2)
test$best_model_fit$opt
test$ace
### Sub-function to handle continuous data ####

prepare_trait_data_for_continuous_data <- function (
    tip_data,
    trait_data_type,
    phylo,
    evolutionary_models = "BM", # Default = "BM" for continuous data
    res = 100,
    plot_map = TRUE,
    PDF_file_path = NULL,
    return_ace = TRUE,
    return_best_model_fit = FALSE,
    return_model_selection_df = FALSE,
    verbose = TRUE,
    ...)
{

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

  ### Get ACE (if return_ace)

  if (verbose) { cat(paste0("\n", Sys.time(), " - Infer Ancestral Character Estimates from the best fitting model: ",best_model_name,".\n")) }

  ## Rescale the tree such as the relative branch length account for trait rates in the best fitting model
  # ?geiger::rescale.phylo

  # Adjust model name
  if (best_model_name == "rate_trend") { best_model_name <- "trend" }
  # Create function to rescale phylo
  rescaled_phylo_fn <- rescale.phylo(phy = phylo, model = best_model_name)
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

    dev.off()

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
  # If is.null(evolutionary_models), use "ARD" as default
  # If NULL, send a message to inform of the model used as default

  #	Run all models and store their results in a list

  #	select_best_trait_model() # Use it to compare model fits with AICc applying on the list of model results
  # Generate a df for model comparison

  # Get ACE (if return_ace)

  # Run stochastic mapping = obtain simmaps

  # Create densityMaps from simmaps

  # Plot densityMaps (if plot_map)
    # If plot_overlay, plot the overlay of densityMaps with alpha
    # If !plot_overlay,  plot one densityMap per state

  # Export densityMaps in PDF (if !is.null(PDF_file_path))

  # Build output
  # Include contMap by default
  # Include trait_data_type <- "categorical" by default
  # Include ACE if return_ace
  # Include simmaps if return_simmaps
  # Include output of best model if best_model_fit
  # Include df for model comparison if return_model_selection_df

  # Return output
}


### Sub-function to handle biogeographic data ####

prepare_trait_data_for_biogeographic_data <- function (
    tip_data,
    phylo,
    evolutionary_models = "DEC+J", # Default = "DEC+J" for biogeographic
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

### Helper function to compare model fits from geiger::fitContinuous with AICc and Akaike's weights ####
# Input = list of models fit with geiger::fitContinuous

#' @title Compare model fits with AICc and Akaike's weights
#'
#' @description Compare models fit with [geiger::fitContinuous()] using AICc and Akaike's weights.
#'   Generate a data.frame summarizing information. Identify the best model and extract its results.
#'
#' @param list_model_fits Named list with the results of a model fit with [geiger::fitContinuous()] in each element.
#'
#' @return The function returns a list with three elements.
#' * `$model_comparison_df` Data.frame summarizing information to compare model fits. It includes the model name (`$model`),
#'   the log-likelihood (`$logLik`), the number of free-parameters (`$k`), the corrected AIC (`$AICc`),
#'   the Akaike weights (`$Akaike_weights`), and their rank based on AICc (`$rank`).
#' * `$best_model_name` Character string. Name of the best model.
#' * `$best_model_fit` List containing the output of [geiger::fitContinuous()] for the model with the best fit.
#'
#' @author Maël Doré
#'
#' @seealso [geiger::fitContinuous()]
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
                                     logL = sapply(X = list_model_fits, FUN = logLik),
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


### Utility functions for tree transformation ####

# Wrapper function for tree transformation from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
rescale.phylo <- function(phy,
                          model = c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda", "kappa", "delta", "white", "depth"),
                          ...)
{

  model = match.arg(model, c("BM", "OU", "EB", "nrate", "lrate", "trend", "lambda",
                             "kappa", "delta", "white", "depth"))

  if (!("phylo" %in% class(phy))) stop("supply 'phy' as a 'phylo' object")

  FUN = switch(model,
               BM = .bm.phylo(phy),
               OU = .ou.phylo(phy),
               EB = .eb.phylo(phy),
               nrate = .nrate.phylo(phy),
               lrate = .lrate.phylo(phy),
               trend = .trend.phylo(phy),
               lambda = .lambda.phylo(phy),
               kappa = .kappa.phylo(phy),
               delta = .delta.phylo(phy),
               white = .white.phylo(phy),
               depth = .depth.phylo(phy)
  )
  class(FUN) = c("transformer", "function")
  dots = list(...)
  if(!missing(...))
  {
    if(!all(names(dots) %in% argn(FUN)))
      stop(paste("The following parameters are expected:\n\t", paste(argn(FUN),
                                                                     collapse = "\n\t", sep = ""), sep = ""))
    return(FUN(...))
  } else {
    return(FUN)
  }
}

# Tree transformation for BM from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.bm.phylo = function(phy)
{
  el = phy$edge.length
  z = function(sigsq)
  {
    phy$edge.length = el * sigsq
    phy
  }
  attr(z, "argn") = "sigsq"
  return(z)
}

# Tree transformation for OU from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.ou.phylo = function(phy)
{
  ht = heights.phylo(phy)
  N = Ntip(phy)
  Tmax = ht$start[N+1]
  mm = match(1:nrow(ht), phy$edge[,2])
  ht$t1 = Tmax - ht$end[phy$edge[mm,1]]
  ht$t2 = ht$start - ht$end + ht$t1
  z = function(alpha, sigsq=1)
  {
    if(alpha<0) stop("'alpha' must be positive valued")
    bl = (1/(2*alpha)) * exp(-2*alpha * (Tmax-ht$t2)) * (1 - exp(-2 * alpha * ht$t2)) - (1/(2*alpha))*exp(-2*alpha * (Tmax-ht$t1)) * (1 - exp(-2 * alpha * ht$t1))
    phy$edge.length = bl[phy$edge[,2]]
    phy$edge.length = phy$edge.length * sigsq
    phy
  }
  attr(z, "argn") = c("alpha", "sigsq")
  return(z)
}

# Tree transformation for EB from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.eb.phylo = function(phy)
{
  ht = heights.phylo(phy)
  N = Ntip(phy)
  Tmax = ht$start[N+1]
  mm = match(1:nrow(ht), phy$edge[,2])
  ht$t1 = Tmax - ht$end[phy$edge[mm,1]]
  ht$t2 = ht$start - ht$end + ht$t1

  z = function(a, sigsq = 1)
  {
    if(a == 0) return(phy)
    bl = (exp(a*ht$t2) - exp(a*ht$t1)) / (a)
    phy$edge.length = bl[phy$edge[,2]]
    phy$edge.length = phy$edge.length * sigsq
    phy
  }
  attr(z, "argn") = c("a", "sigsq")
  return(z)
}

# Tree transformation for multi-regime BM with rates split by time from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.nrate.phylo = function(phy)
{
  ht = heights.phylo(phy)
  N = Ntip(phy)
  Tmax = ht$start[N+1]
  mm = match(1:nrow(ht), phy$edge[,2])
  ht$t = Tmax - ht$end
  ht$e = ht$start - ht$end
  ht$a = ht$t - ht$e
  ht$rS = ht$a/Tmax
  ht$rE = ht$t/Tmax

  dd = phy$edge[,2]

  relscale.brlen = function(start, end, len, dat)
  {
    ss = start <= dat[,"time"]
    strt = min(which(ss))

    ee = dat[,"time"] < end
    etrt = max(which(ee)) + 1

    bl = numeric()
    fragment = numeric()
    marker = start
    for(i in strt:etrt)
    {
      fragment = c(fragment, (nm <- (min(c(end, dat[i, "time"])))) - marker)
      bl = c(bl,dat[i, "rate"])
      marker = nm
    }
    fragment = fragment/(sum(fragment))
    sclbrlen = numeric()
    for (i in 1:length(bl)) sclbrlen = c(sclbrlen, len * fragment[i] * bl[i])
    sc = structure(as.numeric(sclbrlen), names = strt:etrt)
    return(sc)
  }


  z = function(time, rate, sigsq=1)
  {
    if(any(time > 1) | any(time < 0)) stop("supply 'time' as a vector of relative time:\n\tvalues should be in the range 0 (root) to 1 (present)")
    if(any(rate < 0)) stop("'rate' must consist of positive values")
    if(length(time) != length(rate)) stop("'time' and 'rate' must be of equal length")
    phy$edge.length = phy$edge.length * sigsq
    ordx = order(time)
    time = time[ordx]
    rate = rate[ordx]
    dat = cbind(time=c(0,time, 1), rate = (c(1, 1, rate)))
    rs = sapply(dd, function(x) as.numeric(sum(relscale.brlen(ht$rS[x], ht$rE[x], ht$e[x], dat))))
    phy$edge.length = rs
    phy
  }
  attr(z, "argn") = c("time", "rate", "sigsq")
  return(z)
}

# Tree transformation for multi-regime BM with rates split by clades from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.lrate.phylo = function(phy)
{
  N = Ntip(phy)
  n = Nnode(phy)
  cache = .cache.tree(phy)
  cache$phy = phy
  vv = rep(1, N+n-1)

  z = function(node, rate, sigsq = 1)
  {
    shifts = c(sort(node[node > N]), node[node <= N])
    mm = match(shifts, node)
    rates = rate[mm]

    phy$edge.length = phy$edge.length * sigsq
    for(i in 1:length(shifts)) vv = .assigndescendants(vv, shifts[i], rates[i], exclude = shifts, cache = cache)
    phy$edge.length = vv * phy$edge.length
    phy
  }

  attr(z, "argn") = c("node", "rate", "sigsq")
  return(z)
}

# Tree transformation for "rate_trend" model from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.trend.phylo = function(phy)
{
  ht = heights.phylo(phy)
  N = Ntip(phy)
  Tmax = ht$start[N+1]
  mm = match(1:nrow(ht), phy$edge[,2])
  ht$head = Tmax - ht$end[phy$edge[mm,1]] # age
  ht$tail = ht$head + (ht$start - ht$end)

  z = function(slope, sigsq = 1)
  {
    # begin (age): head
    # end: tail
    ht$br = 1 + ht$head * slope
    ht$er = 1 + ht$tail * slope
    scl = sapply(1:nrow(ht),
                 function(idx)
                 {
                   if(idx == N+1) return(NA)
                   if(ht$br[idx] > 0 & ht$er[idx] > 0)
                   {
                     return((ht$br[idx]+ht$er[idx])/2)
                   } else if (ht$br[idx] < 0 & ht$er[idx] < 0) {
                     return(0)
                   } else {
                     si = -1/slope
                     return(ht$br[idx] * (si-ht$head[idx]) / (2 * (ht$tail[idx] - ht$head[idx])))
                   }
                 })

    phy$edge.length = phy$edge.length * scl[phy$edge[,2]]
    phy$edge.length = phy$edge.length * sigsq
    phy
  }
  attr(z, "argn") = c("slope", "sigsq")
  return(z)
}

# Tree transformation for Pagel's lambda model from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.lambda.phylo = function(phy)
{
  ht = heights.phylo(phy)
  N = Ntip(phy)
  Tmax = ht$start[N+1]
  mm = match(1:N, phy$edge[,2])
  ht$e = ht$start - ht$end
  paths = .paths.phylo(phy)

  z = function(lambda, sigsq = 1)
  {
    if(lambda < 0) stop("'lambda' must be positive valued")

    bl = phy$edge.length * lambda
    bl[mm] = bl[mm] + (paths - (paths * lambda))
    phy$edge.length = bl
    if(any(phy$edge.length < 0))
    {
      warning("negative branch lengths encountered:\n\tlambda may be too large")
    }
    phy$edge.length = phy$edge.length * sigsq
    phy
  }
  attr(z, "argn") = c("lambda", "sigsq")
  return(z)
}

# Tree transformation for Pagel's kappa model from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.kappa.phylo = function(phy)
{
  z = function(kappa, sigsq = 1)
  {
    if(kappa < 0) stop("'kappa' must be positive valued")

    phy$edge.length = phy$edge.length^kappa
    phy$edge.length = phy$edge.length * sigsq
    phy
  }
  attr(z, "argn") = c("kappa", "sigsq")
  return(z)
}

# Tree transformation for Pagel's delta model from R package geiger
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.delta.phylo = function(phy)
{
  ht = heights.phylo(phy)
  N = Ntip(phy)
  Tmax = ht$start[N+1]
  mm = match(1:nrow(ht), phy$edge[,2])
  ht$t = Tmax - ht$end
  ht$e = ht$start - ht$end
  ht$a = ht$t - ht$e
  if(sum(ht$a < -0.1) > 0) stop("Calculation error; contact developers (LJ Harmon and JM Eastman).")
  ht$a[ht$a < 0] <- 0

  z = function(delta, sigsq = 1, rescale=TRUE)
  {
    if(delta < 0) stop("'delta' must be positive valued")
    bl = (ht$a + ht$e)^delta - ht$a^delta
    phy$edge.length = bl[phy$edge[,2]]

    if(rescale)
    {
      scl = Tmax^delta
      phy$edge.length = (phy$edge.length/scl) * Tmax
    }
    phy$edge.length = phy$edge.length * sigsq
    phy
  }
  attr(z, "argn") = c("delta", "sigsq")
  return(z)
}

# Tree transformation for a non-phylogenetic "white" model from R package geiger
# that assumes there is no covariance structure among species
# Internal branches are set to a length of zero
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.white.phylo = function(phy)
{
  N = Ntip(phy)
  phy$edge.length[] = 0
  phy$edge.length[phy$edge[,2] <= N] = 1
  el = phy$edge.length
  z = function(sigsq = 1)
  {
    phy$edge.length = el * sigsq
    phy
  }
  attr(z, "argn") = "sigsq"
  return(z)
}

# Tree transformation "depth" model from R package geiger
# Adjusts the total depth of the tree
# which in return increase/decrease rates of evolution
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.depth.phylo = function(phy)
{
  orig = max(heights.phylo(phy))
  z = function(depth)
  {
    phy$edge.length <- (phy$edge.length/orig) * depth
    if(!is.null(phy$root.edge)) phy$root.edge = (phy$root.edge/orig) * depth
    phy
  }
  attr(z, "argn") = "depth"
  z
}

# Internal function from R package geiger
# Compute path length from root to tip
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
.paths.phylo = function(phy, ...)
{

  ## much from ape:::vcv.phylo()
  phy <- reorder(phy, "postorder")

  FUN = function(vcv = FALSE)
  {
    n <- length(phy$tip.label)
    pp <- .cache.descendants(phy)$tips
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length
    xx <- numeric(n + phy$Nnode)
    if(vcv) vmat = matrix(0, n, n)
    for (i in length(e1):1)
    {
      var.cur.node <- xx[e1[i]]
      xx[e2[i]] <- var.cur.node + EL[i]
      if(vcv)
      {
        j <- i - 1L
        while (e1[j] == e1[i] && j > 0)
        {
          left = pp[[e2[j]]]
          right = pp[[e2[i]]]
          vmat[left, right] <- vmat[right, left] <- var.cur.node
          j <- j - 1L
        }
      }
    }
    if(vcv)
    {
      diags <- 1 + 0:(n - 1) * (n + 1)
      vmat[diags] <- xx[1:n]
      colnames(vmat) <- rownames(vmat) <- phy$tip.label
      return(vmat)
    } else {
      return(xx[1:n])
    }
  }

  FUN(...)
}

## Internal function from R package geiger
# Associated C scripts in /src/utilities.cpp
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman

# Declare the use of compiled C code in the package
# Add an import in NAMESPACE
# Need to use Rcpp::Rcpp.package.skeleton() to add the necessary changes to the package to include Rcpp
# Need to add CallEntries in init.c such as R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
#' @useDynLib deepSTRAPP

.cache.descendants = function(phy)
{
  # Fetches all tips subtended by each internal node

  N = as.integer(Ntip(phy))
  n = as.integer(Nnode(phy))

  phy = reorder(phy, "postorder")

  zz = list(N=N,
            MAXNODE=N+n,
            ANC=as.integer(phy$edge[,1]),
            DES=as.integer(phy$edge[,2])
  )

  res = .Call("cache_descendants", phy = zz, PACKAGE = "geiger")
  return(res)
}


# Internal function from R package geiger
# Compute heights of nodes in phylo
# Source: geiger/R/utilities-phylo.R
# Authors: LJ Harmon and JM Eastman
heights.phylo = function(x)
{
  phy = x
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)
  n.node <- phy$Nnode
  xx <- numeric(n + n.node) # ending times
  for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
  root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
  labs = c(phy$tip.label, phy$node.label)
  depth = max(xx)
  tt = depth - xx # time to 'present day' of branch starts
  idx = 1:length(tt)
  #dd = phy$edge.length[idx]
  mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
  dd = c(phy$edge.length, root)[mm] # reordered bls
  ss = tt + dd
  res = cbind(ss, tt)
  rownames(res) = idx
  colnames(res) = c("start", "end")
  res = data.frame(res)
  res
}
