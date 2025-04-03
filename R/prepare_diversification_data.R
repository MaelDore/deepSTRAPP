


# Run a BAMM (Bayesian Analysis of Macroevolutionary Mixtures)
# http://bamm-project.org/

## Explain that BAMM need to be installed first
# Put BAMM in a hidden directory (not included in the package build)
# BAMM is only available on Windows and OS X

# Explain that other diversification models with different assumptions can be used as long as they model regime shifts. Provide other examples. Explain that their output should then be formatted as in bammdata objects
# Structure of STRAPP tests relies on posterior samples so need a Bayesian approach to estimate the model parameters

# Use BAMMtools::setBAMMpriors to define priors from the phylogeny

# Ask the path to the BAMM.exe as argument
# Ask path for directory used to store input/output files generated
# Ask if the directory should be kept or erased

# See my BAMM script (Script 15)!
  # See if subsections can/should be split in sub-functions

## Sub-functions:
  # set_BAMM() => set default config files
  # run_BAMM() => run rjMCMC and deal with outputs
  # evaluate_BAMM() => produce traces for convergence checks and ESS; ask if exported in PDF?
  # import_BAMM() => including burn-in and selection of posterior samples to reach the desire numbers
  # clean_BAMM_files() => remove all files generated
  # plot_BAMM_rates() => See if simply using BAMMtools::plot.bammdata()

# Expected number of shifts based on ???
  # Empirical rule derived from number of tips. Explain the rule and warn that best practice is to try several runs and compare posterior samples to the expected distribution and select the value of expected number of shifts that is the closest match.
  # See if the BAMM manual provide any sort of empirical rules or if BAMMtools allows a default selection.

## Best practice = run multiple runs and check for convergence of the MCMC traces.
 # Here, a single run. Inspect chain stability beyond the burn-in.

## Find how to handle the call to the consoles from different systems!

## Output = the BAMM object used in run_deepSTRAPP, update_rates_and_regimes_for_focal_time and extract_diversification_data_melted_df_for_focal_time, and that can be plotted with BAMMtools::plot.bammdata()

### Outputs from BAMM (default names)
# run_info.txt containing a summary of your parameters/settings
# mcmc_log.txt containing raw MCMC information useful in diagnosing convergence
# event_data.txt containing all evolutionary rate parameters and their topological mappings
# chain_swap.txt containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).
# priors.txt containing default priors used to set up the run, as computed from BAMMtools::setBAMMpriors()

# acceptance_info.txt containing the history of acceptance/proposal of MCMC steps (If additional parameter 'outputAcceptanceInfo' is set to 1)

# Optional outputs
  # BAMM regime tables = $EventData (already included?)



# Include [deepSTRAPP::prepare_diversification_data()] in related functions
# http://bamm-project.org/
# R package BAMMtools = companion package of BAMM for post-processing of BAMM outputs
# Ref to the Rabosky paper

# Examples cannot be ran because they involve the BAMM software that is not installed on CRAN machines
# Still provide examples (but with not run balises) + something to run like plots of the output already loaded

# ape::write.tree
# stringr::str_detect stringr::str_remove
# utils::read.csv write.csv
# @importFrom ggplot2 ggplot geom_line geom_point geom_vline aes labs ggtitle theme element_line element_rect element_text unit margin
# coda::effectiveSize
# cowplot::save_plot
# BAMMtools::plotPrior
# @importFrom grDevices pdf dev.off

# Try on a small example.

# BAMM_install_directory_path <- "./software/bamm-2.5.0/"

# Check if BAMM works with non ultrametric trees (fossils). If not working add validity checks and info in doc.

# library(phytools)
# data(eel.tree)
#
# BAMM_object <- prepare_diversification_data(
#    BAMM_install_directory_path = "./software/bamm-2.5.0/",
#    phylo = eel.tree,
#    prefix_for_files = "eel",
#    numberOfGenerations = 10000,
#    expectedNumberOfShifts = 1,
#    burn_in = 0.5,
#    nb_posterior_samples = 10,
#    additional_BAMM_settings = list(updateLambdaShiftScale = 0.5),
#    skip_evaluations = TRUE,
#    plot_evaluations = TRUE)
#
# plot.bammdata(BAMM_object, labels = TRUE)

# Try to run quick example to see if it works once built
 # With and without keeping files
 # Add the option to remove the directory if empty

prepare_diversification_data <- function (BAMM_install_directory_path, # Ask the path to directory where is the BAMM 'executable'. Use '/' to separate directory and subdirectorys. It must end with '/'.
                                          phylo, # Phylogeny. Object of class phylo. Must be rooted and fully resolved.
                                          prefix_for_files = NULL, # To provide the prefix to add to all BAMM files stored in the 'BAMM_output_directory_path' and kept if 'keep_BAMM_outputs = TRUE'.
                                          # Files will exported such as 'prefix_*' with an underscore separating the prefix and the file name.
                                          seed = NULL, # Set for reproducibility
                                          numberOfGenerations = 10^7, # Number of steps in the MCMC run. Should be set high enough to reach the equilibrium distribution, and allows posterior samples to be decorrelated (check the Effective Sample Size of parameters with coda::effectiveSize() in the Evaluation step)
                                          globalSamplingFraction = 1.0, # Global sampling fraction representing the overall proportion of terminals in the phylogeny compared to the estimated overall richness in the clade. It acts as a multipliers on the rates needed to achieve such extant diversity.
                                          sampleProbsFilename = NULL, # Provide path to clade-specific sampling fractions. See ?BAMMtools::samplingProbs() to generate the file. It must be a '.txt' file. If provided, 'globalSamplingFraction' is ignored.
                                          expectedNumberOfShifts = NULL, # Set the expected number of regime shifts used to set the exponential hyperprior used to modulate reversible jumps across model configuration in the rjMCMC run.
                                                                         # If set to NULL (default), will use an empirical rule to define this value such as we expected a regime shift for every 100 tips in the phylogeny, with a minimum of 1.
                                                                         # Best practice is to try runs with several value and inspect the convergence of the posterior distribution of the regime shift parameter with the prior distribution defined by this hyperprior parameter.
                                          eventDataWriteFreq = NULL, # Set the frequency in which to write the event data to the output file = the sampling frequency of posterior samples.
                                                                     # Aim for 500-5000 posterior samples ideally as some samples will be removed to account for burn-in.
                                                                     # If set to NULL (default), will set frequency such as 2000 posterior samples are recorded.
                                          burn_in = 0.25, # Proportion of posterior samples removed to ensure the remaining samples where drawn once the equilibrium distribution was reached. This can be evaluated looking at the MCMC trace (see Evaluation step). Default is '0.25'.
                                          nb_posterior_samples = 1000,
                                          additional_BAMM_settings = list(), # Additional settings options for BAMM provided as a list of named argument. Ex: list(lambdaInit0 = 0.5, muInit0 = 0). See details in the template file provided with the package files 'BAMM_template_diversification.txt'.
                                          BAMM_output_directory_path = "./BAMM_outputs/", # Ask path for directory used to store input/output files generated. Use '/' to separate directory and subdirectories. It must end with '/'.
                                          keep_BAMM_outputs = TRUE, # Ask if the directory should be kept or erased. If 'BAMM_output_directory' is empty, it will be removed too.
                                          skip_evaluations = FALSE, # To skip the evaluation step (MCMC trace, ESS, and prior/posterior comparisons for LAMBDA = parameter controlling the expected nb of shifts)
                                          plot_evaluations = FALSE, # To display evaluation plots (MCMC trace, and prior/posterior comparisons of the expected nb of shifts)
                                          save_evaluations = TRUE) # To save outputs of evaluations: PDFs and table. 'MCMC_trace_logLik.pdf'. ESS with coda::effectiveSize() => 'ESS_df.csv'. Prior/posterior comparisons of the expected nb of shifts with BAMMtools::plotPrior() => 'PP_lambda_plot.pdf'

{
  ### Check input validity
  {
    ## BAMM_install_directory_path
    # BAMM_install_directory_path must be a directory, so it must end with '/'
    if (!stringr::str_detect(BAMM_install_directory_path, pattern = "/$"))
    {
      stop(paste0("'BAMM_install_directory_path' must end with '/'"))
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

    ## prefix_for_files
    # If provided, prefix_for_files should be character string
    if (!is.null(prefix_for_files))
    {
      if (!is.character(prefix_for_files))
      {
        stop(paste0("'prefix_for_files' must be a character string.\n",
                    "# Files will exported such as 'prefix_*' with an underscore separating the prefix and the file name."))
      }
    }

    ## seed
    if (!is.null(seed))
    {
      if (!is.numeric(seed))
      {
        stop(paste0("'seed' must be an interger."))
      }
    }

    ## numberOfGenerations
    # numberOfGenerations must be a positive integer
    if ((numberOfGenerations != abs(numberOfGenerations)) | (numberOfGenerations != round(numberOfGenerations)))
    {
      stop(paste0("'numberOfGenerations' must be a positive integer defining the number of steps in the MCMC run.\n",
                  "It should be set high enough to reach the equilibrium distribution, and allows posterior samples to be decorrelated.\n",
                  "Check the Effective Sample Size of parameters with coda::effectiveSize() in the Evaluation step.\n",
                  "Default is '10^7'."))
    }

    ## globalSamplingFraction
    # globalSamplingFraction should be a numerical between 0 and 1
    if ((globalSamplingFraction < 0) | (globalSamplingFraction > 1))
    {
      stop(paste0("'globalSamplingFraction' represents the overall proportion of terminals in the phylogeny compared to the estimated overall richness in the clade.\n",
                  "It acts as a multipliers on the rates needed to achieve such extant diversity. It must be between 0 and 1.\n",
                  "Current value of 'globalSamplingFraction' is ",globalSamplingFraction,"."))
    }
    if (globalSamplingFraction == 0)
    {
      stop(paste0("'globalSamplingFraction' represents the overall proportion of terminals in the phylogeny compared to the estimated overall richness in the clade.\n",
                  "It acts as a multipliers on the rates needed to achieve such extant diversity.\n",
                  "It must be between 0 and 1, but cannot be 0.\n"))
    }

    ## sampleProbsFilename
    # If provided, sampleProbsFilename must be a directory, so it must end with '/'
    if (!is.null(sampleProbsFilename))
    {
      if (!stringr::str_detect(sampleProbsFilename, pattern = "\\.txt$"))
      {
        stop(paste0("'sampleProbsFilename' must end with '.txt'"))
      }
    }

    ## expectedNumberOfShifts
    # If provided, expectedNumberOfShifts must be a positive integer
    if (!is.null(expectedNumberOfShifts))
    {
      if ((expectedNumberOfShifts != abs(expectedNumberOfShifts)) | (expectedNumberOfShifts != round(expectedNumberOfShifts)))
      {
        stop(paste0("'expectedNumberOfShifts' must be a positive integer defining the expected number of diversification regime shifts in the phylogeny.\n"),
             "This value is used to set the hyperprior from which the number of shifts is derived.\n",
             "Current value of 'expectedNumberOfShifts' is ",expectedNumberOfShifts,".")
      }
    }

    ## eventDataWriteFreq
    # If provided, eventDataWriteFreq must be a positive integer
    if (!is.null(eventDataWriteFreq))
    {
      if ((eventDataWriteFreq != abs(eventDataWriteFreq)) | (eventDataWriteFreq != round(eventDataWriteFreq)))
      {
        stop(paste0("'eventDataWriteFreq' must be a positive integer defining the frequency at which parameters are sampled in the MCMC run.\n",
                    "It is defined in number of steps between each sampling event.\n",
                    "Current value of 'eventDataWriteFreq' is ",eventDataWriteFreq,"."))
      }
    }

    ## burn_in
    # Burn-in between 0 and 1
    if ((burn_in < 0) | (burn_in > 1))
    {
      stop(paste0("'burn_in' represents the proportion of posterior samples removed to ensure the remaining samples where drawn once the equilibrium distribution was reached.\n",
                  "This can be evaluated looking at the MCMC trace (see Evaluation step). It must be between 0 and 1. Default is '0.25'.\n",
                  "Current value of 'burn_in' is ",burn_in,"."))
    }

    ## nb_posterior_samples
    # nb_posterior_samples must be a positive integer
    if ((nb_posterior_samples != abs(nb_posterior_samples)) | (nb_posterior_samples != round(nb_posterior_samples)))
    {
      stop(paste0("'nb_posterior_samples' must be a positive integer defining the number of posterior samples retainted in the 'BAMM_object' output used for downstream analyses.\n",
                  "Default is '1000'. Current value of 'nb_posterior_samples' is ",nb_posterior_samples,"."))
    }

    ## additional_BAMM_settings
    available_BAMM_settings <- c("runInfoFilename", "sampleFromPriorOnly", "runMCMC", "simulatePriorShifts", "loadEventData", "eventDataInfile",
                                 "initializeModel", "overwrite", "validateEventConfiguration", "lambdaInitPrior", "lambdaShiftPrior", "muInitPrior",
                                 "lambdaIsTimeVariablePrior", "mcmcOutfile", "mcmcWriteFreq", "eventDataOutfile", "printFreq", "outputAcceptanceInfo",
                                 "acceptanceInfoFileName", "acceptanceResetFreq", "updateLambdaInitScale", "updateLambdaShiftScale", "updateMuInitScale",
                                 "updateEventLocationScale", "updateEventRateScale", "updateRateEventNumber", "updateRateEventPosition", "updateRateEventRate",
                                 "updateRateLambda0", "updateRateLambdaShift", "updateRateMu0", "updateRateLambdaTimeMode", "localGlobalMoveRatio",
                                 "lambdaInit0", "lambdaShift0", "muInit0", "initialNumberEvents", "numberOfChains", "deltaT", "swapPeriod", "chainSwapFileName",
                                 "minCladeSizeForShift", "segLength")
    # Check if they match the list of available parameters.
    if (!all(names(additional_BAMM_settings) %in% available_BAMM_settings))
    {
      # Extract settings that do not match with available settings
      error_settings <- names(additional_BAMM_settings)[!(names(additional_BAMM_settings) %in% available_BAMM_settings)]
      initial_warning_length_options <- options()$warning.length
      options(warning.length = 2000L)
      stop(paste0("Names in 'additional_BAMM_settings' do not match with available settings: ",paste(error_settings, collapse = ", "),".\n",
                  "Available setting parameters are: ",paste(available_BAMM_settings, collapse = ", "),".\n",
                  "See details in the template file provided with the package files 'BAMM_template_diversification.txt'."))
      options(warning.length = initial_warning_length_options)
    }

    ## Incompatibility across parameters
    # Check combination of numberOfGenerations, eventDataWriteFreq, burn_in, and nb_posterior_samples to see if enough posterior samples remains after burn-in
    if (is.null(eventDataWriteFreq))
    {
      eventData_freq <- round(as.numeric(numberOfGenerations) / 2000) # Default = write 2000 posterior samples before burn-in
    } else {
      eventData_freq <- eventDataWriteFreq
    }
    nb_posterior_after_burn_in <- (as.numeric(numberOfGenerations)/eventData_freq) * (1-burn_in)
    if (nb_posterior_after_burn_in < nb_posterior_samples)
    {
      stop(paste0("'nb_posterior_samples' requested in not compatible with the number of posterior samples retainted after burn-in.\n",
                  "The number of remaining samples = 'numberOfGenerations'/'eventData_freq' x (1 - 'burn_in'). Here it is ",nb_posterior_after_burn_in,".\n",
                  "'nb_posterior_samples' requested is set to ", nb_posterior_samples,".\n",
                  "Please adjust 'numberOfGenerations', 'eventData_freq', 'burn_in', or 'nb_posterior_samples' to obtain compatible values."))

    }

    ## BAMM_output_directory_path
    # BAMM_output_directory_path must be a directory, so it must end with '/'
    if (!stringr::str_detect(BAMM_output_directory_path, pattern = "/$"))
    {
      stop(paste0("'BAMM_output_directory_path' must end with '/'"))
    }

    ## skip_evaluations & (plot_evaluations | save_evaluations)
    if (skip_evaluations & plot_evaluations)
    {
      cat(paste0("WARNING: 'plot_evaluations' is set to 'TRUE', but 'skip_evaluations' is set to TRUE.\n",
                 "Evaluations will not be plotted as they are skipped.\n\n"))
      warning(paste0("'plot_evaluations' was set to 'TRUE', but 'skip_evaluations' was set to TRUE.\n",
                     "Evaluations were not be plotted as they were skipped."))
    }
    if (skip_evaluations & save_evaluations)
    {
      cat(paste0("WARNING: 'save_evaluations' is set to 'TRUE', but 'skip_evaluations' is set to TRUE.\n",
                 "Evaluations will not be saved as they are skipped.\n\n"))
      warning(paste0("'save_evaluations' is set to 'TRUE', but 'skip_evaluations' is set to TRUE.\n",
                     "Evaluations were not be saved as they were skipped."))
    }

    ## keep_BAMM_outputs & save_evaluations
    if (!keep_BAMM_outputs & !skip_evaluations & save_evaluations)
    {
      warning(paste0("BAMM outputs files were removed ('keep_BAMM_outputs = FALSE'), but evaluation files were preserved in ",BAMM_output_directory_path," as 'save_evaluations' was set to 'TRUE'."))
    }
  }

  ##### set_BAMM() #####

  ## Get names of the additional BAMM settings
  add_settings_names <- names(additional_BAMM_settings)

  ## Create output directory if missing
  if (!dir.exists(paths = file.path(BAMM_output_directory_path)))
  {
    dir.create(path = file.path(BAMM_output_directory_path))
  }

  ## Export the phylogeny in a .tree file

  # Build path
  if (is.null(prefix_for_files))
  {
    phy_path <- file.path(paste0(BAMM_output_directory_path, "phylogeny.tree"))
  } else {
    phy_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files,"_phylogeny.tree"))
  }
  # Export tree
  ape::write.tree(phy = phylo, file = phy_path)

  ## Load control file template for diversification analyses

  # Load it from the root as in binary/installed version of the package
  BAMM_config_file <- tryCatch({
    readLines(con = file.path("./BAMM_template_diversification.txt"))
    }, warning = function(w) { }, # Do nothing
    error = function(e) { }, # Do nothing
    finally = { } # Do nothing
    )
  # If failed, load it from the /inst/ directory as in source/bundled version of the package
  if (is.null(BAMM_config_file))
  {
    BAMM_config_file <- readLines(con = file.path("./inst/BAMM_template_diversification.txt"))
  }

  # Initiate new config file for this analysis
  my_config_file <- BAMM_config_file

  ### 1/ Set general settings and data input ####

  # Set prefix to add to all output files (separated with "_")
  if (is.null(prefix_for_files))
  {
    # Do nothing = Add no prefix (line is commented out in the template)
  } else {
    outName_line <- which(stringr::str_detect(string = my_config_file, pattern = "outName = "))[1]
    my_config_file[outName_line] <- paste0("outName = ", prefix_for_files)
  }

  # Set path to the phylogenetic tree file
  phy_path_line <- which(stringr::str_detect(string = my_config_file, pattern = "treefile ="))
  my_config_file[phy_path_line] <- paste0("treefile = ", phy_path)

  # Set path to the file storing all information on this run
  if ("runInfoFilename" %in% add_settings_names)
  {
    runInfoFilename <- file.path(additional_BAMM_settings$runInfoFilename)
  } else {
    runInfoFilename <- "run_info.txt"
  }
  runInfoFilename_line <- which(stringr::str_detect(string = my_config_file, pattern = "runInfoFilename = "))
  my_config_file[runInfoFilename_line] <- paste0("runInfoFilename = ", runInfoFilename)

  # Should the run sample from prior only?
  # To check that the prior settings are fine
  # To compare with the posterior, to check the influence of the prior on the posterior
  if ("sampleFromPriorOnly" %in% add_settings_names)
  {
    sampleFromPriorOnly <- additional_BAMM_settings$sampleFromPriorOnly
  } else {
    sampleFromPriorOnly <- 0 # 0 = No. 1 = Yes
  }
  sampleFromPriorOnly_line <- which(stringr::str_detect(string = my_config_file, pattern = "sampleFromPriorOnly = "))
  my_config_file[sampleFromPriorOnly_line] <- paste0("sampleFromPriorOnly = ", sampleFromPriorOnly)

  # Should the MCMC simulation be performed?
  # If runMCMC = 0, the program will only check whether the data file can be read and the initial likelihood computed
  if ("runMCMC" %in% add_settings_names)
  {
    runMCMC <- additional_BAMM_settings$runMCMC
  } else {
    runMCMC <- 1  # 0 = No. 1 = Yes
  }
  runMCMC_line <- which(stringr::str_detect(string = my_config_file, pattern = "runMCMC = "))[1]
  my_config_file[runMCMC_line] <- paste0("runMCMC = ", runMCMC)

  # Should the prior distribution of the number of shift events, given the hyperprior on the Poisson rate parameter, be simulated?
  # This was necessary to compute Bayes factors
  # But is now disabled as the exact (analytical) prior is now implemented in BAMMtools.
  if ("simulatePriorShifts" %in% add_settings_names)
  {
    simulatePriorShifts <- additional_BAMM_settings$simulatePriorShifts
  } else {
    simulatePriorShifts <- 0  # 0 = No. 1 = Yes
  }
  simulatePriorShifts_line <- which(stringr::str_detect(string = my_config_file, pattern = "simulatePriorShifts = "))[1]
  my_config_file[simulatePriorShifts_line] <- paste0("simulatePriorShifts = ", simulatePriorShifts)

  # Whether to load a previous event data file
  if ("loadEventData" %in% add_settings_names)
  {
    loadEventData <- additional_BAMM_settings$loadEventData
  } else {
    loadEventData <- 0  # 0 = No. 1 = Yes
  }
  loadEventData_line <- which(stringr::str_detect(string = my_config_file, pattern = "loadEventData = "))[1]
  my_config_file[loadEventData_line] <- paste0("loadEventData = ", loadEventData)

  # Provides file name of the event data file to load, used only if loadEventData = 1
  if ("eventDataInfile" %in% add_settings_names)
  {
    eventDataInfile <- file.path(additional_BAMM_settings$eventDataInfile)
  } else {
    eventDataInfile <- "event_data_in.txt"
  }
  eventDataInfile_line <- which(stringr::str_detect(string = my_config_file, pattern = "eventDataInfile = "))[1]
  my_config_file[eventDataInfile_line] <- paste0("eventDataInfile = ", eventDataInfile)

  # Whether to initialize (but not run) the MCMC.
  # If initializeModel = 0, the program will only ensure that the data files (e.g., treefile) can be read
  if ("initializeModel" %in% add_settings_names)
  {
    initializeModel <- additional_BAMM_settings$initializeModel
  } else {
    initializeModel <- "event_data_in.txt"
  }
  initializeModel <- 1  # 0 = No. 1 = Yes
  initializeModel_line <- which(stringr::str_detect(string = my_config_file, pattern = "initializeModel = "))[1]
  my_config_file[initializeModel_line] <- paste0("initializeModel = ", initializeModel)

  # Whether to use a "global" sampling probability to assign the proportion of terminal represented in the tree.
  # If False (0), expects a file path for clade-specific sampling probabilities (see sampleProbsFilename)
  if (!is.null(sampleProbsFilename))
  {
    useGlobalSamplingProbability = 0  # 0 = No. 1 = Yes
  } else {
    useGlobalSamplingProbability = 1  # 0 = No. 1 = Yes
  }
  useGlobalSamplingProbability_line <- which(stringr::str_detect(string = my_config_file, pattern = "useGlobalSamplingProbability = "))[1]
  my_config_file[useGlobalSamplingProbability_line] <- paste0("useGlobalSamplingProbability = ", useGlobalSamplingProbability)

  # Provides the global sampling fraction
  # If useGlobalSamplingProbability = 0, this is ignored and BAMM looks for a file path to clade-specific sampling fractions
  globalSamplingFraction_line <- which(stringr::str_detect(string = my_config_file, pattern = "globalSamplingFraction = "))[1]
  my_config_file[globalSamplingFraction_line] <- paste0("globalSamplingFraction = ", globalSamplingFraction)

  # Provides path to the file containing clade-specific sampling fractions
  if (is.null(sampleProbsFilename))
  {
    sampleProbsFilename <- "taxa_sampling_probs.txt"
  } else {
    sampleProbsFilename <- file.path(sampleProbsFilename)
  }
  sampleProbsFilename_line <- which(stringr::str_detect(string = my_config_file, pattern = "sampleProbsFilename = "))[1]
  my_config_file[sampleProbsFilename_line] <- paste0("sampleProbsFilename = ", sampleProbsFilename)

  # Set the seed for the random number generator.
  # If not specified (or is -1), a seed is obtained from the system clock
  if (is.null(seed))
  {
    seed <- -1
  }
  seed_line <- which(stringr::str_detect(string = my_config_file, pattern = "seed = "))[1]
  my_config_file[seed_line] <- paste0("seed = ", seed)

  # Should the output files be overwritten?
  # If True (1), the program will overwrite any output files in the current directory (if present)
  if ("overwrite" %in% add_settings_names)
  {
    overwrite <- additional_BAMM_settings$overwrite
  } else {
    overwrite <- 1  # 0 = No. 1 = Yes
  }
  overwrite_line <- which(stringr::str_detect(string = my_config_file, pattern = "overwrite = "))[1]
  my_config_file[overwrite_line] <- paste0("overwrite = ", overwrite)

  # Set limits to valid configurations
  # If 1, rejects proposals that cause a branch and both of its direct descendants to have at least one event.
  # Such an event configuration may cause the parameters of the parent event to change to unrealistic values.
  # If 0, no such proposals are immediately rejected. The default value is 1.
  # validateEventConfiguration <- 0
  if ("validateEventConfiguration" %in% add_settings_names)
  {
    validateEventConfiguration <- additional_BAMM_settings$validateEventConfiguration
  } else {
    validateEventConfiguration <- 1  # 0 = Do not reject. 1 = Reject
  }
  validateEventConfiguration_line <- which(stringr::str_detect(string = my_config_file, pattern = "validateEventConfiguration = "))[1]
  my_config_file[validateEventConfiguration_line] <- paste0("validateEventConfiguration = ", validateEventConfiguration)

  ### 2/ Set (hyper)prior settings ####

  ## Can use this help function to automatically tune prior adapted to your data by scaling the prior distributions based on the age (root depth) of your tree
  # In practice, setBAMMpriors first estimates the rate of speciation for your full tree under a pure birth model of diversification.
  # Then assume, arbitrarily, that a reasonable prior distribution for the initial lambda0/mu0 rate parameters is an exponential distribution with a mean five times greater than this pure birth value.
  # Rationale = having a weakly informative prior that is still in the order of magnitude of the true rate
  # For the shift parameter (alpha), the sd of the normal prior is set such as mean +/- 2s gives an alpha parameter that results in
  # either a 90% decline in the evolutionary rate or a 190% increase in rate on the interval of time from the root to the tips of the tree.

  # Build path to priors file
  if (is.null(prefix_for_files))
  {
    priors_path <- file.path(paste0(BAMM_output_directory_path, "priors.txt"))
  } else {
    priors_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files,"_priors.txt"))
  }
  # Generate priors file
  invisible(capture.output(BAMMtools::setBAMMpriors(phy = phylo, outfile = priors_path, suppressWarning = TRUE)))
  # Read prior file
  default_tuned_priors <- readLines(con = priors_path)

  # Set the expected number of shifts used to set the exponential hyperprior for nb of rate shifts (from which the Λ is drawn)
  # Suggested values (if set to NULL): 1 regime shift for 100 tips.
  # Best practice is to try runs with several value and inspect the convergence of the posterior distribution of the regime shift parameter with the prior distribution defined by this hyperprior parameter.
  if (is.null(expectedNumberOfShifts))
  {
    expectedNumberOfShifts <- max(1, round(length(phylo$tip.label)/100)) # 1 expected regime shift for 100 tips
  }
  expectedNumberOfShifts_line <- which(stringr::str_detect(string = my_config_file, pattern = "expectedNumberOfShifts = "))[1]
  my_config_file[expectedNumberOfShifts_line] <- paste0("expectedNumberOfShifts = ", expectedNumberOfShifts)

  # Set the rate parameter of the exponential prior(s) of initial lambda parameters (lambda0) of speciation rate regimes
    # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
  if ("lambdaInitPrior" %in% add_settings_names)
  {
    # If provided as 'additional_BAMM_settings', use this value
    lambdaInitPrior <- additional_BAMM_settings$lambdaInitPrior
  } else {
    # If not provided, use the default value suggested by BAMMtools::setBAMMpriors
    lambdaInitPrior_default_line <- which(stringr::str_detect(string = default_tuned_priors, pattern = "lambdaInitPrior = "))[1]
    lambdaInitPrior <- as.numeric(stringr::str_remove(string = default_tuned_priors[lambdaInitPrior_default_line], pattern = "lambdaInitPrior = "))

  }
  lambdaInitPrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaInitPrior = "))[1]
  my_config_file[lambdaInitPrior_line] <- paste0("lambdaInitPrior = ", lambdaInitPrior)

  # Set the standard deviation of the normal distribution prior(s) of rate variation parameters (alpha) of speciation rate regimes
    # alpha in lambda(t) = lamba0 x exp(alpha*t)
    # Mean of this prior(s) are fixed to zero such as a constant rate diversification process is the most probable a priori
  if ("lambdaShiftPrior" %in% add_settings_names)
  {
    # If provided as 'additional_BAMM_settings', use this value
    lambdaShiftPrior <- additional_BAMM_settings$lambdaShiftPrior
  } else {
    # If not provided, use the default value suggested by BAMMtools::setBAMMpriors
    lambdaShiftPrior_default_line <- which(stringr::str_detect(string = default_tuned_priors, pattern = "lambdaShiftPrior = "))[1]
    lambdaShiftPrior <- as.numeric(stringr::str_remove(string = default_tuned_priors[lambdaShiftPrior_default_line], pattern = "lambdaShiftPrior = "))

  }
  lambdaShiftPrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaShiftPrior = "))[1]
  my_config_file[lambdaShiftPrior_line] <- paste0("lambdaShiftPrior = ", lambdaShiftPrior)

  # Set the rate parameter of the exponential prior(s) of initial lambda parameters (mu0) of extinction rate regimes
    # mu0 in mu(t) = mu0 x exp(alpha*t)
    # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
  if ("muInitPrior" %in% add_settings_names)
  {
    # If provided as 'additional_BAMM_settings', use this value
    muInitPrior <- additional_BAMM_settings$muInitPrior
  } else {
    # If not provided, use the default value suggested by BAMMtools::setBAMMpriors
    muInitPrior_default_line <- which(stringr::str_detect(string = default_tuned_priors, pattern = "muInitPrior = "))[1]
    muInitPrior <- as.numeric(stringr::str_remove(string = default_tuned_priors[muInitPrior_default_line], pattern = "muInitPrior = "))

  }
  muInitPrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "muInitPrior = "))[1]
  my_config_file[muInitPrior_line] <- paste0("muInitPrior = ", muInitPrior)

  # No prior for rate variation parameters (alpha) of extinction rate regimes as they are assumed to follow constant rates

  # Set the prior (probability) of the time mode (of speciation?) being time-variable (vs. time-constant)
    # By default, allows all regimes to be time-variable, as their rate can still be estimated as constant with alpha = 0
  if ("lambdaIsTimeVariablePrior" %in% add_settings_names)
  {
    lambdaIsTimeVariablePrior <- additional_BAMM_settings$lambdaIsTimeVariablePrior
  } else {
    lambdaIsTimeVariablePrior <- 1
  }
  lambdaIsTimeVariablePrior <- 1
  lambdaIsTimeVariablePrior_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaIsTimeVariablePrior = "))[1]
  my_config_file[lambdaIsTimeVariablePrior_line] <- paste0("lambdaIsTimeVariablePrior = ", lambdaIsTimeVariablePrior)

  ### 3/ Set the MCMC simulation settings, MCMC logs and output options ####

  # Set the number of generations to perform MCMC simulation
  numberOfGenerations <- format(numberOfGenerations, scientific = F) # 10^7
  numberOfGenerations_line <- which(stringr::str_detect(string = my_config_file, pattern = "numberOfGenerations = "))[1]
  my_config_file[numberOfGenerations_line] <- paste0("numberOfGenerations = ", numberOfGenerations)

  # Set the path to the MCMC output file
  # Includes only summary information about MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)
  if ("mcmcOutfile" %in% add_settings_names)
  {
    mcmcOutfile <- file.path(additional_BAMM_settings$mcmcOutfile)
  } else {
    mcmcOutfile <- "mcmc_log.txt"
  }
  mcmcOutfile_line <- which(stringr::str_detect(string = my_config_file, pattern = "mcmcOutfile = "))[1]
  my_config_file[mcmcOutfile_line] <- paste0("mcmcOutfile = ", mcmcOutfile)

  # Set the frequency in which to write the MCMC output to the log file
    # Aim for 500-5000 posterior samples ideally
  if ("mcmcWriteFreq" %in% add_settings_names)
  {
    mcmcWriteFreq <- additional_BAMM_settings$mcmcWriteFreq
  } else {
    mcmcWriteFreq <- round(as.numeric(numberOfGenerations) / 2000) # Default = record 2000 MCMC generations
  }
  mcmcWriteFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "mcmcWriteFreq = "))[1]
  my_config_file[mcmcWriteFreq_line] <- paste0("mcmcWriteFreq = ", mcmcWriteFreq)

  # Set the path to the main output file
    # The raw event data. ALL of the results are contained in this file,
    # All branch-specific speciation rates, shift positions, marginal distributions, etc. can be reconstructed from this output.
    # See ?BAMMtools::getEventData to import this output in R
  if ("eventDataOutfile" %in% add_settings_names)
  {
    eventDataOutfile <- file.path(additional_BAMM_settings$eventDataOutfile)
  } else {
    eventDataOutfile <- "event_data.txt"
  }
  eventDataOutfile_line <- which(stringr::str_detect(string = my_config_file, pattern = "eventDataOutfile = "))[1]
  my_config_file[eventDataOutfile_line] <- paste0("eventDataOutfile = ", eventDataOutfile)

  # Set frequency in which to write the event data to the output file = the sampling frequency of posterior samples
    # Aim for 500-5000 posterior samples ideally
    # Will need to remove some to account for the burn-in
  if (is.null(eventDataWriteFreq))
  {
    eventDataWriteFreq <- round(as.numeric(numberOfGenerations) / 2000) # Default = write 2000 posterior samples before burn-in
  }
  eventDataWriteFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "eventDataWriteFreq = "))[1]
  my_config_file[eventDataWriteFreq_line] <- paste0("eventDataWriteFreq = ", eventDataWriteFreq)

  # Set frequency in which to print MCMC status to the screen
  if ("printFreq" %in% add_settings_names)
  {
    printFreq <- additional_BAMM_settings$printFreq
  } else {
    if (numberOfGenerations >= 10^6)
    {
      printFreq <- 10000 # Print status every 10^4 generations for long runs
    } else {
      printFreq <- 1000 # Print status every 10^3 generations for short runs
    }
  }
  printFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "printFreq = "))[1]
  my_config_file[printFreq_line] <- paste0("printFreq = ", printFreq)

  # Whether acceptance/proposal history should be saved.
    # If 1, outputs whether each proposal was accepted. The number identifying the proposal matches the one in the code.
    # The default value is 0 (i.e., do not output this information).
  if ("outputAcceptanceInfo" %in% add_settings_names)
  {
    outputAcceptanceInfo <- additional_BAMM_settings$outputAcceptanceInfo
  } else {
    outputAcceptanceInfo <- 0 # Do not save acceptance/proposal history
  }
  outputAcceptanceInfo_line <- which(stringr::str_detect(string = my_config_file, pattern = "outputAcceptanceInfo = "))[1]
  my_config_file[outputAcceptanceInfo_line] <- paste0("outputAcceptanceInfo = ", outputAcceptanceInfo)

  # Set the path to the acceptance/proposal history file
    # The path of the file to which to write whether each proposal was accepted.
    # outputAcceptanceInfo must be set to 1 for this information to be written.
  if ("acceptanceInfoFileName" %in% add_settings_names)
  {
    acceptanceInfoFileName <- file.path(additional_BAMM_settings$acceptanceInfoFileName)
  } else {
    acceptanceInfoFileName <- "acceptance_info.txt"
  }
  acceptanceInfoFileName_line <- which(stringr::str_detect(string = my_config_file, pattern = "acceptanceInfoFileName = "))[1]
  my_config_file[acceptanceInfoFileName_line] <- paste0("acceptanceInfoFileName = ", acceptanceInfoFileName)

  # Set frequency in which to update the acceptance rate calculation
    # Acceptance rate = how often new proposal are accepted as the next step in the chain
    # Give information on how efficient is the movement of the chain in the parameter space
    # The acceptance rate is output to both the MCMC data file and print to the screen
  if ("acceptanceResetFreq" %in% add_settings_names)
  {
    acceptanceResetFreq <- additional_BAMM_settings$acceptanceResetFreq
  } else {
    if (numberOfGenerations >= 10^6)
    {
      acceptanceResetFreq <- 10000/2 # Update acceptance rate every 5000 generations for long runs
    } else {
      acceptanceResetFreq <- 1000/2 # Print status every 500 generations for short runs
    }
  }
  acceptanceResetFreq_line <- which(stringr::str_detect(string = my_config_file, pattern = "acceptanceResetFreq = "))[1]
  my_config_file[acceptanceResetFreq_line] <- paste0("acceptanceResetFreq = ", acceptanceResetFreq)

  ### 4/ Set the scaling operators = temperatures, to propose new values for sampled parameters ####

  # The highest scaling operators = temperatures = the bigger changes can be implemented
    # Advantages = allows to escape suboptimum
    # Cons = may be unstable / harder to reach convergence

  # Set scale parameter used for updating the initial speciation rate (lambda0) for each regime/process
    # Updated as lambda0_new ~ lambda0_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
  if ("updateLambdaInitScale" %in% add_settings_names)
  {
    updateLambdaInitScale <- additional_BAMM_settings$updateLambdaInitScale
  } else {
    updateLambdaInitScale <- 2.0
  }
  updateLambdaInitScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateLambdaInitScale = "))[1]
  my_config_file[updateLambdaInitScale_line] <- paste0("updateLambdaInitScale = ", updateLambdaInitScale)

  # Set window size parameter used for updating the rate variation parameter (alpha) of speciation rates for each regime/process
    # Updated as alpha_new ~ alpha_old + U with U a uniform distribution ranging between - window_size_par and + window_size_par
  if ("updateLambdaShiftScale" %in% add_settings_names)
  {
    updateLambdaShiftScale <- additional_BAMM_settings$updateLambdaShiftScale
  } else {
    updateLambdaShiftScale <- 0.1
  }
  updateLambdaShiftScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateLambdaShiftScale = "))[1]
  my_config_file[updateLambdaShiftScale_line] <- paste0("updateLambdaShiftScale = ", updateLambdaShiftScale)

  # Set scale parameter used for updating the initial extinction rate (mu0) for each regime/process
    # Updated as mu0_new ~ mu0_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
  if ("updateMuInitScale" %in% add_settings_names)
  {
    updateMuInitScale <- additional_BAMM_settings$updateMuInitScale
  } else {
    updateMuInitScale <- 2.0
  }
  updateMuInitScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateMuInitScale = "))[1]
  my_config_file[updateMuInitScale_line] <- paste0("updateMuInitScale = ", updateMuInitScale)

  # Set window size parameter used for updating LOCAL moves of the position of shifts on the tree
    # Updated as position_new ~ position_old + U with U a uniform distribution ranging between - window_size_par and + window_size_par
    # Unit = fraction of root_to_tip length. May lead to jump of the shift position to a new branch
    # Ex: For a tree of 100My, with parameter set to 0.05, the proposal window for local position change is +/- 5My around the previous value
  if ("updateEventLocationScale" %in% add_settings_names)
  {
    updateEventLocationScale <- additional_BAMM_settings$updateEventLocationScale
  } else {
    updateEventLocationScale <- 0.05
  }
  updateEventLocationScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateEventLocationScale = "))[1]
  my_config_file[updateEventLocationScale_line] <- paste0("updateEventLocationScale = ", updateEventLocationScale)

  # Set scale parameter used for updating the LAMBDA rate parameter of the Poisson process controlling the number of shifts in the submodels
    # Updated as LAMBDA_new ~ LAMBDA_old x exp(scaling_par x (U - 0.5)) with U a uniform distribution ranging between 0 and 1
  if ("updateEventRateScale" %in% add_settings_names)
  {
    updateEventRateScale <- additional_BAMM_settings$updateEventRateScale
  } else {
    updateEventRateScale <- 4.0
  }
  updateEventRateScale_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateEventRateScale = "))[1]
  my_config_file[updateEventRateScale_line] <- paste0("updateEventRateScale = ", updateEventRateScale)

  ### 5/ Set the relative frequencies of operator uses (frequency of parameter updates) at each generation ####

  # Set the relative frequency of MCMC moves that change the number of events in the submodel (shift from Mk to Mk+1 or Mk-1)
  if ("updateRateEventNumber" %in% add_settings_names)
  {
    updateRateEventNumber <- additional_BAMM_settings$updateRateEventNumber
  } else {
    updateRateEventNumber <- 0.1 # 1/10 generations
  }
  updateRateEventNumber_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateEventNumber = "))[1]
  my_config_file[updateRateEventNumber_line] <- paste0("updateRateEventNumber = ", updateRateEventNumber)

  # Set the relative frequency of MCMC moves that change the location of an event on the tree (update position parameters)
  if ("updateRateEventPosition" %in% add_settings_names)
  {
    updateRateEventPosition <- additional_BAMM_settings$updateRateEventPosition
  } else {
    updateRateEventPosition <- 1 # Every generation !
  }
  updateRateEventPosition_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateEventPosition = "))[1]
  my_config_file[updateRateEventPosition_line] <- paste0("updateRateEventPosition = ", updateRateEventPosition)

  # Set the relative frequency of MCMC moves that change the rate at which events occur (update LAMBDA parameter)
  if ("updateRateEventRate" %in% add_settings_names)
  {
    updateRateEventRate <- additional_BAMM_settings$updateRateEventRate
  } else {
    updateRateEventRate <- 1 # Every generation !
  }
  updateRateEventRate_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateEventRate = "))[1]
  my_config_file[updateRateEventRate_line] <- paste0("updateRateEventRate = ", updateRateEventRate)

  # Set the relative frequency of MCMC moves that change the initial speciation rates (lambda0) associated with a regime
    # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
  if ("updateRateLambda0" %in% add_settings_names)
  {
    updateRateLambda0 <- additional_BAMM_settings$updateRateLambda0
  } else {
    updateRateLambda0 <- 1 # Every generation !
  }
  updateRateLambda0_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateLambda0 = "))[1]
  my_config_file[updateRateLambda0_line] <- paste0("updateRateLambda0 = ", updateRateLambda0)

  # Set the relative frequency of MCMC moves that change the exponential shift parameter (alpha) of the speciation rate associated with a regime
    # alpha in lambda(t) = lamba0 x exp(alpha*t)
  if ("updateRateLambdaShift" %in% add_settings_names)
  {
    updateRateLambdaShift <- additional_BAMM_settings$updateRateLambdaShift
  } else {
    updateRateLambdaShift <- 1 # Every generation !
  }
  updateRateLambdaShift_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateLambdaShift = "))[1]
  my_config_file[updateRateLambdaShift_line] <- paste0("updateRateLambdaShift = ", updateRateLambdaShift)

  # Set the relative frequency of MCMC moves that change the (initial) extinction rate associated with a regime
    # mu0 in mu(t) = mu0 x exp(alpha*t)
    # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
  if ("updateRateMu0" %in% add_settings_names)
  {
    updateRateMu0 <- additional_BAMM_settings$updateRateMu0
  } else {
    updateRateMu0 <- 1 # Every generation !
  }
  updateRateMu0_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateMu0 = "))[1]
  my_config_file[updateRateMu0_line] <- paste0("updateRateMu0 = ", updateRateMu0)

  # Set the relative frequency of MCMC moves that flip the time mode (time-constant <=> time-variable)
  # By default, only use time-variable mode, so the frequency is set to 0.
  if ("updateRateLambdaTimeMode" %in% add_settings_names)
  {
    updateRateLambdaTimeMode <- additional_BAMM_settings$updateRateLambdaTimeMode
  } else {
    updateRateLambdaTimeMode <- 0 # Never shift to time-constant mode
  }
  updateRateLambdaTimeMode_line <- which(stringr::str_detect(string = my_config_file, pattern = "updateRateLambdaTimeMode = "))[1]
  my_config_file[updateRateLambdaTimeMode_line] <- paste0("updateRateLambdaTimeMode = ", updateRateLambdaTimeMode)

  # Set the ratio of local to global moves used to propose new location of events on the tree (update position parameters)
  if ("localGlobalMoveRatio" %in% add_settings_names)
  {
    localGlobalMoveRatio <- additional_BAMM_settings$localGlobalMoveRatio
  } else {
    localGlobalMoveRatio <- 10.0 # Ten times more local changes than global changes
  }
  localGlobalMoveRatio_line <- which(stringr::str_detect(string = my_config_file, pattern = "localGlobalMoveRatio = "))[1]
  my_config_file[localGlobalMoveRatio_line] <- paste0("localGlobalMoveRatio = ", localGlobalMoveRatio)

  ### 6/ Set the initial parameter values to start the MCMC chain(s) ####

  # The MCMC chain start with a model with no shift (M0 submodel). So the initial parameter values are for this unique regime

  # Run a BD model to obtain credible starting values for the root process
  BD_fit <- suppressWarnings(phytools::fit.bd(tree = phylo))

  # Set the initial speciation rate (lambda0) for the first regime starting at the root of the tree (regime 0)
    # lambda0 in lambda(t) = lamba0 x exp(alpha*t)
  if ("lambdaInit0" %in% add_settings_names)
  {
    lambdaInit0 <- additional_BAMM_settings$lambdaInit0
  } else {
    lambdaInit0 <- BD_fit$b # Provide estimates from BD model
  }
  lambdaInit0_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaInit0 = "))[1]
  my_config_file[lambdaInit0_line] <- paste0("lambdaInit0 = ", lambdaInit0)

  # Set the initial shift parameter (alpha) for the root process (regime 0)
    # alpha in lambda(t) = lamba0 x exp(alpha*t)
  if ("lambdaShift0" %in% add_settings_names)
  {
    lambdaShift0 <- additional_BAMM_settings$lambdaShift0
  } else {
    lambdaShift0 <- 0 # Initial value set to 0 such as the process is a constant rate
  }
  lambdaShift0_line <- which(stringr::str_detect(string = my_config_file, pattern = "lambdaShift0 = "))[1]
  my_config_file[lambdaShift0_line] <- paste0("lambdaShift0 = ", lambdaShift0)

  # Set the intial extinction rate (mu0) for the first regime starting at the root of the tree (regime 0)
    # mu0 in mu(t) = mu0 x exp(alpha*t)
    # As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
  if ("muInit0" %in% add_settings_names)
  {
    muInit0 <- additional_BAMM_settings$muInit0
  } else {
    muInit0 <- BD_fit$d # Provide estimates from BD model
  }
  muInit0_line <- which(stringr::str_detect(string = my_config_file, pattern = "muInit0 = "))[1]
  my_config_file[muInit0_line] <- paste0("muInit0 = ", muInit0)

  # Set the initial number of non-root processes/shifts = M0 submodel
  if ("initialNumberEvents" %in% add_settings_names)
  {
    initialNumberEvents <- additional_BAMM_settings$initialNumberEvents
  } else {
    initialNumberEvents <- 0 # Start with a single regime = M0 submodel
  }
  initialNumberEvents_line <- which(stringr::str_detect(string = my_config_file, pattern = "initialNumberEvents = "))[1]
  my_config_file[initialNumberEvents_line] <- paste0("initialNumberEvents = ", muInit0)

  ### 7/ Set the MCMC chain behavior ####

  # Set the number of Markov chains to run
    # Each chain will have a different temperature to favor different exploration behavior of the parameter space
  if ("numberOfChains" %in% add_settings_names)
  {
    numberOfChains <- additional_BAMM_settings$numberOfChains
  } else {
    numberOfChains <- 4 # Default = using 4 coupled-chains to explore pararameter space for a single run.
  }
  numberOfChains_line <- which(stringr::str_detect(string = my_config_file, pattern = "numberOfChains = "))[1]
  my_config_file[numberOfChains_line] <- paste0("numberOfChains = ", numberOfChains)

  # Set the temperature increment parameter that control the difference of temperatures between the chains.
    # This value should be > 0
    # The temperature for the i-th chain is computed as 1 / [1 + deltaT * (i - 1)]
    # Chain 1 is the coldest. Highest chain is the hottest
  if ("deltaT" %in% add_settings_names)
  {
    deltaT <- additional_BAMM_settings$deltaT
  } else {
    deltaT <- 0.01
  }
  deltaT_line <- which(stringr::str_detect(string = my_config_file, pattern = "deltaT = "))[1]
  my_config_file[deltaT_line] <- paste0("deltaT = ", deltaT)

  # Set the frequency of generations at which to propose a chain swap
    # The coupled-MCMC algorithm will check chain state and eventually swap for the one having reach the highest likelihood
  if ("swapPeriod" %in% add_settings_names)
  {
    swapPeriod <- additional_BAMM_settings$swapPeriod
  } else {
    swapPeriod <- 1000 # Check swapping every 10^3 generations
  }
  swapPeriod_line <- which(stringr::str_detect(string = my_config_file, pattern = "swapPeriod = "))[1]
  my_config_file[swapPeriod_line] <- paste0("swapPeriod = ", swapPeriod)

  # Set the path to the file where to store information about each chain swap proposal
    # The format of each line is [generation],[rank_1],[rank_2],[swap_accepted]
    # where [generation] is the generation in which the swap proposal was made,
    # [rank_1] and [rank_2] are the chains that were chosen, and [swap_accepted] is
    # whether the swap was made. The cold chain has a rank of 1.
  if ("chainSwapFileName" %in% add_settings_names)
  {
    chainSwapFileName <- file.path(additional_BAMM_settings$chainSwapFileName)
  } else {
    chainSwapFileName <- "chain_swap_log.txt"
  }
  chainSwapFileName_line <- which(stringr::str_detect(string = my_config_file, pattern = "chainSwapFileName = "))[1]
  my_config_file[chainSwapFileName_line] <- paste0("chainSwapFileName = ", chainSwapFileName)

  ### 8/ Set other parameters ####

  # Set the minimum size of a clade to allow a shift to occur
    # Constrain location of possible rate-change events to occur only on branches with at least this many descendant tips.
    # The default value of 1 allows shifts to occur on all branches.
  if ("minCladeSizeForShift" %in% add_settings_names)
  {
    minCladeSizeForShift <- additional_BAMM_settings$minCladeSizeForShift
  } else {
    minCladeSizeForShift <- 1 # Shift can occur on all branches, even the terminal ones.
  }
  minCladeSizeForShift_line <- which(stringr::str_detect(string = my_config_file, pattern = "minCladeSizeForShift = "))[1]
  my_config_file[minCladeSizeForShift_line] <- paste0("minCladeSizeForShift = ", minCladeSizeForShift)

  # Set the "grain" at which time-continuous calculations are discretized
    # The continuous-time change in diversification rates are approximated by breaking each branch into constant-rate diversification segments
    # with each segment given a length determined by the segLength parameter.
    # segLength is in fraction of the root-to-tip distance of the tree.
    # Ex: For an ultrametric tree of 100My, a segLength of 0.02 lead to a step size of 2My
    # If the value is greater than a given branch length BAMM will not break the branch into segments but use the mean rate across the entire branch.
  if ("segLength" %in% add_settings_names)
  {
    segLength <- additional_BAMM_settings$segLength
  } else {
    segLength = 0.01 # Use segments of length = 1% of tree depth to discretize calculations of time-continuous rates
  }
  segLength_line <- which(stringr::str_detect(string = my_config_file, pattern = "segLength = "))[1]
  my_config_file[segLength_line] <- paste0("segLength = ", segLength)

  ### 9/ Export the updated custom config file ####

  ## Export the config file with the prefix
  # Build path
  if (is.null(prefix_for_files))
  {
    config_file_path <- file.path(paste0(BAMM_output_directory_path, "config_file.txt"))
  } else {
    config_file_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files,"_config_file.txt"))
  }
  # Export my_config_file
  writeLines(text = my_config_file, con = config_file_path)

  ##### run_BAMM() #####

  ## run BAMM
  BAMM_path <- file.path(paste0(BAMM_install_directory_path, "bamm"))
  system(paste0(BAMM_path, " -c ", config_file_path))

  ### Outputs
  # run_info.txt file containing a summary of your parameters/settings
  # mcmc_log.txt file containing raw MCMC information useful in diagnosing convergence
  # event_data.txt file containing all evolutionary rate parameters and their topological mappings
  # chain_swap_log.txt file containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).
  # acceptance_info.txt containing the history of acceptance/proposal of MCMC steps (If additional parameter 'outputAcceptanceInfo' is set to 1)

  ## Clean outputs = move files to the dedicated directory
  # Detect output files
  if (is.null(prefix_for_files))
  {
    file.rename(from = runInfoFilename, to = file.path(paste0(BAMM_output_directory_path, runInfoFilename)))
    file.rename(from = mcmcOutfile, to = file.path(paste0(BAMM_output_directory_path, mcmcOutfile)))
    file.rename(from = eventDataOutfile, to = file.path(paste0(BAMM_output_directory_path, eventDataOutfile)))
    file.rename(from = chainSwapFileName, to = file.path(paste0(BAMM_output_directory_path, chainSwapFileName)))
    if (outputAcceptanceInfo == 1)
    {
      file.rename(from = acceptanceInfoFileName, to = file.path(paste0(BAMM_output_directory_path, acceptanceInfoFileName)))
    }
  } else {
    file.rename(from = file.path(paste0(prefix_for_files, "_", runInfoFilename)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", runInfoFilename)))
    file.rename(from = file.path(paste0(prefix_for_files, "_", mcmcOutfile)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", mcmcOutfile)))
    file.rename(from = file.path(paste0(prefix_for_files, "_", eventDataOutfile)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", eventDataOutfile)))
    file.rename(from = file.path(paste0(prefix_for_files, "_", chainSwapFileName)), to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", chainSwapFileName)))
    if (outputAcceptanceInfo == 1)
    {
      file.rename(from = acceptanceInfoFileName, to = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", acceptanceInfoFileName)))
    }
  }

  ##### evaluate_BAMM() #####

  if (!skip_evaluations)
  {

    ## Plot MCMC trace for logLik

    # Load the MCMC log file
    if (is.null(prefix_for_files))
    {
      MCMC_log <- utils::read.csv(file.path(paste0(BAMM_output_directory_path, mcmcOutfile)), header = T)
    } else {
      MCMC_log <- utils::read.csv(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", mcmcOutfile)), header = T)
    }

    # Find generations used to cut-off the burn-in
    burn_in_threshold <- ceiling(burn_in * MCMC_log$generation[nrow(MCMC_log)])

    # Create binding of new variables to avoid Notes
    logLik <- generation <- NULL

    # Create MCMC trace plot for logLik
    MCMC_logLik_ggplot <- ggplot2::ggplot(data = MCMC_log,
                                          mapping = ggplot2::aes(y = logLik, x = generation)) +

      ggplot2::geom_line(linewidth = 1.0, alpha = 1.0) +
      ggplot2::geom_point(size = 1.5, alpha = 0.8) +

      ggplot2::geom_vline(xintercept = burn_in_threshold, linewidth = 1.5, linetype = 2, color = "red") +

      ggplot2::labs(x = "Generations", y = "LogLikelihood") +

      ggplot2::ggtitle("MCMC trace for logLik") +

      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.5),
                     panel.background = ggplot2::element_rect(fill = NA, color = NA),
                     plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                                        margin = ggplot2::margin(b = 15, t = 5)),
                     axis.title = ggplot2::element_text(size = 20, color = "black"),
                     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
                     axis.line = ggplot2::element_line(linewidth = 1.5),
                     axis.text = ggplot2::element_text(size = 18, color = "black"),
                     axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
                     axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

    # Display plot if requested
    if (plot_evaluations)
    {
      print(MCMC_logLik_ggplot)
    }

    # Save plot if requested
    if (save_evaluations)
    {
      if (is.null(prefix_for_files))
      {
        MCMC_logLik_path <- file.path(paste0(BAMM_output_directory_path, "MCMC_trace_logLik.pdf"))
      } else {
        MCMC_logLik_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_MCMC_trace_logLik.pdf"))
      }
      cowplot::save_plot(plot = MCMC_logLik_ggplot,
                         filename = MCMC_logLik_path,
                         base_height = 8, base_width = 10)
    }

    ## Remove burn-in from MCMC log to explore ESS and prior/posterior of LAMBDA (parameter controlling the nb of rates)
    post_burn_MCMC_log <- MCMC_log[MCMC_log$generation >= burn_in_threshold, ]

    ## Explore effective sample sizes

    # Store ESS info
    ESS_df <- data.frame(ESS_logLik = coda::effectiveSize(post_burn_MCMC_log$logLik),
                         ESS_logPrior = coda::effectiveSize(post_burn_MCMC_log$logPrior),
                         ESS_N_shifts = coda::effectiveSize(post_burn_MCMC_log$N_shifts),
                         ESS_eventRate = coda::effectiveSize(post_burn_MCMC_log$eventRate),
                         ESS_acceptRate = coda::effectiveSize(post_burn_MCMC_log$acceptRate))
    row.names(ESS_df) <- NULL

    cat("\nEffective sample sizes recorded in the MCMC log after removing burn-in:\n")
    print(ESS_df)
    cat("Ideally, ESS should be higher than 200. Increase the 'numberOfGenerations' or 'eventDataWriteFreq' if needed.\n\n")

    if (save_evaluations)
    {
      if (is.null(prefix_for_files))
      {
        utils::write.csv(x = ESS_df, file = file.path(paste0(BAMM_output_directory_path, "ESS_df.csv")), row.names = FALSE)
      } else {
        utils::write.csv(x = ESS_df, file = file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_ESS_df.csv")), row.names = FALSE)
      }
    }

    ## Compare prior and posterior distribution of LAMBDA (parameter controlling the nb of rates)
    # Display plot if requested
    if (plot_evaluations)
    {
      BAMMtools::plotPrior(mcmc = MCMC_log,
                           expectedNumberOfShifts = expectedNumberOfShifts,
                           burnin = burn_in,
                           main = paste0("Comparison prior/posterior distributions\n",
                                         "Number of shifts"))
    }
    # Save plot if requested
    if (save_evaluations)
    {
      if (is.null(prefix_for_files))
      {
        PP_lambda_path <- file.path(paste0(BAMM_output_directory_path, "PP_lambda_plot.pdf"))
      } else {
        PP_lambda_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_PP_lambda_plot.pdf"))
      }
      grDevices::pdf(file = file.path(PP_lambda_path),
                     width = 10, height = 8)

      BAMMtools::plotPrior(mcmc = MCMC_log,
                           expectedNumberOfShifts = expectedNumberOfShifts,
                           burnin = burn_in,
                           main = paste0("Comparison prior/posterior distributions\n",
                                         "of the number of shifts"))

      grDevices::dev.off()
    }
  }

  ##### import_BAMM() #####

  ## Build path to eventData file
  if (is.null(prefix_for_files))
  {
    eventData_path <- file.path(paste0(BAMM_output_directory_path, eventDataOutfile))
  } else {
    eventData_path <- file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", eventDataOutfile))
  }

  ## Create the bammdata summarizing BAMM outputs
  BAMM_data_output <- BAMMtools::getEventData(phy = phylo,
        eventdata = eventData_path,
        burnin = burn_in,
        type = "diversification")

  ## Select the subset of posterior samples

  # Get a subset of a selected number of posterior samples
  if (!is.null(seed))
  {
    set.seed(seed = seed)
  }
  sample_indices <- sample(x = 1:length(BAMM_data_output$eventData), size = nb_posterior_samples)
  BAMM_posterior_samples_data <- BAMMtools::subsetEventData(BAMM_data_output, index = sample_indices)

  ##### clean_BAMM_files() #####

  # If requested, remove all files generated
  if(!keep_BAMM_outputs)
  {
    # Remove phylo, config file and default priors file
    file.remove(phy_path)
    file.remove(config_file_path)
    file.remove(priors_path)

    # Remove files generated during the BAMM run
    if (is.null(prefix_for_files))
    {
      file.remove(file.path(paste0(BAMM_output_directory_path, runInfoFilename)))
      file.remove(file.path(paste0(BAMM_output_directory_path, mcmcOutfile)))
      file.remove(file.path(paste0(BAMM_output_directory_path, eventDataOutfile)))
      file.remove(file.path(paste0(BAMM_output_directory_path, chainSwapFileName)))

      if (outputAcceptanceInfo == 1)
      {
        file.remove(file.path(paste0(BAMM_output_directory_path, acceptanceInfoFileName)))
      }
    } else {
      file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", runInfoFilename)))
      file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", mcmcOutfile)))
      file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", eventDataOutfile)))
      file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", chainSwapFileName)))
      if (outputAcceptanceInfo == 1)
      {
        file.remove(file.path(paste0(BAMM_output_directory_path, prefix_for_files, "_", acceptanceInfoFileName)))
      }
    }

    # If empty, remove the BAMM_output_directory
    if (length(list.files(path = file.path(BAMM_output_directory_path))) == 0)
    {
      unlink(x = file.path(BAMM_output_directory_path), force = TRUE)
    }
  }

  ## Export BAMM object with posterior samples data
  return(invisible(BAMM_posterior_samples_data))


  ##### plot_BAMM_rates() #####

  # See if simply using BAMMtools::plot.bammdata() or a custom function to reuse to plot BAMM rates/regimes at every 'focal_time'?

  # See if adding regime shifts. From what?

  # Do not include it in prepare_diversification_data !
  # Just make it a function to use after!

  ### Split into sub-functions !!!
  # Export sub-functions as they may be useful to may people!

}


## Ignore BAMM directorys during compilation
# usethis::use_build_ignore("software/*")
# usethis::use_build_ignore("BAMM_outputs/*")

