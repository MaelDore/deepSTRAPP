## Scripts to prepare datasets made available in the package

### 1/ Generate Ponerinae_trait_data_10My ####

# Load phylogeny
data(Ponerinae_trait_data, package = "deepSTRAPP")
# Load trait df
data(Ponerinae_tree, package = "deepSTRAPP")

# Extract log(head with) data
Ponerinae_data_ln_HW <- setNames(object = Ponerinae_trait_data$sim_ln_HW,
                                 nm = Ponerinae_trait_data$Taxa)

# Get Ancestral Character Estimates based on a Brownian Motion model
# To obtain values at internal nodes
Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree, x = Ponerinae_data_ln_HW)

# Run a Stochastic Mapping based on a Brownian Motion model
# to interpolate values along branches and obtain a "contMap" object
Ponerinae_contMap <- phytools::contMap(Ponerinae_tree, x = Ponerinae_data_ln_HW,
                                       res = 100, # Number of time steps
                                       plot = FALSE)

# Set focal time to 10 Mya
focal_time <- 10

## Extract trait data and update contMap for the given focal_time

# Extract from tip data and ML estimates of ancestral characters (values are true ML)
Ponerinae_trait_data_10My <- extract_most_likely_trait_values_for_focal_time(
  contMap = Ponerinae_contMap,
  ace = Ponerinae_ACE, tip_data = Ponerinae_data_ln_HW,
  focal_time = focal_time,
  update_contMap = TRUE,
  keep_tip_labels = TRUE)

usethis::use_data(Ponerinae_trait_data_10My, overwrite = TRUE)


### 2/ Generate Ponerinae_BAMM_object_10My ####

## Load the BAMM_object summarizing 1000 posterior samples of BAMM.
data(Ponerinae_BAMM_object, package = "deepSTRAPP")

## Set focal-time to 10 My
focal_time = 10

## Update the BAMM object (May take several minutes to run)
Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
  BAMM_object = Ponerinae_BAMM_object,
  focal_time = focal_time,
  update_rates = TRUE, update_regimes = TRUE,
  update_tree = TRUE, update_plot = TRUE,
  update_all_elements = TRUE,
  keep_tip_labels = TRUE,
  verbose = TRUE)

usethis::use_data(Ponerinae_BAMM_object_10My, overwrite = TRUE)


### 3/ Generate eel_BAMM_object ####

library(stringr)
library(BAMMtools)

##### 1/ Load input files ####

### 1.1/ Inform path to BAMM folder ####

BAMM_path <- "./data-raw/bamm-2.5.0/"

### 1.2/ Load phylogeny and check if it is valid

# Load eel phylogeny from the R package phytools
# Source: Collar et al., 2014; DOI: 10.1038/ncomms6505
library(phytools)
data(eel.tree)

data(package = "phytools")

data(whale.tree)

data(butterfly.data)

plot(butterfly.tree)

whale_phylogeny <- whale.tree

# Check if ultrametric
is.ultrametric(whale_phylogeny )
# Check if fully resolved
is.binary(whale_phylogeny)
# Check that all branch have positive length
min(whale_phylogeny $edge.length) > 0 ; min(whale_phylogeny $edge.length)

### 1.3/ Export phylogeny as .tree file for the analyses
write.tree(phy = whale_phylogeny, file = paste0(BAMM_path, "whale_phylogeny.tree"))

phy_path <- paste0(BAMM_path, "whale_phylogeny.tree")

##### 2/ Configure the control file for the run #####

### 2.1/ Load control file template for diversification analyses ####

div_config_file_template <- readLines(con = paste0(BAMM_path, "template_diversification.txt"))
my_div_control_file <- div_config_file_template

### 2.2/ Set general settings and data input ####

# Set path to the phylogenetic tree file
phy_path_line <- which(str_detect(string = my_div_control_file, pattern = "treefile ="))
my_div_control_file[phy_path_line] <- paste0("treefile = ", phy_path)

# Should the output files be overwritten?
# If True (1), the program will overwrite any output files in the current directory (if present)
# overwrite <- 0
overwrite <- 1
overwrite_line <- which(str_detect(string = my_div_control_file, pattern = "overwrite = "))[1]
my_div_control_file[overwrite_line] <- paste0("overwrite = ", overwrite)

# Set limits to valid configurations
# If 1, rejects proposals that cause a branch and both of its direct descendants to have at least one event.
# Such an event configuration may cause the parameters of the parent event to change to unrealistic values.
# If 0, no such proposals are immediately rejected. The default value is 0.
# validateEventConfiguration <- 0
validateEventConfiguration <- 1
validateEventConfiguration_line <- which(str_detect(string = my_div_control_file, pattern = "validateEventConfiguration = "))[1]
my_div_control_file[validateEventConfiguration_line] <- paste0("validateEventConfiguration = ", validateEventConfiguration)


### 2.3/ Set (hyper)prior settings ####

## Can use this help function to automatically tune prior adapted to your data by scaling the prior distributions based on the age (root depth) of your tree
# In practice, setBAMMpriors first estimates the rate of speciation for your full tree under a pure birth model of diversification.
# Then assume, arbitrarily, that a reasonable prior distribution for the initial lambda0/mu0 rate parameters is an exponential distribution with a mean five times greater than this pure birth value.
# Rationale = having a weakly informative prior that is still in the order of magnitude of the true rate
# For the shift parameter (alpha), the sd of the normal prior is set such as mean +/- 2s gives an alpha parameter that results in
# either a 90% decline in the evolutionary rate or a 190% increase in rate on the interval of time from the root to the tips of the tree.
setBAMMpriors(read.tree(phy_path))
default_tuned_priors <- readLines(con = "./myPriors.txt")
file.remove("./myPriors.txt")

# Set the expected number of shifts used to set the exponential hyperprior for nb of rate shifts (from which the Λ is drawn)
# Suggested values:
#  expectedNumberOfShifts = 1.0 for small trees (< 500 tips)
#	 expectedNumberOfShifts = 10 or even 50 for large trees (> 5000 tips)
# Good practice = set several runs with a range of values to be sure it does not affect results
# expectedNumberOfShifts <- 1.0
expectedNumberOfShifts_default_line <- which(str_detect(string = default_tuned_priors, pattern = "expectedNumberOfShifts = "))[1]
expectedNumberOfShifts <- as.numeric(str_remove(string = default_tuned_priors[expectedNumberOfShifts_default_line], pattern = "expectedNumberOfShifts = "))
expectedNumberOfShifts_line <- which(str_detect(string = my_div_control_file, pattern = "expectedNumberOfShifts = "))[1]
my_div_control_file[expectedNumberOfShifts_line] <- paste0("expectedNumberOfShifts = ", expectedNumberOfShifts)

# Set the rate parameter of the exponential prior(s) of initial lambda parameters (lambda0) of speciation rate regimes
# lambda0 in lambda(t) = lamba0 x exp(alpha*t)
# lambdaInitPrior <- 1.0
# lambdaInitPrior <- 5.79332
lambdaInitPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "lambdaInitPrior = "))[1]
lambdaInitPrior <- as.numeric(str_remove(string = default_tuned_priors[lambdaInitPrior_default_line], pattern = "lambdaInitPrior = "))
lambdaInitPrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaInitPrior = "))[1]
my_div_control_file[lambdaInitPrior_line] <- paste0("lambdaInitPrior = ", lambdaInitPrior)

# Set the standard deviation of the normal distribution prior(s) of rate variation parameters (alpha) of speciation rate regimes
# alpha in lambda(t) = lamba0 x exp(alpha*t)
# Mean of this prior(s) are fixed to zero such as a constant rate diversification process is the most probable a priori
# lambdaShiftPrior <- 0.05
# lambdaShiftPrior <- 0.011629
lambdaShiftPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "lambdaShiftPrior = "))[1]
lambdaShiftPrior <- as.numeric(str_remove(string = default_tuned_priors[lambdaShiftPrior_default_line], pattern = "lambdaShiftPrior = "))
lambdaShiftPrior_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaShiftPrior = "))[1]
my_div_control_file[lambdaShiftPrior_line] <- paste0("lambdaShiftPrior = ", lambdaShiftPrior)

# Set the rate parameter of the exponential prior(s) of initial lambda parameters (mu0) of extinction rate regimes
# mu0 in mu(t) = mu0 x exp(alpha*t)
# As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
# muInitPrior <- 1.0
# muInitPrior <- 5.793324
muInitPrior_default_line <- which(str_detect(string = default_tuned_priors, pattern = "muInitPrior = "))[1]
muInitPrior <- as.numeric(str_remove(string = default_tuned_priors[muInitPrior_default_line], pattern = "muInitPrior = "))
muInitPrior_line <- which(str_detect(string = my_div_control_file, pattern = "muInitPrior = "))[1]
my_div_control_file[muInitPrior_line] <- paste0("muInitPrior = ", muInitPrior)


### 2.4/ Set the MCMC simulation settings, MCMC logs and output options ####

# Set the number of generations to perform MCMC simulation
# numberOfGenerations = format(10000000, scientific = F) # 10^7
numberOfGenerations = format(1000000, scientific = F) # For the test run: 10^6
# numberOfGenerations = 1000 # 10^3
numberOfGenerations_line <- which(str_detect(string = my_div_control_file, pattern = "numberOfGenerations = "))[1]
my_div_control_file[numberOfGenerations_line] <- paste0("numberOfGenerations = ", numberOfGenerations)

# Set the frequency in which to write the MCMC output to the log file
# Aim for 500-5000 posterior samples ideally
# Will need to remove some to account for the burn-in
mcmcWriteFreq <- round(as.numeric(numberOfGenerations) / 2000)
# mcmcWriteFreq <- 500
mcmcWriteFreq_line <- which(str_detect(string = my_div_control_file, pattern = "mcmcWriteFreq = "))[1]
my_div_control_file[mcmcWriteFreq_line] <- paste0("mcmcWriteFreq = ", mcmcWriteFreq)

# Set frequency in which to write the event data to the output file = the sampling frequency of posterior samples
# Aim for 500-5000 posterior samples ideally
# Will need to remove some to account for the burn-in
eventDataWriteFreq <- as.numeric(numberOfGenerations) / 2000 # Sample every 500 generations
# eventDataWriteFreq <- 100 # Sample every 100 generations
eventDataWriteFreq_line <- which(str_detect(string = my_div_control_file, pattern = "eventDataWriteFreq = "))[1]
my_div_control_file[eventDataWriteFreq_line] <- paste0("eventDataWriteFreq = ", eventDataWriteFreq)

# Set frequency in which to print MCMC status to the screen
printFreq <- 1000 # Print status every 10^3 generations for short test runs
# printFreq <- 10000 # Print status every 10^4 generations for long runs
# printFreq <- 100 # Print status every 100 generations
printFreq_line <- which(str_detect(string = my_div_control_file, pattern = "printFreq = "))[1]
my_div_control_file[printFreq_line] <- paste0("printFreq = ", printFreq)

# Set prefix to add to all output files (separated with "_")
# Need to be commented out to do not add any prefix
run_prefix_for_output_files <- "BAMM_whale"
outName_line <- which(str_detect(string = my_div_control_file, pattern = "outName = "))[1]
my_div_control_file[outName_line] <- paste0("outName = ", run_prefix_for_output_files)
# my_div_control_file[outName_line] <- paste0("# outName = ", run_prefix_for_output_files) # Commented version to add any prefix


### 2.7/ Set the initial parameter values to start the MCMC chain(s) ####

# The MCMC chain start with a model with no shift (M0 submodel). So the initial parameter values are for this unique regime
# (But probably also for the new regime if added?)

# Run a BD model to obtain credible starting value for the root process
BD_fit <- phytools::fit.bd(tree = whale_phylogeny) ; BD_fit

# Set the initial speciation rate (lambda0) for the first regime starting at the root of the tree (regime 0)
# lambda0 in lambda(t) = lamba0 x exp(alpha*t)
lambdaInit0 <- 0.02974
lambdaInit0 <- BD_fit$b
lambdaInit0_line <- which(str_detect(string = my_div_control_file, pattern = "lambdaInit0 = "))[1]
my_div_control_file[lambdaInit0_line] <- paste0("lambdaInit0 = ", lambdaInit0)

# Set the intial extinction rate (mu0) for the first regime starting at the root of the tree (regime 0)
# mu0 in mu(t) = mu0 x exp(alpha*t)
# As the extinction rates are actually assumed to follow constant rates, alpha is set to 0, thus mu(t) = mu0 and these are constant extinction rates
muInit0 <- 0.005
muInit0 <- max(0.005, BD_fit$d)
muInit0_line <- which(str_detect(string = my_div_control_file, pattern = "muInit0 = "))[1]
my_div_control_file[muInit0_line] <- paste0("muInit0 = ", muInit0)


### 2.10/ Export the updated custom control file ####
writeLines(text = my_div_control_file, con = paste0(BAMM_path, "BAMM_whale_my_div_control_file.txt"))


##### 3/ Run BAMM with calls to command lines from within r using system() #####

### 3.1/ Run BAMM ####

?system

# Version
system(paste0(BAMM_path,"bamm --version"))

# Test run
system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_whale_my_div_control_file.txt"))

# # Full run
# # system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_my_div_control_file.txt"))
# system(paste0(BAMM_path, "bamm -c ",BAMM_path, "BAMM_Ponerinae_MCC_my_div_control_file.txt"))

### Outputs
# The run_info.txt file, containing a summary of your parameters/settings
# An mcmc_log.txt containing raw MCMC information useful in diagnosing convergence
# An event_data.txt file containing all of evolutionary rate parameters and their topological mappings
# A chain_swap.txt file containing data about each chain swap proposal (when a proposal occurred, which chains might be swapped, and whether the swap was accepted).

### 3.2/ Clean outputs (move them to a dedicated folder) ####

BAMM_output_folder_path <- "./data-raw/bamm-2.5.0/whale_outputs/"
# Create the folder if not existing
if (!file.exists(file.path(BAMM_output_folder_path))) { dir.create(file.path(BAMM_output_folder_path)) }

# Detect output files
output_files_path <- list.files(path = "./", pattern = "BAMM_")
# Move output files to dedicated folder
file.rename(from = paste0("./",output_files_path), to = paste0(BAMM_output_folder_path, output_files_path)) # BAMM output files
file.copy(from = phy_path, to = paste0(BAMM_output_folder_path, "whale_phylogeny.tree")) # Phylo file
# file.rename(from = paste0(BAMM_path, "my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "my_div_control_file.txt")) # Control file
file.rename(from = paste0(BAMM_path, "BAMM_whale_my_div_control_file.txt"), to = paste0(BAMM_output_folder_path, "BAMM_whale_my_div_control_file.txt")) # Control file


## Choose the expectedNumberofShifts the closest to the overall median
selected_expectedNumberofShifts <- expectedNumberOfShifts

# Create the bammdata summarizing BAMM outputs
whale_BAMM_object <- getEventData(phy = whale_phylogeny,
                                 eventdata = paste0(BAMM_output_folder_path, run_prefix_for_output_files, "_event_data.txt"),
                                 burnin = 0.1,
                                 type = "diversification")

### 5.2/ Select the subset of posterior samples ####

# Set the targeted number of posterior samples
nb_samples <- 1000

# Get a subset of a selected number of posterior samples
set.seed(seed = 1234)
sample_indices <- sample(x = 1:length(BAMM_data_output$eventData), size = nb_samples)
whale_BAMM_object <- subsetEventData(whale_BAMM_object, index = sample_indices)

str(whale_BAMM_object, 1)

summary(whale_BAMM_object$numberEvents)

usethis::use_data(whale_BAMM_object, overwrite = TRUE)

