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

str(Ponerinae_BAMM_object_10My, 1)

usethis::use_data(Ponerinae_BAMM_object_10My, overwrite = TRUE)


### 3/ Generate Ponerinae_BAMM_object_25My ####

## Load the BAMM_object summarizing 1000 posterior samples of BAMM.
data(Ponerinae_BAMM_object, package = "deepSTRAPP")

## Set focal-time to 25 My
focal_time = 25

## Update the BAMM object (May take several minutes to run)
Ponerinae_BAMM_object_25My <- update_rates_and_regimes_for_focal_time(
  BAMM_object = Ponerinae_BAMM_object,
  focal_time = focal_time,
  update_rates = TRUE, update_regimes = TRUE,
  update_tree = TRUE, update_plot = TRUE,
  update_all_elements = TRUE,
  keep_tip_labels = TRUE,
  verbose = TRUE)

str(Ponerinae_BAMM_object_25My, 1)

usethis::use_data(Ponerinae_BAMM_object_25My, overwrite = TRUE)


### 4/ Generate whale_BAMM_object ####

library(phytools)
data(whale.tree)

whale_BAMM_object <- prepare_diversification_data(
  BAMM_install_directory_path = "./software/bamm-2.5.0/",
  phylo = whale.tree,
  prefix_for_files = "whale",
  numberOfGenerations = 100000) # Set low for the example

# Load directly the result
data(whale_BAMM_object)

# Explore output
str(whale_BAMM_object, 1)
summary(whale_BAMM_object$numberEvents)

# Plot mean net diversification rates on the phylogeny
plot.bammdata(whale_BAMM_object, labels = TRUE)

usethis::use_data(whale_BAMM_object, overwrite = TRUE)

# ### 5/ Generate whale_BAMM_object_5My ####
#
# ## Load the BAMM_object summarizing 1000 posterior samples of BAMM.
# data(whale_BAMM_object, package = "deepSTRAPP")
#
# ## Set focal-time to 5 My
# focal_time = 5
#
# ## Update the BAMM object (May take several minutes to run)
# whale_BAMM_object_5My <- update_rates_and_regimes_for_focal_time(
#   BAMM_object = whale_BAMM_object,
#   focal_time = focal_time,
#   update_rates = TRUE, update_regimes = TRUE,
#   update_tree = TRUE, update_plot = TRUE,
#   update_all_elements = TRUE,
#   keep_tip_labels = TRUE,
#   verbose = TRUE)
#
# str(whale_BAMM_object_5My, 1)
#
# usethis::use_data(whale_BAMM_object_5My, overwrite = TRUE)


### 6/ Generate Ponerinae_deepSTRAPP_0_40 ####

## deepSTRAPP output for Ponerinae over time from 0 to 40My (steps = 10 My)

## Load trait df
data(Ponerinae_trait_data, package = "deepSTRAPP")

dim(Ponerinae_trait_data)
View(Ponerinae_trait_data)

# Extract continuous trait data as a named vector
Ponerinae_data_ln_HW <- setNames(object = Ponerinae_trait_data$sim_ln_HW,
                                 nm = Ponerinae_trait_data$Taxa)

## Load phylogeny
data(Ponerinae_tree, package = "deepSTRAPP")

## Run prepare_trait_data with default options
# For continuous trait, a BM model is assumed by default.
Ponerinae_trait_object <- prepare_trait_data(tip_data = Ponerinae_data_ln_HW,
                                             trait_data_type = "continuous",
                                             phylo = Ponerinae_tree)

# Explore output
str(Ponerinae_trait_object, 1)

# Extract the contMap representing continuous trait evolution on the phylogeny
Ponerinae_contMap <- Ponerinae_trait_object$contMap
plot(Ponerinae_contMap)

# Extract the Ancestral Character Estimates (ACE) = trait values at nodes
Ponerinae_ACE <- Ponerinae_trait_object$ace
head(Ponerinae_ACE)

## Set for five time steps of 10 My. Will generate deepSTRAPP workflows for 0, 10, 20, 30, and 40 Mya.
nb_time_steps <- 5
time_step_duration <- 10

# Run deepSTRAPP on net diversification rates
## This step is time-consuming. You can skip it and load directly the result if needed
Ponerinae_deepSTRAPP_0_40 <- run_deepSTRAPP_over_time(
  contMap = Ponerinae_contMap,
  ace = Ponerinae_ACE,
  tip_data = Ponerinae_data_ln_HW,
  trait_data_type = "continuous",
  BAMM_object = Ponerinae_BAMM_object,
  nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_trait_data_with_contMap = TRUE,
  return_updated_BAMM_object = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Ponerinae_deepSTRAPP_0_40, max.level = 1)

usethis::use_data(Ponerinae_deepSTRAPP_0_40, overwrite = TRUE)


### 7/ Save diversification template file as .rds #####

# Load it from the /inst/ directory (works only in source/bundled version of the package)
BAMM_template_diversification <- readLines(con = file.path("./inst/BAMM_template_diversification.txt"))

usethis::use_data(BAMM_template_diversification, overwrite = TRUE)

