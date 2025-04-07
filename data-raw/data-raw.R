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


### 3/ Generate Ponerinae_BAMM_object_25My ####

## Load the BAMM_object summarizing 1000 posterior samples of BAMM.
data(Ponerinae_BAMM_object, package = "deepSTRAPP")

## Set focal-time to 25 My
focal_time = 25

## Update the BAMM object (May take several minutes to run)
Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
  BAMM_object = Ponerinae_BAMM_object,
  focal_time = focal_time,
  update_rates = TRUE, update_regimes = TRUE,
  update_tree = TRUE, update_plot = TRUE,
  update_all_elements = TRUE,
  keep_tip_labels = TRUE,
  verbose = TRUE)

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


##### 5/ Save diversification template file as .rds #####

# Load it from the /inst/ directory (works only in source/bundled version of the package)
BAMM_template_diversification <- readLines(con = file.path("./inst/BAMM_template_diversification.txt"))

usethis::use_data(BAMM_template_diversification, overwrite = TRUE)

