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


### 8/ Generate categorical trait evolution data for eel using 3-level factor #####

# Load phylogeny and tip data
library(phytools)
data(eel.tree)
data(eel.data)

# Transform feeding mode data into a 3-level factor
eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
eel_data <- as.character(eel_data)
eel_data[c(1, 5, 6, 7, 10, 11, 15, 16, 17, 24, 25, 28, 30, 51, 52, 53, 55, 58, 60)] <- "kiss"
eel_data <- stats::setNames(eel_data, rownames(eel.data))
table(eel_data)

# Manually define a Q_matrix for rate classes of state transition to use in the 'matrix' model
# Does not allow transitions from state 1 ("bite") to state 2 ("kiss") or state 3 ("suction")
# Does not allow transitions from state 3 ("suction") to state 1 ("bite")
# Set symmetrical rates between state 2 ("kiss") and state 3 ("suction")
Q_matrix = rbind(c(NA, 0, 0), c(1, NA, 2), c(0, 2, NA))

# Set colors per states
colors_per_states <- c("limegreen", "orange", "dodgerblue")
names(colors_per_states) <- c("bite", "kiss", "suction")

eel_cat_data <- prepare_trait_data(tip_data = eel_data, phylo = eel.tree,
                                        trait_data_type = "categorical",
                                        colors_per_states = colors_per_states,
                                        evolutionary_models = c("ER", "SYM", "ARD", "meristic", "matrix"),
                                        Q_matrix = Q_matrix,
                                        nb_simulations = 1000, # Set to 10 to save time.
                                        # But recommended value = 1000.
                                        plot_map = TRUE,
                                        plot_overlay = TRUE,
                                        return_best_model_fit = TRUE,
                                        return_model_selection_df = TRUE)

# Explore output
plot(eel_cat_data$densityMaps[[1]]) # densityMap for state n°1 ("bite")
eel_cat_data$model_selection_df # Summary of model selection
# Parameter estimates and optimization summary of the best model
# (Here, the best model is ER)
print(eel_cat_data$best_model_fit)$ # Summary of the best evolutionary model
eel_cat_data$ace # Posterior probabilities of each state (= ACE) at internal nodes

usethis::use_data(eel_cat_data, overwrite = TRUE)


### 9/ Generate categorical trait evolution data for Ponerinae ants using 3-level factor #####

## Load data

# Load phylogeny
data(Ponerinae_trait_data, package = "deepSTRAPP")
# Load trait df
data(Ponerinae_tree, package = "deepSTRAPP")
# Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(Ponerinae_BAMM_object, package = "deepSTRAPP")

## Prepare trait data

# Extract log(head with) data
Ponerinae_data_ln_HW <- setNames(object = Ponerinae_trait_data$sim_ln_HW,
                                 nm = Ponerinae_trait_data$Taxa)
# Convert to three-factor categorical traits
Ponerinae_data <- Ponerinae_data_ln_HW
Ponerinae_data[seq_along(Ponerinae_data)] <- "small"
Ponerinae_data[Ponerinae_data_ln_HW > -1] <- "medium"
Ponerinae_data[Ponerinae_data_ln_HW > 0] <- "large"
table(Ponerinae_data)

# Select color scheme for states
colors_per_states <- c("darkblue", "dodgerblue", "lightblue")
names(colors_per_states) <- c("large", "medium", "small")

## Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
Ponerinae_cat_data <- prepare_trait_data(tip_data = Ponerinae_data, phylo = Ponerinae_tree,
                                         trait_data_type = "categorical",
                                         colors_per_states = colors_per_states,
                                         evolutionary_models = c("ER", "SYM", "ARD", "meristic"),
                                         nb_simulations = 1000,
                                         return_best_model_fit = TRUE,
                                         return_model_selection_df = TRUE,
                                         plot_map = FALSE)

# Explore output
plot(Ponerinae_cat_data$densityMaps[[1]]) # densityMap for state n°1 ("large")
Ponerinae_cat_data$model_selection_df # Summary of model selection
# Parameter estimates and optimization summary of the best model
# (Here, the best model is ARD)
print(Ponerinae_cat_data$best_model_fit) # Summary of the best evolutionary model
Ponerinae_cat_data$ace # Posterior probabilities of each state (= ACE) at internal nodes

usethis::use_data(Ponerinae_cat_data, overwrite = TRUE)

