#
#
# # See also BAMMtools::traitDependentBAMM
#
# # State in the description that this is a wrapped up of the original function traitDependentBAMM written by Dan Rabosky & Huateng Huang, 2014 in the BAMMtools package
# # Describe changes from the original
#   # Additions
#     # Use the names of the tipStates, tipLambda, and tipMu to check compatibility between tips rather than the tip.label in the phylogeny
#     # Allows to use the function to compute STRAPP test for any given time
#     # Allows to choose if random sampling of posterior configurations must be done with replacement or not.
#     # Default to using one permutation per posterior sample if not requested otherwise.
#     # Provide outputs tailored for histogram plots [ref_fun] and p-value time-series plots [ref_fun]
#     # Add verbose parameter to provide insights on progress
#   # Changes
#     # Function split in multiple sub-functions according to $data_type for clarity
#     # Prevent using Pearson's correlation tests and applying log-transformation as their is no rationale
#     # to assume that tip rates are distributed normally or log-normally. Favor Spearman's rank correlation tests.
#
#
#
# ### Inputs
# # From extract_most_likely_trait_values_for_focal_time() => trait_data named integer (may also be a list with two items, trait_data + a contMap)
# # From update_rates_and_regimes_for_focal_time() => Updated bammdata list with $tipStates, $tipLambda and $tipMu holding data from all posterior samples
#
# # Make a master function run_STRAPP_test()
#   # Prepare the input data and select the proper sub-function to be ran
# # Make sub-function for each type of data
#   # run_STRAPP_test_for_continuous_data()
#   # run_STRAPP_test_for_binary_data()
#      # Special case of categorical/biogeographic data with only two states
#   # run_STRAPP_test_for_categorical_data()
#      # Case of categorical/biogeographic data with more than two states
#      # Include biogeographic data as the tests are the same
#      # Can define the pairwise tests to run/plot
#
# # Make a wrapped-up to be run over time-series
#   # run_STRAPP_test_over_time()
#
# # Make a unique help documentation for all functions?
#
# ### Outputs
# # focal_time, P-value, Q5% Δ stat
# # Info needed to plot histogram
#   # Focal-time
#   # Melted STRAPP data for this focal_time
#      # BAMM ID, permut stat, observed stat, Δ stat
#   # Will be used to plot a single histogram for a given focal_time
#     # Will compute summary stats from melted df
#     # Include option of a PDF output
#   # Will be merged in a melted data.frame encompassing all focal_time to produce sets of histogram
#     # Include option of a merged PDF output with one histogram per page
#
#
# # ## Load the BAMM_object summarizing 1000 posterior samples of BAMM
# # data(Ponerinae_BAMM_object, package = "deepSTRAPP")
# #
# # ## Set focal-time to 10 My
# # focal_time = 10
# #
# # ## Update the BAMM object (May take several minutes to run)
# # Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
# #   BAMM_object = Ponerinae_BAMM_object,
# #   focal_time = focal_time,
# #   update_rates = TRUE, update_regimes = TRUE,
# #   update_tree = TRUE, update_plot = TRUE,
# #   update_all_elements = TRUE,
# #   keep_tip_labels = TRUE,
# #   verbose = TRUE)
#
# Ponerinae_BAMM_object_10My <- readRDS(file = "./src/Ponerinae_BAMM_object_10My.rds")
#
# BAMM_object <- Ponerinae_BAMM_object_10My
#
# run_STRAPP_test <- function (BAMM_object, trait_data_list,
#                              nb_permutations = NULL, # If NULL use all posterior samples by default
#                              rate = "net_diversification",# Can use speciation, extinction and net diversification
#                              return_test_data = FALSE, # return.full in initial function # To provide the Melted STRAPP data for this focal_time needed to plot an histogram
#                              two.tailed = TRUE,
#                              replace = FALSE, # To allow or not multiple uses of posterior samples in random permutations
#                              traitorder = NA, # For one-tailed tests. To select the hypothesis (negative/positive correlation or which state with greater rates). Inherited from BAMMtools::traitDependentBAMM
#                              nthreads = 1, # Number of threads to use for parallelization of the function. The R package parallel must be loaded for nthreads > 1.
#                              verbose = TRUE)
# {
#
#   ### Check input validity
#
#   # BAMM_object must be a 'bammdata' object
#   # Number of posterior sample data must be equal between $tipStates, $tipLambda and $tipMu
#
#   # trait_data_list must be a list with $trait_data and $data_type
#   # $trait_data type must match with $data_type
#      # Numerical if $data_type = "continuous"
#      # Character string if $data_type = "categorical" or "biogeographic"
#
#   # Length of $trait_data should match length of $tipStates, $tipLambda and $tipMu (for each posterior sample)
#
#   # Names of $trait_data should match names in $tipStates, $tipLambda and $tipMu (for each posterior sample)
#     # If need to be reordered, do it, but send a warning
#
#   # If nb_permutations is higher than number of posterior samples (length of $tipStates, $tipLambda and $tipMu) AND replace = FALSE,
#   # Send an error to say that replace should be set to TRUE to allow multiple samplings of posterior in order to reach the requested number of permutations
#
#   ## See other validation checks from BAMMtools::traitDependentBAMM
#
#   ## Extract BAMM rates data
#
#   ### Use an exported function to build melted data.frame from BAMM object with diversification data
#   extract_diversification_data_melted_df_for_focal_time(BAMM_object, verbose = TRUE)
#     # Tip ID, BAMM ID, regime, rate
#
#
#   # Prepare the input data and select the proper sub-function to be ran
#
#   ##
# }
#
# ## 12.3.1/ Get STRAPP data ####
#
# # Get tip regimes and tip rates for focal time i
# BAMM_data_updated_i <- update_tipStates_and_tipRates_for_focal_time(BAMM_object = BAMM_posterior_samples_data, time = time_i, update_rates = T, verbose = T)
#
# # Extract only $tipStates, $tipLambda and $tipMu in the loop per time to save place
# BAMM_data_updated_i <- BAMM_data_updated_i[c("tipStates", "tipLambda", "tipMu", "type", "tip.label")]
#
# # Get tip states for focal time i
# bioregion_data_OW_vs_NW_i <- get_most_likely_binary_states_for_focal_time(density_map = DEC_J_density_map_OW_vs_NW, time = time_i)
#
# ## Store input data for STRAPP test
# STRAPP_data_OW_vs_NW_all_time[[i]] <- list(BAMM_data = BAMM_data_updated_i,
#                                            Bioregion_data_OW_vs_NW = bioregion_data_OW_vs_NW_i)
#
# # Save input data for STRAPP test
# # saveRDS(STRAPP_data_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_data_OW_vs_NW_all_time.rds")
# saveRDS(STRAPP_data_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_data_OW_vs_NW_all_time.rds")
#
# ## 12.3.2/ Run STRAPP test ####
#
# if ((is.null(bioregion_data_OW_vs_NW_i) | (any(!c("Old World", "New World") %in% bioregion_data_OW_vs_NW_i))))
# { # If not edge data, provide NA
#
#   STRAPP_test_OW_vs_NW_i <- list()
#   STRAPP_test_OW_vs_NW_i$estimate <- NA
#   STRAPP_test_OW_vs_NW_i$p.value <- NA
#   STRAPP_test_OW_vs_NW_i$obs.corr <- NA
#   STRAPP_test_OW_vs_NW_i$gen <- NA
#   STRAPP_test_OW_vs_NW_i$null <- NA
#   STRAPP_test_OW_vs_NW_i$test_stats <- NA
#   STRAPP_test_OW_vs_NW_i$stat_median <- NA
#   STRAPP_test_OW_vs_NW_i$stat_Q5 <- NA
#
# } else {
#
#   # Run STRAPP test if edge data is present
#   STRAPP_test_OW_vs_NW_i <- run_STRAPP_test(BAMM_data = BAMM_data_updated_i,
#                                             trait_data = bioregion_data_OW_vs_NW_i,
#                                             reps = 1000, rate = "net diversification",
#                                             return.full = T,
#                                             method = "mann-whitney", # For categorical binomial data (G = 2)
#                                             logrates = F, # Do not use log for Mann-Whitney as it is a rank test, so there is no need, and it prevents removing the negative rates
#                                             two.tailed = T,
#                                             replace = FALSE,
#                                             nthreads = 1)
#
#
#
#
#   # # Explore output
#   # STRAPP_test_OW_vs_NW_i$estimate # Mean tip rates per categories
#   # STRAPP_test_OW_vs_NW_i$p.value
#   # STRAPP_test_OW_vs_NW_i$obs.corr # Observed statistic for each posterior sample
#   # STRAPP_test_OW_vs_NW_i$gen # Generation ID of the selected posterior sample
#   # STRAPP_test_OW_vs_NW_i$null # Null statistic for each posterior sample
#
#   # Compute info for histogram
#   STRAPP_test_OW_vs_NW_i$test_stats <- STRAPP_test_OW_vs_NW_i$null - STRAPP_test_OW_vs_NW_i$obs.corr
#   # summary(STRAPP_test_OW_vs_NW_i$test_stats)
#   # table(STRAPP_test_OW_vs_NW_i$test_stats > 0)
#   STRAPP_test_OW_vs_NW_i$stat_median <- median(STRAPP_test_OW_vs_NW_i$test_stats)
#   STRAPP_test_OW_vs_NW_i$stat_Q5 <- quantile(STRAPP_test_OW_vs_NW_i$test_stats, p = 0.05)
# }
#
# # Store focal time information
# STRAPP_test_OW_vs_NW_i$focal_time <- time_i
#
# ## Store test output
# STRAPP_tests_OW_vs_NW_all_time[[i]] <- STRAPP_test_OW_vs_NW_i
#
# # Save STRAPP test outputs
# # saveRDS(STRAPP_tests_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_rough_phylogeny_1534t/STRAPP_tests_OW_vs_NW_all_time.rds")
# saveRDS(STRAPP_tests_OW_vs_NW_all_time, file = "./outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/STRAPP_tests_OW_vs_NW_all_time.rds")
