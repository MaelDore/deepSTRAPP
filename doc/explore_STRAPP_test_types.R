## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  eval = FALSE, # Chunks of codes will not be evaluated by default
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, fig.height = 5   # Set device size at rendering time (when plots are generated)
)

## ----setup, eval = TRUE, include = FALSE--------------------------------------
library(deepSTRAPP)

is_dev_version <- function (pkg = "deepSTRAPP")
{
  # # Check if ran on CRAN
  # not_cran <- identical(Sys.getenv("NOT_CRAN"), "true") # || interactive()

  # Version number check
  version <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "")
  dev_version <- grepl("\\.9000", version)

  # not_cran || dev_version
  
  return(dev_version)
}


## ----STRAPP_tests_cont_two_tailed---------------------------------------------
# # ------ Example 1: Continuous trait data ------ #
# 
# #### Step 1: Prepare data ####
# 
# ## Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# ## Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# 
# ## Load BAMM diversification outputs
# data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")
# # This dataset is only available in development versions installed from GitHub.
# # It is not available in CRAN versions.
# # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
# 
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# 
# # This not valid biological data. For the sake of example, we will assume this is size data.
# 
# # Reorder tip data as in phylogeny
# Ponerinae_cont_tip_data <- Ponerinae_cont_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# # Run prepare_trait_data with default options
# # For continuous trait, a BM model is assumed by default.
# Ponerinae_trait_object <- prepare_trait_data(tip_data = Ponerinae_cont_tip_data,
#                                              trait_data_type = "continuous",
#                                              phylo = Ponerinae_tree_old_calib,
#                                              evolutionary_model = "BM",
#                                              seed = 1234) # Set seed for reproducibility
# 
# 
# #### Step 2: Run a two-tailed test ####
# 
# # No assumption is made on the direction of the effect tested.
# # The null hypothesis is that no correlation can be found between rates and trait values.
# # The alternative hypothesis is that there is a correlation between rates and trait values,
# # in any direction (i.e., "positive" or "negative").
# 
# # Here, we run a two-tailed STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_two_tailed <- run_deepSTRAPP_for_focal_time(
#    # Provide ML trait estimates as a contMap
#    contMap = Ponerinae_trait_object$contMap,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cont_tip_data,
#    trait_data_type = "continuous",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    # Deal with uncertainty in estimates by combining trait ML estimates
#    # with all BAMM posterior samples
#    uncertainty_strategy = "rates_only",
#    two_tailed = TRUE, # Select a two-tailed test (Default)
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    # Needed to access trait and rate data to evaluate the direction of the correlation
#    extract_trait_data_melted_df = TRUE,
#    extract_diversification_data_melted_df = TRUE)
# 
# ## Explore output
# str(deepSTRAPP_output_two_tailed, max.level = 1)
# 
# # Access deepSTRAPP results
# deepSTRAPP_output_two_tailed$STRAPP_results[1:5]
# 
# # We obtain a p-value = 0.018 leading us to discard the null hypothesis (no correlation).
# # The estimate is the quantile 5% of the distribution of differences in
# # absolute Spearman rho-stats between observed and permuted data.
# # The test is significant as the estimate is higher than zero.
# # However, without plotting the rates vs. trait data, we do not know the direction of this correlation.
# 
# # Plot histogram of Spearman sum-rank test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_two_tailed)
# 
# # The test is significant as the red line (quantile 5% = 0.066) is above the null expectation
# # that Δ abs(Spearman rho-stats) = 0.

## ----STRAPP_tests_cont_two_tailed_eval_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# ## Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# ## Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# ## Load BAMM diversification outputs
# data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# # Reorder tip data as in phylogeny
# Ponerinae_cont_tip_data <- Ponerinae_cont_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# # Run prepare_trait_data with default options
# Ponerinae_trait_object <- suppressWarnings(prepare_trait_data(
#   tip_data = Ponerinae_cont_tip_data,
#   trait_data_type = "continuous",
#   phylo = Ponerinae_tree_old_calib,
#   evolutionary_model = "BM",
#   seed = 1234,# Set seed for reproducibility
#   plot_map = FALSE,
#   verbose = FALSE))
# 
# # Here, we run a two-tailed STRAPP test for t = 10 Mya.
# deepSTRAPP_output_two_tailed <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    contMap = Ponerinae_trait_object$contMap,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cont_tip_data,
#    trait_data_type = "continuous",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    two_tailed = TRUE, # Select a two-tailed test (Default)
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    # Needed to access trait and rate data to evaluate the direction of the correlation
#    extract_trait_data_melted_df = TRUE,
#    extract_diversification_data_melted_df = TRUE,
#    verbose = FALSE))
# 
# # Plot histogram of Spearman sum-rank test stats
# plot_histogram_STRAPP_test_for_focal_time(
#    deepSTRAPP_outputs = deepSTRAPP_output_two_tailed)

## ----STRAPP_tests_cont_two_tailed_eval_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_1.2_Example_continuous_two_tailed_1.PNG")


## ----STRAPP_tests_cont_two_tailed_2-------------------------------------------
# ## Plot rates vs. trait values
# 
# # We can plot the mean rates vs. trait values inferred by deepSTRAPP for T = 10Mya
# # to unravel the direction of the correlation detected.
# 
# ?plot_rates_vs_trait_data_for_focal_time()
# 
# # Select a color scheme from lowest to highest values (i.e., small ants to large ants)
# color_scale = c("darkgreen", "limegreen", "orange", "red")
# 
# # Plot mean rates vs. traits aggregated across all BAMM posterior samples
# plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_outputs = deepSTRAPP_output_two_tailed,
#                                         color_scale = color_scale)
# 
# # The plot indicates that mean rates are globally negatively correlated with mean trait values as
# # "small" ants exhibit higher mean rates than "large" ants for T = 10 Mya.
# # The plot also displays summary statistics for the STRAPP test associated with the data shown:
# #   * An observed statistic computed across the mean rates and trait values shown on the plot.
# #     Here, rho obs = -0.86, indicates a negative correlation.
# #     This is not the statistic of the STRAPP test itself, which is conducted across all BAMM posterior samples.
# #   * The quantile of null statistic distribution at the significant threshold used to define test significance.
# #     The test will be considered significant (i.e., the null hypothesis is rejected)
# #     if this value is higher than zero, as shown on the histogram above.
# #     Here, Q5% = 0.066, so the test is significant (according to this significance threshold).
# #   * The p-value of the associated STRAPP test. Here, p = 0.018.
# 
# # We could have directly tested for hypothesis with a one-tailed test (See Step 3 below).
# 

## ----STRAPP_tests_cont_two_tailed_eval_2_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# 
# ## Plot rates vs. trait values
# 
# # Select a color scheme from lowest to highest values (i.e., small ants to large ants)
# color_scale = c("darkgreen", "limegreen", "orange", "red")
# 
# # Plot mean rates vs. traits aggregated across all BAMM posterior samples
# plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_outputs = deepSTRAPP_output_two_tailed,
#                                         color_scale = color_scale)
# 

## ----STRAPP_tests_cont_two_tailed_eval_2_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_1.2_Example_continuous_two_tailed_2.PNG")


## ----STRAPP_tests_cont_one_tailed---------------------------------------------
# #### Step 3: Run a one-tailed test ####
# 
# # Here, we make an assumption on the direction of the effect tested.
# # The null hypothesis is that there is a positive or no correlation between rates and trait values.
# # The alternative hypothesis is that there is a negative correlation between rates and trait values.
# 
# # We run a one-tailed STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_one_tailed <- run_deepSTRAPP_for_focal_time(
#    contMap = Ponerinae_trait_object$contMap,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cont_tip_data,
#    trait_data_type = "continuous",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    uncertainty_strategy = "rates_only",
#    two_tailed = FALSE, # Select a one-tailed test
#    one_tailed_hypothesis = "negative", # We select a "negative" correlation as the alternative hypothesis
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE) # Keep permutation data to be able to plot histograms of tests
# 
# ## Explore output
# str(deepSTRAPP_output_one_tailed, max.level = 1)
# 
# # Access deepSTRAPP results
# deepSTRAPP_output_one_tailed$STRAPP_results[1:5]
# 
# # We obtain a p-value = 0.009 leading us to discard the null hypothesis
# # (no correlation or positive correlation).
# # Thus, the alternative hypothesis (negative correlation) is supported by the test.
# # The estimate is the quantile 95% of the distribution of differences in Spearman rho-stats
# # between observed and permuted data.
# # The test is significant as the estimate is lower than zero.
# 
# # Plot histogram of Spearman sum-rank test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_one_tailed)
# 
# # The test is significant as the red line (quantile 95% = -0.147) is below the null expectation
# # that Δ Spearman rho stat = 0.
# 

## ----STRAPP_tests_cont_one_tailed_eval_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # We run a one-tailed STRAPP test for t = 10 Mya.
# deepSTRAPP_output_one_tailed <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    contMap = Ponerinae_trait_object$contMap,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cont_tip_data,
#    trait_data_type = "continuous",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    two_tailed = FALSE, # Select a one-tailed test
#    one_tailed_hypothesis = "negative", # We select a "negative" correlation as the alternative hypothesis
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    verbose = FALSE))
# 
# # Plot histogram of Spearman sum-rank test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_one_tailed)

## ----STRAPP_tests_cont_one_tailed_eval_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_1.3_Example_continuous_one_tailed.PNG")


## ----STRAPP_tests_cat_2lvl_two_tailed-----------------------------------------
# # ------ Example 2: Categorical binary trait data ------ #
# 
# #### Step 1: Prepare data ####
# 
# ## Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# ## Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# 
# ## Load BAMM diversification outputs
# data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")
# # This dataset is only available in development versions installed from GitHub.
# # It is not available in CRAN versions.
# # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
# 
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cat_2lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_2lvl_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# # Reorder tip data as in phylogeny
# Ponerinae_cat_2lvl_tip_data <- Ponerinae_cat_2lvl_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# # Select color scheme for states
# colors_per_states <- c("darkblue", "lightblue")
# names(colors_per_states) <- c("large", "small")
# 
# # Run prepare_trait_data with default options
# Ponerinae_trait_object <- prepare_trait_data(tip_data = Ponerinae_cat_2lvl_tip_data,
#                                              trait_data_type = "categorical",
#                                              phylo = Ponerinae_tree_old_calib,
#                                              evolutionary_model = "ARD", # Select an All-Rates-Different model
#                                              colors_per_states = colors_per_states,
#                                              nb_simulations = 100,
#                                              seed = 1234) # Set seed for reproducibility
# # Load directly output to save time
# data(Ponerinae_cat_2lvl_data_old_calib, package = "deepSTRAPP")
# Ponerinae_trait_object <- Ponerinae_cat_2lvl_data_old_calib
# 
# 
# #### Step 2: Run a two-tailed test ####
# 
# # No assumption is made about which states exhibit the highest rates.
# # The null hypothesis is that there is no difference in rates between the two states.
# # The alternative hypothesis is that there is a difference in rates between the two states,
# # in any direction (i.e., "small > large" or "large > small").
# 
# # Here, we run a two-tailed STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_two_tailed <- run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_2lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    two_tailed = TRUE, # Select a two-tailed test (Default)
#    rate_type = "net_diversification",
#    uncertainty_strategy = "rates_only",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    # Needed to access traits and rates data to evaluate the direction of the correlation
#    extract_trait_data_melted_df = TRUE,
#    extract_diversification_data_melted_df = TRUE)
# 
# ## Explore output
# str(deepSTRAPP_output_two_tailed, max.level = 1)
# 
# # Access deepSTRAPP results
# deepSTRAPP_output_two_tailed$STRAPP_results[1:5]
# 
# # We obtain a p-value = 0.047 leading us to (barely) discard the null hypothesis (no differences).
# # The estimate is the quantile 5% of the distribution of differences
# # in absolute Mann-Whitney-Wilcoxon U-stats between observed and permuted data.
# # The test is significant as the estimate is higher than zero.
# # However, without plotting the rates vs. states, we do not know the direction of this difference.
# 
# # Plot histogram of Mann-Whitney-Wilcoxon test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_two_tailed)
# 
# # The test is significant as the red line (quantile 5% = 176.45) is above the null expectation
# # that Δ abs(Mann-Whitney-Wilcoxon U-stats) = 0.

## ----STRAPP_tests_cat_2lvl_two_tailed_eval_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# # Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# ## Load BAMM diversification outputs
# data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cat_2lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_2lvl_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# # Reorder tip data as in phylogeny
# Ponerinae_cat_2lvl_tip_data <- Ponerinae_cat_2lvl_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# # Load directly deepSTRAPP output to save time
# data(Ponerinae_cat_2lvl_data_old_calib, package = "deepSTRAPP")
# Ponerinae_trait_object <- Ponerinae_cat_2lvl_data_old_calib
# 
# # Select color scheme for states
# colors_per_states <- c("darkblue", "lightblue")
# names(colors_per_states) <- c("large", "small")
# 
# # Here, we run a two-tailed STRAPP test for t = 10 Mya.
# deepSTRAPP_output_two_tailed <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_2lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    two_tailed = TRUE, # Select a two-tailed test (Default)
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    # Needed to access traits and rates data to evaluate the direction of the correlation
#    extract_trait_data_melted_df = TRUE,
#    extract_diversification_data_melted_df = TRUE,
#    verbose = FALSE))
# 
# # Plot histogram of Mann-Whitney-Wilcoxon test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_two_tailed)

## ----STRAPP_tests_cat_2lvl_two_tailed_eval_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_2.2_Example_cat_2lvl_two_tailed.PNG")


## ----STRAPP_tests_cat_2lvl_two_tailed_2---------------------------------------
# ## Plot rates vs. trait values
# 
# # We can plot the mean rates vs. states inferred by deepSTRAPP for T = 10Mya
# # to unravel the direction of the difference detected.
# 
# ?plot_rates_vs_trait_data_for_focal_time()
# 
# # Plot mean rates vs. trait states aggregated across all BAMM posterior samples
# plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_outputs = deepSTRAPP_output_two_tailed,
#                                         colors_per_levels = colors_per_states)
# 
# # The plot highlights that "small" ants exhibited higher rates than "large" ants for T = 10 Mya.
# # The plot also displays summary statistics for the STRAPP test associated with the data shown:
# #   * An observed statistic computed across the mean rates and trait states shown on the plot.
# #     Here, the very low U-stats obs = -38713, indicates the first group "large" ants have lower rates than "small" ants.
# #     This is not the statistic of the STRAPP test itself, which is conducted across all BAMM posterior samples.
# #   * The quantile of null statistic distribution at the significant threshold used to define test significance.
# #     The test will be considered significant (i.e., the null hypothesis is rejected)
# #     if this value is higher than zero, as shown on the histogram above.
# #     Here, Q5% = 176.45, so the test is significant (according to this significance threshold).
# #   * The p-value of the associated STRAPP test. Here, p = 0.047.
# 
# # We could have directly tested for this hypothesis (i.e., "small" > "large") with a one-tailed test (See Step 3 below).
# 

## ----STRAPP_tests_cat_2lvl_two_tailed_eval_2_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# ## Plot rates vs. trait values
# 
# # Plot mean rates vs. trait states aggregated across all BAMM posterior samples
# plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_outputs = deepSTRAPP_output_two_tailed,
#                                         colors_per_levels = colors_per_states)
# 

## ----STRAPP_tests_cat_2lvl_two_tailed_eval_2_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_2.2_Example_cat_2lvl_two_tailed_2.PNG")


## ----STRAPP_tests_cat_2lvl_one_tailed-----------------------------------------
# #### Step 3: Run a one-tailed test ####
# 
# # Here, we make an assumption on the direction of the effect tested.
# # The null hypothesis is that diversification rates are lower or equal
# # for "small" ants compared to "large" ants.
# # The alternative hypothesis is that diversification rates are higher
# # for "small" ants than "large" ants.
# 
# # We run a one-tailed STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_one_tailed <- run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_2lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    uncertainty_strategy = "rates_only",
#    two_tailed = FALSE, # Select a one-tailed test
#    one_tailed_hypothesis = "small > large", # We select "small > large" as the alternative hypothesis
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE) # Keep permutation data to be able to plot histograms of tests
# 
# ## Explore output
# str(deepSTRAPP_output_one_tailed, max.level = 1)
# 
# # Access deepSTRAPP results
# deepSTRAPP_output_one_tailed$STRAPP_results[1:5]
# 
# # We obtain a p-value = 0.020 leading us to discard the null hypothesis ("small <= large").
# # Thus, the alternative hypothesis ("small > large") is supported by the test.
# # The estimate is the quantile 5% of the distribution of differences in Mann-Whitney-Wilcoxon U-stats
# # between observed and permuted data
# # The test is significant as the estimate is higher than zero.
# 
# # Plot histogram of Mann-Whitney-Wilcoxon test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_one_tailed)
# 
# # The test is significant as the red line (quantile 5% = 3335.4) is above the null expectation
# # that Δ Mann-Whitney-Wilcoxon U-stats = 0.
# 

## ----STRAPP_tests_cat_2lvl_one_tailed_eval_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # We run a one-tailed STRAPP test for t = 10 Mya.
# deepSTRAPP_output_one_tailed <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_2lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    two_tailed = FALSE, # Select a one-tailed test
#    one_tailed_hypothesis = "small > large", # We select "small > large" as the alternative hypothesis
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    verbose = FALSE))
# 
# # Plot histogram of Mann-Whitney-Wilcoxon test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_one_tailed)

## ----STRAPP_tests_cat_2lvl_one_tailed_eval_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_2.3_Example_cat_2lvl_one_tailed.PNG")


## ----STRAPP_tests_cat_3lvl_overall--------------------------------------------
# # ------ Example 3: Categorical multinominal (3-levels) trait data ------ #
# 
# #### Step 1: Prepare data ####
# 
# ## Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# ## Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# 
# ## Load BAMM diversification outputs
# data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")
# # This dataset is only available in development versions installed from GitHub.
# # It is not available in CRAN versions.
# # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
# 
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# # Reorder tip data as in phylogeny
# Ponerinae_cat_3lvl_tip_data <- Ponerinae_cat_3lvl_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# ## Select color scheme for states (i.e., habitats)
# colors_per_states <- c("forestgreen", "sienna", "goldenrod")
# names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
# 
# # Run prepare_trait_data with default options
# Ponerinae_trait_object <- prepare_trait_data(tip_data = Ponerinae_cat_3lvl_tip_data,
#                                              trait_data_type = "categorical",
#                                              phylo = Ponerinae_tree_old_calib,
#                                              evolutionary_model = "ARD", # Select an All-Rates-Different model
#                                              colors_per_states = colors_per_states,
#                                              nb_simulations = 100,
#                                              seed = 1234) # Set seed for reproducibility
# # Load directly output to save time
# data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
# Ponerinae_trait_object <- Ponerinae_cat_3lvl_data_old_calib
# 
# #### Step 2: Run the overall Kruskal-Wallis test ####
# 
# # For multinominal data (i.e., with more than two states/ranges), the overall test
# # is by definition a one-tailed Kruskal-Wallis test for rate differences across all states/ranges.
# # There is no two-tailed version of this test, thus deepSTRAPP assumes a one-tailed test by default.
# # The null hypothesis is that there is no difference in rates across all states/ranges.
# # The alternative hypothesis is that there is at least one pairs of states/ranges
# # with differences in diversification rates.
# 
# # To test for differences between pairs of states/ranges, you need to run post hoc pairwise tests.
# # These Dunn's test can be one-tailed or two-tailed depending on the null hypotheses stated.
# # See the next Steps 3 & 4 for post hoc pairwise tests.
# 
# # Here, we run an overall STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_overall <- run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_3lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    uncertainty_strategy = "rates_only",
#    alpha = 0.10, # Adjust the significance threshold
#    posthoc_pairwise_tests = FALSE, # To run the overall
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    # Needed to access traits and rates data to evaluate the direction of the correlation
#    extract_trait_data_melted_df = TRUE,
#    extract_diversification_data_melted_df = TRUE)
# 
# ## Explore output
# str(deepSTRAPP_output_overall, max.level = 1)
# 
# # Access deepSTRAPP results
# deepSTRAPP_output_overall$STRAPP_results[1:5]
# 
# # We obtain a p-value = 0.083 leading us to discard the null hypothesis (no differences).
# # The estimate is the quantile 10% of the distribution of differences
# # in Kruskal-Wallis H-stats between observed and permuted data.
# # The test is significant as the estimate is higher than zero.
# # However, without plotting the rates vs. states, or running post hoc pairwise tests,
# # we do not know which pairs of states exhibit differences.
# 
# # Plot histogram of Kruskall-Wallis one-way ANOVA test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_overall)
# 
# # The test is significant as the red line (quantile 10% = 5.606) is above the null expectation
# # that Δ Kruskal-Wallis H-stats = 0.

## ----STRAPP_tests_cat_3lvl_overall_eval_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# # Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# ## Load BAMM diversification outputs
# data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# # Reorder tip data as in phylogeny
# Ponerinae_cat_3lvl_tip_data <- Ponerinae_cat_3lvl_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# # Load directly deepSTRAPP output to save time
# data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
# Ponerinae_trait_object <- Ponerinae_cat_3lvl_data_old_calib
# 
# ## Select color scheme for states
# colors_per_states <- c("forestgreen", "sienna", "goldenrod")
# names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
# 
# # Here, we run an overall STRAPP test for t = 10 Mya.
# deepSTRAPP_output_overall <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_3lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    alpha = 0.10, # Adjust the significance threshold
#    posthoc_pairwise_tests = FALSE, # To run the overall
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    # Needed to access traits and rates data to evaluate the direction of the correlation
#    extract_trait_data_melted_df = TRUE,
#    extract_diversification_data_melted_df = TRUE,
#    verbose = FALSE))
# 
# # Plot histogram of Kruskall-Wallis one-way ANOVA test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_overall)

## ----STRAPP_tests_cat_3lvl_overall_eval_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_3.2_Example_cat_3lvl_overall.PNG")


## ----STRAPP_tests_cat_3lvl_overall_2------------------------------------------
# ## Plot rates vs. trait values
# 
# # We can plot the mean rates vs. states inferred by deepSTRAPP for T = 10Mya
# # to unravel the direction of the differences detected.
# 
# ?plot_rates_vs_trait_data_for_focal_time()
# 
# # Plot mean rates vs. trait states aggregated across all BAMM posterior samples
# plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_outputs = deepSTRAPP_output_overall,
#                                         colors_per_levels = colors_per_states)
# 
# # The plot indicates that net diversification rates are higher in "terricolous" ants,
# # than "subterranean" ants, and lower in "arboreal" ants, for T = 10 Mya.
# # The plot also displays summary statistics for the STRAPP test associated with the data shown:
# #   * An observed statistic computed across the mean rates and trait states (i.e., habitats) shown on the plot.
# #     Here, the high H-stat obs = 374.82 indicates there is a difference in rates between states.
# #     Please note that this is not the statistic of the STRAPP test itself,
# #     which is conducted across all BAMM posterior samples.
# #   * The quantile of null statistic distribution at the significant threshold used to define test significance.
# #     The test will be considered significant (i.e., the null hypothesis is rejected)
# #     if this value is higher than zero, as shown on the histogram above.
# #     Here, Q10% = 5.606, so the test is significant (according to this significance threshold),
# #     meaning we detected differences between habitats overall.
# #   * The p-value of the associated STRAPP test. Here, p = 0.083.
# 
# # We can explore those potential differences with post hoc pairwise tests,
# # with both one-tailed hypothesis tests or two-tailed tests (See Steps 3 & 4 below).
# 

## ----STRAPP_tests_cat_3lvl_overall_eval_2_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# ## Plot rates vs. trait values
# 
# # Plot mean rates vs. trait states aggregated across all BAMM posterior samples
# plot_rates_vs_trait_data_for_focal_time(deepSTRAPP_outputs = deepSTRAPP_output_overall,
#                                         colors_per_levels = colors_per_states)
# 

## ----STRAPP_tests_cat_3lvl_overall_eval_2_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_3.2_Example_cat_3lvl_overall_2.PNG")


## ----STRAPP_tests_cat_3lvl_two_tailed_eval2-----------------------------------
# #### Step 3: Run all post hoc Dunn's two-tailed tests ####
# 
# # Post hoc tests are pairwise tests for differences between pairs of states/ranges.
# # They can be two-tailed, without assumption regarding which states exhibit higher rates than the others.
# # or they can be one-tailed, with assumption on which states exhibit higher rates within each pair.
# 
# # deepSTRAPP runs post hoc pairwise tests across all pairs of states.
# # You can choose which type of tests (i.e., "two-tailed" or "one-tailed") to run.
# # Here, the two-tailed version is presented. See Step 4 for the one-tailed version.
# 
# # In the two-tailed tests, the null hypothesis for each pair is that there is
# # no difference in rates between the two states.
# # The alternative hypothesis is that there is a difference in rates between the two states,
# # in any direction (e.g., "arboreal > terricolous" or "terricolous > arboreal").
# 
# # You can also select a method to adjust p-values for multiple testing.
# # See stats::p.adjust() for the available methods.
# # Here, we select none.
# 
# # We run a two-tailed STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_two_tailed <- run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_3lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    uncertainty_strategy = "rates_only",
#    two_tailed = TRUE, # Select two-tailed tests for the post hoc tests (Default)
#    posthoc_pairwise_tests = TRUE, # Ask for post hoc tests to be carried out
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE) # Keep permutation data to be able to plot histograms of tests
# 
# ## Explore output
# str(deepSTRAPP_output_two_tailed, max.level = 2)
# 
# # Access deepSTRAPP results for post hoc tests
# deepSTRAPP_output_two_tailed$STRAPP_results$posthoc_pairwise_tests$summary_df
# 
# # We obtain a p-value for each pair of states. Here there are 3 possible pairs.
# # Only the pair "arboreal != terricolous" yields to a significant result.
# # For this pair, we obtain a p-value = 0.025 leading us to discard the null hypothesis (no differences).
# # The estimate is the quantile 5% of the distribution of differences
# # in absolute Dunn's Z-stats between observed and permuted data.
# # The test is significant as the estimate is higher than zero (Q5% = 0.759).
# # However, without plotting the rates vs. states, we would have not know the direction of this difference.
# 
# # Plot histogram of all post hoc two-tailed test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_two_tailed,
#   plot_posthoc_tests = TRUE)
# 
# # Only the test for the pair "arboreal =! terricolous" is significant as the red line
# # (quantile 5% = 0.532) is above the null expectation that Δ Dunn's Z-stats = 0.
# # The conclusion is that it is the difference in rates between "arboreal" ants
# # and "terricolous" ants that drives the differences in rates recorded overall
# # with the Kruskal-Wallis test.
# # Yet, with this test only, we do not know which states has the higher rates
# # We can look at the plot of rates vs. states above (See Step 2),
# # but we can also test directly for the two hypotheses (i.e., ("arboreal" > "terricolous")
# # or ("terricolous" > "arboreal") using one-tailed tests.
# # See Step 4 below for the one-tailed tests.

## ----STRAPP_tests_cat_3lvl_two_tailed_eval_dev, fig.width = 8.5, fig.height = 7, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # Load phylogeny with old time-calibration
# data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
# # Load trait df
# data("Ponerinae_trait_tip_data", package = "deepSTRAPP")
# 
# # Extract continuous trait data as a named vector
# Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#                                     nm = Ponerinae_trait_tip_data$Taxa)
# # Reorder tip data as in phylogeny
# Ponerinae_cat_3lvl_tip_data <- Ponerinae_cat_3lvl_tip_data[Ponerinae_tree_old_calib$tip.label]
# 
# # Load directly deepSTRAPP output to save time
# data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
# Ponerinae_trait_object <- Ponerinae_cat_3lvl_data_old_calib
# 
# ## Select color scheme for states
# colors_per_states <- c("forestgreen", "sienna", "goldenrod")
# names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
# 
# # Here, we run a two-tailed STRAPP test for t = 10 Mya.
# deepSTRAPP_output_two_tailed <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_3lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    two_tailed = TRUE, # Select two-tailed tests for the post hoc tests (Default)
#    posthoc_pairwise_tests = TRUE, # Ask for post hoc tests to be carried out
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE,
#    verbose = FALSE)) # Keep permutation data to be able to plot histograms of tests
# 
# # Plot histogram of all post hoc two-tailed test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_two_tailed,
#   plot_posthoc_tests = TRUE)

## ----STRAPP_tests_cat_3lvl_two_tailed_eval_CRAN, fig.width = 8.5, fig.height = 7, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_3.3_Example_cat_3lvl_two_tailed.PNG")


## ----STRAPP_tests_cat_3lvl_one_tailed-----------------------------------------
# #### Step 4: Run all post hoc Dunn's one-tailed tests ####
# 
# # Here, we make assumptions on the direction of the effects tested.
# # We conduct all possible one-tailed tests across pairs of states.
# # Thus, for three pairs, we conduct six tests, one for each possible alternative hypothesis
# # (i.e., "A" > "B" and "B" > "A").
# 
# # We run a one-tailed STRAPP test for t = 10 Mya.
# focal_time <- 10
# 
# deepSTRAPP_output_one_tailed <- run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_3lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = focal_time,
#    uncertainty_strategy = "rates_only",
#    two_tailed = FALSE, # Select one-tailed tests for the post hoc tests
#    posthoc_pairwise_tests = TRUE, # Ask for post hoc tests to be carried out
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE) # Keep permutation data to be able to plot histograms of tests
# 
# ## Explore output
# str(deepSTRAPP_output_one_tailed, max.level = 2)
# 
# # Access deepSTRAPP results for post hoc tests
# deepSTRAPP_output_one_tailed$STRAPP_results$posthoc_pairwise_tests$summary_df
# 
# # We obtain a p-value for each pair of states in each direction.
# # Here there are 3 possible pairs = 6 one-tailed tests.
# # Only the alternative hypothesis "arboreal < terricolous" yields to a significant result.
# # For this test, we obtain a p-value = 0.014 leading us to discard the null hypothesis that
# # "arboreal" ants have equivalent or higher rates than "terricolous" ants.
# # The estimate is the quantile 5% of the distribution of differences
# # in Dunn's Z-stats between observed and permuted data.
# # The test is significant as the estimate is higher than zero (Q5% = 1.713).
# # We can conclude that our data support the idea that "arboreal" ants
# # have lower diversification rates than "terricolous" ants.
# 
# # Plot histogram of all post hoc one-tailed test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_one_tailed,
#   plot_posthoc_tests = TRUE)
# 
# # Only the test for the alternative hypothesis "arboreal < terricolous" is significant as the red line
# # (quantile 5% = 1.713) is above the null expectation that Δ Dunn's Z-stats = 0.
# # The conclusion is that it is the higher rates in "terricolous" ants vs. lower rates in "arboreal" ants
# # that drives the differences in rates recorded overall with the Kruskal-Wallis test.

## ----STRAPP_tests_cat_3lvl_one_tailed_eval_dev, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # We run a one-tailed STRAPP test for t = 10 Mya.
# deepSTRAPP_output_one_tailed <- suppressWarnings(run_deepSTRAPP_for_focal_time(
#    densityMaps = Ponerinae_trait_object$densityMaps,
#    ace = Ponerinae_trait_object$ace,
#    tip_data = Ponerinae_cat_3lvl_tip_data,
#    trait_data_type = "categorical",
#    BAMM_object = Ponerinae_BAMM_object_old_calib,
#    focal_time = 10,
#    uncertainty_strategy = "rates_only",
#    two_tailed = FALSE, # Select one-tailed tests for the post hoc tests
#    posthoc_pairwise_tests = TRUE, # Ask for post hoc tests to be carried out
#    rate_type = "net_diversification",
#    seed = 1234, # Set seed for reproducibility
#    return_perm_data = TRUE, # Keep permutation data to be able to plot histograms of tests
#    verbose = FALSE))
# 
# # Plot histogram of all post hoc one-tailed test stats
# plot_histogram_STRAPP_test_for_focal_time(
#   deepSTRAPP_outputs = deepSTRAPP_output_one_tailed,
#   plot_posthoc_tests = TRUE)

## ----STRAPP_tests_cat_3lvl_one_tailed_eval_CRAN, fig.width = 8.5, fig.height = 6, out.width = "100%", eval = !is_dev_version(), echo = FALSE----

# Plot pre-rendered graph
knitr::include_graphics("figures/4_Explore_STRAPP_hypotheses_3.4_Example_cat_3lvl_one_tailed.PNG")


