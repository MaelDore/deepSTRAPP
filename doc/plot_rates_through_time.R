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


## ----adjust_dpi_CRAN, include = FALSE, eval = !is_dev_version()---------------
knitr::opts_chunk$set(
  dpi = 50   # Lower DPI to save space
)

## ----adjust_dpi_dev, include = FALSE, eval = is_dev_version()-----------------
# knitr::opts_chunk$set(
#   dpi = 72   # Default DPI for the dev version
# )

## ----plot_RTT_cont------------------------------------------------------------
# # ------ Example 1: Continuous trait data ------ #
# 
# #### Step 1: Prepare data ####
# 
# # Load results of a STRAPP test workflow run on continuous trait data
# data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
# # This dataset is only available in development versions installed from GitHub.
# # It is not available in CRAN versions.
# # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
# 
# # Visualize trait data
# hist(Ponerinae_deepSTRAPP_cont_old_calib_0_40$trait_data_df_over_time$trait_value)
# # Visualize diversification rates data
# hist(Ponerinae_deepSTRAPP_cont_old_calib_0_40$diversification_data_df_over_time$rates)
# 
# # Select a color scheme from lowest to highest values
# color_scale = c("darkgreen", "limegreen", "orange", "red")
# 
# ### Step 2: Plot RTT ####
# 
# # For continuous trait data, the trait data are divided into groups based on 'quantile_ranges'.
# # They form groups of trait values used to visualize the existence of a correlation with rates.
# # If groups with higher quartile ranges exhibit higher rates over time, and conversely, there is a positive correlation.
# # If groups with higher quartile ranges exhibit lower rates over time, and conversely, there is a negative correlation.
# 
# ## Generate default plot
# plot_RTT_continuous <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    color_scale = color_scale)
#    # PDF_file_path = "./plot_RTT_continuous.pdf")
# 
# # Here, we observed a negative correlation as quartiles of lower trait values are consistently
# # associated with higher rates, and conversely.
# 
# ## Adjust aesthetics of plot a posteriori
# plot_RTT_continuous_adj <- plot_RTT_continuous$rates_TT_ggplot +
#     # Change size and color of the title
#     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 18))
# print(plot_RTT_continuous_adj)
# 
# ## Adjust quartile ranges
# 
# # You can define different trait quantile groups based on 'quantile_ranges'.
# # The default is c(0, 0.25, 0.5, 0.75, 1.0) which yields for equivalent groups
# # with each spanning 25% of the trait data.
# 
# # You can try a different option with for instance five groups:
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    quantile_ranges = c(0, 0.20, 0.40, 0.60, 0.80, 1.0),
#    color_scale = color_scale)
#    # PDF_file_path = "./plot_RTT_continuous_5groups.pdf")
# 
# ## Plot Confidence Intervals (CI)
# 
# # You can include confidence intervals on the plot.
# # They are build from the BAMM posterior samples, which each yields different diversification rates,
# # that can be plotted to show uncertainty in the evaluation.
# 
# # The function offers two types of CI (CI_type):
#   # "fuzzy": draws a line for each posterior sample with high transparency to highlight
#            # redundancy across BAMM posterior samples.
#   # "quantiles_rect": draws plain polygons whose range can be controlled with 'CI_quantiles'.
#                     # By default, 95% CI are shown (CI_quantiles = 0.95).
# 
# # Plot "fuzzy" CI
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "fuzzy") # Select type of CI
# 
# # Plot "quantiles_rect" CI
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "quantiles_rect", # Select type of CI
#    CI_quantiles = 0.9) # Adjust range of CI
# 
# # The lack of overlap between the two most extreme trait quantile groups tends to indicate
# # that the correlation is significant. However, a proper STRAPP test is needed to confirm what is observed.
# 
# ## Plot different types of rates
# 
# # deepSTRAPP also let you plot different types of rates.
# # Even if you carried out the analysis on "net_diversification" rates (as the default)
# # you can still plot the evolution of "speciation" and "extinction" rates
# 
# # Plot "speciation" rates
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    rate_type = "speciation",
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "quantiles_rect", # Select type of CI
#    CI_quantiles = 0.9) # Adjust range of CI
# 
# # Plot "extinction" rates
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    rate_type = "extinction",
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "quantiles_rect", # Select type of CI
#    CI_quantiles = 0.9) # Adjust range of CI
# 

## ----plot_RTT_cont_eval_dev, fig.width = 16, fig.height = 15, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # Load results of a STRAPP test workflow run on continuous trait data
# data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
# 
# # Select a color scheme from lowest to highest values
# color_scale = c("darkgreen", "limegreen", "orange", "red")
# 
# ## Generate default plot
# plot_RTT_continuous_1 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    color_scale = color_scale, display_plot = FALSE, verbose = FALSE)
# plot_RTT_continuous_1 <- plot_RTT_continuous_1$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Default plot")
# 
# # Plot with five trait quantile groups
# plot_RTT_continuous_2 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    quantile_ranges = c(0, 0.20, 0.40, 0.60, 0.80, 1.0),
#    color_scale = color_scale, display_plot = FALSE, verbose = FALSE)
# plot_RTT_continuous_2 <- plot_RTT_continuous_2$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Five groups")
# 
# # Plot "fuzzy" CI
# plot_RTT_continuous_3 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "fuzzy", # Select type of CI
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_continuous_3 <- plot_RTT_continuous_3$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Fuzzy CI")
# 
# # Plot "quantiles_rect" CI
# plot_RTT_continuous_4 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "quantiles_rect", # Select type of CI
#    CI_quantiles = 0.9, # Adjust range of CI
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_continuous_4 <- plot_RTT_continuous_4$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Rectangular CI")
# 
# # Plot "speciation" rates
# plot_RTT_continuous_5 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    rate_type = "speciation",
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "quantiles_rect", # Select type of CI
#    CI_quantiles = 0.9, # Adjust range of CI
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_continuous_5 <- plot_RTT_continuous_5$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Speciation rates")
# 
# # Plot "extinction" rates
# plot_RTT_continuous_6 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#    rate_type = "extinction",
#    color_scale = color_scale,
#    plot_CI = TRUE, # To add CI on the plot
#    CI_type = "quantiles_rect", # Select type of CI
#    CI_quantiles = 0.9, # Adjust range of CI
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_continuous_6 <- plot_RTT_continuous_6$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Extinction rates")
# 
# cowplot::plot_grid(plotlist = list(plot_RTT_continuous_1, plot_RTT_continuous_2,
#                                    plot_RTT_continuous_3, plot_RTT_continuous_4,
#                                    plot_RTT_continuous_5, plot_RTT_continuous_6),
#                    ncol = 2, nrow = 3)
# 

## ----plot_RTT_cont_eval_CRAN, eval = !is_dev_version(), echo = FALSE, out.width = "100%"----

# Plot pre-rendered graph
knitr::include_graphics("figures/5_Explore_plot_RTT_1_Example_continuous.PNG")


## ----plot_RTT_cat_3lvl--------------------------------------------------------
# # ------ Example 2: Categorical multinominal trait data ------ #
# 
# #### Step 1: Prepare data ####
# 
# # Load results of a STRAPP test workflow run on continuous trait data
# data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, package = "deepSTRAPP")
# # This dataset is only available in development versions installed from GitHub.
# # It is not available in CRAN versions.
# # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
# 
# # Visualize trait data
# table(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40$trait_data_df_over_time$trait_value)
# # Visualize diversification rates data
# hist(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40$diversification_data_df_over_time$rates)
# 
# ## Select color scheme for states
# colors_per_states <- c("forestgreen", "sienna", "goldenrod")
# names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
# 
# ### Step 2: Plot RTT ####
# 
# # For categorical trait data, the trait data are states.
# # We can directly plot the evolution of rates between states.
# # In the case of multinominal data (with more than two states/ranges),
# # the function allows you to select the states you wish to plot, if not all.
# 
# ## Generate default plot with all states
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    colors_per_levels = colors_per_states)
#    # PDF_file_path = "./plot_RTT_categorical.pdf")
# 
# # Plot with "fuzzy" CI
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    colors_per_levels = colors_per_states,
#    plot_CI = TRUE,
#    CI_type = "fuzzy")
# 
# # Plot with "quantiles_rect" CI
# plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    colors_per_levels = colors_per_states,
#    plot_CI = TRUE,
#    CI_type = "quantiles_rect")
# 
# ## Subset to plot only "arboreal" and "terricolous" states
# plot_RTT_categorical <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    select_trait_levels = c("arboreal", "terricolous"), # List the states to plot
#    colors_per_levels = colors_per_states[c("arboreal", "terricolous")], # Subset colors
#    plot_CI = TRUE,
#    CI_type = "quantiles_rect")
# 

## ----plot_RTT_cat_3lvl_eval_dev, fig.width = 14, fig.height = 10, out.width = "100%", eval = is_dev_version(), echo = FALSE----
# # Load results of a STRAPP test workflow run on continuous trait data
# data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, package = "deepSTRAPP")
# 
# ## Select color scheme for states
# colors_per_states <- c("forestgreen", "sienna", "goldenrod")
# names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
# 
# ## Generate default plot
# plot_RTT_categorical_1 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    colors_per_levels = colors_per_states,
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_categorical_1 <- plot_RTT_categorical_1$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Default plot")
# 
# # Plot "fuzzy" CI
# plot_RTT_categorical_2 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    colors_per_levels = colors_per_states,
#    plot_CI = TRUE,
#    CI_type = "fuzzy",
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_categorical_2 <- plot_RTT_categorical_2$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Fuzzy CI")
# 
# # Plot "quantiles_rect" CI
# plot_RTT_categorical_3 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    colors_per_levels = colors_per_states,
#    plot_CI = TRUE,
#    CI_type = "quantiles_rect",
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_categorical_3 <- plot_RTT_categorical_3$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Rectangular CI")
# 
# ## Subset to plot only "arboreal" and "terricolous" states
# plot_RTT_categorical_4 <- plot_rates_through_time(
#    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#    select_trait_levels = c("arboreal", "terricolous"), # List the states to plot
#    colors_per_levels = colors_per_states[c("arboreal", "terricolous")], # Subset colors
#    plot_CI = TRUE,
#    CI_type = "quantiles_rect",
#    display_plot = FALSE, verbose = FALSE)
# plot_RTT_categorical_4 <- plot_RTT_categorical_4$rates_TT_ggplot +
#   ggplot2::ggtitle(label = "Subset two states")
# 
# cowplot::plot_grid(plotlist = list(plot_RTT_categorical_1, plot_RTT_categorical_2,
#                                    plot_RTT_categorical_3, plot_RTT_categorical_4),
#                    ncol = 2, nrow = 2)
# 

## ----plot_RTT_cat_3lvl_eval_CRAN, eval = !is_dev_version(), echo = FALSE, out.width = "100%"----

# Plot pre-rendered graph
knitr::include_graphics("figures/5_Explore_plot_RTT_2_Example_multinominal.PNG")


