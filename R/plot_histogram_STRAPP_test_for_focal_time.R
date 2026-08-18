
#' @title Plot histogram of STRAPP test statistics to assess results
#'
#' @description Plot an histogram of the distribution of the test statistics
#'   obtained from a STRAPP test carried out for a unique `focal_time`.
#'   (See [deepSTRAPP::compute_STRAPP_test_for_focal_time()] and
#'   [deepSTRAPP::run_deepSTRAPP_for_focal_time()]).
#'
#'   Returns a single histogram for overall tests.
#'   If `plot_posthoc_tests = TRUE`, it will return a faceted plot with an histogram per post hoc tests.
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plot will be saved directly in a PDF file.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_for_focal_time()],
#'   that summarize the results of a STRAPP test for a specific time in the past (i.e. the `focal_time`).
#'   `deepSTRAPP_outputs` can also be extracted from the output of [deepSTRAPP::run_deepSTRAPP_over_time()] that
#'   run the whole deepSTRAPP workflow and store results inside `$STRAPP_results_over_time`.
#' @param focal_time Numerical. (Optional) If `deepSTRAPP_outputs` comprises results over multiple time-steps
#'   (i.e., output of [deepSTRAPP::run_deepSTRAPP_over_time()], this is the time of the STRAPP test targeted for plotting.
#' @param display_plot Logical. Whether to display the histogram(s) generated in the R console. Default is `TRUE`.
#' @param plot_posthoc_tests Logical. For multinominal data only. Whether to plot the histogram for the overall Kruskal-Wallis test across all states (`plot_posthoc_tests = FALSE`),
#'   or plot the histograms for all the pairwise post hoc Dunn's test across pairs of states (`plot_posthoc_tests = TRUE`). Default is `FALSE`.
#' @param PDF_file_path Character string. If provided, the plot will be saved in a PDF file following the path provided here. The path must end with ".pdf".
#'
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes geom_vline labs ggtitle theme element_line element_rect element_text unit margin annotate annotation_custom
#' @importFrom grid gpar textGrob
#' @importFrom cowplot plot_grid save_plot
#'
#' @details The main input `deepSTRAPP_outputs` is the typical output of [deepSTRAPP::run_deepSTRAPP_for_focal_time()].
#'   It provides information on results of a STRAPP test performed at a given `focal_time`.
#'
#'   Histograms are built based on the distribution of the test statistics.
#'   Such distributions are recorded in the outputs of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_for_focal_time()]
#'   when `return_perm_data = TRUE` so that the distributions of test stats computed across posterior samples are returned
#'   among the outputs under `$STRAPP_results$perm_data_df`.
#'
#'   For multinominal data (categorical or biogeographic data with more than 2 states), it is possible to plot the histograms of post hoc pairwise tests.
#'   Set `plot_posthoc_tests = TRUE` to generate histograms for all the pairwise post hoc Dunn's test across pairs of states.
#'   To achieve this, the `deepSTRAPP_outputs` input object must contain a `$STRAPP_results$posthoc_pairwise_tests$perm_data_array` element
#'   that summarizes test statistics computed across posterior samples for all pairwise post hoc tests.
#'   This is obtained from [deepSTRAPP::run_deepSTRAPP_for_focal_time()] when setting both `posthoc_pairwise_tests = TRUE` to carry out post hoc tests,
#'   and `return_perm_data = TRUE` to record distributions of test statistics.
#'
#'   Alternatively, the main input `deepSTRAPP_outputs` can be the output of [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   providing results of STRAPP tests over multiple time-steps. In this case, you must provide a `focal_time` to select the
#'   unique time-step used for plotting.
#'   * `return_perm_data` must be set to `TRUE` so that the permutation data used to compute the tests are returned among the outputs
#'     under `$STRAPP_results_over_time[[i]]$perm_data_df`.
#'   * `posthoc_pairwise_tests` must be set to `TRUE` so that the permutation data used to performed the post hoc tests are also
#'     returned among the outputs under `$STRAPP_results_over_time[[i]]$posthoc_pairwise_tests$perm_data_array`.
#'
#'  For plotting all time-steps at once, see [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()].
#'
#' @return By default, the function returns a list of classes `gg` and `ggplot`.
#'   This object is a ggplot that can be displayed on the console with `print(output)`.
#'   It corresponds to the histogram being displayed on the console when the function is run, if `display_plot = TRUE`, and can be further
#'   modify for aesthetics using the ggplot2 grammar.
#'
#'   If using multinominal data and set `plot_posthoc_tests = TRUE`, the function will return a list of objects.
#'   Each object is the ggplot associated with a pairwise post hoc test.
#'   To plot each histogram `i` individually, use `print(output_list[[i]])`.
#'   To plot all histograms at once in a multifaceted plot, as displayed on the console if `display_plot = TRUE`, use `cowplot::plot_grid(plotlist = output_list)`.
#'
#'   Each plot also displays summary statistics for the STRAPP test associated with the data displayed.
#'   * QX% = The quantile of null statistic distribution at the significant threshold used to define test significance. This is the value found on the red dashed line.
#'     The test will be considered significant (i.e., the null hypothesis is rejected) if this value is higher than zero (the black dashed line).
#'   * P-value =  The p-value of the STRAPP test which correspond the proportion of cases in which the statistics was lower than expected under the null hypothesis
#'     (i.e., the proportion of the histogram found below / on the left-side of the black dashed line).
#'   * N = The number of test statistics in the null distribution.
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plot. For post hoc tests, this will save the multifaceted plot.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::run_deepSTRAPP_for_focal_time()] [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()]
#'
#' @examples
#' if (deepSTRAPP::is_dev_version())
#' {
#'  # ----- Example 1: Continuous trait ----- #
#'
#'  # Load fake trait df
#'  data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#'  # Load phylogeny with old calibration
#'  data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
#'
#'  # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#'  data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Prepare trait data
#'
#'  # Extract continuous trait data as a named vector
#'  Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
#'                                      nm = Ponerinae_trait_tip_data$Taxa)
#'
#'  # Select a color scheme from lowest to highest values
#'  color_scale = c("darkgreen", "limegreen", "orange", "red")
#'
#'  # Get Ancestral Character Estimates based on a Brownian Motion model
#'  # To obtain values at internal nodes
#'  Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)
#'
#'  \donttest{ # (May take several minutes to run)
#'  # Run a Stochastic Mapping based on a Brownian Motion model
#'  # to interpolate values along branches and obtain a "contMap" object
#'  Ponerinae_contMap <- phytools::contMap(Ponerinae_tree, x = Ponerinae_cont_tip_data,
#'                                         res = 100, # Number of time steps
#'                                         plot = FALSE)
#'  # Plot contMap = stochastic mapping of continuous trait
#'  plot_contMap(contMap = Ponerinae_contMap,
#'               color_scale = color_scale)
#'
#'  ## Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#'  Ponerinae_deepSTRAPP_cont_old_calib_10My <- run_deepSTRAPP_for_focal_time(
#'    contMap = Ponerinae_contMap,
#'    ace = Ponerinae_ACE,
#'    tip_data = Ponerinae_cont_tip_data,
#'    trait_data_type = "continuous",
#'    BAMM_object = Ponerinae_BAMM_object_old_calib,
#'    focal_time = focal_time,
#'    rate_type = "net_diversification",
#'    uncertainty_strategy = "rates_only",
#'    return_perm_data = TRUE)
#'
#'  ## Explore output
#'  str(Ponerinae_deepSTRAPP_cont_old_calib_10My, max.level = 1)
#'
#'  # ----- Plot histogram of STRAPP overall test results from run_deepSTRAPP_for_focal_time() ----- #
#'
#'  histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_10My,
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_STRAPP_histogram_10My.pdf",
#'     plot_posthoc_tests = FALSE)
#'  # Adjust aesthetics a posteriori
#'  histogram_ggplot_adj <- histogram_ggplot +
#'     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#'  print(histogram_ggplot_adj) }
#'
#'  # ----- Plot histogram of STRAPP overall test results from run_deepSTRAPP_over_time() ----- #
#'
#'  ## Load directly outputs from run_deepSTRAPP_over_time()
#'  data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  # Select focal_time = 10My
#'  focal_time <- 10
#'
#'  ## Plot histogram for overall test
#'  histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     focal_time = focal_time,
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_STRAPP_histogram_10My.pdf",
#'     plot_posthoc_tests = FALSE)
#'  # Adjust aesthetics a posteriori
#'  histogram_ggplot_adj <- histogram_ggplot +
#'    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#'  print(histogram_ggplot_adj)
#'
#'  # ----- Example 2: Categorical trait ----- #
#'
#'  ## Load data
#'
#'  # Load phylogeny
#'  data(Ponerinae_tree, package = "deepSTRAPP")
#'  # Load trait df
#'  data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#'
#'  # Load the BAMM_object summarizing 1000 posterior samples of BAMM
#'  data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Prepare trait data
#'
#'  # Extract categorical data with 3-levels
#'  Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
#'                                          nm = Ponerinae_trait_tip_data$Taxa)
#'  table(Ponerinae_cat_3lvl_tip_data)
#'
#'  # Select color scheme for states
#'  colors_per_states <- c("forestgreen", "sienna", "goldenrod")
#'  names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Produce densityMaps using stochastic character mapping based on an ARD Mk model
#'  Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     phylo = Ponerinae_tree_old_calib,
#'     trait_data_type = "categorical",
#'     colors_per_states = colors_per_states,
#'     evolutionary_models = "ARD", # Use default ARD model
#'     nb_simulations = 100, # Reduce number of simulations to save time
#'     seed = 1234, # Seet seed for reproducibility
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE,
#'     plot_map = FALSE) }
#'
#'  # Load directly output
#'  data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
#'
#'  ## Set focal time to 10 Mya
#'  focal_time <- 10
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.
#'
#'  Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My <- run_deepSTRAPP_for_focal_time(
#'     densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
#'     nb_simulations = 100,
#'     ace = Ponerinae_cat_3lvl_data_old_calib$ace,
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     trait_data_type = "categorical",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     focal_time = focal_time,
#'     rate_type = "net_diversification",
#'     uncertainty_strategy = "paired",
#'     posthoc_pairwise_tests = TRUE,
#'     return_perm_data = TRUE,
#'     return_updated_Maps = TRUE,
#'     return_updated_BAMM_object = TRUE)
#'
#'  ## Explore output
#'  str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My, max.level = 1)
#'
#'  # ------ Plot histogram of STRAPP overall test results ------ #
#'
#'  histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My,
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_STRAPP_histogram_overall_test.pdf",
#'     plot_posthoc_tests = FALSE)
#'  # Adjust aesthetics a posteriori
#'  histogram_ggplot_adj <- histogram_ggplot +
#'     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#'  print(histogram_ggplot_adj)
#'
#'  # ------ Plot histograms of STRAPP post hoc test results ------ #
#'
#'  histograms_ggplot_list <- plot_histogram_STRAPP_test_for_focal_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My,
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_STRAPP_histograms_posthoc_tests.pdf",
#'     plot_posthoc_tests = TRUE)
#'  # Plot all histograms one by one
#'  print(histograms_ggplot_list)
#'  # Plot all histograms on one faceted plot
#'  cowplot::plot_grid(plotlist = histograms_ggplot_list) }
#' }
#'


plot_histogram_STRAPP_test_for_focal_time <- function (deepSTRAPP_outputs,
                                                       focal_time = NULL,
                                                       display_plot = TRUE,
                                                       plot_posthoc_tests = FALSE,
                                                       PDF_file_path = NULL)

{
  ### Check input validity
  {
    ## deepSTRAPP_outputs
    # Check presence of STRAPP_results
    if (is.null(deepSTRAPP_outputs$STRAPP_results) & is.null(deepSTRAPP_outputs$STRAPP_results_over_time))
    {
      stop(paste0("`deepSTRAPP_outputs` must have a `$STRAPP_results` or `$STRAPP_results_over_time` element.\n",
                  "Be sure to set `return_perm_data = TRUE` in [deepSTRAPP::run_deepSTRAPP_for_focal_time] or [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                  "This element is needed to plot the histogram of STRAPP test stats."))
    }

    ## Identify the type of inputs
    if (is.null(deepSTRAPP_outputs$STRAPP_results_over_time))
    {
      inputs_over_time <- FALSE
    } else {
      inputs_over_time <- TRUE
    }

    ## focal_time
    # Ensure a focal_time is provided if output is from [deepSTRAPP::run_deepSTRAPP_over_time()]
    if (inputs_over_time)
    {
      if (is.null(focal_time))
      {
        stop(paste0("You provided as input a `deepSTRAPP_outputs` object with multiple time-steps resulting from [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                    "You must provide a `focal_time` to select the appropriate time-step to be plotted.\n",
                    "For plotting all time-steps at once, please see [deepSTRAPP::plot_histogram_STRAPP_tests_over_time]."))
      }
      # Ensure focal_time match (any) time in deepSTRAPP_outputs
      if (!(focal_time %in% deepSTRAPP_outputs$time_steps))
      {
        stop(paste0("You provided as input a `deepSTRAPP_outputs` object with multiple time-steps resulting from [deepSTRAPP::run_deepSTRAPP_over_time].\n",
                    "You must provide a `focal_time` that matches with the `$time_steps` recorded in the `deepSTRAPP_outputs` object."))
      }
    } else {
      # Ensure focal_time match with the focal_time recorded in deepSTRAPP_outputs
      if (!is.null(focal_time))
      {
        if (!(focal_time %in% deepSTRAPP_outputs$focal_time))
        {
          stop(paste0("You provided as input a `deepSTRAPP_outputs` object with a unique time-step resulting from [deepSTRAPP::run_deepSTRAPP_for_focal_time].\n",
                      "However, the `focal_time` you provided does not that match with the `$focal_time` recorded in the `deepSTRAPP_outputs` object.\n",
                      "focal_time provided: ",focal_time,".\n",
                      "focal_time recorded in `deepSTRAPP_outputs`: ",deepSTRAPP_outputs$focal_time,"."))
        }
      }
    }

    ## Extract STRAPP_results
    if (!inputs_over_time)
    {
      # For outputs from run_deepSTRAPP_for_focal_time
      STRAPP_results <- deepSTRAPP_outputs$STRAPP_results
    } else {
      # For outputs from run_deepSTRAPP_over_time
      focal_time_ID <- which(deepSTRAPP_outputs$time_steps == focal_time)
      STRAPP_results <- deepSTRAPP_outputs$STRAPP_results_over_time[[focal_time_ID]]
    }

    ## STRAPP_results

    # STRAPP_results must have recorded a STRAPP test.
    if (STRAPP_results$trait_data_type_for_stats == "none")
    {
      if (!inputs_over_time)
      {
        stop(paste0("STRAPP test results are missing from 'deepSTRAPP_outputs$STRAPP_results'.\n",
                    "A unique ML state/range was inferred across branches for 'focal_time' = ",STRAPP_results$focal_time,".\n",
                    "No STRAPP test for differences in rates between states/ranges can be computed with a unique state/range.\n",
                    "Therefore, no histogram of test stats can be plotted for this 'focal_time'."))
      } else {
        stop(paste0("STRAPP test results are missing from 'deepSTRAPP_outputs$STRAPP_results_over_time' for the given 'focal_time'.\n",
                    "A unique ML state/range was inferred across branches for 'focal_time' = ",STRAPP_results$focal_time,".\n",
                    "No STRAPP test for differences in rates between states/ranges can be computed with a unique state/range.\n",
                    "Therefore, no histogram of test stats can be plotted for this 'focal_time'."))
      }
    }

    # STRAPP_results must have the needed elements: $focal_time, $estimate, $p_value, $method, $trait_data_type_for_stats, $perm_data_df
    # if (is.null(STRAPP_results$focal_time))
    # {
    #   stop(paste0("'$focal_time' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
    #               "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects."))
    # }
    # if (is.null(STRAPP_results$estimate))
    # {
    #   stop(paste0("'$estimate' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
    #               "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects."))
    # }
    # if (is.null(STRAPP_results$p_value))
    # {
    #   stop(paste0("'$p_value' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
    #               "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects."))
    # }
    # if (is.null(STRAPP_results$method))
    # {
    #   stop(paste0("'$method' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
    #               "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects."))
    # }
    # if (is.null(STRAPP_results$trait_data_type_for_stats))
    # {
    #   stop(paste0("'$trait_data_type_for_stats' is missing from 'STRAPP_results'. You can inspect the structure of the input object with 'str(STRAPP_results, 2)'.\n",
    #               "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects."))
    # }

    if (is.null(STRAPP_results$perm_data_df))
    {
      if (!inputs_over_time)
      {
        stop(paste0("'$perm_data_df' is missing from 'deepSTRAPP_outputs$STRAPP_results'.\n",
                    "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                    "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects.\n",
                    "Especially, check if you used 'return_perm_data = TRUE' to save the permutated data needed for the histogram plot."))
      } else {
        stop(paste0("'$perm_data_df' is missing from 'deepSTRAPP_outputs$STRAPP_results_over_time'.\n",
                    "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs$STRAPP_results_over_time, 2)'.\n",
                    "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                    "Especially, check if you used 'return_perm_data = TRUE' to save the permutated data needed for the histogram plot."))
      }

    }

    ## plot_posthoc_tests
    if (plot_posthoc_tests)
    {
      # plot_posthoc_tests = TRUE only for "multinominal" data
      if (STRAPP_results$trait_data_type_for_stats != "multinominal")
      {
        stop(paste0("'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                    "Set 'posthoc_pairwise_tests = FALSE', or provide 'deepSTRAPP_outputs' for a trait with more than two states/ranges for the given 'focal_time'."))
      }
      # Check if STRAPP_results$posthoc_pairwise_tests$summary_df is present.
      if (is.null(STRAPP_results$posthoc_pairwise_tests$summary_df))
      {
        if (!inputs_over_time)
        {
          stop(paste0("'$posthoc_pairwise_tests$summary_df' is missing from 'deepSTRAPP_outputs$STRAPP_results'.\n",
                      "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                      "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects.\n",
                      "Especially, check if you used 'posthoc_pairwise_tests = TRUE' to compute pairwise tests.\n",
                      "In addition, you must also select 'return_perm_data = TRUE' to save the permutated data needed for the histogram plots."))
        } else {
          stop(paste0("'$posthoc_pairwise_tests$summary_df' is missing from 'deepSTRAPP_outputs$STRAPP_results_over_time'.\n",
                      "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs$STRAPP_results_over_time, 2)'.\n",
                      "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                      "Especially, check if you used 'posthoc_pairwise_tests = TRUE' to compute pairwise tests.\n",
                      "In addition, you must also select 'return_perm_data = TRUE' to save the permutated data needed for the histogram plots."))
        }

      }
      # Check if $posthoc_pairwise_tests$perm_data_array is present.
      if (is.null(STRAPP_results$posthoc_pairwise_tests$perm_data_array))
      {
        if (!inputs_over_time)
        {
          stop(paste0("'$posthoc_pairwise_tests$perm_data_array' is missing from 'deepSTRAPP_outputs$STRAPP_results'.\n",
                      "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                      "See ?deepSTRAPP::run_deepSTRAPP_for_focal_time() to learn how to generate those objects.\n",
                      "Especially, check if you used 'posthoc_pairwise_tests = TRUE' AND 'return_perm_data = TRUE' to compute pairwise tests ",
                      "and save the permutated data needed for the histogram plots."))
        } else {
          stop(paste0("'$posthoc_pairwise_tests$perm_data_array' is missing from 'deepSTRAPP_outputs$STRAPP_results_over_time'.\n",
                      "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs$STRAPP_results_over_time, 2)'.\n",
                      "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                      "Especially, check if you used 'posthoc_pairwise_tests = TRUE' AND 'return_perm_data = TRUE' to compute pairwise tests ",
                      "and save the permutated data needed for the histogram plots."))
        }
      }
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
  }

  ## Save initial par() and reassign them on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (!plot_posthoc_tests)
  {
    ## Case for overall test plot

    # Extract name of the statistic
    stat_name <- names(STRAPP_results$perm_data_df)[ncol(STRAPP_results$perm_data_df)]
    # Extract null distribution of the statistic
    stat_null_distri <- STRAPP_results$perm_data_df[,stat_name]
    # Extract quantile of the critical threshold
    estimate_quantile <- names(STRAPP_results$estimate)
    # Extract the number of stats in the null distribution
    nb_stats <- nrow(STRAPP_results$perm_data_df)

    # Build ggplot object
    ggplot_histo <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

      # Generate histogram
      ggplot2::geom_histogram(mapping = ggplot2::aes(x = stat_null_distri),
                              bins = 30, # Set to avoid displaying automatic message from ggplot
                              fill = "grey80", color = "grey50", alpha = 0.8) +

      # Add vline for x = 0 (evaluation threshold)
      ggplot2::geom_vline(xintercept = 0, color = "black",
                          linetype = "dashed", linewidth = 1.0) +

      # Add vline for x = estimate (significance threshold)
      ggplot2::geom_vline(xintercept = STRAPP_results$estimate, color = "red",
                          linetype = "dashed", linewidth = 1.0) +

      # Add test summary
      # alpha value : Estimate,  p-value
      annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                   label = paste0("Q", estimate_quantile, " = ", round(STRAPP_results$estimate, digits = 3), "\n",
                                  "P-value = ", round(STRAPP_results$p_value, digits = 3), "\n",
                                  "N = ", nb_stats)) +

      # Adjust axis labels
      ggplot2::labs(x = paste0("Distribution of the test statistics"),
                    y = "Counts across posterior samples") +

      # Add title
      ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results$method," test\n",
                              "Focal time = ", STRAPP_results$focal_time)) +

      # Adjust aesthetics
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3),
                     panel.background = ggplot2::element_rect(fill = NA, color = NA),
                     plot.margin = ggplot2::unit(c(0.3, 0.5, 0.5, 0.5), "inches"),
                     plot.title = ggplot2::element_text(size = 20, hjust = 0.5, color = "black",
                                                        margin = ggplot2::margin(b = 15, t = 5)),
                     axis.title = ggplot2::element_text(size = 20, color = "black"),
                     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
                     axis.line = ggplot2::element_line(linewidth = 1.5),
                     axis.text = ggplot2::element_text(size = 18, color = "black"),
                     axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
                     axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)))

    ## Display plot if requested
    if (display_plot)
    {
      print(ggplot_histo)
    }

    ## Export plot if requested
    if (!is.null(PDF_file_path))
    {
      cowplot::save_plot(plot = ggplot_histo,
                         filename = PDF_file_path,
                         base_height = 8, base_width = 10)
    }

    ## Return ggplot
    return(invisible(ggplot_histo))

   } else {

    ## Case for post hoc tests plot

    # Create list of pairwise plots
    ggplots_histo_list <- list()

    # Extract pairs
    all_pairs <- dimnames(STRAPP_results$posthoc_pairwise_tests$perm_data_array)$pairs
    # Compute size factor to adjust size of everything
    size_factor <- sqrt(length(all_pairs))

    for (i in seq_along(all_pairs))
    {
      # i <- 1

      # Extract perm data for pair i
      perm_data_df_i <- as.data.frame(STRAPP_results$posthoc_pairwise_tests$perm_data_array[i, , ])

      # Extract pair names
      pair_i <- all_pairs[i]
      # # Replace "!=" with unicode
      # pair_i <- gsub(pattern = "!=", replacement = "\u2260", x = pair_i )

      # Extract name of the statistic
      stat_name <- names(perm_data_df_i)[ncol(perm_data_df_i)]
      # Extract null distribution of the statistic
      stat_null_distri <- unlist(perm_data_df_i[,stat_name])
      # Extract quantile of the critical threshold (same the one used for the main test)
      estimate_quantile <- names(STRAPP_results$estimate)
      # Extract estimate of the critical threshold
      estimate_value <- STRAPP_results$posthoc_pairwise_tests$summary_df$estimates[i]
      # Extract p-value
      p_value <- STRAPP_results$posthoc_pairwise_tests$summary_df$p_values[i]
      # Extract nb of stats in the null distribution
      nb_stats <- nrow(perm_data_df_i)

      # Detect if adjusted p-values differ from p-value
      p_value_adj_to_plot <- any(STRAPP_results$posthoc_pairwise_tests$summary_df$p_values != STRAPP_results$posthoc_pairwise_tests$summary_df$p_values_adjusted)
      # Extract p-value adjusted if different
      if (p_value_adj_to_plot)
      {
        p_value_adj <- STRAPP_results$posthoc_pairwise_tests$summary_df$p_values_adjusted[i]
      }

      # Build ggplot object
      ggplot_histo_i <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

        # Generate histogram
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = stat_null_distri),
                                bins = 30, # Set to avoid displaying automatic message from ggplot
                                fill = "grey80", color = "grey50", alpha = 0.8) +

        # Add vline for x = 0 (evaluation threshold)
        ggplot2::geom_vline(xintercept = 0, color = "black",
                            linetype = "dashed", linewidth = 1.0/size_factor) +

        # Add vline for x = estimate (significance threshold)
        ggplot2::geom_vline(xintercept = estimate_value, color = "red",
                            linetype = "dashed", linewidth = 1.0/size_factor) +

        # Add test summary
        # alpha value : Estimate,  p-value
        annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18/size_factor),
                     label = paste0("Q", estimate_quantile, " = ", round(estimate_value, digits = 3), "\n",
                                    "P-value = ", round(p_value, digits = 3), "\n",
                                    "N = ", nb_stats)) +

        # Adjust axis labels
        ggplot2::labs(x = paste0("Distribution of the test statistics"),
                      y = "Counts across posterior samples") +

        # Add title
        ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results$posthoc_pairwise_tests$method," test\n",
                                "Focal time = ", STRAPP_results$focal_time, "\n",
                                "Hypothesis = ",pair_i, "\n")) +

        # Adjust aesthetics
        ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey70", linetype = "dashed", linewidth = 0.3/size_factor),
                       panel.background = ggplot2::element_rect(fill = NA, color = NA),
                       plot.margin = ggplot2::unit(c(0.3/(size_factor*2), 0.5/(size_factor*2), 0.5/(size_factor*2), 0.5/(size_factor*2)), "inches"),
                       plot.title = ggplot2::element_text(size = 20/size_factor, hjust = 0.5, color = "black",
                                                          margin = ggplot2::margin(b = 15/size_factor, t = 5/size_factor)),
                       axis.title = ggplot2::element_text(size = 20/size_factor, color = "black"),
                       axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10/size_factor)),
                       axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12/size_factor)),
                       axis.line = ggplot2::element_line(linewidth = 1.5/size_factor),
                       axis.text = ggplot2::element_text(size = 18/size_factor, color = "black"),
                       axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5/size_factor)),
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5/size_factor)))

      # Add adjusted p-value if present
      if (p_value_adj_to_plot)
      {
        ggplot_histo_i <- ggplot_histo_i +
          # Add test summary including adjusted p-value
          annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18/size_factor),
                       label = paste0("Q", estimate_quantile, " = ", round(estimate_value, digits = 3), "\n",
                                      "P-value = ", round(p_value, digits = 3), "\n",
                                      "P-value adj = ", round(p_value_adj, digits = 3)))
      }

      # Store ggplot for pairs i
      ggplots_histo_list[[i]] <- ggplot_histo_i

    }

    ## Display plot if requested
    if (display_plot)
    {
      cow_plot_to_print <- cowplot::plot_grid(plotlist = ggplots_histo_list)
      print(cow_plot_to_print)
    }

    ## Export plot if requested
    if (!is.null(PDF_file_path))
    {
      cowplot::save_plot(plot = cowplot::plot_grid(plotlist = ggplots_histo_list),
                         filename = PDF_file_path,
                         base_height = 8, base_width = 10)
    }

    ## Return list of ggplots
    return(invisible(ggplots_histo_list))
  }
}



### Helper function to enable the use of "npc" units in ggplot2::annotate()

#' @noRd

annotate_npc <- function(label, x, y, ...)
{
  ggplot2::annotation_custom(
    grob = grid::textGrob(x = ggplot2::unit(x, "npc"),
                          y = ggplot2::unit(y, "npc"),
                          label = label, ...))
}
