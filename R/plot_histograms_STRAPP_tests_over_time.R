
#' @title Plot multiple histograms of STRAPP test statistics over time-steps
#'
#' @description Plot an histogram of the distribution of the test statistics
#'   obtained from a deepSTRAPP workflow carried out for each focal time in `$time_steps`.
#'   Main input = output of a deepSTRAPP run over time using [deepSTRAPP::run_deepSTRAPP_over_time()]).
#'
#'   Returns one histogram for overall tests for each focal time in `$time_steps`.
#'   If `plot_posthoc_tests = TRUE`, it will return one faceted plot with an histogram
#'   per post hoc tests for each focal time in `$time_steps`.
#'
#'   If a PDF file path is provided in `PDF_file_path`, the plots will be saved directly in a PDF file,
#'   with one page per focal time in `$time_steps`.
#'
#' @param deepSTRAPP_outputs List of elements generated with [deepSTRAPP::run_deepSTRAPP_over_time()],
#'   that summarize the results of multiple deepSTRAPP across `$time_steps`. It needs to include the `$STRAPP_results_over_time`
#'   element with `$perm_data_df` obtained when setting both `return_STRAPP_results = TRUE` and `return_perm_data = TRUE`.
#' @param display_plots Logical. Whether to display the histograms generated in the R console. Default is `TRUE`.
#' @param plot_posthoc_tests Logical. For multinominal data only. Whether to plot the histograms for the overall Kruskal-Wallis test across all states (`plot_posthoc_tests = FALSE`),
#'   or plot the histograms for all the pairwise post hoc Dunn's tests across pairs of states (`plot_posthoc_tests = TRUE`).
#'   Time-steps at which the data does not yield more than two states/ranges will show a warning and generate no plot. Default is `FALSE`.
#' @param PDF_file_path Character string. If provided, the plots will be saved in a unique PDF file following the path provided here. The path must end with ".pdf".
#'   Each page of the PDF corresponds to a focal time in `$time_steps`.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_histogram aes geom_vline labs ggtitle theme element_line element_rect element_text unit margin annotate annotation_custom
#' @importFrom grid gpar
#' @importFrom cowplot plot_grid save_plot
#' @importFrom grDevices pdf dev.off
#'
#' @details The main input `deepSTRAPP_outputs` is the typical output of [deepSTRAPP::run_deepSTRAPP_over_time()].
#'   It provides information on results of a STRAPP tests performed over multiple time-steps.
#'
#'   Histograms are built based on the distribution of the test statistics.
#'   Such distributions are recorded in the outputs of a deepSTRAPP run carried out with [deepSTRAPP::run_deepSTRAPP_over_time()]
#'   when `return_STRAPP_results = TRUE` AND `return_perm_data = TRUE`. The `$STRAPP_results_over_time` objects provided within the input are lists that must contain
#'   a `$perm_data_df` element that summarizes test statistics computed across posterior samples.
#'
#'   For multinominal data (categorical or biogeographic data with more than 2 states), it is possible to plot the histograms of post hoc pairwise tests.
#'   Set `plot_posthoc_tests = TRUE` to generate histograms for all the pairwise post hoc Dunn's test across pairs of states.
#'   To achieve this, the `$STRAPP_results_over_time` objects must contain a `$posthoc_pairwise_tests$perm_data_array` element that summarizes test statistics
#'   computed across posterior samples for all pairwise post hoc tests. This is obtained from [deepSTRAPP::run_deepSTRAPP_over_time()] when setting
#'   `return_STRAPP_results = TRUE` to return the STRAPP results, `posthoc_pairwise_tests = TRUE` to carry out post hoc tests,
#'   and `return_perm_data = TRUE` to record distributions of test statistics.
#'   Time-steps for which the data do not yield more than two states/ranges will show a warning and generate no plot.
#'
#' @return By default, the function returns a list of sub-lists of classes `gg` and `ggplot` ordered as in `$time_steps`.
#'   Each sub-list corresponds to a ggplot for a given `focal_time` `i` that can be displayed on the console with `print(output[[i]])`.
#'   If `display_plots = TRUE`, the histograms are being displayed on the console one by one while generated.
#'
#'   If using multinominal data and set `plot_posthoc_tests = TRUE`, the function will return a list of sub-lists of objects ordered as in `$time_steps`.
#'   Each sub-list is a list of the ggplots associated with pairwise post hoc tests carried out for this a given `focal_time`.
#'   For a given `focal_time` i, to plot each histogram j individually, use `print(output_list[[i]][[j]])`.
#'   To plot all histograms of a given `focal_time` `i` at once in a multifaceted plot, as displayed sequentially on the console if `display_plots = TRUE`,
#'   use `cowplot::plot_grid(plotlist = output_list[[i]])`.
#'
#'   Each plot also displays summary statistics for the STRAPP test associated with the data displayed.
#'   * The quantile of null statistic distribution at the significant threshold used to define test significance. This is the value found on the red dashed line.
#'     The test will be considered significant (i.e., the null hypothesis is rejected) if this value is higher than zero (the black dashed line).
#'   * The p-value of the STRAPP test which correspond the proportion of cases in which the statistics was lower than expected under the null hypothesis
#'     (i.e., the proportion of the histogram found below / on the left-side of the black dashed line).
#'
#'   If a `PDF_file_path` is provided, the function will also generate a PDF file of the plots with one page per `$time_steps`.
#'   For post hoc tests, this will save the multifaceted plots.
#'
#' @author Maël Doré
#'
#' @seealso Associated functions in deepSTRAPP: [deepSTRAPP::run_deepSTRAPP_over_time()] [deepSTRAPP::plot_histogram_STRAPP_test_for_focal_time()]
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
#'  ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
#'  # nb_time_steps <- 5
#'  time_step_duration <- 5
#'  time_range <- c(0, 40)
#'
#'  ## Run deepSTRAPP on net diversification rates
#'  Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     contMap = Ponerinae_contMap,
#'     ace = Ponerinae_ACE,
#'     tip_data = Ponerinae_cont_tip_data,
#'     trait_data_type = "continuous",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     # nb_time_steps = nb_time_steps,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     return_perm_data = TRUE,
#'     extract_trait_data_melted_df = TRUE,
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_trait_data_with_Map = TRUE,
#'     return_updated_BAMM_object = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly trait data output
#'  data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Plot histograms of STRAPP overall test results
#'  # Tests are Spearman's rank correlation tests
#'
#'  # Plot all histograms
#'  histogram_ggplots <- plot_histograms_STRAPP_tests_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_STRAPP_histogram_overall_test.pdf",
#'     plot_posthoc_tests = FALSE)
#'
#'  # Print histogram for time step 1 = 0 My
#'  print(histogram_ggplots[[1]])
#'  # Adjust aesthetics of plot for time step 1 a posteriori
#'  histogram_ggplot_adj <- histogram_ggplots[[1]] +
#'     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#'  print(histogram_ggplot_adj)
#'
#'  # ----- Example 2: Categorical data ----- #
#'
#'  ## Load data
#'
#'  # Load trait df
#'  data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
#'  # Load phylogeny
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
#'     colors_per_levels = colors_per_states,
#'     evolutionary_models = "ARD",
#'     nb_simulations = 100,
#'     return_best_model_fit = TRUE,
#'     return_model_selection_df = TRUE,
#'     plot_map = FALSE) }
#'
#'  # Load directly trait data output
#'  data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")
#'
#'  ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
#'  # nb_time_steps <- 5
#'  time_step_duration <- 5
#'  time_range <- c(0, 40)
#'
#'  \donttest{ # (May take several minutes to run)
#'  ## Run deepSTRAPP on net diversification rates across time-steps.
#'  Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- run_deepSTRAPP_over_time(
#'     densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
#'     ace = Ponerinae_cat_3lvl_data_old_calib$ace,
#'     tip_data = Ponerinae_cat_3lvl_tip_data,
#'     trait_data_type = "categorical",
#'     BAMM_object = Ponerinae_BAMM_object_old_calib,
#'     # nb_time_steps = nb_time_steps,
#'     time_range = time_range,
#'     time_step_duration = time_step_duration,
#'     rate_type = "net_diversification",
#'     seed = 1234, # Set for reproducibility
#'     alpha = 0.10, # Select a generous level of significance for the sake of the example
#'     posthoc_pairwise_tests = TRUE,
#'     return_perm_data = TRUE,
#'     extract_trait_data_melted_df = TRUE,
#'     extract_diversification_data_melted_df = TRUE,
#'     return_STRAPP_results = TRUE,
#'     return_updated_trait_data_with_Map = TRUE,
#'     return_updated_BAMM_object = TRUE,
#'     verbose = TRUE,
#'     verbose_extended = TRUE) }
#'
#'  ## Load directly deepSTRAPP output
#'  data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, package = "deepSTRAPP")
#'  ## This dataset is only available in development versions installed from GitHub.
#'  # It is not available in CRAN versions.
#'  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.
#'
#'  ## Explore output
#'  str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, max.level = 1)
#'
#'  ## Plot histograms of STRAPP overall test results #
#'  # Tests are Kruskall-Wallis H tests when more than two states/ranges are present.
#'  # Tests are Mann–Whitney–Wilcoxon rank-sum tests when only two states/ranges are present.
#'
#'  histogram_ggplots <- plot_histograms_STRAPP_tests_over_time(
#'     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#'     display_plot = TRUE,
#'     # PDF_file_path = "./plot_STRAPP_histograms_overall_tests.pdf",
#'     plot_posthoc_tests = FALSE)
#'
#'  # Print histogram for time step 1 = 0 My
#'  print(histogram_ggplots[[1]])
#'  # Adjust aesthetics of plot for time step 1 a posteriori
#'  histogram_ggplot_adj <- histogram_ggplots[[1]] +
#'      ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
#'  print(histogram_ggplot_adj)
#'
#'  ## Plot histograms of STRAPP post hoc test results ------ #
#'  # Tests are Dunn's multiple comparison pairwise post hoc tests possible
#'  # only when more than two states/ranges are present.
#'
#'  histograms_ggplots_list <- plot_histograms_STRAPP_tests_over_time(
#'      deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
#'      display_plot = TRUE,
#'      # PDF_file_path = "./plot_STRAPP_histograms_posthoc_tests.pdf",
#'      plot_posthoc_tests = TRUE)
#'
#'  # Print all histograms for time step 1 (= 0 My) one by one
#'  print(histograms_ggplots_list[[1]])
#'  # Plot all histograms for time step 1 (= 0 My) on one faceted plot
#'  cowplot::plot_grid(plotlist = histograms_ggplots_list[[1]])
#' }
#'


plot_histograms_STRAPP_tests_over_time <- function (
    deepSTRAPP_outputs,
    display_plots = TRUE,
    plot_posthoc_tests = FALSE,
    PDF_file_path = NULL)

{
  ### Check input validity
  {
    ## deepSTRAPP_outputs$STRAPP_results_over_time
    # deepSTRAPP_outputs must have element $STRAPP_results_over_time
    if (is.null(deepSTRAPP_outputs$STRAPP_results_over_time))
    {
      stop(paste0("'$STRAPP_results_over_time' is missing from 'deepSTRAPP_outputs'. You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                  "Especially, check if you used 'return_STRAPP_results = TRUE' to save the STRAPP_results for each time step needed for the histogram plots."))
    }
    ## deepSTRAPP_outputs$time_steps
    # deepSTRAPP_outputs must have element $time_steps
    if (is.null(deepSTRAPP_outputs$time_steps))
    {
      stop(paste0("'$time_steps' is missing from 'deepSTRAPP_outputs'. You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects."))
    }
    # Extract $focal_time from deepSTRAPP_outputs$STRAPP_results_over_time
    focal_time_list <- unlist(lapply(X = deepSTRAPP_outputs$STRAPP_results_over_time, FUN = function (x) { x$focal_time } ))
    # deepSTRAPP_outputs$time_steps must match with focal_time found in deepSTRAPP_outputs$STRAPP_results_over_time
    if (!all(deepSTRAPP_outputs$time_steps == focal_time_list))
    {
      stop(paste0("'deepSTRAPP_outputs$time_steps' does not match with 'focal_time' found in 'deepSTRAPP_outputs$STRAPP_results_over_time'.\n",
                  "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects."))
    }

    ## deepSTRAPP_outputs$STRAPP_results_over_time is a list.
    # All elements are STRAPP_results that should have all the needed sub-elements: $focal_time, $estimate, $p_value, $method, $trait_data_type_for_stats, $perm_data_df

    # Detect time-steps with no STRAPP test results
    trait_data_type_for_stats_list <- unlist(lapply(X = deepSTRAPP_outputs$STRAPP_results_over_time, FUN = function (x) { x$trait_data_type_for_stats }))
    time_steps_with_STRAPP_results_ID <- which(trait_data_type_for_stats_list != "none")
    time_steps_with_no_STRAPP_results <- focal_time_list[which(trait_data_type_for_stats_list == "none")]

    # Send warning for time-steps without STRAPP results
    if (any(trait_data_type_for_stats_list == "none"))
    {
      if (sum(trait_data_type_for_stats_list == "none") == 1)
      {
        warning(paste0("STRAPP test results were missing from 'deepSTRAPP_outputs$STRAPP_results_over_time'.\n",
                       "A unique ML state/range was inferred across branches for this time-step: ",time_steps_with_no_STRAPP_results,".\n",
                       "No STRAPP test for differences in rates between states/ranges can be computed with a unique state/range.\n",
                       "Therefore, no histogram of test stats can be plotted for this time-step."))
      } else {
        warning(paste0("STRAPP test results were missing from 'deepSTRAPP_outputs$STRAPP_results_over_time'.\n",
                       "A unique ML state/range was inferred across branches for these time-steps: ",paste(time_steps_with_no_STRAPP_results, collapse = ", "),".\n",
                       "No STRAPP test for differences in rates between states/ranges can be computed with a unique state/range.\n",
                       "Therefore, no histogram of test stats can be plotted for these time-steps."))
      }
    }

    # Extract deepSTRAPP_outputs$STRAPP_results_over_time only for time-steps with STRAPP results
    STRAPP_results_over_time_with_STRAPP_results <- deepSTRAPP_outputs$STRAPP_results_over_time[time_steps_with_STRAPP_results_ID]

    # Check for presence of $perm_data_df
    perm_data_df_check_list <- unlist(lapply(X = STRAPP_results_over_time_with_STRAPP_results, FUN = function (x) { "perm_data_df" %in% names(x) } ))
    if (!all(perm_data_df_check_list))
    {
      stop(paste0("Elements in 'deepSTRAPP_outputs$STRAPP_results_over_time' do not have a '$perm_data_df'.\n",
                  "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs$STRAPP_results_over_time, 2)'.\n",
                  "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                  "Especially, check if you used 'return_STRAPP_results = TRUE' AND 'return_perm_data = TRUE' to save the permutated data",
                  "needed for the histogram plot to save the permutated data needed for the histogram plot."))
    }

    ## plot_posthoc_tests
    if (plot_posthoc_tests)
    {
      # Extract $trait_data_type_for_stats from deepSTRAPP_outputs$STRAPP_results_over_time
      trait_data_type_for_stats <- deepSTRAPP_outputs$trait_data_type_for_stats
      # plot_posthoc_tests = TRUE only for "multinominal" data
      if (trait_data_type_for_stats != "multinominal")
      {
        stop(paste0("'posthoc_pairwise_tests = TRUE' only makes sense for categorical/biogeographic data with more than two states/ranges.\n",
                    "Set 'posthoc_pairwise_tests = FALSE', or provide 'deepSTRAPP_outputs' for a trait with more than two states/ranges'.\n",
                    "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects."))
      }

      # Check if $posthoc_pairwise_tests is present in elements of deepSTRAPP_outputs$STRAPP_results_over_time.
      posthoc_pairwise_tests_check_list <- unlist(lapply(X = deepSTRAPP_outputs$STRAPP_results_over_time, FUN = function (x) { "posthoc_pairwise_tests" %in% names(x) } ))
      posthoc_pairwise_tests_check_list_with_STRAPP_results <- unlist(lapply(X = STRAPP_results_over_time_with_STRAPP_results, FUN = function (x) { "posthoc_pairwise_tests" %in% names(x) } ))
      if (all(!posthoc_pairwise_tests_check_list_with_STRAPP_results))
      {
        stop(paste0("Elements in 'deepSTRAPP_outputs$STRAPP_results_over_time' do not have a '$posthoc_pairwise_tests'.\n",
                    "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs$STRAPP_results_over_time, 2)'.\n",
                    "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                    "Especially, check if you used 'return_STRAPP_results = TRUE' AND 'posthoc_pairwise_tests = TRUE' to run post hoc tests ",
                    "and save the results in 'deepSTRAPP_outputs'."))
      } else { # Case with '$posthoc_pairwise_tests' missing only for some time-steps
        if (!all(posthoc_pairwise_tests_check_list))
        {
          time_steps_no_posthoc_tests <- focal_time_list[!posthoc_pairwise_tests_check_list]
          warning(paste0("No posthoc pairwise tests recorded for these time-steps: ", paste(time_steps_no_posthoc_tests, collapse = ", "),".\n",
                        "This is likely because two or less states/ranges where present at these time-steps.\n",
                         "No histogram of test stats was generated for these time-steps."))
        }
      }

      # Extract only STRAPP_results for time-step with multinominal data
      STRAPP_results_over_time_for_post_hoc <- deepSTRAPP_outputs$STRAPP_results_over_time[posthoc_pairwise_tests_check_list]

      # Check if $posthoc_pairwise_tests$perm_data_array is present in elements of STRAPP_results_over_time_for_post_hoc.
      posthoc_pairwise_tests_list <- lapply(X = STRAPP_results_over_time_for_post_hoc, FUN = function (x) { x$posthoc_pairwise_tests } )
      perm_data_array_check_list <- unlist(lapply(X = posthoc_pairwise_tests_list, FUN = function (x) { "perm_data_array" %in% names(x) } ))
      if (all(!perm_data_array_check_list))
      {
        warning(paste0("Elements in 'deepSTRAPP_outputs$STRAPP_results_over_time' do not have a '$posthoc_pairwise_tests$perm_data_array' despite recording posthoc test results otherwise.\n",
                       "You can inspect the structure of the input object with 'str(deepSTRAPP_outputs$STRAPP_results_over_time, 3)'.\n",
                       "See ?deepSTRAPP::run_deepSTRAPP_over_time() to learn how to generate those objects.\n",
                       "Especially, check if you used 'return_STRAPP_results = TRUE', 'posthoc_pairwise_tests = TRUE' AND return_perm_data = TRUE to run post hoc tests ",
                       "and save the results and the permutated data needed for the histogram plots in 'deepSTRAPP_outputs'."))
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


  if (!plot_posthoc_tests)
  {
    ## Case for overall test plot

    # Initiate list of ggplots
    ggplot_histo_list <- list()

    # Loop per time steps with recorded STRAPP results
    for (i in time_steps_with_STRAPP_results_ID)
    {
      # i <- 1

      # Extract STRAPP_results
      STRAPP_results_i <- deepSTRAPP_outputs$STRAPP_results_over_time[[i]]

      # Extract name of the statistic
      stat_name <- names(STRAPP_results_i$perm_data_df)[4]
      # Extract null distribution of the statistic
      stat_null_distri <- STRAPP_results_i$perm_data_df[,stat_name]
      # Extract quantile of the critical threshold
      estimate_quantile <- names(STRAPP_results_i$estimate)

      # Build ggplot object
      ggplot_histo_i <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

        # Generate histogram
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = stat_null_distri),
                                bins = 30, # Set to avoid displaying automatic message from ggplot
                                fill = "grey80", color = "grey50", alpha = 0.8) +

        # Add vline for x = 0 (evaluation threshold)
        ggplot2::geom_vline(xintercept = 0, color = "black",
                            linetype = "dashed", linewidth = 1.0) +

        # Add vline for x = estimate (significance threshold)
        ggplot2::geom_vline(xintercept = STRAPP_results_i$estimate, color = "red",
                            linetype = "dashed", linewidth = 1.0) +

        # Add test summary
        # alpha value : Estimate,  p-value
        annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18),
                     label = paste0("Q", estimate_quantile, " = ", round(STRAPP_results_i$estimate, digits = 3), "\n",
                                    "P-value = ", round(STRAPP_results_i$p_value, digits = 3))) +

        # Adjust axis labels
        ggplot2::labs(x = paste0("Distribution of the test statistics"),
                      y = "Counts across posterior samples") +

        # Add title
        ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results_i$method," test\n",
                                "Focal time = ", STRAPP_results_i$focal_time)) +

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

      ## Display plot at each time step if requested
      if (display_plots)
      {
        print(ggplot_histo_i)
      }

      ## Store in list of ggplots
      ggplot_histo_list[[i]] <- ggplot_histo_i
    }

    ## Export all plots in one PDF if requested
    if (!is.null(PDF_file_path))
    {
      grDevices::pdf(file = file.path(PDF_file_path), height = 8, width = 10, onefile = TRUE)
      for (i in seq_along(ggplot_histo_list))
      {
        print(ggplot_histo_list[[i]])
      }
      grDevices::dev.off()
    }

    ## Return list of ggplots
    return(invisible(ggplot_histo_list))

  } else {

    ## Case for post hoc tests plot

    # Initiate list of lists of ggplots
    ggplot_histo_list_of_lists <- list()

    # Detect ID of time-steps with recorded posthoc tests
    time_steps_ID_with_posthoc_tests <- which(posthoc_pairwise_tests_check_list)

    # Loop per time steps
    # for (i in seq_along(deepSTRAPP_outputs$time_steps))
    for (i in seq_along(time_steps_ID_with_posthoc_tests)) # Run only across time-steps with recorded posthoc tests
    {
      # i <- 1

      # Extract ID of the time-step
      time_steps_ID_i <- time_steps_ID_with_posthoc_tests[i]

      # Extract STRAPP_results
      STRAPP_results_i <- deepSTRAPP_outputs$STRAPP_results_over_time[[time_steps_ID_i]]

      # Create list of pairwise plots for time step i
      ggplots_histo_list_i <- list()

      # Extract pairs
      all_pairs <- dimnames(STRAPP_results_i$posthoc_pairwise_tests$perm_data_array)$pairs
      # Compute size factor to adjust size of everything
      size_factor <- sqrt(length(all_pairs))

      for (j in seq_along(all_pairs))
      {
        # j <- 1

        # Extract perm data for pair j
        perm_data_df_j <- as.data.frame(STRAPP_results_i$posthoc_pairwise_tests$perm_data_array[j, , ])

        # Extract pair names
        pair_j <- all_pairs[j]
        # # Replace "!=" with unicode
        # pair_j <- gsub(pattern = "!=", replacement = "\u2260", x = pair_j )

        # Extract name of the statistic
        stat_name <- names(perm_data_df_j)[3]
        # Extract null distribution of the statistic
        stat_null_distri <- perm_data_df_j[,stat_name]
        # Extract quantile of the critical threshold (same the one used for the main test)
        estimate_quantile <- names(STRAPP_results_i$estimate)
        # Extract estimate of the critical threshold
        estimate_value <- STRAPP_results_i$posthoc_pairwise_tests$summary_df$estimates[j]
        # Extract p-value
        p_value <- STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values[j]

        # Detect if adjusted p-values differ from p-value
        p_value_adj_to_plot <- any(STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values != STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values_adjusted)
        # Extract p-value adjusted if different
        if (p_value_adj_to_plot)
        {
          p_value_adj <- STRAPP_results_i$posthoc_pairwise_tests$summary_df$p_values_adjusted[j]
        }

        # Build ggplot object for pair j
        ggplot_histo_j <- ggplot2::ggplot(data = as.data.frame(stat_null_distri)) +

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
                                      "P-value = ", round(p_value, digits = 3))) +

          # Adjust axis labels
          ggplot2::labs(x = paste0("Distribution of the test statistics"),
                        y = "Counts across posterior samples") +

          # Add title
          ggplot2::ggtitle(paste0("STRAPP based on ",STRAPP_results_i$posthoc_pairwise_tests$method," test\n",
                                  "Focal time = ", STRAPP_results_i$focal_time, "\n",
                                  "Hypothesis = ",pair_j, "\n")) +

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
          ggplot_histo_j <- ggplot_histo_j +
            # Add test summary including adjusted p-value
            annotate_npc(x = 0.05, y = 0.95, hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 18/size_factor),
                         label = paste0("Q", estimate_quantile, " = ", round(estimate_value, digits = 3), "\n",
                                        "P-value = ", round(p_value, digits = 3), "\n",
                                        "P-value adj = ", round(p_value_adj, digits = 3)))
        }

        # Store ggplot for pairs j
        ggplots_histo_list_i[[j]] <- ggplot_histo_j

      }

      ## Display plot at each time step if requested
      if (display_plots)
      {
        cow_plot_to_print_i <- cowplot::plot_grid(plotlist = ggplots_histo_list_i)
        print(cow_plot_to_print_i)
      }

      ## Store in list of lists of ggplots
      ggplot_histo_list_of_lists[[i]] <- ggplots_histo_list_i
    }

    ## Export all plots in one PDF if requested
    if (!is.null(PDF_file_path))
    {
      # # Find the number of rows and columns
      # nb_plots <- length(cow_plot_to_print_i$layers)
      # nb_cols <- ceiling(sqrt(nb_plots))
      # nb_rows <- ceiling(nb_plots/nb_cols)

      # Export PDF witn one page per plot
      grDevices::pdf(file = file.path(PDF_file_path), height = 8, width = 10, onefile = TRUE)
      # pdf(file = PDF_file_path, height = 8*nb_rows, width = 10*nb_cols, onefile = TRUE)
      for (i in seq_along(ggplot_histo_list_of_lists))
      {
        cow_plot_to_print_i <- cowplot::plot_grid(plotlist = ggplot_histo_list_of_lists[[i]])
        print(cow_plot_to_print_i)
      }
      grDevices::dev.off()
    }

    ## Return list of lists of ggplots
    return(invisible(ggplot_histo_list_of_lists))
  }
}

