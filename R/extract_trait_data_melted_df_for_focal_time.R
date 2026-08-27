#' @title Extract trait data from a trait_data object
#'
#' @description Extract trait data from a trait_data_list object as produced within the deepSTRAPP workflow
#'   for a specific time in the past (i.e. the `focal_time`) and produce a melted dataset recording trait values
#'   across (stochastic maps x) edges for the associated `focal_time`.
#'
#' @param trait_data_list Object summarizing trait data extracted for a given `focal_time`,
#'   including a `$trait_data` element storing trait data recorded across (stochastic maps x) edges for the associated `$focal_time`.
#'   Use [deepSTRAPP::extract_all_trait_values_for_focal_time()] to obtain trait data across multiple stochastic maps.
#'   Use [deepSTRAPP::extract_most_likely_trait_values_for_focal_time()] to obtain only the most likely trait data values.
#'
#' @export
#' @importFrom tidyr pivot_longer last_col
#'
#' @return Returns a data.frame with five columns.
#'
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted. Should be equal for all rows
#'   as a unique BAMM_object updated for a unique `focal_time` is being extracted.
#'   * `$Map_ID` Character string. ID of the stochastic from which the trait data are extracted.
#'     If using 'rate_only' strategy to account for uncertainty in ancestral estimate, this is fixed to "Map_ML".
#'     If using 'paired' or 'full' strategies, this records either true stochastic maps "Map_X", or dummy maps "Dummy_map_X",
#'     depending if you provided respectively stochastic maps (`contMaps`/`simmaps`) or simply `densityMaps` as inputs.
#'   * `$tip_ID` Character string. Tip labels of the branches cut-off at `focal_time`.
#'     + If `keep_tip_labels = TRUE` was used in [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'     cut-off branches with a single descendant tip retain their initial `tip.label`.
#'     + If `keep_tip_labels = FALSE` was used in [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'     all cut-off branches are labeled using their tipward node ID.
#'   * `$trait_data_type` Character string. Records the type of traits: "continuous", "categorical", or "biogeographic".
#'   * `$trait_value` Numerical or Character string.
#'     + For "continuous" traits: numerical values.
#'     + For "categorical" traits: states.
#'     + For "biogeographic" data: ranges.
#'
#' @author Maël Doré
#'
#' @examples
#' # ----- Example 1: Extract ML estimates data ----- #
#'
#' ## Load categorical trait data mapped on a phylogeny
#' data(eel_cat_3lvl_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_cat_3lvl_data, 1)
#' eel_cat_3lvl_data$densityMaps # Three density maps: one per state
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' ## Extract ML estimates of ancestral states for the given focal_time
#'
#' # Extract from the densityMaps
#' eel_cat_3lvl_data_10My <- extract_most_likely_trait_values_for_focal_time(
#'    densityMaps = eel_cat_3lvl_data$densityMaps,
#'    trait_data_type = "categorical",
#'    focal_time = focal_time)
#'
#' ## Format ML states data as a melted df
#'
#' melted_df <- extract_trait_data_melted_df_for_focal_time(trait_data_list = eel_cat_3lvl_data_10My)
#' head(melted_df)
#'
#'
#' # ----- Example 2: Extract data across multiple stochastic maps ----- #
#'
#' ## Load biogeographic range data mapped on a phylogeny
#' data(eel_biogeo_data, package = "deepSTRAPP")
#'
#' # Explore data
#' str(eel_biogeo_data, 1)
#' length(eel_biogeo_data$simmaps) # 100 stochastic maps: one per simulated biogeographic history
#'
#' # Set focal time to 10 Mya
#' focal_time <- 10
#'
#' ## Extract range data across stochastic maps for the given focal_time
#'
#' # Extract from the simmaps
#' eel_biogeo_data_10My <- extract_all_trait_values_for_focal_time(
#'    simmaps = eel_biogeo_data$simmaps,
#'    trait_data_type = "biogeographic",
#'    focal_time = focal_time)
#'
#' ## Format range data as a melted df
#'
#' melted_df <- extract_trait_data_melted_df_for_focal_time(trait_data_list = eel_biogeo_data_10My)
#' head(melted_df)
#'

extract_trait_data_melted_df_for_focal_time <- function (trait_data_list)
{
  ### Check input validity
  {
    ## trait_data_list
    # Must have a $trait_data element. '$trait_data_type' elements
    if (is.null(trait_data_list$trait_data))
    {
      stop(paste0("'trait_data_list' must have a $trait_data element recording trait data extracted for the associated focal_time.\n",
                  "See ?deepSTRAPP::extract_all_trait_values_for_focal_time() and ?deepSTRAPP::extract_most_likely_trait_values_for_focal_time()\n",
                  "to learn how to generate a valid 'trait_data_list' to use as input."))
    }

    # Must have a $focal_time element. If not, send a message saying it is likely not an output of update_rates_and_regimes_for_focal_time(). Must go through that function, even if local_time is set to 0 My (current time).
    if (is.null(trait_data_list$focal_time))
    {
      stop(paste0("'trait_data_list' must have a $focal_time element indicating for which 'focal_time' the trait data have been extracted.\n",
                  "This also applies if investigating the present ('focal_time' = 0).\n",
                  "See ?deepSTRAPP::extract_all_trait_values_for_focal_time() and ?deepSTRAPP::extract_most_likely_trait_values_for_focal_time()\n",
                  "to learn how to generate a valid 'trait_data_list' to use as input."))
    }

    # Must have a $trait_data_type element.
    if (is.null(trait_data_list$trait_data_type))
    {
      stop(paste0("'trait_data_list' must have a $trait_data_type element indicating which type of data have been extracted.\n",
                  "See ?deepSTRAPP::extract_all_trait_values_for_focal_time() and ?deepSTRAPP::extract_most_likely_trait_values_for_focal_time()\n",
                  "to learn how to generate a valid 'trait_data_list' to use as input."))
    }
  }

  ## Extract elements
  trait_data <- trait_data_list$trait_data
  focal_time <- trait_data_list$focal_time
  trait_data_type <- trait_data_list$trait_data_type

  ## Check if trait data is recorded for ML estimates (as a named vector) or across multiple stochastic maps (as a list of named vectors)
  if (is.list(trait_data_list$trait_data))
  {
    data_is_ML_estimates <- FALSE
  } else {
    data_is_ML_estimates <- TRUE
  }

  ## Extract data for ML estimates
  if (data_is_ML_estimates)
  {
    melted_df <- data.frame(focal_time = focal_time, Map_ID = "Map_ML", tip_ID = names(trait_data), trait_data_type = trait_data_type, trait_value = trait_data)
    row.names(melted_df) <- NULL
  }

  ## Extract data across multiple stochastic maps
  if (!data_is_ML_estimates)
  {
    # Aggregate data as matrix of Map x edge
    unmelted_df <- as.data.frame(do.call(rbind, trait_data))
    unmelted_df$Map_ID <- row.names(unmelted_df)
    # Melt data.frame in a tidy way so that one trait value = one row
    melted_df <- unmelted_df |>
      tidyr::pivot_longer(cols = -tidyr::last_col(),
                          names_to = "tip_ID", values_to = "trait_value") |>
      as.data.frame()
    # Add missing columns
    melted_df$focal_time <- focal_time
    melted_df$trait_data_type <- trait_data_type
    # Reorganize columns
    melted_df <- melted_df[, c("focal_time", "Map_ID", "tip_ID", "trait_data_type", "trait_value")]
  }

  ## Return melted df
  return(melted_df)
}

