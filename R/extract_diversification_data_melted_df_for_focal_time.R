#' @title Extracts diversification data from a BAMM_object
#'
#' @description Extracts regimes ID and tip rates from a `BAMM_object` that have been
#'   updated to provide diversification data for a specific time in the past (i.e. the `focal_time`).
#'   Use [deepSTRAPP::update_rates_and_regimes_for_focal_time()] to obtain
#'   a `BAMM_object` updated for a given `focal_time`.
#'
#' @param BAMM_object Object of class `"bammdata"`, typically generated with
#'   [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'   that contains a phylogenetic tree and associated diversification rate mapping
#'   across selected posterior samples.
#' @param verbose Logical. Should progression be displayed? A message will be printed for every batch of
#'   100 BAMM posterior samples extracted. Default is `TRUE`.
#'
#' @export
#' @importFrom tidyr pivot_longer
#'
#' @return Returns a data.frame with five columns.
#'
#'   * `$focal_time` Integer. The time, in terms of time distance from the present, at which the trait data were extracted. Should be equal for all rows
#'   as a unique BAMM_object updated for a unique `focal_time` is being extracted.
#'   * `$tip_ID` Character string. Tip labels of the branches cut-off at `focal_time`.
#'     + If `keep_tip_labels = TRUE` was used in [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'     cut-off branches with a single descendant tip retain their initial `tip.label`.
#'     + If `keep_tip_labels = FALSE` was used in [deepSTRAPP::update_rates_and_regimes_for_focal_time()],
#'     all cut-off branches are labeled using their tipward node ID.
#'   * `$regime_ID` Integer. The regime ID on tips at `focal_time`. IDs are integer. The root process equals "1", then they are incremented from oldest to youngest.
#'   Regime IDs are independent across posterior samples.
#'   * `$rate_type` Character string. Type of rates: "lambda" for speciation rates, "mu" for extinction rates,
#'   and "net_diversification" for net diversification rates (lambda - mu).
#'   * `$rates` Numerical. Rates in \\[number of events / branch / evolutionary time\\].
#'
#' @author Maël Doré
#'
#' @examples
#' # Load the BAMM_object summarizing 1000 posterior samples of BAMM updated for a focal_time of 10 My
#' data(Ponerinae_BAMM_object_10My)
#'
#' # Extract diversification data
#' diversification_data_df <- extract_diversification_data_melted_df_for_focal_time(
#'    BAMM_object = Ponerinae_BAMM_object_10My,
#'    verbose = TRUE)
#' # Print output
#' head(diversification_data_df)
#'

extract_diversification_data_melted_df_for_focal_time <- function (BAMM_object, verbose = TRUE)
{
  # BAMM_object must be a 'bammdata' object
  # Number of posterior sample data must be equal between $tipStates, $tipLambda and $tipMu
  # Must have a $focal_time element. If not, send a warning message saying it is likely not an output of update_rates_and_regimes_for_focal_time(). Must go through that function, even if local_time is set to 0 My (current time).

  ## Extract focal_time
  focal_time <- BAMM_object$focal_time

  ## Extract tip.labels
  tip_labels <- names(BAMM_object$tipStates[[1]])

  # Initiate unmelted data.frame
  unmelted_df <- data.frame(focal_time = NULL, tip_ID = NULL, regime_ID = NULL, lambda = NULL, mu = NULL)

  ## Loop per posterior samples to extract tip data
  for (i in seq_along(BAMM_object$tipStates))
  {
    # i <- 1

    # Extract diversification data for posterior sample i
    unmelted_df_i <- data.frame(focal_time = focal_time, tip_ID = tip_labels,
                                regime_ID = BAMM_object$tipStates[[i]],
                                lambda = BAMM_object$tipLambda[[i]], mu = BAMM_object$tipMu[[i]])

    # Merge data
    unmelted_df <- rbind(unmelted_df, unmelted_df_i)

    if (verbose & (i %% 100 == 0))
    {
      cat(paste0(Sys.time(), " - Tip states/rates extracted for BAMM posterior sample n\u00B0", i, "/", length(BAMM_object$tipStates),"\n"))
    }
  }

  ## Compute net_diversification
  unmelted_df$net_diversification <- unmelted_df$lambda - unmelted_df$mu

  ## Melt data.frame in a tidy way so that one rate value = one row
  melted_df <- unmelted_df |>
    tidyr::pivot_longer(cols = c("lambda", "mu", "net_diversification"),
                        names_to = "rate_type", values_to = "rates")

  ## Return melted df
  return(melted_df)
}

