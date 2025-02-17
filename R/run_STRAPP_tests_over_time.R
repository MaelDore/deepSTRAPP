
# Make a wrapped-up to be run over time-series
# run_STRAPP_tests_over_time()

### Outputs

# Main output: melted data.frame for p-value time-series plot
    # Focal_time, P-value, Qalpha% Δ stat (estimate)
# Optional output (raw data): List of outputs from the compute_STRAPP_test_for_focal_time() function which is already wrapped in run_STRAPP_test_for_focal_time
# Optional output (for downstream analyses):
  # Melted data.frame needed to plot histogram with plot_histogram_STRAPP_tests_over_time()
     # Focal-time, BAMM ID, Δ stat distribution
     # Will be used to plot an histogram for each  focal_time
     # Include option of a PDF output
     # Include option of a merged PDF output with one histogram per page
  # Melted data.frame with rates/regimes data
     # Include a run of extract_diversification_data_for_focal_time() within each loop of run_STRAPP_test_for_focal_time()
     # Will be used by plot_rates_through_time().
  # Melted data.frame of trait data
     # Focal-time, branch ID, trait value/state
     # Will be used by plot_rates_through_time().
# Optional output (raw data):
  # List of updated BAMM from update_rates_and_regimes_for_focal_time() => make it clear it is heavy. Can be used to plot rates on cut phylo with plot.bammdata()
  # List of updated contMap/simmaps from extract_most_likely_trait/states_for_focal_time() => make it clear it is heavy. Can be used to plot traits on cut phylo with plot.contMap()/plot.simmaps()

