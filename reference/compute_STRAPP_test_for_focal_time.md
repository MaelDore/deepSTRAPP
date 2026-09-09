# Compute STRAPP to test for a relationship between diversification rates and trait data

Carries out the appropriate statistical method to test for a
relationship between diversification rates and trait data for a given
point in the past (i.e. the `focal_time`). Tests are based on
block-permutations: rates data are randomized across tips following
blocks defined by the diversification regimes identified on each tip
(typically from a BAMM).

Such tests are called STructured RAte Permutations on Phylogenies
(STRAPP) as described in Rabosky, D. L., & Huang, H. (2016). A robust
semi-parametric test for detecting trait-dependent diversification.
Systematic biology, 65(2), 181-193.
[doi:10.1093/sysbio/syv066](https://doi.org/10.1093/sysbio/syv066) .

The function is an extension of the original
[`BAMMtools::traitDependentBAMM()`](https://rdrr.io/pkg/BAMMtools/man/traitDependentBAMM.html)
function used to carry out STRAPP test on extant time-calibrated
phylogenies.

Tests can be carried out on speciation, extinction and net
diversification rates.

`deepSTRAPP::compute_STRAPP_test_for_focal_time()` can handle three
types of statistical tests depending on the type of trait data provided:

### Continuous trait data

Tests for correlations between trait and rates carried out with
`deepSTRAPP::compute_STRAPP_test_for_continuous_data()`. The associated
test is the Spearman's rank correlation test (See
[stats::cor.test](https://rdrr.io/r/stats/cor.test.html)).

### Binary trait data

For categorical and biogeographic trait data that have only two states
(ex: 'Nearctic' vs. 'Neotropics'). Tests for differences in rates
between states are carried out with
`deepSTRAPP::compute_STRAPP_test_for_binary_data()`. The associated test
is the Mann-Whitney-Wilcoxon rank-sum test (See
[stats::wilcox.test](https://rdrr.io/r/stats/wilcox.test.html)).

### Multinominal trait data

For categorical and biogeographic trait data with more than two states
(ex: 'No leg' vs. 'Two legs' vs. 'Four legs'). Tests for differences in
rates between states are carried out with
`deepSTRAPP::compute_STRAPP_test_for_multinomial_data()`. The associated
test for all states is the Kruskal-Wallis H test (See
[stats::kruskal.test](https://rdrr.io/r/stats/kruskal.test.html)). If
`posthoc_pairwise_tests = TRUE`, post hoc pairwise tests between pairs
of states will be carried out too. The associated test for post hoc
pairwise tests is Dunn's post hoc pairwise rank-sum test (See
[dunn.test::dunn.test](https://rdrr.io/pkg/dunn.test/man/dunn.test.html)).

## Usage

``` r
compute_STRAPP_test_for_focal_time(
  BAMM_object,
  trait_data_list,
  rate_type = "net_diversification",
  uncertainty_strategy = "paired",
  trait_maps_vs_BAMM_samples_list = NULL,
  seed = NULL,
  nb_permutations = NULL,
  alpha = 0.05,
  two_tailed = TRUE,
  one_tailed_hypothesis = NULL,
  posthoc_pairwise_tests = FALSE,
  p.adjust_method = "none",
  return_perm_data = FALSE,
  nthreads = 1,
  print_hypothesis = TRUE
)
```

## Arguments

- BAMM_object:

  Object of class `"bammdata"`, typically generated with
  [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
  that contains a phylogenetic tree and associated diversification rates
  across selected posterior samples updated to a specific time in the
  past (i.e. the `focal_time`).

- trait_data_list:

  List obtained from
  [`extract_most_likely_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_most_likely_trait_values_for_focal_time.md)
  or
  [`extract_all_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_all_trait_values_for_focal_time.md)
  that contains at least a `$trait_data` element, a `$focal_time`
  element, and a `$trait_data_type`. `$trait_data` is a named vector
  with the trait data found on the phylogeny at `focal_time`.
  `$focal_time` informs on the time in the past at which the trait and
  rates data will be tested. `$trait_data_type` informs on the type of
  trait data: 'continuous', 'categorical', or 'biogeographic'.

- rate_type:

  A character string specifying the type of diversification rates to
  use. Must be one of 'speciation', 'extinction' or
  'net_diversification' (default).

- uncertainty_strategy:

  Character string. To select the strategy used to account for
  uncertainty in estimates.

  - `"rates_only"`: Only accounts for diversification-rate uncertainty
    across BAMM posterior samples. Uses ML estimates for continuous
    traits and the most frequent state/range observed across stochastic
    maps for categorical and biogeographic data.

  - `"paired"`: Default option. Accounts for both diversification-rate
    and ancestral trait/range reconstruction uncertainty by pairing BAMM
    posterior samples with stochastic maps. When the number of BAMM
    samples and stochastic maps differ, random pairing with replacement
    from the smaller set is used so that all posterior samples and
    stochastic maps contribute to the analysis.

  - `"full"`: Exhaustive option that accounts for trait/range- and rate-
    uncertainty by crossing all BAMM posterior samples with all
    stochastic maps. Accounts for both diversification-rate and
    ancestral reconstruction uncertainty by evaluating every combination
    of BAMM posterior sample and stochastic map. WARNING: This
    exhaustive approach can substantially increase computation time and
    memory requirements and is therefore recommended only for
    moderate-sized analyses.

- trait_maps_vs_BAMM_samples_list:

  (Optional) Data.frame of two variables manually providing the names to
  associate stochastic maps (`$trait_map_ID`) with BAMM samples
  (`$BAMM_posterior_sample_ID`). This is typically used to ensure the
  same stochastic maps and BAMM samples are used to test across multiple
  time-steps. Values are the names of the objects such as "Map_X" and
  "BAMM_X". Default = `NULL`. \* For uncertainty_strategy ==
  'rates_only', the `$trait_map_ID` must be "Map_ML" as only the ML
  estimates of trait values/states/ranges are used. \* For
  uncertainty_strategy == 'paired', each pair of stochastic map and BAMM
  sample will be used once. \* For uncertainty_strategy == 'full', all
  stochastic maps will be matched with all BAMM samples. Those may
  partly differ from the actual maps and BAMM samples used for the test,
  as recorded in `$perm_data_df`, because invalid maps with not enough
  states/ranges are discarded.

- seed:

  Integer. Set the seed to ensure reproducibility. Default is `NULL` (a
  random seed is used).

- nb_permutations:

  Integer. To select the number of random permutations to perform during
  the tests. If NULL (default), all BAMM posterior samples will be used
  once.

- alpha:

  Numeric. Significance level to use to compute the `estimate`
  corresponding to the values of the test statistic used to assess
  significance of the test. This does NOT affect p-values. Default is
  `0.05`.

- two_tailed:

  Logical. To define the type of tests. If `TRUE` (default), tests for
  correlations/differences in rates will be carried out with a null
  hypothesis that rates are not correlated with trait values (continuous
  data) or equal between trait states (categorical and biogeographic
  data). If `FALSE`, one-tailed tests are carried out.

  - For continuous data, it involves defining a `one_tailed_hypothesis`
    testing for either a "positive" or "negative" correlation under the
    alternative hypothesis.

  - For binary data (two states), it involves defining a
    `one_tailed_hypothesis` indicating which states have higher rates
    under the alternative hypothesis.

  - For multinomial data (more than two states), it defines the type of
    post hoc pairwise tests to carry out between pairs of states. If
    `posthoc_pairwise_tests = TRUE`, all two-tailed (if
    `two_tailed = TRUE`) or one-tailed (if `two_tailed = FALSE`) tests
    are automatically carried out.

- one_tailed_hypothesis:

  A character string specifying the alternative hypothesis in the
  one-tailed test. For continuous data, it is either "negative" or
  "positive" correlation. For binary data, it lists the trait states
  with states ordered in increasing rates under the alternative
  hypothesis, separated by a greater-than such as c('A \> B').

- posthoc_pairwise_tests:

  Logical. Only for multinomial data (with more than two states). If
  `TRUE`, all possible post hoc pairwise (Dunn) tests will be computed
  across all pairs of states. This is a way to detect which pairs of
  states have significant differences in rates if the overall test
  (Kruskal-Wallis) is significant. Default is `FALSE`.

- p.adjust_method:

  A character string. Only for multinomial data (with more than two
  states). It specifies the type of correction to apply to the p-values
  in the post hoc pairwise tests to account for multiple comparisons.
  See [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  the available methods. Default is `none`.

- return_perm_data:

  Logical. Whether to return the stats data computed from the posterior
  samples for observed and permuted data in the output. This is needed
  to plot the histogram of the null distribution used to assess
  significance of the test with
  [`plot_histogram_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histogram_STRAPP_test_for_focal_time.md).
  Default is `FALSE`.

- nthreads:

  Integer. Number of threads to use for parallel computing of the tests
  across the permutations. The R package `parallel` must be loaded for
  `nthreads > 1`. Default is `1`.

- print_hypothesis:

  Logical. Whether to print information on what test is carried out,
  detailing the null and alternative hypotheses, what significance level
  is used to reject or not the null hypothesis, and how uncertainty in
  trait estimates is handled. Default is `TRUE`.

## Value

The function returns a list with at least eleven elements.

Summary elements for the main test:

- `$estimate` Named numeric. Value of the test statistic used to assess
  significance of the test according to the significance level provided
  (`alpha`). The test is significant if `$estimate` is higher than zero.

- `$stats_median` Numeric. Median value of the distribution of test
  statistics.

- `$nb_test_stats` Integer. Number of test stats in the distribution.

- `$p-value` Numeric. P-value of the test. The test is considered
  significant if `$p-value` is lower than `alpha`.

- `$method` Character string. The statistical method used to carry out
  the test.

- `$rate_type` Character string. The type of diversification rates
  tested. One of 'speciation', 'extinction' or 'net_diversification'.

- `$trait_data_type` Character string. The type of trait data as found
  in 'trait_data_list\$trait_data_type'. One of 'continuous',
  'categorical', or 'biogeographic'.

- `$trait_data_type_for_stats` Character string. The type of trait data
  used to select statistical method. One of 'continuous', 'binary', or
  'multinomial'.

- `$uncertainty_strategy` Character string. The strategy used to account
  for uncertainty in estimates.

- `$trait_maps_vs_BAMM_samples_list` List of two elements recording the
  stochastic maps (`$trait_map_ID`) and BAMM samples
  (`$BAMM_posterior_sample_ID`) chosen for testing. Those may partly
  differ from the actual maps and BAMM samples used for the test as
  recorded in `$perm_data_df` because invalid maps with not enough
  states/ranges are discarded.

- `$focal_time` The time in the past at which the trait and rates data
  were tested.

If using continuous or binary data:

- `$two-tailed` Logical. Record the type of test used: two-tailed if
  `TRUE`, one-tailed if `FALSE`. If `one_tailed_hypothesis` is provided
  (only for continuous and binary trait data):

- `$one_tailed_hypothesis` Character string. Record of the alternative
  hypothesis used for the one-tailed tests.

If `posthoc_pairwise_tests = TRUE` (only for multinomial trait data):

- `$posthoc_pairwise_tests` List of at least 3 sub-elements:

  - `$summary_df` Data.frame of six variables providing the summary
    results of post hoc pairwise tests

  - `$method` Character string. The statistical method used to carry out
    the test. Here, "Dunn".

  - `$two-tailed` Logical. Record the type of post hoc pairwise tests
    used: two-tailed if `TRUE`, one-tailed if `FALSE`.

If `return_perm_data = TRUE`, the stats data computed from the posterior
samples for observed and permuted data are provided. This is needed to
plot the histogram of the null distribution used to assess significance
of the test with
[`plot_histogram_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histogram_STRAPP_test_for_focal_time.md).

- `$perm_data_df` A data.frame with four variables summarizing the data
  generated during the STRAPP test:

  - `$trait_map_ID` Integer. ID of the stochastic map from which trait
    data were extracted and used for the STRAPP test.

  - `$BAMM_posterior_sample_ID` Integer. ID of the posterior samples
    randomly drawn and used for the STRAPP test.

  - `$*_obs` Numeric. Test stats computed from the observed data in the
    posterior samples. Name depends on the test used.

  - `$*_perm` Numeric. Test stats computed from the permuted data in the
    posterior samples. Name depends on the test used.

  - `$delta_*` OR `$abs_delta_*` Numeric. Test stats computed for the
    STRAPP test comparing observed stats and permuted stats. Name
    depends on the test used and the type of tests (two-tailed compare
    absolute values; one-tailed compare raw values). Combined with
    `posthoc_pairwise_tests = TRUE`, the stats data are also provided
    for the post hoc pairwise tests:

- `$posthoc_pairwise_tests$perm_data_array` A 3D array containing stats
  data for all post hoc pairwise tests in a similar format to
  `$perm_data_df`.

If no STRAPP test was performed in the case of categorical/biogeographic
data with a single state/range at `focal_time`, only the
`$trait_data_type`, `$trait_data_type_for_stats` = "none", and
`$focal_time` are returned.

## Details

This set of functions carries out the STructured RAte Permutations on
Phylogenies (STRAPP) test as defined in Rabosky, D. L., & Huang, H.
(2016). A robust semi-parametric test for detecting trait-dependent
diversification. Systematic biology, 65(2), 181-193.

It is an extension of the original
[`BAMMtools::traitDependentBAMM()`](https://rdrr.io/pkg/BAMMtools/man/traitDependentBAMM.html)
function used to carry out STRAPP test on extant time-calibrated
phylogenies, but allowing here to test for differences/correlations at
any point in the past (i.e. the `focal_time`).

It takes an object of class `"bammdata"` (`BAMM_object`) that was
updated such that its diversification rates (`$tipLambda` and `$tipMu`)
and regimes (`$tipStates`) are reflecting values observed at a specific
time in the past (i.e. the `$focal_time`). Similarly, it takes a list
(`trait_data_list`) that provides `$trait_data` as observed on branches
at the same `focal_time` as the diversification rates and regimes.

A STRAPP test is carried out by drawing a random set of posterior
samples from the `BAMM_object`, then randomly permuting rates across
blocks of tips defined by the macroevolutionary regimes. Test statistics
are then computed across the initial observed data and the permuted data
for each sample. In a two-tailed test, the p-value is the proportion of
posterior samples in which the test stats is more extreme in the
permuted than in the observed data. In a one-tailed test, the p-value is
the proportion of posterior samples in which the test stats is higher in
the permuted than in the observed data.

———- Major changes compared to
[`BAMMtools::traitDependentBAMM()`](https://rdrr.io/pkg/BAMMtools/man/traitDependentBAMM.html)
———-

- Allow to account for uncertainty in trait estimates by pairing/mapping
  multiple trait data extracted from stochastic maps across the BAMM
  samples. See the `uncertainty_strategy` argument for details.

- Add post hoc pairwise tests (Dunn test) for multinomial data. Use
  `posthoc_pairwise_tests = TRUE`.

- Provide outputs tailored for histogram plots
  [`plot_histogram_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histogram_STRAPP_test_for_focal_time.md)
  and p-value time-series plots
  [`plot_STRAPP_pvalues_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_STRAPP_pvalues_over_time.md).

- Add prints detailing what test is carried out, what are the null and
  alternative hypotheses, and what significance level is used to reject
  or not the null hypothesis. (Enabled with `print_hypothesis = TRUE`).

- Split the function in multiple sub-functions according to the type of
  data (`$trait_data_type`).

- Prevent using Pearson's correlation tests and applying
  log-transformation for continuous data. The rationale is that there is
  no reason to assume that tip rates are distributed normally or
  log-normally. Thus, a Spearman's rank correlation test is favored.

## References

For STRAPP: Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric
test for detecting trait-dependent diversification. Systematic biology,
65(2), 181-193.
[doi:10.1093/sysbio/syv066](https://doi.org/10.1093/sysbio/syv066) .

For STRAPP in deep times: Doré, M., Borowiec, M. L., Branstetter, M. G.,
Camacho, G. P., Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B.
B., (2025), Evolutionary history of ponerine ants highlights how the
timing of dispersal events shapes modern biodiversity, Nature
Communications.
[doi:10.1038/s41467-025-63709-3](https://doi.org/10.1038/s41467-025-63709-3)

## See also

Associated functions in deepSTRAPP:
[`extract_most_likely_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_most_likely_trait_values_for_focal_time.md)
[`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md)

Original function in BAMMtools:
[`BAMMtools::traitDependentBAMM()`](https://rdrr.io/pkg/BAMMtools/man/traitDependentBAMM.html)

Statistical tests:
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)
[`dunn.test::dunn.test()`](https://rdrr.io/pkg/dunn.test/man/dunn.test.html)

For a guided tutorial, see this vignette:
[`vignette("explore_STRAPP_test_types", package = "deepSTRAPP")`](https://maeldore.github.io/deepSTRAPP/articles/explore_STRAPP_test_types.md)

## Author

Maël Doré

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
 # ------ Prepare data ------ #

 ## Load the BAMM_object summarizing 1000 posterior samples of BAMM with diversification rates
 # for ponerine ants extracted for 10My ago.
 data(Ponerinae_BAMM_object_10My, package = "deepSTRAPP")
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
# Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 # Plot the associated phylogeny with mapped rates
 plot_BAMM_rates(Ponerinae_BAMM_object_10My)

 ## Load the object containing head width trait data for ponerine ants extracted for 10My ago.
 data(Ponerinae_trait_cont_tip_data_10My, package = "deepSTRAPP")

 # Plot the associated contMap (continuous trait stochastic map)
 plot_contMap(Ponerinae_trait_cont_tip_data_10My$contMap)

 # Check that objects are ordered in the same fashion
 identical(names(Ponerinae_BAMM_object_10My$tipStates[[1]]),
           names(Ponerinae_trait_cont_tip_data_10My$trait_data))

 # Save continuous data
 trait_data_continuous <- Ponerinae_trait_cont_tip_data_10My

 ## Transform trait data into binary and multinomial data

 # Binarize data into two states
 trait_data_binary <- trait_data_continuous
 trait_data_binary$trait_data[trait_data_continuous$trait_data < 0.5] <- "state_A"
 trait_data_binary$trait_data[trait_data_continuous$trait_data >= 0.5] <- "state_B"
 trait_data_binary$trait_data_type <- "categorical"

 table(trait_data_binary$trait_data)

 # Categorize data into three states
 trait_data_multinomial <- trait_data_continuous
 trait_data_multinomial$trait_data[trait_data_continuous$trait_data < 0.6] <- "state_B"
 trait_data_multinomial$trait_data[trait_data_continuous$trait_data < 0.4] <- "state_A"
 trait_data_multinomial$trait_data[trait_data_continuous$trait_data >= 0.6] <- "state_C"
 trait_data_multinomial$trait_data_type <- "categorical"

 table(trait_data_multinomial$trait_data)

 ## Duplicate trait data as if extracted from multiple stochastic maps

 trait_data_continuous_multimaps <- Ponerinae_trait_cont_tip_data_10My
 trait_data_initial <- trait_data_continuous_multimaps$trait_data

 trait_data_multimaps <- list()
 for (i in 1:10)
 {
   # Add a bit of randomness across all stochastic maps data
   trait_data_multimaps[[i]] <- trait_data_initial +
       rnorm(n = length(trait_data_initial), mean = 0, sd = sd(trait_data_initial) / 10)
 }
 names(trait_data_multimaps) <- paste0("Map_", 1:10)
 trait_data_continuous_multimaps$trait_data <- trait_data_multimaps


  # (May take several minutes to run)
 # ------ Compute STRAPP test for continuous data ------ #

 plot(x = trait_data_continuous$trait_data, y = Ponerinae_BAMM_object_10My$tipLambda[[1]])

 # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
 # between "net_diversification" rates and trait data
 STRAPP_results <- compute_STRAPP_test_for_focal_time(
    BAMM_object = Ponerinae_BAMM_object_10My,
    trait_data_list = trait_data_continuous,
    uncertainty_strategy = "rates_only",
    two_tailed = FALSE,
    one_tailed_hypothesis = "negative",
    return_perm_data = TRUE)
 str(STRAPP_results, max.level = 2)
 # Data from the posterior samples is available in STRAPP_results$perm_data_df
 head(STRAPP_results$perm_data_df)

 # ------ Compute STRAPP test for binary data ------ #

 # Compute STRAPP test under the alternative hypothesis that "state_A" is associated
 # with higher "net_diversification" that "state_B"
 STRAPP_results <- compute_STRAPP_test_for_focal_time(
    BAMM_object = Ponerinae_BAMM_object_10My,
    trait_data_list = trait_data_binary,
    uncertainty_strategy = "rates_only",
    two_tailed = FALSE,
    one_tailed_hypothesis = c("state_A > state_B"))
 str(STRAPP_results, max.level = 1)

 # Compute STRAPP test under the alternative hypothesis that "state_B" is associated
 # with higher "net_diversification" that "state_A"
 STRAPP_results <- compute_STRAPP_test_for_focal_time(BAMM_object = Ponerinae_BAMM_object_10My,
    trait_data_list = trait_data_binary,
    uncertainty_strategy = "rates_only",
    two_tailed = FALSE,
    one_tailed_hypothesis = c("state_B > state_A"))
 str(STRAPP_results, max.level = 1)

 # ------ Compute STRAPP test for multinomial data ------ #

 # Compute STRAPP test between all three states, and compute post hoc tests
 # for differences in rates between all possible pairs of states
 # with a p-value adjusted for multiple comparison using Bonferroni's correction
 STRAPP_results <- compute_STRAPP_test_for_focal_time(
    BAMM_object = Ponerinae_BAMM_object_10My,
    trait_data_list = trait_data_multinomial,
    uncertainty_strategy = "rates_only",
    posthoc_pairwise_tests = TRUE,
    two_tailed = TRUE,
    p.adjust_method = "bonferroni",
    return_perm_data = TRUE)
 str(STRAPP_results, max.level = 2)
 # All post hoc pairwise test summaries are available in $summary_df
 STRAPP_results$posthoc_pairwise_tests$summary_df

 # ------ Compute STRAPP test with the 'paired' strategy ------ #

 # Account for uncertainty in trait estimates
 # by pairing stochastic maps with BAMM posterior samples

 # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
 # between "net_diversification" rates and trait data
 STRAPP_results <- compute_STRAPP_test_for_focal_time(
   BAMM_object = Ponerinae_BAMM_object_10My,
   trait_data_list = trait_data_continuous_multimaps,
   uncertainty_strategy = "paired",
   two_tailed = FALSE,
   one_tailed_hypothesis = "negative",
   return_perm_data = TRUE)
 str(STRAPP_results, max.level = 2)
 # Data from the paired stochastic maps and BAMM posterior samples
 # is available in STRAPP_results$perm_data_df
 head(STRAPP_results$perm_data_df)
 # Tests were performed across 1000 iterations since the stochastic maps
 # have been paired with the 1000 BAMM samples
 nrow(STRAPP_results$perm_data_df)
 # Each of the 10 maps was paired ca. 100 times with a unique BAMM sample
 table(STRAPP_results$perm_data_df$trait_map_ID)

 # ------ Compute STRAPP test with the 'full' strategy ------ #

 # Account for uncertainty in trait estimates
 # by crossing stochastic maps with BAMM posterior samples

 # Compute STRAPP test under the alternative hypothesis of a "negative" correlation
 # between "net_diversification" rates and trait data
 STRAPP_results <- compute_STRAPP_test_for_focal_time(
   BAMM_object = Ponerinae_BAMM_object_10My,
   trait_data_list = trait_data_continuous_multimaps,
   uncertainty_strategy = "full",
   two_tailed = FALSE,
   one_tailed_hypothesis = "negative",
   return_perm_data = TRUE)
 str(STRAPP_results, max.level = 2)
 # Data from the combined stochastic maps X BAMM posterior samples
 # is available in STRAPP_results$perm_data_df
 head(STRAPP_results$perm_data_df)
 # Tests were performed across 10000 iterations since the 10 stochastic maps
 # have been combined with the 1000 BAMM samples
 nrow(STRAPP_results$perm_data_df)
 # Each of the 10 maps was paired with all 1000 BAMM samples
 table(STRAPP_results$perm_data_df$trait_map_ID)
 
}
```
