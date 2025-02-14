

# See also BAMMtools::traitDependentBAMM stats::wilcox.test stats::cor.test stats::kruskal.test dunn.test::dunn.test stats::p.adjust
# Dependency too

### Add parallel package as Suggested in DESCRIPTION?


# State in the description that this is a wrapped up of the original function traitDependentBAMM written by Dan Rabosky & Huateng Huang, 2014 in the BAMMtools package
# Describe changes from the original
  # Additions
    # Use the names of the tipStates, tipLambda, and tipMu to check compatibility between tips rather than the tip.label in the phylogeny
    # Allows to use the function to compute STRAPP test for any given time
    # Allows to choose if random sampling of posterior configurations must be done with replacement or not.
    # Default to using one permutation per posterior sample if not requested otherwise.
    # Add post hoc pairwise tests (Dunn test) for multinominal data
    # Provide outputs tailored for histogram plots [ref_fun] and p-value time-series plots [ref_fun]
    # Add print_hypothesis parameter to display or not the hypothesis tested.
  # Changes
    # Function split in multiple sub-functions according to $data_type for clarity
    # Prevent using Pearson's correlation tests and applying log-transformation as their is no rationale
    # to assume that tip rates are distributed normally or log-normally. Favor Spearman's rank correlation tests.



### Inputs
# From extract_most_likely_trait_values_for_focal_time() => trait_data named integer (may also be a list with two items, trait_data + a contMap)
# From update_rates_and_regimes_for_focal_time() => Updated bammdata list with $tipStates, $tipLambda and $tipMu holding data from all posterior samples

# Make a master function run_STRAPP_test()
  # Prepare the input data and select the proper sub-function to be ran
# Make sub-function for each type of data
  # run_STRAPP_test_for_continuous_data()
  # run_STRAPP_test_for_binary_data()
     # Special case of categorical/biogeographic data with only two states
  # run_STRAPP_test_for_multinominal_data()
     # Case of categorical/biogeographic data with more than two states
     # Include biogeographic data as the tests are the same
     # Can define the pairwise tests to run/plot

# Make a wrapped-up to be run over time-series
  # run_STRAPP_test_over_time()

# Make a unique help documentation for all functions?

### Outputs
# focal_time, P-value, Qalpha% Δ stat
# Info needed to plot histogram
  # Focal-time
  # Melted STRAPP data for this focal_time
     # Time-step, BAMM ID, Δ stat distribution
  # Will be used to plot a single histogram for a given focal_time
    # Will compute summary stats from melted df
    # Include option of a PDF output
  # Will be merged in a melted data.frame encompassing all focal_time to produce sets of histogram
    # Include option of a merged PDF output with one histogram per page


# ## Load the BAMM_object summarizing 1000 posterior samples of BAMM
# data(Ponerinae_BAMM_object, package = "deepSTRAPP")
#
# ## Set focal-time to 10 My
# focal_time = 10
#
# ## Update the BAMM object (May take several minutes to run)
# Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
#   BAMM_object = Ponerinae_BAMM_object,
#   focal_time = focal_time,
#   update_rates = TRUE, update_regimes = TRUE,
#   update_tree = TRUE, update_plot = TRUE,
#   update_all_elements = TRUE,
#   keep_tip_labels = TRUE,
#   verbose = TRUE)

Ponerinae_BAMM_object_10My <- readRDS(file = "./src/Ponerinae_BAMM_object_10My.rds")

Ponerinae_data <- readRDS("../Ponerinae_Historical_Biogeography/input_data/Traits_data/Ponerinae_data.rds")
Ponerinae_tree <- readRDS("../Ponerinae_Historical_Biogeography/outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")
identical(Ponerinae_tree$tip.label, Ponerinae_data$Taxa)

attr(x = Ponerinae_BAMM_object_10My, which = "order") # Cladewise
attr(x = Ponerinae_tree, which = "order") # Postorder

# Reorder edge topology as in "cladewise order"
Ponerinae_tree <- ape::reorder.phylo(Ponerinae_tree, order = "cladewise")
identical(Ponerinae_tree$tip.label, Ponerinae_data$Taxa) # Does not affect tip order

Ponerinae_data_ln_HW <- setNames(object = Ponerinae_data$sim_ln_HW,
                                 nm = Ponerinae_data$Taxa)

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

Ponerinae_trait_data_10My$trait_data

identical(names(Ponerinae_BAMM_object_10My$tipStates[[1]]), names(Ponerinae_trait_data_10My$trait_data))

# Plot node labels on initial stochastic map with cut-off
plot(Ponerinae_contMap, lwd = 2)
ape::nodelabels()
abline(v = max(phytools::nodeHeights(Ponerinae_contMap$tree)[,2]) - focal_time,
       col = "red", lty = 2, lwd = 2)

# Plot initial node labels on cut stochastic map
plot(Ponerinae_trait_data_10My$contMap)
nodelabels(text = Ponerinae_trait_data_10My$contMap$tree$initial_nodes_ID)

BAMM_object <- Ponerinae_BAMM_object_10My
trait_data_list <- Ponerinae_trait_data_10My

## Example for continuous data

identical(names(BAMM_data$tipStates[[1]]), names(trait_data))

run_STRAPP_test_for_continuous_data(BAMM_data = BAMM_data, trait_data = trait_data_continuous, two_tailed = FALSE, one_tailed_hypothesis = "positive", return_perm_data = FALSE)


## Example for binary data

trait_data_continuous <- trait_data
trait_data_binary <- trait_data_continuous
trait_data_binary[trait_data_continuous < 0] <- "state_A"
trait_data_binary[trait_data_continuous >= 0] <- "state_B"

table(trait_data_binary)

run_STRAPP_test_for_binary_data(BAMM_data = BAMM_data, trait_data = trait_data_binary, two_tailed = TRUE)
run_STRAPP_test_for_binary_data(BAMM_data = BAMM_data, trait_data = trait_data_binary, two_tailed = FALSE, one_tailed_hypothesis = c("state_A > state_B"))
run_STRAPP_test_for_binary_data(BAMM_data = BAMM_data, trait_data = trait_data_binary, two_tailed = FALSE, one_tailed_hypothesis = c("state_B > state_A"))


## Example for multinominal data

# trait_data_continuous <- trait_data
trait_data_multinominal <- trait_data_continuous
trait_data_multinominal[trait_data_continuous < 0] <- "state_B"
trait_data_multinominal[trait_data_continuous < -1] <- "state_A"
trait_data_multinominal[trait_data_continuous >= 0] <- "state_C"

hist(trait_data_continuous)

table(trait_data_multinominal)

# trait_data_continuous <- trait_data
trait_data_multinominal_reverse <- trait_data_continuous
trait_data_multinominal_reverse[trait_data_continuous < 0] <- "state_B"
trait_data_multinominal_reverse[trait_data_continuous < -1] <- "state_C"
trait_data_multinominal_reverse[trait_data_continuous >= 0] <- "state_A"


run_STRAPP_test_for_multinominal_data(BAMM_data = BAMM_data, trait_data = trait_data_multinominal, alpha = 0.23, return_perm_data = TRUE)
test_output <- run_STRAPP_test_for_multinominal_data(BAMM_data = BAMM_data, trait_data = trait_data_multinominal, posthoc_pairwise_tests = TRUE, two_tailed = FALSE, alpha = 0.23)
test_output <- run_STRAPP_test_for_multinominal_data(BAMM_data = BAMM_data, trait_data = trait_data_multinominal, posthoc_pairwise_tests = TRUE, two_tailed = TRUE, return_perm_data = TRUE)

str(test_output$posthoc_pairwise_tests$perm_data_array, max.level = 2)
test_output$posthoc_pairwise_tests$summary_df

##



run_STRAPP_test <- function (BAMM_object, trait_data_list,
                             nb_permutations = NULL, # If NULL use all posterior samples by default
                             rate_type = "net_diversification", # Can use speciation, extinction and net diversification
                             two_tailed = TRUE, # Only for posthoc tests if applied on multinominal data.
                             alpha = 0.05, # Provide significant level use for estimates = significant threshold for the test statistics
                             replace = FALSE, # To allow or not multiple uses of posterior samples in random permutations
                             one_tailed_hypothesis = NULL, # For one-tailed tests. # Only for continuous and binary data. To select the hypothesis (negative/positive correlation or which state with greater rates). Inherited of traitorder from BAMMtools::traitDependentBAMM
                             posthoc_pairwise_tests = FALSE, # Only for categorical and biogeographic data with more than two states. If TRUE, run all possible pairs. Can also provide a list of pairs to run.
                             p.adjust_method = "none", # To adjust p-values in posthoc_pairwise_tests depending on the number of unique pairs of states
                             return_perm_data = FALSE, # return.full in initial function # To provide the data for this focal_time needed to plot an histogram
                             nthreads = 1, # Number of threads to use for parallelization of the function. The R package parallel must be loaded for nthreads > 1.
                             print_hypothesis = TRUE) # For printing what is tested
{

  ### Check input validity

  # BAMM_object must be a 'bammdata' object
  # Number of posterior sample data must be equal between $tipStates, $tipLambda and $tipMu

  # trait_data_list must be a list with $trait_data and $data_type
  # $trait_data must be a named vector (can be numerical or character string)
  # $data_type can only be "continuous", categorical" or "biogeographic"
  # $trait_data type must match with $data_type
     # Numerical if $data_type = "continuous"
     # Character string if $data_type = "categorical" or "biogeographic"

  # Length of $trait_data should match length of $tipStates, $tipLambda and $tipMu (for each posterior sample)

  # Names of $trait_data should match names in $tipStates, $tipLambda and $tipMu (for each posterior sample)
    # If need to be reordered, do it, but send a warning explaing differences may be due the fact the initial phylogenies used to model trait evolution and diversification dynamics may have not been ordered in the same fashion
    # Do NOT use $tip.label to check that as tip.label main contains labels for fossils older than focal_time

  # If nb_permutations is higher than number of posterior samples (length of $tipStates, $tipLambda and $tipMu) AND replace = FALSE,
  # Send an error to say that replace should be set to TRUE to allow multiple samplings of posterior in order to reach the requested number of permutations

  # rate_type must be either "speciation", "extinction" or "net diversification"

  # If two_tailed = FALSE, one_tailed_hypothesis must not be FALSE and contain credible hypotheses according to the type of data and state names
    # Provide error-specific warnings in case of mismatch (see warnings in traitDependentBAMM)

  ## See other validation checks from BAMMtools::traitDependentBAMM
  if (nthreads > 1) {
    if (!"package:parallel" %in% search()) {
      stop("Please load package 'parallel' for using the multi-thread option\n")
    }
  }


  ## Extract BAMM rates and regimes data
  BAMM_data <- list(tipStates = BAMM_object$tipStates, tipLambda = BAMM_object$tipLambda, tipMu = BAMM_object$tipMu)

  ## Extract trait data
  trait_data <- trait_data_list$trait_data

  ## Extract type of data
  data_type <- trait_data_list$data_type

  # If data_type is "categorical" or "biogeographic", reclassify according to the number of states
  if (data_type %in% c("categorical", "biogeographic"))
  {
    nb_levels <- nlevels(as.factor(trait_data))
    if (nb_levels == 2) # Case with two states
    {
      data_type <- "binary"
    } else { # Case with more than two states
      data_type <- "multinominal"
    }
  }

  ## Run the appropriate internal function depending on the type of data

  switch(EXPR = data_type,
    continuous =   { # Case for continuous data
                     # Stat test = Spearman's rank Rho test
                     STRAPP_test_output <- run_STRAPP_test_for_continuous_data(
                       BAMM_data = BAMM_data,
                       trait_data = trait_data,
                       nb_permutations = nb_permutations,
                       rate_type = rate_type,
                       two_tailed = two_tailed,
                       alpha = alpha,
                       replace = replace,
                       one_tailed_hypothesis = one_tailed_hypothesis,
                       return_perm_data = return_perm_data,
                       nthreads = nthreads,
                       print_hypothesis = print_hypothesis)
                   },
    binary =       { # Case for binary data (Special case of categorical/biogeographic data with only two states)
                     # Stat test = Mann-Whitney U test
                     STRAPP_test_output <- run_STRAPP_test_for_binary_data(
                       BAMM_data = BAMM_data,
                       trait_data = trait_data,
                       nb_permutations = nb_permutations,
                       rate_type = rate_type,
                       two_tailed = two_tailed,
                       alpha = alpha,
                       replace = replace,
                       one_tailed_hypothesis = one_tailed_hypothesis,
                       return_perm_data = return_perm_data,
                       nthreads = nthreads,
                       print_hypothesis = print_hypothesis)
                   },
    multinominal = { # Case for multinominal data (Case of categorical/biogeographic data with more than two states)
                     # Stat test = Kruskal-Wallis H test
                     # Can define the post-hoc pairwise tests to run
                     STRAPP_test_output <- run_STRAPP_test_for_multinominal_data(
                       BAMM_data = BAMM_data,
                       trait_data = trait_data,
                       nb_permutations = nb_permutations,
                       rate_type = rate_type,
                       alpha = alpha,
                       replace = replace,
                       posthoc_pairwise_tests = posthoc_pairwise_tests, # See if I implement that for pairwise posthoc tests. Need to provide list of pairs with hypotheses
                       two_tailed = two_tailed,
                       p.adjust_method = p.adjust_method,
                       return_perm_data = return_perm_data,
                       nthreads = nthreads,
                       print_hypothesis = print_hypothesis)
                   }
  )

  ## Include focal_time in the output
  STRAPP_test_output$focal_time <- BAMM_object$focal_time

  ## Export the STRAPP test output
  return(STRAPP_test_output)
}



### Needs doc. See if making a unique doc for the whole family of functions

run_STRAPP_test_for_continuous_data <- function (
    BAMM_data, trait_data,
    nb_permutations = NULL,
    rate_type = "net_diversification",
    two_tailed = TRUE,
    alpha = 0.05,
    replace = FALSE,
    one_tailed_hypothesis = NULL,
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{
  ### Check input validity

  # Apply the same as in the master function in case this function is used independently
  # Also check that the data_type is the good one

  ## See other validation checks from BAMMtools::traitDependentBAMM
  if (nthreads > 1)
  {
    if (!"package:parallel" %in% search())
    {
      stop("Please load package 'parallel' for using the multi-thread option\n")
    }
  }

  ## Extract rates data
  if (rate_type == "speciation")
  {
    rates_data <- BAMM_data$tipLambda
  }
  else if (rate_type == "extinction")
  {
    rates_data <- BAMM_data$tipMu
  }
  else if (rate_type == "net_diversification")
  {
    rates_data <- lapply(X = 1:length(BAMM_data$tipLambda), FUN = function (i) { BAMM_data$tipLambda[[i]] - BAMM_data$tipMu[[i]] })
  }

  ## Extract regime data
  regimes_data <- BAMM_data$tipStates

  ## Set number of permutations
  if (is.null(nb_permutations))
  {
    # If NULL, set to the number of posterior samples
    nb_permutations <- length(rates_data)
  }

  # Randomly sample posteriors to use for each permutation
  posterior_samples_random_ID <- sample(x = 1:length(rates_data), size = nb_permutations, replace = replace)

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:length(posterior_samples_random_ID))
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[posterior_samples_random_ID[l]]],
                                                           regimes = regimes_data[[posterior_samples_random_ID[l]]], stringsAsFactors = FALSE)
  }

  # Permute tip rates on tips using blocks defined by regime membership
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run permutation in parallel
    posterior_samples_permuted_rates_data <- parallel::parLapply(
      cl = cl, X = posterior_samples_random_rates_data, fun = block_permute_rates_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    posterior_samples_permuted_rates_data <- lapply(posterior_samples_random_rates_data, block_permute_rates_data)
  }

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })

  ## Print what is tested for One-tailed test
  if (print_hypothesis)
  {
    if (!two_tailed)
    {
      cat(paste0("Selected two-tailed Spearman's rank correlation test:\n\n",
                 "Null hypothesis: no correlation between trait data and diversification rates.\n\n",
                 "Alternative hypothesis: negative or positive correlation between trait data diversification rates.\n\n",
                 "'Estimate' stats is the ",alpha*100,"% quantile of abolute rho differences between observed and permuted data.\n",
                 "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))

      if (one_tailed_hypothesis == "positive")
      {
        cat(paste0("Selected one-tailed positive Spearman's rank correlation test:\n\n",
                   "Null hypothesis: negative or no correlation between trait data and diversification rates.\n\n",
                   "Alternative hypothesis: positive correlation between trait data diversification rates.\n",
                   "Low trait values associated with low diversification rates, and conversely.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of rho differences between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      } else {
        cat(paste0("Selected one-tailed negative Spearman's rank correlation test:\n\n",
                   "Null hypothesis: positive or no correlation between trait data and diversification rates.\n\n",
                   "Alternative hypothesis: negative correlation between trait data diversification rates.\n",
                   "Low trait values associated with high diversification rates, and conversely.\n\n",
                   "'Estimate' stats is the ",(1-alpha)*100,"% quantile of rho differences between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is lower than zero / p-value lower than ",alpha,".\n\n"))
      }
    }
  }

  ## Wrapped-up function to extract rho stats from Spearman's correlation test
  spearman_test <- function(rates, trait_data)
  {
    if (sd(rates, na.rm = TRUE) == 0)
    { # Case with no variance in rates. Rho = 0.
      return(0)
    } else { # Default case
      test_output <- stats::cor.test(rates, trait_data, method = "spearman", exact = FALSE)
      return(test_output$estimate)
    }
  }

  ## Run correlation test on each permutation. For observed data and permuted data.
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run test on observed data
    rho_obs <- parallel::parLapply(cl = cl,
                               X = posterior_samples_obs_rates_data,
                               fun = spearman_test,
                               trait_data = trait_data)
    # Run test on permuted data
    rho_perm <- parallel::parLapply(cl = cl,
                                    X = posterior_samples_permuted_rates_data,
                                    fun = spearman_test,
                                    trait_data = trait_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    rho_obs <- lapply(X = posterior_samples_obs_rates_data,
                      FUN = spearman_test,
                      trait_data = trait_data)
    rho_perm <- lapply(X = posterior_samples_permuted_rates_data,
                       FUN = spearman_test,
                       trait_data = trait_data)
  }

  ## Unlist outputs
  rho_obs <- unlist(rho_obs)
  rho_perm <- unlist(rho_perm)

  ## Compute p-value for two-tailed test
  if (two_tailed)
  {
    # Ho: correlation is equal in observed data than in permuted data
    # Ha: correlation is lower or higher in observed data than in permuted data
    # P-value = frequency of cases where observed stats is less extreme
    # (closer from null hypothesis) than the permuted stats
    p_value <- sum(abs(rho_obs) <= abs(rho_perm)) / length(rho_perm)
  } else {

    ## Compute p-value for one-tailed tests

    if (one_tailed_hypothesis == "positive")
    { # Test for positive correlation
      # Ho: correlation is lower in observed data than in permuted data
      # Ha: correlation is higher in observed data than in permuted data
      # P-value = frequency of cases where observed stats is lower than the permuted stats
      p_value <- sum(rho_obs <= rho_perm) / length(rho_perm)
    } else { # Test for negative correlation
      # Ho: correlation is higher in observed data than in permuted data
      # Ha: correlation is lower in observed data than in permuted data
      # P-value = frequency of cases where observed stats is higher than the permuted stats
      p_value <- sum(rho_obs >= rho_perm) / length(rho_perm)
    }
  }

  ## Save test stats
  if (two_tailed)
  {
    # If two-tailed test, need to compare the abs_delta_rho with alpha % quantile to see if higher than zero.
    STRAPP_test_output <- list(
      estimate = quantile(abs(rho_obs) - abs(rho_perm), p = alpha))
  } else {

    if (one_tailed_hypothesis == "positive")
    {
      # If one-tailed test for positive correlation, need to compare the delta_rho with alpha % quantile to see if higher than zero.
      STRAPP_test_output <- list(
        estimate = quantile(as.numeric(rho_obs) - as.numeric(rho_perm), p = alpha))
    } else {
      # If one-tailed test for negative correlation, need to compare the delta_rho with (1-alpha) % quantile to see if lower than zero.
      STRAPP_test_output <- list(
        estimate = quantile(as.numeric(rho_obs) - as.numeric(rho_perm), p = 1 - alpha))
    }
  }

  ## Save test summary results
  if (two_tailed)
  {
    # For two-tailed test, distribution based on difference in absolute correlations
    STRAPP_test_output$stats_median <- median(abs(rho_obs) - abs(rho_perm))
  } else {
    # For one-tailed test, distribution based on difference in correlations
    STRAPP_test_output$stats_median <- median(as.numeric(rho_obs) - as.numeric(rho_perm))
  }
  STRAPP_test_output$p_value <- p_value # P-value of the test
  STRAPP_test_output$method <- "Spearman" # Stats method
  STRAPP_test_output$two_tailed <- two_tailed # Type of test: two-tailed or not
  STRAPP_test_output$one_tailed_hypothesis <- one_tailed_hypothesis # Type of hypothesis if one-tailed test
  STRAPP_test_output$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    perm_data_df <- data.frame(posterior_samples_random_ID = posterior_samples_random_ID,
                               rho_obs = as.numeric(rho_obs),
                               rho_perm = as.numeric(rho_perm))
    if (two_tailed)
    { # For two-tailed test, distribution based on difference in absolute correlations
      perm_data_df$abs_delta_rho <- abs(rho_obs) - abs(rho_perm)
    } else { # For one-tailed test, distribution based on difference in correlations
      perm_data_df$delta_rho <- as.numeric(rho_obs) - as.numeric(rho_perm)
    }
    STRAPP_test_output$perm_data_df <- perm_data_df
  }

  ## Export output
  return(STRAPP_test_output)
}



### Needs doc. See if making a unique doc for the whole family of functions

run_STRAPP_test_for_binary_data <- function (
    BAMM_data, trait_data,
    nb_permutations = NULL,
    rate_type = "net_diversification",
    two_tailed = TRUE,
    alpha = 0.05,
    replace = FALSE,
    one_tailed_hypothesis = NULL,
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{

  ### Check input validity

  # Apply the same as in the master function in case this function is used independently
  # Also check that the data_type is the good one

  if (!two_tailed & is.null(one_tailed_hypothesis))
  {
    stop(paste0("You selected a one-tailed test ('two_tailed' = FALSE), but 'one_tailed_hypothesis' is not specified.\n",
                "You must specify the hypothesis by providing a character string vector with states ordered in increasing rates under the alternative hypothesis, separated by a greater-than such as c('A > B').\n"))
  }

  if (two_tailed & !is.null(one_tailed_hypothesis))
  {
    stop(paste0("You selected a two-tailed test ('two_tailed' = TRUE), but also specified a 'one_tailed_hypothesis': '",one_tailed_hypothesis,"'.\n",
                "If you want to test that hypothesis, please select a one-tailed test ('two_tailed' = FALSE).\n",
                "If you want to run a two-tailed test, replace the 'one_tailed_hypothesis' with 'NULL'.\n"))
  }

  ## See other validation checks from BAMMtools::traitDependentBAMM
  if (nthreads > 1)
  {
    if (!"package:parallel" %in% search())
    {
      stop("Please load package 'parallel' for using the multi-thread option\n")
    }
  }

  ## Extract rates data
  if (rate_type == "speciation")
  {
    rates_data <- BAMM_data$tipLambda
  }
  else if (rate_type == "extinction")
  {
    rates_data <- BAMM_data$tipMu
  }
  else if (rate_type == "net_diversification")
  {
    rates_data <- lapply(X = 1:length(BAMM_data$tipLambda), FUN = function (i) { BAMM_data$tipLambda[[i]] - BAMM_data$tipMu[[i]] })
  }

  ## Extract regime data
  regimes_data <- BAMM_data$tipStates

  ## Set number of permutations
  if (is.null(nb_permutations))
  {
    # If NULL, set to the number of posterior samples
    nb_permutations <- length(rates_data)
  }

  # Randomly sample posteriors to use for each permutation
  posterior_samples_random_ID <- sample(x = 1:length(rates_data), size = nb_permutations, replace = replace)

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:length(posterior_samples_random_ID))
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[posterior_samples_random_ID[l]]],
                                                           regimes = regimes_data[[posterior_samples_random_ID[l]]], stringsAsFactors = FALSE)
  }

  # Permute tip rates on tips using blocks defined by regime membership
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run permutation in parallel
    posterior_samples_permuted_rates_data <- parallel::parLapply(
      cl = cl, X = posterior_samples_random_rates_data, fun = block_permute_rates_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    posterior_samples_permuted_rates_data <- lapply(posterior_samples_random_rates_data, block_permute_rates_data)
  }

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })

  ## Prepare trait data

  # one_tailed_hypothesis <- c("state_A > state_B")

  obs_trait_states <- unique(trait_data)[order(unique(trait_data))]
  trait_states <- NA
  if (two_tailed)
  { # Case for two-tailed test

    if (print_hypothesis)
    {
      cat(paste0("Selected two-tailed Mann-Whitney-Wilcoxon rank-sum test:\n\n",
                 "Null hypothesis: taxa with trait '",
                 obs_trait_states[1], "' have equal ",rate_type," rates than those with trait '",
                 obs_trait_states[2], "'.\n",
                 "Alternative hypothesis: taxa with trait '",
                 obs_trait_states[1], "' have higher or lower ",rate_type," rates than those with trait '",
                 obs_trait_states[2], "'.\n\n",
                 "'Estimate' stats is the ",alpha*100,"% quantile of absolute differences in U-stats between observed and permuted data.\n",
                 "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
    }

  } else { # Case for one-tailed test

    # Parse one_tailed_hypothesis
    one_tailed_hypothesis_parsed <- gsub(pattern = " ", replacement = "", x = one_tailed_hypothesis)
    trait_states <- as.character(unlist(strsplit(x = one_tailed_hypothesis_parsed, split = ">")))

    if (!identical(trait_states, obs_trait_states) & !identical(rev(trait_states), obs_trait_states))
    {
      stop(paste0("States specified in the 'one_tailed_hypothesis' do not match with observed states in 'trait_data'.\n",
                  "States in 'one_tailed_hypothesis' = ", paste(trait_states, collapse = ", "), ". States in 'trait_data' = ", paste(obs_trait_states, collapse = ", "),".\n",
                  "Please use this format to provide the 'one_tailed_hypothesis': Two states separated by greater-than sign, with the state that is expected to have higher ",rate_type," rates in first.\n",
                  "Example: 'A > B'.\n"))
    } else {
      if (print_hypothesis)
      {
        cat(paste0("Selected one-tailed Mann-Whitney-Wilcoxon rank-sum test:\n\n",
                   "Null hypothesis: taxa with trait '",
                   trait_states[1], "' have lower or equal ",rate_type," rates than those with trait '",
                   trait_states[2], "'.\n",
                   "Alternative hypothesis: taxa with trait '",
                   trait_states[1], "' have higher ",rate_type," rates than those with trait '",
                   trait_states[2],"'.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of differences in U-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      }
    }
  }

  ## Wrapped-up function to extract U-stats from Mann-Whitney-Wilcoxon's rank-sum test
  mann_whitney_wilcoxon_test <- function(rates, trait_data, two_tailed, trait_states)
  {
    if (two_tailed)
    { # Case for two-tailed test
      test_output <- stats::wilcox.test(formula = rates ~ trait_data, exact = FALSE)
    } else { # Case for one-tailed test
      test_output <- stats::wilcox.test(x = rates[which(trait_data == trait_states[1])], # State with the higher ranked rates in Ha
                                        y = rates[which(trait_data == trait_states[2])], # State with the lower ranked rates in Ha
                                        exact = FALSE)
    }
    return(test_output$statistic)
  }

  ## Run MWW test on each permutation. For observed data and permuted data.
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run test on observed data
    U_obs <- parallel::parLapply(cl = cl,
                                 X = posterior_samples_obs_rates_data,
                                 fun = mann_whitney_wilcoxon_test,
                                 trait_data = trait_data,
                                 two_tailed = two_tailed,
                                 trait_states = trait_states)
    # Run test on permuted data
    U_perm <- parallel::parLapply(cl = cl,
                                  X = posterior_samples_permuted_rates_data,
                                  fun = mann_whitney_wilcoxon_test,
                                  trait_data = trait_data,
                                  two_tailed = two_tailed,
                                  trait_states = trait_states)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    U_obs <- lapply(X = posterior_samples_obs_rates_data,
                    FUN = mann_whitney_wilcoxon_test,
                    trait_data = trait_data,
                    two_tailed = two_tailed,
                    trait_states = trait_states)
    U_perm <- lapply(X = posterior_samples_permuted_rates_data,
                     FUN = mann_whitney_wilcoxon_test,
                     trait_data = trait_data,
                     two_tailed = two_tailed,
                     trait_states = trait_states)
  }

  ## Unlist outputs
  U_obs <- unlist(U_obs)
  U_perm <- unlist(U_perm)

  ## Center stats around location shift of the null hypothesis (mu)
  # Null hypothesis is that ranks of the values of the two groups are random
  # Compute location shift (mu) from state frequencies as average of the products of frequencies
  trait_data_counts <- table(trait_data)
  trait_data_counts <- trait_data_counts[!is.na(names(trait_data_counts))] # Remove NA
  stat_mu <- prod(trait_data_counts)/2
  # Center U-stats to get an estimate of how greater/lower (far away) than the null hypothesis (mu) are the calculated U-stats
  U_obs <- U_obs - stat_mu
  U_perm <- U_perm - stat_mu

  ## Compute p-value for two-tailed test
  if (two_tailed)
  {
    # Ho: rate differences in ranks between states are equal in observed data and permuted data
    # Ha: rate differences in ranks between states are lower or higher in observed data than in permuted data
    # P-value = frequency of cases where observed stats is less extreme
    # (closer from null hypothesis) than the permuted stats
    p_value <- sum(abs(U_obs) <= abs(U_perm)) / length(U_perm)

  } else {

    ## Compute p-value for one-tailed tests

    # Ho: rate differences in ranks between states are lower or equal in observed data than in permuted data
    # Ha: rate differences in ranks between states are higher in observed data than in permuted data
    # P-value = frequency of cases where observed stats is lower than the permuted stats
    p_value <- sum(U_obs <= U_perm) / length(U_perm)
  }

  ## Save test stats
  if (two_tailed)
  {
    # If two-tailed test, need to compare the abs_delta_U with alpha % quantile to see if higher than zero.
    STRAPP_test_output <- list(
      estimate = quantile(abs(U_obs) - abs(U_perm), p = alpha))
  } else {

    # If one-tailed test, need to compare the delta_U with alpha % quantile to see if higher than zero.
    STRAPP_test_output <- list(
      estimate = quantile(as.numeric(U_obs) - as.numeric(U_perm), p = alpha))
  }

  ## Save test summary results
  if (two_tailed)
  {
    # For two-tailed test, distribution based on difference in absolute U-stats
    STRAPP_test_output$stats_median <- median(abs(U_obs) - abs(U_perm))
  } else {
    # For one-tailed test, distribution based on difference in U-stats
    STRAPP_test_output$stats_median <- median(as.numeric(U_obs) - as.numeric(U_perm))
  }
  STRAPP_test_output$p_value <- p_value # P-value of the test
  STRAPP_test_output$method <- "Mann-Whitney-Wilcoxon" # Stats method
  STRAPP_test_output$two_tailed <- two_tailed # Type of test: two-tailed or not
  STRAPP_test_output$one_tailed_hypothesis <- one_tailed_hypothesis # Type of hypothesis if one-tailed test
  STRAPP_test_output$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification

  ## Save permutation results in a data.frame
  if (return_perm_data)
  {
    perm_data_df <- data.frame(posterior_samples_random_ID = posterior_samples_random_ID,
                               U_obs = as.numeric(U_obs),
                               U_perm = as.numeric(U_perm))
    if (two_tailed)
    { # For two-tailed test, distribution based on difference in absolute U-stats
      perm_data_df$abs_delta_U <- abs(U_obs) - abs(U_perm)
    } else { # For one-tailed test, distribution based on difference in U-stats
      perm_data_df$delta_U <- as.numeric(U_obs) - as.numeric(U_perm)
    }
    STRAPP_test_output$perm_data_df <- perm_data_df
  }

  ## Export output
  return(STRAPP_test_output)
}


### Needs doc. See if making a unique doc for the whole family of functions

run_STRAPP_test_for_multinominal_data <- function (
    BAMM_data, trait_data,
    nb_permutations = NULL,
    rate_type = "net_diversification",
    alpha = 0.05,
    replace = FALSE,
    posthoc_pairwise_tests = FALSE,
    two_tailed = TRUE,
    p.adjust_method = "none",
    return_perm_data = FALSE,
    nthreads = 1,
    print_hypothesis = TRUE)
{

  ### Check input validity

  # Apply the same as in the master function in case this function is used independently
  # Also check that the data_type is the good one

  ## See other validation checks from BAMMtools::traitDependentBAMM
  if (nthreads > 1)
  {
    if (!"package:parallel" %in% search())
    {
      stop("Please load package 'parallel' for using the multi-thread option\n")
    }
  }

  ## Extract rates data
  if (rate_type == "speciation")
  {
    rates_data <- BAMM_data$tipLambda
  }
  else if (rate_type == "extinction")
  {
    rates_data <- BAMM_data$tipMu
  }
  else if (rate_type == "net_diversification")
  {
    rates_data <- lapply(X = 1:length(BAMM_data$tipLambda), FUN = function (i) { BAMM_data$tipLambda[[i]] - BAMM_data$tipMu[[i]] })
  }

  ## Extract regime data
  regimes_data <- BAMM_data$tipStates

  ## Set number of permutations
  if (is.null(nb_permutations))
  {
    # If NULL, set to the number of posterior samples
    nb_permutations <- length(rates_data)
  }

  # Randomly sample posteriors to use for each permutation
  posterior_samples_random_ID <- sample(x = 1:length(rates_data), size = nb_permutations, replace = replace)

  # Build list of data.frame with rates and regimes ID data for each permutation
  posterior_samples_random_rates_data <- list()
  for (l in 1:length(posterior_samples_random_ID))
  {
    posterior_samples_random_rates_data[[l]] <- data.frame(rates = rates_data[[posterior_samples_random_ID[l]]],
                                                           regimes = regimes_data[[posterior_samples_random_ID[l]]], stringsAsFactors = FALSE)
  }

  # Permute tip rates on tips using blocks defined by regime membership
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run permutation in parallel
    posterior_samples_permuted_rates_data <- parallel::parLapply(
      cl = cl, X = posterior_samples_random_rates_data, fun = block_permute_rates_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    posterior_samples_permuted_rates_data <- lapply(posterior_samples_random_rates_data, block_permute_rates_data)
  }

  # Extract initial observed tip rates
  posterior_samples_obs_rates_data <- lapply(X = posterior_samples_random_rates_data, FUN = function (x) { x$rates })

  if (print_hypothesis)
  {
    cat(paste0("Selected Kruskal-Wallis's one-way ANOVA on ranks test:\n\n",
               "Null hypothesis: taxa have equal ",rate_type," rates independent from states.\n",
               "Alternative hypothesis: taxa have different ",rate_type," rates between states.\n\n",
               "'Estimate' stats is the ",alpha*100,"% quantile of differences in H-stats between observed and permuted data.\n",
               "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
  }

  ## Wrapped-up function to extract H-stats from Kruskal-Wallis's one-way ANOVA on ranks test
  kruskal_wallis_test <- function(rates, trait_data)
  {
    # Run the Kruskal-Wallis test
    test_output <- stats::kruskal.test(rates ~ trait_data)

    # If the test failed to provide a statistic because the value is reaching the ceiling for computation,
    # use the Khi-squared approximation by setting an extremely high p-value
    if (is.na(test_output$statistic))
    {
      H_approximation <- qchisq(p = 0.999, df = test_output$parameter)
      return(H_approximation)
    } else { # Otherwise, provide the computed H-stats
      return(test_output$statistic)
    }
  }

  ## Run Kruskal-Wallis test on each permutation. For observed data and permuted data.
  if (nthreads > 1) # In parallel
  {
    # Open cluster
    cl <- parallel::makePSOCKcluster(nthreads)
    # Run test on observed data
    H_obs <- parallel::parLapply(cl = cl,
                                 X = posterior_samples_obs_rates_data,
                                 fun = kruskal_wallis_test,
                                 trait_data = trait_data)
    # Run test on permuted data
    H_perm <- parallel::parLapply(cl = cl,
                                  X = posterior_samples_permuted_rates_data,
                                  fun = kruskal_wallis_test,
                                  trait_data = trait_data)
    # Close cluster
    parallel::stopCluster(cl)
  } else { # In series
    H_obs <- lapply(X = posterior_samples_obs_rates_data,
                    FUN = kruskal_wallis_test,
                    trait_data = trait_data)
    H_perm <- lapply(X = posterior_samples_permuted_rates_data,
                     FUN = kruskal_wallis_test,
                     trait_data = trait_data)
  }

  ## Unlist outputs
  H_obs <- unlist(H_obs)
  H_perm <- unlist(H_perm)

  ## Compute p-value

  # Ho: rate differences in ranks between states are equal in observed data and permuted data
  # Ha: rate differences in ranks between states are more extreme in observed data than in permuted data
  # Because the H-stat follows a Khi-squared distribution, thus is strictly positive and increases when
  # ranks between states are biased indifferently towards lower or higher ranks,
  # the comparison is made on H-stat differences.
  # P-value = frequency of cases where observed stats is higher than the permuted stats
  p_value <- sum(H_obs <= H_perm) / length(H_perm)

  ## Save test stats

  # Need to compare the delta_H with alpha % quantile to see if higher than zero.
  STRAPP_test_output <- list(
    estimate = quantile(as.numeric(H_obs) - as.numeric(H_perm), p = alpha))

  ## Save test summary results
  STRAPP_test_output$stats_median <- median(as.numeric(H_obs) - as.numeric(H_perm))
  STRAPP_test_output$p_value <- p_value # P-value of the test
  STRAPP_test_output$method <- "Kruskal-Wallis" # Stats method
  STRAPP_test_output$rate_type <- rate_type # Type of rates: speciation, extinction, or net diversification

  ## Save permutation results
  if (return_perm_data)
  {
    perm_data_df <- data.frame(posterior_samples_random_ID = posterior_samples_random_ID,
                               H_obs = as.numeric(H_obs),
                               H_perm = as.numeric(H_perm),
                               delta_U = as.numeric(H_obs) - as.numeric(H_perm))
    STRAPP_test_output$perm_data_df <- perm_data_df
  }

  ## Deal with post hoc pairwise tests
  if (posthoc_pairwise_tests)
  {
    ## Initiate output elements for post hoc pairwise tests
    STRAPP_test_output$posthoc_pairwise_tests <- list()

    ## Print hypothesis if requested
    if (print_hypothesis)
    {
      if (two_tailed)
      {
        cat(paste0("# --------- Post hoc pairwise tests --------- #\n\n",
                   "Selected two-tailed Dunn's post hoc pairwise rank-sum test:\n",
                   "Tests will be ran across all possible unique pairs of states.\n\n",
                   "Null hypothesis: taxa have equal ",rate_type," rates independent from states.\n",
                   "Alternative hypothesis: taxa have different ",rate_type," rates between states.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of absolute differences in Z-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      } else {
        cat(paste0("# --------- Post hoc pairwise tests --------- #\n\n",
                   "Selected one-tailed Dunn's post hoc pairwise rank-sum test:\n",
                   "Tests will be ran across all possible asymmetric pairs of states.\n\n",
                   "Null hypothesis: taxa with trait in the first state have lower or equal ",rate_type," rates than taxa with the second state in 'pairs'.\n",
                   "Alternative hypothesis: taxa with trait the first state have higher ",rate_type," rates than taxa with the second state in 'pairs'.\n\n",
                   "'Estimate' stats is the ",alpha*100,"% quantile of differences in Z-stats between observed and permuted data.\n",
                   "Null hypothesis is rejected if 'estimate' is higher than zero / p-value lower than ",alpha,".\n\n"))
      }
    }

    ## Wrapped-up function to extract Z-stats from Dunn's post hoc pairwise rank-sum tests
    dunn_test <- function(rates, trait_data, two_tailed)
    {
      # Run the Dunn test for all possible unique pairs of states
      invisible(capture.output(test_output <- dunn.test::dunn.test(x = rates, g = trait_data)))

      # Reformat test output
      if (two_tailed) # For two-tailed tests
      {
        test_output_df <- data.frame(pairs = gsub(pattern = " - ", replacement = " != ", x = test_output$comparisons),
                                     Z_stats = test_output$Z)
      } else { # For one-tailed tests
        test_output_df <- data.frame(pairs = c(gsub(pattern = " - ", replacement = " > ", x = test_output$comparisons), gsub(pattern = " - ", replacement = " < ", x = test_output$comparisons)),
                                     Z_stats = c(test_output$Z, -test_output$Z))
      }

      # If the test failed to provide a statistic because the value is reaching the ceiling for computation,
      # use the normal distribution and set an extremely high p-value to approximate a value
      if (any(is.na(test_output_df$Z_stats)))
      {
        Z_approximation <- qnorm(p = 0.999)
        test_output_df$Z_stats[is.na(test_output_df$Z_stats)] <- Z_approximation

      }

      # Export df
      return(test_output_df)
    }

    ## Run Dunn test on each permutation. For observed data and permuted data.
    if (nthreads > 1) # In parallel
    {
      # Open cluster
      cl <- parallel::makePSOCKcluster(nthreads)
      # Run test on observed data
      Dunn_obs <- parallel::parLapply(cl = cl,
                                      X = posterior_samples_obs_rates_data,
                                      fun = dunn_test,
                                      trait_data = trait_data,
                                      two_tailed = two_tailed)
      # Run test on permuted data
      Dunn_perm <- parallel::parLapply(cl = cl,
                                       X = posterior_samples_permuted_rates_data,
                                       fun = dunn_test,
                                       trait_data = trait_data,
                                       two_tailed = two_tailed)
      # Close cluster
      parallel::stopCluster(cl)
    } else { # In series
      Dunn_obs <- lapply(X = posterior_samples_obs_rates_data,
                         FUN = dunn_test,
                         trait_data = trait_data,
                         two_tailed = two_tailed)
      Dunn_perm <- lapply(X = posterior_samples_permuted_rates_data,
                           FUN = dunn_test,
                           trait_data = trait_data,
                           two_tailed = two_tailed)
    }

    ## Extract Z-scores from outputs

    # Get list of pairs
    pairs_list <- Dunn_obs[[1]]$pairs
    # Initiate list for Z-scores
    Z_obs <- list()
    Z_perm <- list()

    # Loop per pairs of states tested
    for (i in seq_along(pairs_list))
    {
      # i <- 1
      Z_obs[[i]] <- unlist(lapply(X = Dunn_obs, FUN = function (x) { x$Z_stats[i] }))
      Z_perm[[i]] <- unlist(lapply(X = Dunn_perm, FUN = function (x) { x$Z_stats[i] }))
    }
    names(Z_obs) <- names(Z_perm) <- pairs_list

    ## Compute p-value
    p_values <- list()

    # Loop per pairs of states tested
    for (i in seq_along(Z_obs))
    {
      # i <- 1

      if (two_tailed)
      {
        # For two-tailed, compare differences of absolute values of Z-scores. alpha % should be higher than 0.
        # Ho: rate differences in ranks between states are equal in observed data and permuted data
        # Ha: rate differences in ranks between states are more extreme in observed data than in permuted data
        # P-value = frequency of cases where observed stats is less extreme
        # (closer from null hypothesis) than the permuted stats

        p_values[[i]] <- sum(abs(Z_obs[[i]]) <= abs(Z_perm[[i]])) / length(Z_perm[[i]])
      } else {
        # For one-tailed, compare differences of values of Z-scores. alpha % should be higher than 0
        # Ho: rate differences in ranks between states are lower or equal in observed data than in permuted data
        # Ha: rate differences in ranks between states are higher in observed data than in permuted data
        # P-value = frequency of cases where observed stats is lower than the permuted stats
        p_values[[i]] <- sum(Z_obs[[i]] <= Z_perm[[i]]) / length(Z_perm[[i]])
      }
    }
    p_values <- unlist(p_values)

    ## Adjust p-values for multiple comparisons
    if (two_tailed)
    {
      # For two-tailed, use the number of unique pairs for adjustment
      p_values <- stats::p.adjust(p = p_values, method = p.adjust_method)
    } else {
      # For one-tailed, use the number of unique pairs by splitting p-values in two blocks
      # because reciprocal tests (A > B vs. B > A) are not independent tests.
      n_pairs <- length(p_values)
      p_values_forward <- p_values[1:(n_pairs/2)]
      p_values_backward <- p_values[((n_pairs/2)+1):n_pairs]
      p_values_forward <- stats::p.adjust(p = p_values_forward, method = p.adjust_method)
      p_values_backward <- stats::p.adjust(p = p_values_backward, method = p.adjust_method)
      p_values <- c(p_values_forward, p_values_backward)
    }

    ## Save test stats
    estimates <- list()
    stats_median <- list()

    # Loop per pairs of states tested
    for (i in seq_along(Z_obs))
    {
      # i <- 1

      if (two_tailed)
      {
        # If two-tailed test, need to compare the abs_delta_Z with alpha % quantile to see if higher than zero.
        estimates[[i]] <- quantile(abs(Z_obs[[i]]) - abs(Z_perm[[i]]), p = alpha)
        stats_median[[i]] <- median(abs(Z_obs[[i]]) - abs(Z_perm[[i]]))
      } else {
        # If one-tailed test, need to compare the delta_U with alpha % quantile to see if higher than zero.
        estimates[[i]] <- quantile(as.numeric(Z_obs[[i]]) - as.numeric(Z_perm[[i]]), p = alpha)
        stats_median[[i]] <- median(as.numeric(Z_obs[[i]]) - as.numeric(Z_perm[[i]]))
      }
    }

    ## Build summary df
    summary_df <- data.frame(pairs = pairs_list, estimates = unlist(estimates), stats_median = unlist(stats_median), p_values = p_values)

    ## Save test summary results
    STRAPP_test_output$posthoc_pairwise_tests$summary_df <- summary_df # Tests per pairs: estimates and p-values
    STRAPP_test_output$posthoc_pairwise_tests$method <- "Dunn" # Stats method
    STRAPP_test_output$posthoc_pairwise_tests$two_tailed <- two_tailed # Type of test: two-tailed or not

    ## Save permutation results in an array
    if (return_perm_data)
    {
      ## 3D-array to save permutation data
      # 1D = pairs
      # 2D = Posterior samples
      # 3D = Stats: Z_obs, Z_perm, Z_delta/Z_abs_delta

      if (two_tailed)
      { # For two-tailed test, distribution based on difference in absolute Z-scores
        perm_data_array <- array(data = NA,
                                 dim = c(pairs = length(pairs_list), posterior_samples = length(posterior_samples_random_ID), stats = 3),
                                 dimnames = list(pairs = pairs_list, posterior_samples = as.character(posterior_samples_random_ID), stats = c("Z_obs", "Z_perm", "abs_delta_Z")))
      } else { # For one-tailed test, distribution based on difference in Z-scores
        perm_data_array <- array(data = NA,
                                 dim = c(pairs = length(pairs_list), posterior_samples = length(posterior_samples_random_ID), stats = 3),
                                 dimnames = list(pairs = pairs_list, posterior_samples = as.character(posterior_samples_random_ID), stats = c("Z_obs", "Z_perm", "delta_Z")))
      }

      # Extract Z_obs across pairs
      Z_obs_df <- do.call(rbind.data.frame, Z_obs)
      perm_data_array[,,"Z_obs"] <- as.matrix(Z_obs_df)

      # Extract Z_perm across pairs
      Z_perm_df <- do.call(rbind.data.frame, Z_perm)
      perm_data_array[,,"Z_perm"] <- as.matrix(Z_perm_df)

      # Extract delta-stats distribution
      if (two_tailed)
      { # For two-tailed test, distribution based on difference in absolute Z-scores
        perm_data_array[,,"abs_delta_Z"] <- as.matrix(abs(Z_obs_df) - abs(Z_perm_df))
      } else { # For one-tailed test, distribution based on difference in Z-scores
        perm_data_array[,,"delta_Z"] <- as.matrix((Z_obs_df - Z_perm_df))
      }

      # Store array
      STRAPP_test_output$posthoc_pairwise_tests$perm_data_array <- perm_data_array
    }
  }

  ## Export output
  return(STRAPP_test_output)
}



### Helper function to permute rates on tips using regime-blocks
# Input = data.frame with $regimes and $rates

#' @title Permutes rates on tips using regime-blocks
#'
#' @description Permutes rates on tips using regime membership to define
#'   blocks used for permutation as in a STRAPP test.
#'   Each block of tips assigned the a regime has the same rates at a given point in time.
#'   During permutation, each block of tips is assigned a unique rate from any regime drawn randomly.
#'
#' @param rates_regimes_df Data.frame with `$regimes` (integer) and `$rates` (numerical)
#'   found in each tip (one row per tip).
#'
#' @return The function returns a numerical vector with permuted rates across tips.
#'
#' @author Maël Doré, Dan Rabosky, Huateng Huang
#'
#' @references Rabosky, D. L. and Huang, H., 2015. A Robust Semi-Parametric Test for Detecting Trait-Dependent Diversification. Systematic Biology 65: 181-193.
#'
#' @noRd
#'

block_permute_rates_data <- function (rates_regimes_df)
{
  # Extract regimes and rates data
  regimes <- rates_regimes_df$regimes
  rates <- rates_regimes_df$rates

  # Extract ID of unique regimes
  regimes_ID <- unique(regimes)

  # Extract rates for each regime
  # Rates are all equals in a given regime since they were extracted at the same focal_time
  rates_per_regimes <- numeric(length(regimes_ID))
  for (k in 1:length(regimes_ID))
  {
    rates_per_regimes[k] <- rates[regimes == regimes_ID[k]][1]
  }
  # Randomly rearrange regimes ID
  new_regimes_ID <- sample(regimes_ID, size = length(regimes_ID), replace = FALSE)

  # Replace rates with rates of newly assigned regimes = block-permutation based on regimes
  new_rates <- rep(0, length(regimes))
  for (k in 1:length(regimes_ID))
  {
    new_rates[which(regimes == regimes_ID[k])] <- rates_per_regimes[which(regimes_ID == new_regimes_ID[k])]
  }
  new_rates
}

