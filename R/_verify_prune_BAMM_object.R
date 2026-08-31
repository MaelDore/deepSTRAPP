###########################################################################################
##                                                                                       ##
##   Verification tests for deepSTRAPP::prune_BAMM_object()                              ##
##                                                                                       ##
##   Author: Maël Doré                                                                   ##
##                                                                                       ##
##   These tests check that a pruned BAMM_object is internally consistent, that it        ##
##   describes exactly the same diversification dynamics as the initial BAMM_object on     ##
##   the retained branches, and that it remains usable by downstream deepSTRAPP and        ##
##   BAMMtools functions.                                                                  ##
##                                                                                        ##
##   Run this script from the root of the deepSTRAPP package directory.                    ##
##   It is a development script: it is not run by R CMD check, and it can be turned        ##
##   into {testthat} tests once the function is stabilized.                                ##
##                                                                                        ##
###########################################################################################


### 0/ Setup ##############################################################################

## Load the package, including its internal helper functions
# 'compute_rate_at_time()' and 'getRecursiveSequence()' are internal, and are needed below
if (!exists(".deepSTRAPP_already_loaded"))
{
  devtools::load_all(".")
}

library(ape)
library(BAMMtools)

## Load the test dataset
data(whale_BAMM_object, package = "deepSTRAPP")

## Number of posterior samples used for the tests
# Set to NULL to use all 1000 posterior samples. A subset makes the tests much faster,
# but note that consecutive MCMC samples are strongly autocorrelated.
nb_test_samples <- 50

subset_posterior_samples <- function (BAMM_object, sample_indices)
{
  for (element in c("eventData", "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs"))
  {
    BAMM_object[[element]] <- BAMM_object[[element]][sample_indices]
  }
  BAMM_object$numberEvents <- BAMM_object$numberEvents[sample_indices]

  return(BAMM_object)
}

if (is.null(nb_test_samples))
{
  BAMM_object <- whale_BAMM_object
} else {
  BAMM_object <- subset_posterior_samples(whale_BAMM_object, 1:nb_test_samples)
}

initial_Ntip <- length(BAMM_object$tip.label)


### 1/ Test helpers #######################################################################

## Collect the outcome of every test, to print a summary at the end
test_results <- data.frame(section = character(0), test = character(0),
                           passed = logical(0), value = character(0),
                           stringsAsFactors = FALSE)

record_test <- function (section, test, passed, value = "")
{
  test_results <<- rbind(test_results,
                         data.frame(section = section, test = test,
                                    passed = isTRUE(passed), value = as.character(value),
                                    stringsAsFactors = FALSE))
  cat(sprintf("  [%s] %-42s %s\n", ifelse(isTRUE(passed), "OK  ", "FAIL"), test, value))

  invisible(isTRUE(passed))
}

## Extract a plain "phylo" object from any "bammdata" object
get_phylo <- function (x)
{
  phylo <- list(edge = x$edge, Nnode = x$Nnode, tip.label = x$tip.label, edge.length = x$edge.length)
  class(phylo) <- "phylo"
  attr(x = phylo, which = "order") <- "cladewise"

  return(phylo)
}

## Check the internal consistency of a "bammdata" object
# Returns the vector of individual checks, and records them all
check_bammdata <- function (x, section, label)
{
  Ntip <- length(x$tip.label)
  nb_samples <- length(x$eventData)
  phylo <- get_phylo(x)
  ok <- c()

  ## --- Phylogeny ---
  ok["binary"] <- ape::is.binary(phylo)
  ok["rooted"] <- ape::is.rooted(phylo)
  ok["Nnode_vs_Ntip"] <- (x$Nnode == Ntip - 1)
  ok["root_numbering"] <- identical(as.integer(setdiff(unique(x$edge[, 1]), x$edge[, 2])), as.integer(Ntip + 1L))
  ok["no_duplicated_tips"] <- (anyDuplicated(x$tip.label) == 0)

  ## --- Branch times ---
  node_heights <- ape::node.depth.edgelength(phylo)
  ok["begin_matches_tree"] <- max(abs(x$begin - node_heights[x$edge[, 1]])) < 1e-8
  ok["end_matches_tree"] <- max(abs(x$end - node_heights[x$edge[, 2]])) < 1e-8

  ## --- Tree traversal ---
  ok["downseq"] <- identical(x$downseq, getRecursiveSequence(phylo)$downseq)
  ok["lastvisit"] <- identical(x$lastvisit, getRecursiveSequence(phylo)$lastvisit)

  ## --- eventData ---
  ok["numberEvents"] <- all(x$numberEvents == sapply(x$eventData, nrow))
  ok["root_regime_is_first"] <- all(sapply(x$eventData, function (e) {
    (e$index[1] == 1) & (abs(e$time[1]) < 1e-8) & (e$node[1] == Ntip + 1) }))
  ok["index_is_sequence"] <- all(sapply(x$eventData, function (e) identical(e$index, seq_len(nrow(e)))))
  ok["events_sorted_by_time"] <- all(sapply(x$eventData, function (e) !is.unsorted(e$time)))
  ok["shifts_within_branches"] <- all(sapply(x$eventData, function (e) {
    shifts <- e[e$index != 1L, , drop = FALSE]
    if (nrow(shifts) == 0) return(TRUE)
    focal_edges <- match(shifts$node, x$edge[, 2])
    all(!is.na(focal_edges)) &
      all(shifts$time >= x$begin[focal_edges] - 1e-8) &
      all(shifts$time <= x$end[focal_edges] + 1e-8) }))

  ## --- Regime assignment ---
  ok["tipStates_valid"] <- all(mapply(function (ts, e) all(ts %in% e$index), x$tipStates, x$eventData))
  ok["eventVectors_valid"] <- all(mapply(function (ev, e) all(ev %in% e$index), x$eventVectors, x$eventData))
  ok["eventVectors_length"] <- all(sapply(x$eventVectors, length) == nrow(x$edge))
  ok["tipStates_length"] <- all(sapply(x$tipStates, length) == Ntip)
  ok["tip_rates_length"] <- all(sapply(x$tipLambda, length) == Ntip) & all(sapply(x$tipMu, length) == Ntip)
  # Tip vectors are named only when the initial BAMM_object named them:
  # BAMMtools::getEventData() leaves them unnamed, deepSTRAPP::update_rates_and_regimes_for_focal_time() names them.
  # Whichever convention applies, the names must be correct and consistent across elements.
  ok["tip_names"] <- all(sapply(x$tipStates, function (ts) is.null(names(ts)) | identical(names(ts), x$tip.label))) &
    all(sapply(x$tipLambda, function (tl) is.null(names(tl)) | identical(names(tl), x$tip.label))) &
    all(sapply(x$tipMu, function (tm) is.null(names(tm)) | identical(names(tm), x$tip.label))) &
    (length(unique(sapply(c(x$tipStates, x$tipLambda, x$tipMu), function (v) is.null(names(v))))) == 1)
  ok["tipStates_vs_eventVectors"] <- all(sapply(seq_len(nb_samples), function (i) {
    all(x$tipStates[[i]] == x$eventVectors[[i]][match(seq_len(Ntip), x$edge[, 2])]) }))

  ## --- Branch segments ---
  ok["nb_segments"] <- all(sapply(seq_len(nb_samples), function (i) {
    nrow(x$eventBranchSegs[[i]]) == nrow(x$edge) + x$numberEvents[i] - 1 }))
  ok["segments_tile_branches"] <- all(sapply(x$eventBranchSegs, function (m) {
    branch_totals <- tapply(m[, 3] - m[, 2], m[, 1], sum)
    branch_lengths <- x$edge.length[match(as.integer(names(branch_totals)), x$edge[, 2])]
    max(abs(as.numeric(branch_totals) - branch_lengths)) < 1e-8 }))
  ok["segments_monotonic"] <- all(sapply(x$eventBranchSegs, function (m) all(m[, 3] >= m[, 2] - 1e-10)))
  ok["segments_vs_eventVectors"] <- all(sapply(seq_len(nb_samples), function (i) {
    m <- x$eventBranchSegs[[i]]
    m <- m[order(m[, 1], m[, 2]), , drop = FALSE]
    last_segment <- !duplicated(m[, 1], fromLast = TRUE)
    all(m[last_segment, 4] == x$eventVectors[[i]][match(m[last_segment, 1], x$edge[, 2])]) }))

  cat("\n ", label, "\n", sep = "")
  for (test_name in names(ok))
  {
    record_test(section, paste0(label, " | ", test_name), ok[test_name])
  }

  invisible(ok)
}

## Recompute the tip rates from the macroevolutionary regimes alone, and compare them to
## the stored $tipLambda and $tipMu.
# This is the end-to-end check of the re-anchoring of the background/root regime: if the rate
# function of a regime were altered by the pruning, the recomputed tip rates would drift.
recompute_tip_rates <- function (x)
{
  Ntip <- length(x$tip.label)
  tip_ages <- x$end[match(seq_len(Ntip), x$edge[, 2])]
  lambda_errors <- mu_errors <- numeric(length(x$eventData))

  for (i in seq_along(x$eventData))
  {
    eventData_i <- x$eventData[[i]]
    tipStates_i <- x$tipStates[[i]]

    lambda_i <- mapply(FUN = function (state, tip_age) {
      compute_rate_at_time(eventData_i$lam1[state], eventData_i$lam2[state], tip_age - eventData_i$time[state]) },
      tipStates_i, tip_ages)
    mu_i <- mapply(FUN = function (state, tip_age) {
      compute_rate_at_time(eventData_i$mu1[state], eventData_i$mu2[state], tip_age - eventData_i$time[state]) },
      tipStates_i, tip_ages)

    lambda_errors[i] <- max(abs(lambda_i - x$tipLambda[[i]]))
    mu_errors[i] <- max(abs(mu_i - x$tipMu[[i]]))
  }

  return(c(lambda = max(lambda_errors), mu = max(mu_errors)))
}

## Compare the tip rates of a pruned object with those of the initial object
compare_tip_rates_with_initial <- function (pruned, initial)
{
  matching_tips <- match(pruned$tip.label, initial$tip.label)

  return(c(meanTipLambda = max(abs(pruned$meanTipLambda - initial$meanTipLambda[matching_tips])),
           meanTipMu = max(abs(pruned$meanTipMu - initial$meanTipMu[matching_tips])),
           tipLambda = max(sapply(seq_along(pruned$tipLambda), function (i) {
             max(abs(pruned$tipLambda[[i]] - initial$tipLambda[[i]][matching_tips])) })),
           tipMu = max(sapply(seq_along(pruned$tipMu), function (i) {
             max(abs(pruned$tipMu[[i]] - initial$tipMu[[i]][matching_tips])) }))))
}


### 2/ Baseline: the initial BAMM_object must pass the same checks ########################

cat("\n\n=========== 2/ Baseline on the initial BAMM_object ===========\n")

check_bammdata(BAMM_object, section = "2 - baseline", label = "initial BAMM_object")

initial_rate_error <- recompute_tip_rates(BAMM_object)
record_test("2 - baseline", "tip rates recomputed from eventData",
            max(initial_rate_error) < 1e-8, signif(max(initial_rate_error), 3))
# This residual is BAMM's own numerical precision, and is the reference against which the
# pruned objects are judged below
reference_rate_error <- max(initial_rate_error)


### 3/ Pruning nothing, and pruning a single tip ##########################################

cat("\n\n=========== 3/ No-op pruning and single-tip pruning ===========\n\n")

## Retaining every tip must short-circuit and return the object unchanged
pruned_all <- prune_BAMM_object(BAMM_object, tips_to_keep = BAMM_object$tip.label)
record_test("3 - minimal pruning", "retaining all tips returns the object unchanged",
            identical(pruned_all, BAMM_object))

## Removing a single tip exercises the whole machinery
pruned_one <- prune_BAMM_object(BAMM_object, tips_to_prune = BAMM_object$tip.label[1])
check_bammdata(pruned_one, section = "3 - minimal pruning", label = "one tip removed")

record_test("3 - minimal pruning", "one tip removed from the phylogeny",
            length(pruned_one$tip.label) == initial_Ntip - 1)
record_test("3 - minimal pruning", "removed tip recorded",
            identical(pruned_one$pruned_tip_labels, BAMM_object$tip.label[1]))
## Removing one tip suppresses exactly one internal node, hence exactly one merged branch
record_test("3 - minimal pruning", "exactly one branch results from merging",
            sum(pruned_one$pruning_edges_ID_df$nb_merged_edges > 1) == 2,
            paste0(sum(pruned_one$pruning_edges_ID_df$nb_merged_edges > 1) / 2, " merged branch(es)"))

## Branches that were not merged must carry exactly the same partition of regimes as before
untouched_edges <- pruned_one$pruning_edges_ID_df[pruned_one$pruning_edges_ID_df$nb_merged_edges == 1, ]
record_test("3 - minimal pruning", "regime partition preserved on untouched branches",
            all(sapply(seq_along(pruned_one$eventVectors), function (i) {
              new_regimes <- pruned_one$eventVectors[[i]][untouched_edges$new_edge_ID]
              initial_regimes <- BAMM_object$eventVectors[[i]][untouched_edges$initial_edge_ID]
              identical(as.integer(as.factor(new_regimes)), as.integer(as.factor(initial_regimes))) })))


### 4/ Arbitrary, non-monophyletic set of tips ############################################

cat("\n\n=========== 4/ Arbitrary, non-monophyletic set of tips ===========\n\n")

set.seed(1234)
tips_to_keep <- sample(BAMM_object$tip.label, size = 30, replace = FALSE)

pruned <- prune_BAMM_object(BAMM_object, tips_to_keep = tips_to_keep, verbose = TRUE)

check_bammdata(pruned, section = "4 - arbitrary tips", label = "pruned object")
check_bammdata(pruned$MAP_BAMM_object, section = "4 - arbitrary tips", label = "pruned MAP object")
check_bammdata(pruned$MSC_BAMM_object, section = "4 - arbitrary tips", label = "pruned MSC object")

cat("\n")

## The tips retained must be exactly the ones requested
record_test("4 - arbitrary tips", "retained tips are the ones requested",
            setequal(pruned$tip.label, tips_to_keep) & (length(pruned$tip.label) == 30))

## Providing the tips to prune instead of the tips to keep must give the same object
pruned_bis <- prune_BAMM_object(BAMM_object, tips_to_prune = setdiff(BAMM_object$tip.label, tips_to_keep))
record_test("4 - arbitrary tips", "'tips_to_keep' and 'tips_to_prune' agree",
            isTRUE(all.equal(pruned, pruned_bis)))

## Pruning must not touch the rates estimated at the retained tips
tip_rate_differences <- compare_tip_rates_with_initial(pruned, BAMM_object)
record_test("4 - arbitrary tips", "tip rates identical to the initial object",
            max(tip_rate_differences) == 0, signif(max(tip_rate_differences), 3))

## Merged branches must have a length equal to the sum of the initial branches they merge
merged_branch_lengths <- tapply(BAMM_object$edge.length[pruned$pruning_edges_ID_df$initial_edge_ID],
                                pruned$pruning_edges_ID_df$new_edge_ID, sum)
record_test("4 - arbitrary tips", "merged branch lengths add up",
            max(abs(as.numeric(merged_branch_lengths) - pruned$edge.length)) < 1e-10,
            signif(max(abs(as.numeric(merged_branch_lengths) - pruned$edge.length)), 3))

## The pruned topology must be the one produced by a plain {ape} pruning
ape_pruned_tree <- ape::keep.tip(get_phylo(BAMM_object), tips_to_keep)
record_test("4 - arbitrary tips", "topology matches ape::keep.tip()",
            ape::all.equal.phylo(ape_pruned_tree, get_phylo(pruned), use.edge.length = TRUE))

## The node conversion table must be consistent with the retained tips
record_test("4 - arbitrary tips", "node conversion table matches tip labels",
            identical(pruned$pruning_nodes_ID_df$initial_node_ID[pruned$pruning_nodes_ID_df$node_type == "tip"],
                      match(pruned$tip.label, BAMM_object$tip.label)))


### 5/ Subclade extraction with MRCA_node #################################################

cat("\n\n=========== 5/ Subclade extraction with MRCA_node ===========\n\n")

## Pick an internal node subtending a reasonably large subclade, that is not the root
initial_phylo <- get_phylo(BAMM_object)
candidate_nodes <- setdiff(unique(BAMM_object$edge[, 1]), initial_Ntip + 1L)
clade_sizes <- sapply(candidate_nodes, function (focal_node) {
  length(ape::extract.clade(initial_phylo, focal_node)$tip.label) })
MRCA_node <- candidate_nodes[which(clade_sizes >= 8)[1]]

cat("  Focal MRCA node:", MRCA_node, "subtending", clade_sizes[which(clade_sizes >= 8)[1]], "tips\n")

pruned_clade <- prune_BAMM_object(BAMM_object, MRCA_node = MRCA_node)

check_bammdata(pruned_clade, section = "5 - subclade", label = "pruned subclade")

cat("\n")

## The retained tips must be the descendants of the focal node
record_test("5 - subclade", "retained tips are the descendants of MRCA_node",
            setequal(pruned_clade$tip.label, ape::extract.clade(initial_phylo, MRCA_node)$tip.label))

## The root shift must equal the difference in root age
# Pruning never changes the distance of a retained node to the present
root_age_difference <- max(BAMM_object$end) - max(pruned_clade$end)
record_test("5 - subclade", "root shift equals the drop in root age",
            abs(pruned_clade$pruning_root_shift - root_age_difference) < 1e-8,
            paste0("shift = ", signif(pruned_clade$pruning_root_shift, 5)))
record_test("5 - subclade", "root shift is strictly positive",
            pruned_clade$pruning_root_shift > 0)

## Even with a shifted root, tip rates must be untouched
clade_rate_differences <- compare_tip_rates_with_initial(pruned_clade, BAMM_object)
record_test("5 - subclade", "tip rates identical to the initial object",
            max(clade_rate_differences) == 0, signif(max(clade_rate_differences), 3))

## Compare with BAMMtools::subtreeBAMM(), which is reliable for subclade extraction
subtree_reference <- BAMMtools::subtreeBAMM(BAMM_object, node = MRCA_node)
record_test("5 - subclade", "same tips as BAMMtools::subtreeBAMM()",
            setequal(pruned_clade$tip.label, subtree_reference$tip.label))
record_test("5 - subclade", "same topology as BAMMtools::subtreeBAMM()",
            ape::all.equal.phylo(get_phylo(pruned_clade), get_phylo(subtree_reference), use.edge.length = TRUE))
## Regime IDs are arbitrary labels, so compare the partition of tips into regimes.
# Note that BAMMtools::subtreeBAMM() returns UNNAMED $tipStates (its match() call drops the names),
# so tips must be aligned through $tip.label rather than through the names of the vectors.
reference_tip_order <- match(pruned_clade$tip.label, subtree_reference$tip.label)
record_test("5 - subclade", "same regime partition as BAMMtools::subtreeBAMM()",
            all(sapply(seq_along(pruned_clade$tipStates), function (i) {
              new_states <- as.integer(pruned_clade$tipStates[[i]])
              reference_states <- as.integer(subtree_reference$tipStates[[i]])[reference_tip_order]
              identical(as.integer(as.factor(new_states)), as.integer(as.factor(reference_states))) })))
record_test("5 - subclade", "same tip rates as BAMMtools::subtreeBAMM()",
            max(sapply(seq_along(pruned_clade$tipLambda), function (i) {
              max(abs(as.numeric(pruned_clade$tipLambda[[i]]) -
                        as.numeric(subtree_reference$tipLambda[[i]])[reference_tip_order])) })) < 1e-12)
## The re-anchored background regime must match the one computed by BAMMtools::subtreeBAMM()
record_test("5 - subclade", "same background regime as BAMMtools::subtreeBAMM()",
            max(sapply(seq_along(pruned_clade$eventData), function (i) {
              max(abs(pruned_clade$eventData[[i]]$lam1[1] - subtree_reference$eventData[[i]]$lam1[1]),
                  abs(pruned_clade$eventData[[i]]$mu1[1] - subtree_reference$eventData[[i]]$mu1[1])) })) < 1e-12)


### 6/ Rate functions are preserved on the retained branches ##############################

cat("\n\n=========== 6/ Rate functions preserved on the retained branches ===========\n\n")

## Recomputing the tip rates from the pruned $eventData alone must reproduce $tipLambda/$tipMu.
# This is the decisive test of the re-anchoring of the background/root regime: $lam1 and $mu1 of
# that regime do change (they now report the rates reached at the new root age), but the rate
# function they parametrize must be strictly identical on the retained branches.
for (test_case in list(list(label = "pruned, root unchanged", object = pruned),
                       list(label = "pruned subclade, root shifted", object = pruned_clade)))
{
  rate_error <- recompute_tip_rates(test_case$object)
  record_test("6 - rate functions", paste0(test_case$label, " | tip rates recomputed from eventData"),
              max(rate_error) <= max(reference_rate_error, 1e-8), signif(max(rate_error), 3))
}

## Report the re-anchoring of the background regime, for the record
cat("\n  Background regime of posterior sample 1, subclade case:\n")
cat("    initial: time =", signif(BAMM_object$eventData[[1]]$time[1], 5),
    "| lam1 =", signif(BAMM_object$eventData[[1]]$lam1[1], 6), "\n")
cat("    pruned : time =", signif(pruned_clade$eventData[[1]]$time[1], 5),
    "| lam1 =", signif(pruned_clade$eventData[[1]]$lam1[1], 6), "\n")
cat("    expected lam1 =", signif(compute_rate_at_time(BAMM_object$eventData[[1]]$lam1[1],
                                                       BAMM_object$eventData[[1]]$lam2[1],
                                                       pruned_clade$pruning_root_shift - BAMM_object$eventData[[1]]$time[1]), 6), "\n")
record_test("6 - rate functions", "background regime re-anchored as expected",
            abs(pruned_clade$eventData[[1]]$lam1[1] -
                  compute_rate_at_time(BAMM_object$eventData[[1]]$lam1[1], BAMM_object$eventData[[1]]$lam2[1],
                                       pruned_clade$pruning_root_shift - BAMM_object$eventData[[1]]$time[1])) < 1e-12)


### 7/ Marginal Shift Probability tree ####################################################

cat("\n\n=========== 7/ Marginal Shift Probability tree ===========\n\n")

## The stored $MSP_tree must match a recomputation on the pruned posterior samples
MSP_reference <- BAMMtools::marginalShiftProbsTree(pruned)
record_test("7 - MSP tree", "matches BAMMtools::marginalShiftProbsTree()",
            max(abs(pruned$MSP_tree$edge.length - MSP_reference$edge.length)) < 1e-12)

## The MSP of a merged branch must equal the joint probability of at least one shift
## anywhere along the initial branches it merges
merged_branch_IDs <- unique(pruned$pruning_edges_ID_df$new_edge_ID[pruned$pruning_edges_ID_df$nb_merged_edges > 1])
joint_probabilities <- sapply(merged_branch_IDs, function (focal_edge) {
  initial_edges <- pruned$pruning_edges_ID_df$initial_edge_ID[pruned$pruning_edges_ID_df$new_edge_ID == focal_edge]
  initial_nodes <- BAMM_object$edge[initial_edges, 2]
  mean(sapply(BAMM_object$eventData, function (e) any(initial_nodes %in% e$node))) })
record_test("7 - MSP tree", "merged branches hold the joint shift probability",
            max(abs(pruned$MSP_tree$edge.length[merged_branch_IDs] - joint_probabilities)) < 1e-12,
            paste0(length(merged_branch_IDs), " merged branches checked"))

## Marginal probabilities must remain valid probabilities
record_test("7 - MSP tree", "probabilities within [0, 1]",
            all(pruned$MSP_tree$edge.length >= 0) & all(pruned$MSP_tree$edge.length <= 1))
record_test("7 - MSP tree", "MSP tree topology matches the pruned tree",
            identical(pruned$MSP_tree$edge, pruned$edge) &
              identical(pruned$MSP_tree$tip.label, pruned$tip.label))


### 8/ MAP and MSC configurations #########################################################

cat("\n\n=========== 8/ MAP and MSC configurations ===========\n\n")

## By default, the posterior sample indices identified on the full phylogeny are kept
record_test("8 - MAP/MSC", "MAP indices unchanged by default",
            identical(pruned$MAP_indices, BAMM_object$MAP_indices))
record_test("8 - MAP/MSC", "MSC indices unchanged by default",
            identical(pruned$MSC_indices, BAMM_object$MSC_indices))
record_test("8 - MAP/MSC", "expectedNumberOfShifts carried over",
            identical(pruned$expectedNumberOfShifts, BAMM_object$expectedNumberOfShifts))

## The MAP/MSC sub-objects must describe the same phylogeny as the main pruned object
record_test("8 - MAP/MSC", "MAP object shares the pruned topology",
            identical(pruned$MAP_BAMM_object$edge, pruned$edge) &
              identical(pruned$MAP_BAMM_object$tip.label, pruned$tip.label))
record_test("8 - MAP/MSC", "MSC object shares the pruned topology",
            identical(pruned$MSC_BAMM_object$edge, pruned$edge) &
              identical(pruned$MSC_BAMM_object$tip.label, pruned$tip.label))

## Detecting the configurations again on the pruned phylogeny
pruned_recomputed <- prune_BAMM_object(BAMM_object, tips_to_keep = tips_to_keep,
                                       recompute_shift_configurations = TRUE, verbose = TRUE)

check_bammdata(pruned_recomputed$MAP_BAMM_object, section = "8 - MAP/MSC", label = "recomputed MAP object")
check_bammdata(pruned_recomputed$MSC_BAMM_object, section = "8 - MAP/MSC", label = "recomputed MSC object")

cat("\n")
cat("  MAP samples on the full phylogeny :", length(BAMM_object$MAP_indices), "\n")
cat("  MAP samples on the pruned phylogeny:", length(pruned_recomputed$MAP_indices), "\n")
record_test("8 - MAP/MSC", "recomputed MAP indices are valid sample indices",
            all(pruned_recomputed$MAP_indices %in% seq_along(pruned$eventData)))
record_test("8 - MAP/MSC", "recomputed MSC indices are valid sample indices",
            all(pruned_recomputed$MSC_indices %in% seq_along(pruned$eventData)))
## Only the shift configurations may differ between the two policies
record_test("8 - MAP/MSC", "posterior samples identical under both policies",
            isTRUE(all.equal(pruned_recomputed$eventData, pruned$eventData)))


### 9/ Errors and edge cases ##############################################################

cat("\n\n=========== 9/ Errors and edge cases ===========\n\n")

## Every call below must fail, with an informative message
expect_error <- function (label, expr)
{
  outcome <- tryCatch({ expr ; NULL }, error = function (e) e)
  record_test("9 - errors", label, !is.null(outcome))
  if (!is.null(outcome))
  {
    cat("        ", gsub("\n", "\n         ", conditionMessage(outcome)), "\n")
  }
}

expect_error("no tip selection provided",
             prune_BAMM_object(BAMM_object))
expect_error("two tip selections provided",
             prune_BAMM_object(BAMM_object, tips_to_keep = tips_to_keep, tips_to_prune = BAMM_object$tip.label[1]))
expect_error("tips provided as indices",
             prune_BAMM_object(BAMM_object, tips_to_keep = 1:5))
expect_error("unknown tip label",
             prune_BAMM_object(BAMM_object, tips_to_keep = c("not_a_tip", tips_to_keep[1:5])))
expect_error("a single tip retained",
             prune_BAMM_object(BAMM_object, tips_to_keep = tips_to_keep[1]))
expect_error("no tip retained",
             prune_BAMM_object(BAMM_object, tips_to_prune = BAMM_object$tip.label))
expect_error("MRCA_node is a tip",
             prune_BAMM_object(BAMM_object, MRCA_node = 1))
expect_error("MRCA_node holds several nodes",
             prune_BAMM_object(BAMM_object, MRCA_node = c(89, 90)))
expect_error("BAMM_object is not a bammdata object",
             prune_BAMM_object(get_phylo(BAMM_object), tips_to_keep = tips_to_keep))

## Retaining a few tips must warn about the loss of statistical power, but still work
small_pruning <- withCallingHandlers({ prune_BAMM_object(BAMM_object, tips_to_keep = tips_to_keep[1:5]) },
                                     warning = function (w) {
                                       record_test("9 - errors", "warning raised for very small phylogenies", TRUE)
                                       invokeRestart("muffleWarning") })
check_bammdata(small_pruning, section = "9 - errors", label = "five tips retained")

## Two sister tips: the smallest valid output, with the deepest possible root shift
sister_parents <- BAMM_object$edge[BAMM_object$edge[, 2] <= initial_Ntip, 1]
cherry_node <- as.integer(names(which(table(sister_parents) == 2))[1])
sister_tips <- BAMM_object$tip.label[BAMM_object$edge[BAMM_object$edge[, 1] == cherry_node, 2]]

cat("\n  Sister tips:", paste(sister_tips, collapse = " & "), "\n")
two_tips_pruning <- suppressWarnings(prune_BAMM_object(BAMM_object, tips_to_keep = sister_tips))
check_bammdata(two_tips_pruning, section = "9 - errors", label = "two sister tips")

cat("\n")
record_test("9 - errors", "two sister tips | single regime remains",
            all(two_tips_pruning$numberEvents == 1))
two_tips_rate_error <- recompute_tip_rates(two_tips_pruning)
record_test("9 - errors", "two sister tips | tip rates recomputed from eventData",
            max(two_tips_rate_error) <= max(reference_rate_error, 1e-8), signif(max(two_tips_rate_error), 3))


### 10/ Downstream compatibility ##########################################################

cat("\n\n=========== 10/ Downstream compatibility ===========\n\n")

## BAMMtools must accept the pruned object
record_test("10 - downstream", "BAMMtools::dtRates() runs",
            !inherits(tryCatch(BAMMtools::dtRates(pruned, tau = 0.01), error = function (e) e), "error"))
record_test("10 - downstream", "BAMMtools::getTipRates() runs",
            !inherits(tryCatch(BAMMtools::getTipRates(pruned), error = function (e) e), "error"))

## deepSTRAPP plotting must accept the pruned object, including the regime shift overlays
plot_outcome <- tryCatch({
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  plot_BAMM_rates(pruned, add_regime_shifts = TRUE, configuration_type = "MAP",
                  adjust_size_to_prob = TRUE, regimes_size = 3, bg = "black")
  plot_BAMM_rates(pruned, add_regime_shifts = TRUE, configuration_type = "MSC",
                  regimes_size = 3, bg = "black")
  plot_BAMM_rates(pruned_clade, add_regime_shifts = TRUE, configuration_type = "MAP",
                  regimes_size = 3, bg = "black")
  grDevices::dev.off()
  NULL }, error = function (e) { try(grDevices::dev.off(), silent = TRUE) ; e })
record_test("10 - downstream", "plot_BAMM_rates() runs on pruned objects", is.null(plot_outcome),
            ifelse(is.null(plot_outcome), "", conditionMessage(plot_outcome)))

## A pruned object must remain usable by the deepSTRAPP time-slicing machinery.
## Note that the maximum valid 'focal_time' shrinks by 'pruning_root_shift'.
cut_outcome <- tryCatch({
  update_rates_and_regimes_for_focal_time(BAMM_object = pruned, focal_time = 10,
                                          update_all_elements = TRUE, verbose = FALSE) },
  error = function (e) e)
record_test("10 - downstream", "update_rates_and_regimes_for_focal_time() runs on a pruned object",
            !inherits(cut_outcome, "error"),
            ifelse(inherits(cut_outcome, "error"), conditionMessage(cut_outcome), ""))


### 11/ Summary ###########################################################################

cat("\n\n=========== 11/ Summary ===========\n\n")

summary_per_section <- aggregate(passed ~ section, data = test_results,
                                 FUN = function (x) { paste0(sum(x), "/", length(x)) })
names(summary_per_section) <- c("Section", "Passed")
print(summary_per_section, row.names = FALSE)

cat("\n  TOTAL:", sum(test_results$passed), "/", nrow(test_results), "tests passed\n\n")

if (any(!test_results$passed))
{
  print(test_results[!test_results$passed, c("section", "test", "value")], row.names = FALSE)
  stop("Some verification tests failed for prune_BAMM_object().")
} else {
  cat("  All verification tests passed for prune_BAMM_object().\n\n")
}
