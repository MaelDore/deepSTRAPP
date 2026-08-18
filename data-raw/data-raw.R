##### Scripts to prepare datasets made available in the package #####


##### Trait datasets #####

### 1/ Include motmot::mammals dataset within deepSTRAPP directly to not have to rely on motmot as dependency ####

library(motmot)

data("mammals")
force(mammals)

?mammals

# Export in deepSTRAPP
usethis::use_data(mammals, overwrite = TRUE)

# ### 2/ Simulate Ponerinae HW data ####
#
# library(phytools)
#
# # Load phylogeny
# Ponerinae_tree <- readRDS(file = "../Ponerinae_Historical_Biogeography/outputs/Grafting_missing_taxa/Ponerinae_MCC_phylogeny_1534t_short_names.rds")
#
# # Load length data for Ponerinae (continuous trait)
# Ponerinae_data <- readRDS(file = "../Ponerinae_Historical_Biogeography/input_data/Traits_data/Trait_database.rds")
# Ponerinae_data <- Ponerinae_data |>
#   dplyr::select(Current_name, HW) |>
#   dplyr::distinct() |>
#   dplyr::mutate(ln_HW = log(HW))
#
# # Extract known data
# Ponerinae_data_temp <-as.data.frame(Ponerinae_tree$tip.label)
# names(Ponerinae_data_temp) <- "Taxa"
# Ponerinae_data <- dplyr::left_join(x = Ponerinae_data_temp, y = Ponerinae_data, by = dplyr::join_by("Taxa" == "Current_name"))
# Ponerinae_data <- Ponerinae_data |>
#   dplyr::distinct(Taxa, .keep_all = T)
#
# ## Simulate fake data per Genera using known data
# Ponerinae_data$sim_ln_HW <- Ponerinae_data$ln_HW
#
# # Detect genera
# genera_list <- unique(stringr::str_split(string = Ponerinae_tree$tip.label, pattern = "_", simplify = T)[, 1])
#
# # Loop per genus
# for (i in seq_along(genera_list))
# {
#   # i <- 1
#
#   genus_i <- genera_list[i]
#   # Get list of taxa
#   taxa_list_i <- which(stringr::str_detect(string = Ponerinae_tree$tip.label, pattern = genus_i))
#
#   # Remove Neoponera_bucki
#   taxa_list_i <- setdiff(x = taxa_list_i, y = which(Ponerinae_tree$tip.label %in% "Neoponera_bucki"))
#
#   if (length(taxa_list_i) > 1)
#   {
#     # Get MRCA
#     MRCA_i <- ape::getMRCA(phy = Ponerinae_tree, tip = taxa_list_i)
#     # Get all descending taxa
#     descendents_list_i <- phytools::getDescendants(tree = Ponerinae_tree, node = MRCA_i)
#     descendents_list_i <- descendents_list_i[descendents_list_i <= length(Ponerinae_tree$tip.label)]
#     descendents_list_names_i <- Ponerinae_tree$tip.label[descendents_list_i]
#
#     # Prune tree to MRCA
#     genus_tree_i <- ape::keep.tip(phy = Ponerinae_tree, tip = descendents_list_i)
#
#     # Get mean and var of available data
#     genus_data_i <- Ponerinae_data[Ponerinae_data$Taxa %in% descendents_list_names_i, ]
#     mean_i <- mean(genus_data_i$ln_HW, na.rm = T)
#     var_i <- var(genus_data_i$ln_HW, na.rm = T)
#     if (is.na(var_i)) { var_i <- 0.5 }
#
#     # Run BM to simulate data
#     simBM_i <- phytools::fastBM(tree = genus_tree_i, a = mean_i, sig2 = var_i)
#
#     # Replace missing data with simulated data
#     genus_data_i$sim_ln_HW[is.na(genus_data_i$ln_HW)] <- simBM_i[is.na(genus_data_i$ln_HW)]
#     Ponerinae_data$sim_ln_HW[match(x = genus_data_i$Taxa, table = Ponerinae_data$Taxa)] <- genus_data_i$sim_ln_HW
#   }
# }
#
# # Convert back to natural scale
# Ponerinae_data$sim_HW <- exp(Ponerinae_data$sim_ln_HW)
# table(is.na(Ponerinae_data$sim_HW))
#
# # Save simulated data
# saveRDS(object = Ponerinae_data, file = "../Ponerinae_Historical_Biogeography/input_data/Traits_data/Ponerinae_data.rds")


### 2/ Produce trait dataset for old_calib to include in Ponerinae_trait_tip_data ####

# # Load phylogeny with UCE data
# Ponerinae_phylogeny_789t_calibrated_short_names <- readRDS("../Ponerinae_Historical_Biogeography/input_data/Phylogenies/Ponerinae_phylogeny_789t_calibrated_short_names.rds")
# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
# Load trait df
data(Ponerinae_trait_tip_data, package = "deepSTRAPP")

## 2.1/ Produce fake_cont_tip_data ####

# Ponerinae_phylogeny_789t_calibrated_short_names$tip.label[Ponerinae_phylogeny_789t_calibrated_short_names$tip.label == "NewGenus_bucki"] <- "Neoponera_bucki"

# ## Extract trait data for reduced tree
# ln_HW_789t <- Ponerinae_trait_tip_data$ln_HW[match(x = Ponerinae_phylogeny_789t_calibrated_short_names$tip.label, table = Ponerinae_trait_tip_data$Taxa)]
# # Extract continuous trait data as a named vector
# ln_HW_789t <- setNames(object = ln_HW_789t,
#                     nm = Ponerinae_phylogeny_789t_calibrated_short_names$tip.label)
# names(ln_HW_789t)[is.na(ln_HW_789t)]
#
# ln_HW_789t["Myopias_mj_ohu_1"] <- log(0.750)
# ln_HW_789t["Myopias_my11"] <- log(0.850)
# ln_HW_789t["Hypoponera_cn01"] <- log(0.750)
# ln_HW_789t["Hypoponera_my13"] <- log(0.550)
# ln_HW_789t["Ectomomyrmex_th05"] <- log(1.200)
# ln_HW_789t["Ponera_swezeyi"] <- log(0.400)
# ln_HW_789t["Harpegnathos_my01"] <- log(2.200)
# ln_HW_789t["Thaumatomyrmex_br02"] <- log(0.950)
# ln_HW_789t["Neoponera_bra549385"] <- log(2.000)
#
# ##  Compute ACE on reduced tree
# Ponerinae_trait_789t <- prepare_trait_data(tip_data = ln_HW_789t,
#                                            trait_data_type = "continuous",
#                                            phylo = Ponerinae_phylogeny_789t_calibrated_short_names,
#                                            seed = 1234) # Set seed for reproducibility
#
# # Extract the Ancestral Character Estimates (ACE) = trait values at nodes
# Ponerinae_ACE_789t <- Ponerinae_trait_789t$ace
# head(Ponerinae_ACE_789t)
#
# ## Adjust ACE data to improve significance for the example
#
# pdf(file = "Ponerinae_phylogeny_789t.pdf", width = 30, height = 150)
# plot(Ponerinae_phylogeny_789t_calibrated_short_names)
# nodelabels()
# dev.off()
#
# # Simopelta
# Ponerinae_ACE_789t["1538"] <- 0.7
# # Raspone
# Ponerinae_ACE_789t["1514"] <- 0.9
# # Thaumatomyrmex
# Ponerinae_ACE_789t["1437"] <- 1.0
# # Ectomomyrmex
# Ponerinae_ACE_789t["1380"] <- -0.8
# # Plectroctena
# Ponerinae_ACE_789t["1176"] <- 0.0
# # Myopias
# Ponerinae_ACE_789t["1074"] <- -0.6
# # Leptogenys
# Ponerinae_ACE_789t["916"] <- -0.8

# ### Use BM per clades at 80 My to produce missing data
#
# all_clades <- cut_phylo_for_focal_time(tree = Ponerinae_tree_old_calib,
#                                        focal_time = 80,
#                                        keep_tip_labels = FALSE)$tip.label
#
# sim_ln_HW_old_calib <- setNames(object = rep(NA, length(Ponerinae_tree_old_calib$tip.label)), nm = Ponerinae_tree_old_calib$tip.label)
# for (i in seq_along(all_clades))
# {
#   # i <- 1
#   clade_i <- all_clades[i]
#
#   all_descendents_i <- phytools::getDescendants(tree = Ponerinae_tree_old_calib, node = clade_i)
#
#   if (length(all_descendents_i) > 2)
#   {
#     ## Extract clade tree
#     clade_tree_i <- extract.clade(phy = Ponerinae_tree_old_calib,
#                                   node = clade_i)
#     clade_tree_789t_i <- drop.tip(phy = clade_tree_i, tip = clade_tree_i$tip.label[!(clade_tree_i$tip.label %in% names(ln_HW_789t))])
#
#     if (length(clade_tree_789t_i$tip.label) > 2)
#     {
#       ## Fit BM
#       fit_BM_i <- geiger::fitContinuous(phy = clade_tree_789t_i,
#                                         dat = ln_HW_789t[clade_tree_789t_i$tip.label])
#
#       ## Get matching ancestral value
#       MRCA_i <- ape::getMRCA(phy = Ponerinae_phylogeny_789t_calibrated_short_names, tip = clade_tree_789t_i$tip.label)
#       ACE_i <- Ponerinae_ACE_789t[as.character(MRCA_i)]
#
#       ## Simulate data
#       sim_HW_i <- fastBM(tree = clade_tree_i,
#                          a = ACE_i,
#                          sig2 = fit_BM_i$opt$sigsq/3,
#                          nsim = 1)
#
#       ## Store data
#       sim_ln_HW_old_calib[names(sim_HW_i)] <- sim_HW_i
#     } else {
#       next
#     }
#   } else {
#     next
#   }
# }
#
# # Replace simulated data with real data if available
# # sim_ln_HW_old_calib[names(ln_HW_789t)] <- ln_HW_789t
#
# # Fix last issues
# which(is.na(sim_ln_HW_old_calib))
# sim_ln_HW_old_calib["Platythyrea_viehmeyeri"] <- log(1.400)
# sim_ln_HW_old_calib["Platythyrea_strenua"] <- log(1.500)
# sim_ln_HW_old_calib["Neoponera_bucki"] <- log(1.230)
# sim_ln_HW_old_calib["Neoponera_foetida"] <- log(2.200)
# sim_ln_HW_old_calib["Neoponera_fisheri"] <- log(2.200)
# sim_ln_HW_old_calib["Platythyrea_arnoldi"] <- log(1.500)


# ### Use BM per genus to produce missing data
#
# ## Get Genus MRCA
#
# all_genus_names <- unique(str_remove(string = Ponerinae_tree_old_calib$tip.label, pattern = "_.*"))
# all_genus_names <- all_genus_names[order(all_genus_names)]
#
# all_taxa_per_genus_list <- list()
# for (i in seq_along(all_genus_names))
# {
#   # i <- 1
#
#   # Extract genus
#   genus_i <- all_genus_names[i]
#   # Extract associated tip labels
#   all_taxa_i <- Ponerinae_tree_old_calib$tip.label[str_which(string = Ponerinae_tree_old_calib$tip.label, pattern = genus_i)]
#   # Store taxa
#   all_taxa_per_genus_list[[i]] <- all_taxa_i
# }
# names(all_taxa_per_genus_list) <- all_genus_names
#
# # Adjust outliers
# all_taxa_per_genus_list[["Neoponera"]] <- setdiff(all_taxa_per_genus_list[["Neoponera"]], "Neoponera_bucki")
# all_taxa_per_genus_list[["Pachycondyla"]] <- setdiff(all_taxa_per_genus_list[["Pachycondyla"]], "Pachycondyla_vidua")
# all_taxa_per_genus_list[["Pachycondyla"]] <- setdiff(all_taxa_per_genus_list[["Pachycondyla"]], "Pachycondyla_solitaria")
# all_taxa_per_genus_list[["Pachycondyla"]] <- setdiff(all_taxa_per_genus_list[["Pachycondyla"]], "Pachycondyla_jonesii")
# all_taxa_per_genus_list[["Brachyponera"]] <- setdiff(all_taxa_per_genus_list[["Brachyponera"]], "Brachyponera_sennaarensis")
# all_taxa_per_genus_list[["Anochetus"]] <- setdiff(all_taxa_per_genus_list[["Anochetus"]], "Anochetus_filicornis")
#
# all_descendant_taxa_per_genus_list <- list()
# MRCA_list <- c()
# # Get MRCA
# for (i in seq_along(all_genus_names))
# {
#   # i <- 1
#
#   # Extract list of taxa in genus
#   all_taxa_i <- all_taxa_per_genus_list[[i]]
#
#   # Get MRCA
#   if (length(all_taxa_i) > 1) {
#     MRCA_i <- getMRCA(phy = Ponerinae_tree_old_calib, tip = all_taxa_i)
#   } else {
#     MRCA_i <- which(Ponerinae_tree_old_calib$tip.label == all_taxa_i)
#   }
#
#   # Store MRCA
#   MRCA_list[i] <- MRCA_i
#
#   # Get descendants
#   all_descendant_taxa_i <- getDescendants(tree = Ponerinae_tree_old_calib, node = MRCA_i)
#
#   # Store descendant taxa
#   all_descendant_taxa_per_genus_list[[i]] <- all_taxa_i
# }
# names(all_descendant_taxa_per_genus_list) <- all_genus_names
# names(MRCA_list) <- all_genus_names
#
# # Check MRCA visually
# pdf(file = "./phlyo_old_calib.pdf", width = 30, height = 200)
# plot(Ponerinae_tree_old_calib)
# nodelabels()
# dev.off()
#
# ## Fit BM on reduced tree and simulate data on full tree
#
# set.seed(seed = 654321)
#
# MRCA_789t_list <- rep(NA, times = length(all_descendant_taxa_per_genus_list))
# sim_ln_HW_old_calib <- setNames(object = rep(NA, length(Ponerinae_tree_old_calib$tip.label)), nm = Ponerinae_tree_old_calib$tip.label)
# for (i in seq_along(all_descendant_taxa_per_genus_list))
# {
#   # i <- 15
#   genus_i <- names(all_descendant_taxa_per_genus_list)[i]
#   MRCA_i <- MRCA_list[i]
#   all_descendents_i <- all_descendant_taxa_per_genus_list[[i]]
#
#   if (length(all_descendents_i) > 2)
#   {
#     ## Extract clade tree
#     clade_tree_i <- extract.clade(phy = Ponerinae_tree_old_calib,
#                                   node = MRCA_i)
#     clade_tree_789t_i <- drop.tip(phy = clade_tree_i, tip = clade_tree_i$tip.label[!(clade_tree_i$tip.label %in% names(ln_HW_789t))])
#
#     if (length(clade_tree_789t_i$tip.label) > 2)
#     {
#       ## Fit BM
#       fit_BM_i <- geiger::fitContinuous(phy = clade_tree_789t_i,
#                                         dat = ln_HW_789t[clade_tree_789t_i$tip.label])
#
#       ## Get matching ancestral value
#       MRCA_789t_i <- ape::getMRCA(phy = Ponerinae_phylogeny_789t_calibrated_short_names, tip = clade_tree_789t_i$tip.label)
#       ACE_i <- Ponerinae_ACE_789t[as.character(MRCA_789t_i)]
#
#       ## Simulate data
#       sim_HW_i <- phytools::fastBM(tree = clade_tree_i,
#                          a = ACE_i,
#                          sig2 = fit_BM_i$opt$sigsq/3, # Divide by 3 to reduce variance
#                          nsim = 1)
#
#       ## Store data
#       sim_ln_HW_old_calib[names(sim_HW_i)] <- sim_HW_i
#       MRCA_789t_list[i] <- MRCA_789t_i
#     } else {
#       next
#     }
#   } else {
#     next
#   }
# }
# names(MRCA_789t_list) <- names(all_descendant_taxa_per_genus_list)
#
# view(sim_ln_HW_old_calib)
#
# # # Replace simulated data with real data if available
# # sim_ln_HW_old_calib[names(ln_HW_789t)] <- ln_HW_789t
#
# # Fix last issues
# which(is.na(sim_ln_HW_old_calib))
#
# ## Add to Ponerinae_trait_tip_data
# Ponerinae_trait_tip_data$sim_ln_HW_old_calib <- sim_ln_HW_old_calib[match(x = Ponerinae_trait_tip_data$Taxa, table = names(sim_ln_HW_old_calib))]
# Ponerinae_trait_tip_data$sim_HW_old_calib <- exp(sim_ln_HW_old_calib)[match(x = Ponerinae_trait_tip_data$Taxa, table = names(sim_ln_HW_old_calib))]
#
# View(Ponerinae_trait_tip_data)
#
# ## Check consistency on contMap
# ln_HW_old_calib <- setNames(object = Ponerinae_trait_tip_data$sim_ln_HW_old_calib,
#                          nm = Ponerinae_trait_tip_data$Taxa)
# ln_HW_old_calib <- ln_HW_old_calib[Ponerinae_tree_old_calib$tip.label]
#
# Ponerinae_trait_cont_old_calib <- prepare_trait_data(tip_data = ln_HW_old_calib,
#                                            trait_data_type = "continuous",
#                                            phylo = Ponerinae_tree_old_calib,
#                                            seed = 1234) # Set seed for reproducibility
#
# plot_contMap(Ponerinae_trait_cont_old_calib$contMap, color_scale = c("darkblue", "yellow", "red"))
#
# ## Export the updated Ponerinae_trait_tip_data in deepSTRAPP
# usethis::use_data(Ponerinae_trait_tip_data, overwrite = TRUE)


### Different approach: produce trait data based on rates to ensure significance and eventually add noise!

# Extract rates
Ponerinae_net_div_rates <- Ponerinae_BAMM_object_old_calib$meanTipLambda - Ponerinae_BAMM_object_old_calib$meanTipMu
# Apply continuous transformation
fake_cont_tip_data <- log(Ponerinae_net_div_rates + 0.1) + 2.5
# Scale between 0 and 1.2
hist(fake_cont_tip_data)
names(fake_cont_tip_data) <- Ponerinae_BAMM_object_old_calib$tip.label

# Plot
plot(x = Ponerinae_net_div_rates, y = fake_cont_tip_data)

# Quick test for correlation
quick_test <- BAMMtools::traitDependentBAMM(ephy = Ponerinae_BAMM_object_old_calib, traits = fake_cont_tip_data, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE)
str(quick_test, 2)

## Adjust fake_cont_tip_data for deepSTRAPP

fake_cont_tip_data_for_deepSTRAPP <- fake_cont_tip_data

# Reverse value
fake_cont_tip_data_for_deepSTRAPP <- (-1*fake_cont_tip_data_for_deepSTRAPP)+1.2
hist(fake_cont_tip_data_for_deepSTRAPP)
# Add noise
fake_cont_tip_data_for_deepSTRAPP <- fake_cont_tip_data_for_deepSTRAPP + rnorm(n = length(fake_cont_tip_data_for_deepSTRAPP), mean = 0, sd = sd(fake_cont_tip_data_for_deepSTRAPP)/5)
names(fake_cont_tip_data_for_deepSTRAPP) <- Ponerinae_BAMM_object_old_calib$tip.label

# Adjust Plectroctena and Ectomomyrmex to buffer significance in present
fake_cont_tip_data_for_deepSTRAPP[str_which(string = names(fake_cont_tip_data_for_deepSTRAPP), pattern = "Plectroctena")] <- fake_cont_tip_data_for_deepSTRAPP[str_which(string = names(fake_cont_tip_data_for_deepSTRAPP), pattern = "Plectroctena")] + 0.5
fake_cont_tip_data_for_deepSTRAPP[str_which(string = names(fake_cont_tip_data_for_deepSTRAPP), pattern = "Ectomomyrmex")] <- fake_cont_tip_data_for_deepSTRAPP[str_which(string = names(fake_cont_tip_data_for_deepSTRAPP), pattern = "Ectomomyrmex")] + 0.5
Hypoponera_Madagascar <- c("Hypoponera_mg024", "Hypoponera_sc_akir", "Hypoponera_sc_mano", "Hypoponera_sc_mora", "Hypoponera_mg025", "Hypoponera_sc_befi", "Hypoponera_sc_bina",
                           "Hypoponera_sc_tamp", "Hypoponera_sc_ambo", "Hypoponera_sc_ampa", "Hypoponera_sc_andr", "Hypoponera_sc_ivoh",  "Hypoponera_sc_mah", "Hypoponera_sc_anta",
                           "Hypoponera_sc_bema",  "Hypoponera_sc_befa", "Hypoponera_sc_zaha", "Hypoponera_sb_mano",  "Hypoponera_sc_ando", "Hypoponera_sc_beta", "Hypoponera_sc_ant",
                           "Hypoponera_sc_beka", "Hypoponera_sc_nosy", "Hypoponera_sc_amba", "Hypoponera_sc_mand", "Hypoponera_sc_isa", "Hypoponera_sc_maka")
fake_cont_tip_data_for_deepSTRAPP[names(fake_cont_tip_data_for_deepSTRAPP) %in% Hypoponera_Madagascar] <- fake_cont_tip_data_for_deepSTRAPP[names(fake_cont_tip_data_for_deepSTRAPP) %in% Hypoponera_Madagascar] + 0.8

Hypoponera_NW <- c("Hypoponera_dias12", "Hypoponera_dias20", "Hypoponera_psw_pe01", "Hypoponera_psw_pe03", "Hypoponera_dias25", "Hypoponera_jtl031", "Hypoponera_jtl006", "Hypoponera_menozzii", "Hypoponera_distinguenda",
                   "Hypoponera_collegiana", "Hypoponera_dias19_2", "Hypoponera_trigona", "Hypoponera_vc01", "Hypoponera_idelettae", "Hypoponera_perplexa", "Hypoponera_jtl010", "Hypoponera_jtl035", "Hypoponera_reichenspergeri",
                   "Hypoponera_jtl032", "Hypoponera_jtl034", "Hypoponera_jtl002", "Hypoponera_jtl023", "Hypoponera_clavatula")
fake_cont_tip_data_for_deepSTRAPP[names(fake_cont_tip_data_for_deepSTRAPP) %in% Hypoponera_NW] <- fake_cont_tip_data_for_deepSTRAPP[names(fake_cont_tip_data_for_deepSTRAPP) %in% Hypoponera_NW] + 0.7


# Plot
plot(x = Ponerinae_net_div_rates, y = fake_cont_tip_data_for_deepSTRAPP)

fake_cont_tip_data_for_deepSTRAPP_reordered <- fake_cont_tip_data_for_deepSTRAPP[match(x = Ponerinae_trait_tip_data$Taxa, table = names(fake_cont_tip_data_for_deepSTRAPP))]

# Quick test for correlation
quick_test <- BAMMtools::traitDependentBAMM(ephy = Ponerinae_BAMM_object_old_calib, traits = fake_cont_tip_data_for_deepSTRAPP_reordered, rate = 'net diversification', return.full = TRUE, method = 'spearman', logrates = FALSE)
str(quick_test, 2)

## Add to Ponerinae_trait_tip_data
Ponerinae_trait_tip_data$fake_cont_tip_data <- fake_cont_tip_data_for_deepSTRAPP_reordered
Ponerinae_trait_tip_data

## Export the updated Ponerinae_trait_tip_data in deepSTRAPP
usethis::use_data(Ponerinae_trait_tip_data, overwrite = TRUE)


### 2.2/ Produce fake_cat_3lvl_tip_data ####

## Convert to three-factor categorical traits
fake_cat_3lvl_tip_data <- fake_cont_tip_data
fake_cat_3lvl_tip_data[fake_cont_tip_data > 0] <- "arboreal"
fake_cat_3lvl_tip_data[fake_cont_tip_data > 0.5] <- "subterranean"
fake_cat_3lvl_tip_data[fake_cont_tip_data > 1.1] <- "terricolous"
table(fake_cat_3lvl_tip_data)
names(fake_cat_3lvl_tip_data) <- Ponerinae_BAMM_object_old_calib$tip.label
by(data = Ponerinae_net_div_rates, INDICES = fake_cat_3lvl_tip_data, FUN = mean)

## Select color scheme for states
colors_per_states <- c("forestgreen", "sienna", "goldenrod")
names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")

# # Adjust Plectroctena and Ectomomyrmex to buffer significance in present
# fake_cat_3lvl_tip_data[str_which(string = names(fake_cat_3lvl_tip_data), pattern = "Plectroctena")] <- "arboreal"
# fake_cat_3lvl_tip_data[str_which(string = names(fake_cat_3lvl_tip_data), pattern = "Ectomomyrmex")] <- "arboreal"

## Add to Ponerinae_trait_tip_data
Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data <- fake_cat_3lvl_tip_data[match(x = Ponerinae_trait_tip_data$Taxa, table = names(fake_cat_3lvl_tip_data))]
Ponerinae_trait_tip_data

## Export the updated Ponerinae_trait_tip_data in deepSTRAPP
usethis::use_data(Ponerinae_trait_tip_data, overwrite = TRUE)


# Quick test for differences
quick_test <- BAMMtools::traitDependentBAMM(ephy = Ponerinae_BAMM_object_old_calib, traits = fake_cat_3lvl_tip_data, rate = 'net diversification', return.full = TRUE, method = 'kruskal', logrates = FALSE, reps = 1000)
str(quick_test, 2)


## Plot data on tips
pdf(file = "./Ponerinae_fake_cat_3lvl_data_old_calib_on_phylo.pdf", width = 20, height = 200)

# Set plotting parameters
par(mar = c(0.1,0.1,0.1,0.1), oma = c(0,0,0,0)) # bltr
# Graph presence/absence using plotTree.datamatrix
range_map <- phytools::plotTree.datamatrix(
  tree = Ponerinae_tree_old_calib,
  X = as.data.frame(fake_cat_3lvl_tip_data),
  fsize = 0.7, yexp = 1.1,
  header = TRUE, xexp = 1.25,
  colors = colors_per_states)

# Get plot info in "last_plot.phylo"
plot_info <- get("last_plot.phylo", envir=.PlotPhyloEnv)

# Add time line

# Extract root age
root_age <- max(phytools::nodeHeights(Ponerinae_tree_old_calib))

# Define ticks
# ticks_labels <- seq(from = 0, to = 100, by = 20)
ticks_labels <- seq(from = 0, to = 120, by = 20)
axis(side = 1, pos = 0, at = (-1 * ticks_labels) + root_age, labels = ticks_labels, cex.axis = 1.5)
legend(x = root_age/2,
       y = 0 - 5, adj = 0,
       bty = "n", legend = "", title = "Time  [My]", title.cex = 1.5)

# Add a legend
legend(x = plot_info$x.lim[2] - 10,
       y = mean(plot_info$y.lim),
       # adj = c(0,0),
       # x = "topleft",
       legend = c("Absence", "Presence"),
       pch = 22, pt.bg = c("white","gray30"), pt.cex =  1.8,
       cex = 1.5, bty = "n")

dev.off()


### 2.3/ Produce fake_cat_2lvl_tip_data ####

## Transform into two categories
fake_cat_2lvl_tip_data <- fake_cont_tip_data
fake_cat_2lvl_tip_data[fake_cont_tip_data > 0] <- "large"
fake_cat_2lvl_tip_data[fake_cont_tip_data > 0.65] <- "small"
table(fake_cat_2lvl_tip_data)
names(fake_cat_2lvl_tip_data) <- Ponerinae_BAMM_object_old_calib$tip.label
by(data = Ponerinae_net_div_rates, INDICES = fake_cat_2lvl_tip_data, FUN = mean)

# Adjust Plectroctena and Ectomomyrmex to buffer significance in present
fake_cat_2lvl_tip_data[str_which(string = names(fake_cat_2lvl_tip_data), pattern = "Plectroctena")] <- "large"
fake_cat_2lvl_tip_data[str_which(string = names(fake_cat_2lvl_tip_data), pattern = "Ectomomyrmex")] <- "large"
# fake_cat_2lvl_tip_data[str_which(string = names(fake_cat_2lvl_tip_data), pattern = "Ponera")] <- "small"
Hypoponera_Madagascar <- c("Hypoponera_mg024", "Hypoponera_sc_akir", "Hypoponera_sc_mano", "Hypoponera_sc_mora", "Hypoponera_mg025", "Hypoponera_sc_befi", "Hypoponera_sc_bina",
                           "Hypoponera_sc_tamp", "Hypoponera_sc_ambo", "Hypoponera_sc_ampa", "Hypoponera_sc_andr", "Hypoponera_sc_ivoh",  "Hypoponera_sc_mah", "Hypoponera_sc_anta",
                           "Hypoponera_sc_bema",  "Hypoponera_sc_befa", "Hypoponera_sc_zaha", "Hypoponera_sb_mano",  "Hypoponera_sc_ando", "Hypoponera_sc_beta", "Hypoponera_sc_ant",
                           "Hypoponera_sc_beka", "Hypoponera_sc_nosy", "Hypoponera_sc_amba", "Hypoponera_sc_mand", "Hypoponera_sc_isa", "Hypoponera_sc_maka")
fake_cat_2lvl_tip_data[names(fake_cat_2lvl_tip_data) %in% Hypoponera_Madagascar] <- "large"

# Quick test for differences
quick_test_2lvl <- BAMMtools::traitDependentBAMM(ephy = Ponerinae_BAMM_object_old_calib, traits = fake_cat_2lvl_tip_data, rate = 'net diversification', return.full = TRUE, method = 'mann-whitney', logrates = FALSE, reps = 1000)
str(quick_test_2lvl, 2)

## Add to Ponerinae_trait_tip_data
Ponerinae_trait_tip_data$fake_cat_2lvl_tip_data <- fake_cat_2lvl_tip_data[match(x = Ponerinae_trait_tip_data$Taxa, table = names(fake_cat_2lvl_tip_data))]
Ponerinae_trait_tip_data

## Export the updated Ponerinae_trait_tip_data in deepSTRAPP
usethis::use_data(Ponerinae_trait_tip_data, overwrite = TRUE)


##### Trait evolution data #####

### 3/ Generate categorical trait evolution data for eel using 3-level factor #####

# Load phylogeny and tip data
library(phytools)
data(eel.tree)
data(eel.data)

# Transform feeding mode data into a 3-level factor
eel_tip_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
eel_tip_data <- as.character(eel_tip_data)
eel_tip_data[c(1, 5, 6, 7, 10, 11, 15, 16, 17, 24, 25, 28, 30, 51, 52, 53, 55, 58, 60)] <- "kiss"
eel_tip_data <- stats::setNames(eel_tip_data, rownames(eel.data))
table(eel_tip_data)

# Manually define a Q_matrix for rate classes of state transition to use in the 'matrix' model
# Does not allow transitions from state 1 ("bite") to state 2 ("kiss") or state 3 ("suction")
# Does not allow transitions from state 3 ("suction") to state 1 ("bite")
# Set symmetrical rates between state 2 ("kiss") and state 3 ("suction")
Q_matrix = rbind(c(NA, 0, 0), c(1, NA, 2), c(0, 2, NA))

# Set colors per states
colors_per_states <- c("limegreen", "orange", "dodgerblue")
names(colors_per_states) <- c("bite", "kiss", "suction")

eel_cat_3lvl_data <- prepare_trait_data(
   tip_data = eel_tip_data, phylo = eel.tree,
   trait_data_type = "categorical",
   colors_per_states = colors_per_states,
   evolutionary_models = c("ER", "SYM", "ARD", "meristic", "matrix"),
   Q_matrix = Q_matrix,
   nb_simulations = 1000, # Set to 10 to save time.
   # But recommended value = 1000.
   plot_map = TRUE,
   plot_overlay = TRUE,
   return_best_model_fit = TRUE,
   return_model_selection_df = TRUE)

# Explore output
plot(eel_cat_3lvl_data$densityMaps[[1]]) # densityMap for state n°1 ("bite")
eel_cat_3lvl_data$model_selection_df # Summary of model selection
# Parameter estimates and optimization summary of the best model
# (Here, the best model is ER)
print(eel_cat_3lvl_data$best_model_fit)$ # Summary of the best evolutionary model
  eel_cat_3lvl_data$ace # Posterior probabilities of each state (= ACE) at internal nodes

# Export in deepSTRAPP
usethis::use_data(eel_cat_3lvl_data, overwrite = TRUE)


### 4/ Generate biogeographic evolution data for eel #####

# Load phylogeny and tip data
library(phytools)
data(eel.tree)
data(eel.data)

# Transform feeding mode data into biogeographic data with ranges A, B, and AB.
eel_tip_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
eel_tip_data <- as.character(eel_tip_data)
eel_tip_data[eel_tip_data == "bite"] <- "A"
eel_tip_data[eel_tip_data == "suction"] <- "B"
eel_tip_data[c(5, 6, 7, 15, 25, 32, 33, 34, 50, 52, 57, 58, 59)] <- "AB"
eel_tip_data <- stats::setNames(eel_tip_data, rownames(eel.data))
table(eel_tip_data)

colors_per_levels <- c("dodgerblue3", "gold")
names(colors_per_levels) <- c("A", "B")

eel_biogeo_data <- prepare_trait_data(
  tip_data = eel_tip_data,
  trait_data_type = "biogeographic",
  phylo = eel.tree,
  evolutionary_models = c("BAYAREALIKE", "DIVALIKE", "DEC", "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J"), # Default = "DEC" for biogeographic
  prefix_for_files = "eel",
  max_range_size = 2,
  split_multi_area_ranges = TRUE, # Set to TRUE to display the two outputs
  nb_simulations = 100, # Reduce the number of Stochastic Mapping simulations to save time (Default = '1000')
  colors_per_levels = colors_per_levels,
  return_BSM = TRUE,
  return_simmaps = TRUE,
  return_best_model_fit = TRUE,
  return_model_selection_df = TRUE,
  verbose = TRUE)

# Explore output
str(eel_biogeo_data, 1)
eel_biogeo_data$model_selection_df # Summary of model selection
# Parameter estimates and optimization summary of the best model
# (Here, the best model is DEC+J)

# Posterior probabilities of each state (= ACE) at internal nodes
eel_biogeo_data$ace # Only with unique areas
eel_biogeo_data$ace_all_ranges # Including multi-area ranges (Here, AB)

# Plot densityMaps
plot(eel_biogeo_data$densityMaps[[1]]) # densityMap for range n°1 ("A")
plot_densityMaps_overlay(eel_biogeo_data$densityMaps) # densityMaps with all unique areas overlaid
plot_densityMaps_overlay(eel_biogeo_data$densityMaps_all_ranges) # densityMaps with all ranges (including multi-area ranges) overlaid

# Export in deepSTRAPP
usethis::use_data(eel_biogeo_data, overwrite = TRUE)

## Remove problematic BioGeoBEARS classes so the object can be loaded by CRAN even if BioGeoBEARS is not installed
eel_biogeo_data_for_CRAN <- eel_biogeo_data

class(eel_biogeo_data_for_CRAN$best_model_fit)
best_model_fit_unclassed <- unclass(eel_biogeo_data_for_CRAN$best_model_fit)
class(best_model_fit_unclassed)
eel_biogeo_data_for_CRAN$best_model_fit <- best_model_fit_unclassed

class(eel_biogeo_data_for_CRAN$best_model_fit$inputs$BioGeoBEARS_model_object)
BioGeoBEARS_model_object <- lapply(slotNames(eel_biogeo_data_for_CRAN$best_model_fit$inputs$BioGeoBEARS_model_object), function(s) slot(eel_biogeo_data_for_CRAN$best_model_fit$inputs$BioGeoBEARS_model_object, s))
names(BioGeoBEARS_model_object) <- slotNames(eel_biogeo_data_for_CRAN$best_model_fit$inputs$BioGeoBEARS_model_object)
class(BioGeoBEARS_model_object)
eel_biogeo_data_for_CRAN$best_model_fit$inputs$BioGeoBEARS_model_object <- BioGeoBEARS_model_object

class(eel_biogeo_data_for_CRAN$best_model_fit$outputs)
outputs_list <- lapply(slotNames(eel_biogeo_data_for_CRAN$best_model_fit$outputs), function(s) slot(eel_biogeo_data_for_CRAN$best_model_fit$outputs, s))
names(outputs_list) <- slotNames(eel_biogeo_data_for_CRAN$best_model_fit$outputs)
class(outputs_list)
eel_biogeo_data_for_CRAN$best_model_fit$outputs <- outputs_list

# Export in deepSTRAPP
# usethis::use_data(eel_biogeo_data_for_CRAN, overwrite = TRUE)

# Need to be renamed before export so it loads with the same name
eel_biogeo_data <- eel_biogeo_data_for_CRAN
usethis::use_data(eel_biogeo_data, overwrite = TRUE)

### 5/ Generate continuous trait evolution data for Ponerinae ants_old_calib #####

# Load phylogeny
data("Ponerinae_tree_old_calib", package = "deepSTRAPP")

## Load the BAMM_object summarizing 1000 posterior samples of BAMM
data("Ponerinae_BAMM_object_old_calib", package = "deepSTRAPP")

## Load trait df
data("Ponerinae_trait_tip_data", package = "deepSTRAPP")

dim(Ponerinae_trait_tip_data)
View(Ponerinae_trait_tip_data)

## Prepare trait data

# Extract continuous trait data as a named vector
# Ponerinae_tip_data_ln_HW <- setNames(object = Ponerinae_trait_tip_data$sim_ln_HW,
#                                  nm = Ponerinae_trait_tip_data$Taxa)
Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
                                    nm = Ponerinae_trait_tip_data$Taxa)

## Run prepare_trait_data with default options
# For continuous trait, a BM model is assumed by default.
Ponerinae_cont_data_old_calib <- prepare_trait_data(# tip_data = Ponerinae_data_ln_HW,
  tip_data = Ponerinae_cont_tip_data,
  trait_data_type = "continuous",
  phylo = Ponerinae_tree_old_calib,
  seed = 1234,
  evolutionary_models = "BM",
  plot_map = FALSE,
  run_stochastic_maps = TRUE,
  nb_simulations = 100, # Run 100 simulations of trait evolution
  verbose = TRUE)

# Explore output
str(Ponerinae_cont_data_old_calib, 1)

## Export Ponerinae_cat_2lvl_data_old_calib in deepSTRAPP
usethis::use_data(Ponerinae_cont_data_old_calib, overwrite = TRUE)


### 6/ Generate categorical (2-lvl) trait evolution data for Ponerinae ants_old_calib #####

## Load data

# Load trait df
data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
# Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")

## Prepare trait data

# Extract categorical data with 3-levels
Ponerinae_cat_2lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_2lvl_tip_data,
                                        nm = Ponerinae_trait_tip_data$Taxa)
table(Ponerinae_cat_2lvl_tip_data)

# Select color scheme for states
colors_per_states <- c("darkblue", "lightblue")
names(colors_per_states) <- c("large", "small")

## Produce densityMaps using stochastic character mapping based on an ARD Mk model
Ponerinae_cat_2lvl_data_old_calib <- prepare_trait_data(
  tip_data = Ponerinae_cat_2lvl_tip_data,
  phylo = Ponerinae_tree_old_calib,
  trait_data_type = "categorical",
  seed = 1234,
  colors_per_levels = colors_per_states,
  evolutionary_models = "ARD",
  nb_simulations = 100,
  return_best_model_fit = TRUE,
  return_model_selection_df = TRUE,
  plot_map = FALSE)

# Explore output
plot(Ponerinae_cat_2lvl_data_old_calib$densityMaps[[1]]) # densityMap for state n°1 ("large")
Ponerinae_cat_2lvl_data_old_calib$model_selection_df # Summary of model selection
# Parameter estimates and optimization summary of the best model
# (Here, the best model is ARD)
print(Ponerinae_cat_2lvl_data_old_calib$best_model_fit) # Summary of the best evolutionary model
Ponerinae_cat_2lvl_data_old_calib$ace # Posterior probabilities of each state (= ACE) at internal nodes

# Export in deepSTRAPP
usethis::use_data(Ponerinae_cat_2lvl_data_old_calib, overwrite = TRUE)

### 7/ Generate categorical (3-lvl) trait evolution data for Ponerinae ants_old_calib #####

## Load data

# Load trait df
data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
# Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")

## Prepare trait data

# Extract categorical data with 3-levels
Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
                                        nm = Ponerinae_trait_tip_data$Taxa)
table(Ponerinae_cat_3lvl_tip_data)

# Select color scheme for states
colors_per_states <- c("forestgreen", "sienna", "goldenrod")
names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")

## Produce densityMaps using stochastic character mapping based on an ARD Mk model
Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
   tip_data = Ponerinae_cat_3lvl_tip_data,
   phylo = Ponerinae_tree_old_calib,
   trait_data_type = "categorical",
   seed = 1234,
   colors_per_levels = colors_per_states,
   evolutionary_models = "ARD",
   nb_simulations = 100,
   return_best_model_fit = TRUE,
   return_model_selection_df = TRUE,
   plot_map = FALSE)

# Explore output
plot(Ponerinae_cat_3lvl_data_old_calib$densityMaps[[1]]) # densityMap for state n°1 ("arboreal")
Ponerinae_cat_3lvl_data_old_calib$model_selection_df # Summary of model selection
# Parameter estimates and optimization summary of the best model
# (Here, the best model is ARD)
print(Ponerinae_cat_3lvl_data_old_calib$best_model_fit) # Summary of the best evolutionary model
Ponerinae_cat_3lvl_data_old_calib$ace # Posterior probabilities of each state (= ACE) at internal nodes

# Export in deepSTRAPP
usethis::use_data(Ponerinae_cat_3lvl_data_old_calib, overwrite = TRUE)


### 8/ Generate biogeographic evolution data for Ponerinae_old_calib  #####

### Run Biogeographic inference for Ponerinae across New World and Old World ranges using the old ill-calibrated phylogeny for demonstration purpose

## Load data

# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
# Load range data
data(Ponerinae_binary_range_table, package = "deepSTRAPP")
# Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(Ponerinae_BAMM_object, package = "deepSTRAPP")

## Prepare range data for Old World vs. New World

# No overlap in ranges
table(Ponerinae_binary_range_table$Old_World, Ponerinae_binary_range_table$New_World)

Ponerinae_NO_tip_data <- stats::setNames(object = Ponerinae_binary_range_table$Old_World,
                                         nm = Ponerinae_binary_range_table$Taxa)
Ponerinae_NO_tip_data <- as.character(Ponerinae_ON_tip_data)
Ponerinae_NO_tip_data[Ponerinae_NO_tip_data == "TRUE"] <- "O" # O = Old World
Ponerinae_NO_tip_data[Ponerinae_NO_tip_data == "FALSE"] <- "N" # N = New World
names(Ponerinae_NO_tip_data) <- Ponerinae_binary_range_table$Taxa
table(Ponerinae_NO_tip_data)

colors_per_levels <- c("mediumpurple2", "peachpuff2")
names(colors_per_levels) <- c("N", "O")

## Run evolutionary models for biogeographic inferences
Ponerinae_biogeo_data_old_calib <- prepare_trait_data(
  tip_data = Ponerinae_NO_tip_data,
  trait_data_type = "biogeographic",
  phylo = Ponerinae_tree_old_calib,
  seed = 1234,
  evolutionary_models = "DEC+J", # Default = "DEC" for biogeographic
  prefix_for_files = "Ponerinae_old_calib",
  max_range_size = 2,
  split_multi_area_ranges = TRUE, # Set to TRUE to display the two outputs
  nb_simulations = 100,
  colors_per_levels = colors_per_levels,
  return_model_selection_df = TRUE,
  verbose = TRUE)

## Explore output
str(Ponerinae_biogeo_data_old_calib, 1)

# Posterior probabilities of each state (= ACE) at internal nodes
Ponerinae_biogeo_data_old_calib$ace # Only with unique areas
Ponerinae_biogeo_data_old_calib$ace_all_ranges # Including multi-area ranges (Here, AB)

## Plot densityMaps
# densityMap for range n°1 (N = "New World")
plot(Ponerinae_biogeo_data_old_calib$densityMaps[[1]])
# densityMaps with all unique areas overlaid
plot_densityMaps_overlay(Ponerinae_biogeo_data_old_calib$densityMaps)
# densityMaps with all ranges (including multi-area ranges) overlaid
plot_densityMaps_overlay(Ponerinae_biogeo_data_old_calib$densityMaps_all_ranges)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_biogeo_data_old_calib, overwrite = TRUE)


### 9/ Extract continuous trait data for Ponerinae_trait_data_10My ####

# Load phylogeny
data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
# Load trait df
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")

# Extract log(head with) data
Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
                                     nm = Ponerinae_trait_tip_data$Taxa)

# Get Ancestral Character Estimates based on a Brownian Motion model
# To obtain values at internal nodes
Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)

# Run a Stochastic Mapping based on a Brownian Motion model
# to interpolate values along branches and obtain a "contMap" object
Ponerinae_contMap <- phytools::contMap(Ponerinae_tree, x = Ponerinae_cont_tip_data,
                                       res = 100, # Number of time steps
                                       plot = FALSE)

# Set focal time to 10 Mya
focal_time <- 10

## Extract trait data and update contMap for the given focal_time

# Extract from tip data and ML estimates of ancestral characters (values are true ML)
Ponerinae_trait_cont_tip_data_10My <- extract_most_likely_trait_values_for_focal_time(
  contMap = Ponerinae_contMap,
  ace = Ponerinae_ACE,
  tip_data = Ponerinae_cont_tip_data,
  trait_data_type = "continuous",
  focal_time = focal_time,
  update_Map = TRUE,
  keep_tip_labels = TRUE)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_trait_cont_tip_data_10My, overwrite = TRUE)



##### BAMM data #####

### 10/ Save diversification template file as .rds #####

# Load it from the /inst/ directory (works only in source/bundled version of the package)
BAMM_template_diversification <- readLines(con = file.path("./inst/BAMM_template_diversification.txt"))

# Export in deepSTRAPP
usethis::use_data(BAMM_template_diversification, overwrite = TRUE)


### 11/ Generate whale_BAMM_object ####

library(phytools)
data(whale.tree)

whale_BAMM_object <- prepare_diversification_data(
  BAMM_install_directory_path = "./software/bamm-2.5.0/",
  phylo = whale.tree,
  prefix_for_files = "whale",
  numberOfGenerations = 100000, # Set low for the example
  BAMM_output_directory_path =  "./BAMM_outputs/")

# Load directly the result
data(whale_BAMM_object)

# Explore output
str(whale_BAMM_object, 1)
summary(whale_BAMM_object$numberEvents)

# Plot mean net diversification rates on the phylogeny
plot.bammdata(whale_BAMM_object, labels = TRUE)

# Export in deepSTRAPP
usethis::use_data(whale_BAMM_object, overwrite = TRUE)


### 12/ Include Ponerinae_BAMM_object ####

Ponerinae_BAMM_object <- readRDS(file = "../Ponerinae_Historical_Biogeography/outputs/BAMM/Ponerinae_MCC_phylogeny_1534t/BAMM_posterior_samples_data.rds")

usethis::use_data(Ponerinae_BAMM_object, overwrite = TRUE)


### 13/ Generate Ponerinae_BAMM_object_10My ####

## Load the BAMM_object summarizing 1000 posterior samples of BAMM.
data(Ponerinae_BAMM_object, package = "deepSTRAPP")

## Set focal-time to 10 My
focal_time = 10

## Update the BAMM object (May take several minutes to run)
Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
  BAMM_object = Ponerinae_BAMM_object,
  focal_time = focal_time,
  update_rates = TRUE, update_regimes = TRUE,
  update_tree = TRUE, update_plot = TRUE,
  update_all_elements = TRUE,
  keep_tip_labels = TRUE,
  verbose = TRUE)

str(Ponerinae_BAMM_object_10My, 1)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_BAMM_object_10My, overwrite = TRUE)


### 14/ Generate Ponerinae_BAMM_object_old_calib ####

# Load phylogeny
data("Ponerinae_tree_old_calib", package = "deepSTRAPP")
plot(Ponerinae_tree_old_calib, show.tip.label = FALSE)

# Run BAMM workflow with deepSTRAPP
Ponerinae_BAMM_object_old_calib <- prepare_diversification_data(
  BAMM_install_directory_path = "./software/bamm-2.5.0/",
  phylo = Ponerinae_tree_old_calib,
  prefix_for_files = "Ponerinae_old_calib",
  numberOfGenerations = 10^7, # Set high for optimal run, but will take ages
  BAMM_output_directory_path =  "./BAMM_outputs/")

# Explore output
str(Ponerinae_BAMM_object_old_calib, 1)

# Plot mean net diversification rates and regime shifts on the phylogeny
plot_BAMM_rates(Ponerinae_BAMM_object_old_calib,
                labels = FALSE, legend = TRUE)


# Export in deepSTRAPP
usethis::use_data(Ponerinae_BAMM_object_old_calib, overwrite = TRUE)



##### deepSTRAPP data #####


### 15/ Generate Ponerinae_deepSTRAPP_cont_old_calib_0_40 ####

## deepSTRAPP output for Ponerinae over time from 0 to 40My (steps = 5 My)
# Tests for continuous trait

## Load trait df
data("Ponerinae_trait_tip_data", package = "deepSTRAPP")

# Extract continuous trait data as a named vector
Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
                                    nm = Ponerinae_trait_tip_data$Taxa)

## Load trait evolution data
data("Ponerinae_cont_data_old_calib", package = "deepSTRAPP")

# Extract the contMap representing ML estimates of continuous trait evolution on the phylogeny
Ponerinae_contMap <- Ponerinae_cont_data_old_calib$contMap
plot_contMap(contMap = Ponerinae_contMap,
             color_scale = c("darkgreen", "limegreen", "orange", "red"))

# Extract the contMaps representing all stochastic simulations of continuous trait evolution
Ponerinae_contMaps <- Ponerinae_cont_data_old_calib$contMaps
# Plot simulation n°1
plot_contMap(contMap = Ponerinae_contMaps[[1]],
             color_scale = c("darkgreen", "limegreen", "orange", "red"))

# Extract the Ancestral Character Estimates (ACE) = ML estimates of trait values at nodes
Ponerinae_ACE <- Ponerinae_cont_data_old_calib$ace
head(Ponerinae_ACE)

## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
# nb_time_steps <- 5
time_step_duration <- 5
time_range <- c(0, 40)
# time_range <- c(0, 5)

# Run deepSTRAPP on net diversification rates
## This step is time-consuming. You can skip it and load directly the result if needed
Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- run_deepSTRAPP_over_time(
  contMap = Ponerinae_contMap, # Include contMap to plot ML estimates
  # contMaps = Ponerinae_contMaps, # Include contMaps to extract trait estimates across all simulations
  ace = Ponerinae_ACE,
  # tip_data = Ponerinae_tip_data_ln_HW,
  tip_data = Ponerinae_cont_tip_data,
  trait_data_type = "continuous",
  BAMM_object = Ponerinae_BAMM_object_old_calib,
  # nb_time_steps = nb_time_steps,
  time_range = time_range,
  time_step_duration = time_step_duration,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  seed = 1234,
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Ponerinae_deepSTRAPP_cont_old_calib_0_40, max.level = 1)

# Access p-values df
Ponerinae_deepSTRAPP_cont_old_calib_0_40$pvalues_summary_df

# Plot p-values of Spearman tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
  alpha = 0.05)

# Plot rates through time
plot_rates_through_time(deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
                        plot_CI = TRUE)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_deepSTRAPP_cont_old_calib_0_40, overwrite = TRUE)


### 16/ Generate Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40 ####

## deepSTRAPP output for Ponerinae over time from 0 to 40My (steps = 5 My)
# Tests for 2-factor categorical traits

## Load data

# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
Ponerinae_tree_old_calib$node.label <- NULL
# Load trait df
data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
# Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")

## Extract trait data

# Extract categorical data with 2-levels
Ponerinae_cat_2lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_2lvl_tip_data,
                                    nm = Ponerinae_trait_tip_data$Taxa)
table(Ponerinae_cat_2lvl_tip_data)

# Select color scheme for states
colors_per_states <- c("darkblue", "lightblue")
names(colors_per_states) <- c("large", "small")

## Plot data on tips
pdf(file = "./Ponerinae_cat_2lvl_data_old_calib_on_phylo.pdf", width = 20, height = 200)

# Set plotting parameters
par(mar = c(0.1,0.1,0.1,0.1), oma = c(0,0,0,0)) # bltr
# Graph presence/absence using plotTree.datamatrix
range_map <- phytools::plotTree.datamatrix(
  tree = Ponerinae_tree_old_calib,
  X = as.data.frame(Ponerinae_cat_2lvl_data),
  fsize = 0.7, yexp = 1.1,
  header = TRUE, xexp = 1.25,
  colors = colors_per_states)

# Get plot info in "last_plot.phylo"
plot_info <- get("last_plot.phylo", envir=.PlotPhyloEnv)

# Add time line

# Extract root age
root_age <- max(phytools::nodeHeights(Ponerinae_tree_old_calib))

# Define ticks
# ticks_labels <- seq(from = 0, to = 100, by = 20)
ticks_labels <- seq(from = 0, to = 120, by = 20)
axis(side = 1, pos = 0, at = (-1 * ticks_labels) + root_age, labels = ticks_labels, cex.axis = 1.5)
legend(x = root_age/2,
       y = 0 - 5, adj = 0,
       bty = "n", legend = "", title = "Time  [My]", title.cex = 1.5)

# Add a legend
legend(x = plot_info$x.lim[2] - 10,
       y = mean(plot_info$y.lim),
       # adj = c(0,0),
       # x = "topleft",
       legend = c("Absence", "Presence"),
       pch = 22, pt.bg = c("white","gray30"), pt.cex =  1.8,
       cex = 1.5, bty = "n")

dev.off()

## Check differences in mean rates

# Extract rates
Ponerinae_net_div_rates <- Ponerinae_BAMM_object_old_calib$meanTipLambda - Ponerinae_BAMM_object_old_calib$meanTipMu
names(Ponerinae_net_div_rates) <- Ponerinae_BAMM_object_old_calib$tip.label
# Reorder rates
Ponerinae_net_div_rates <- Ponerinae_net_div_rates[names(Ponerinae_cat_2lvl_tip_data)]
# Compute means
by(data = Ponerinae_net_div_rates, INDICES = Ponerinae_cat_2lvl_tip_data, FUN = summary)

## Run quick test for t = 0
quick_test <- BAMMtools::traitDependentBAMM(ephy = Ponerinae_BAMM_object_old_calib, traits = Ponerinae_cat_2lvl_tip_data, rate = 'net diversification', return.full = TRUE, method = 'kruskal', logrates = FALSE)
str(quick_test, 1)

## Produce densityMaps using stochastic character mapping based on an all-rates-different (ARD) Mk model
Ponerinae_cat_2lvl_data_old_calib <- prepare_trait_data(
  # tip_data = Ponerinae_data,
  # tip_data = fake_cat_tip_data,
  # tip_data = fake_cat_2lvl_tip_data,
  tip_data = Ponerinae_cat_2lvl_tip_data,
  phylo = Ponerinae_tree_old_calib,
  seed = 1234,
  trait_data_type = "categorical",
  colors_per_states = colors_per_states,
  evolutionary_models = "ARD",
  nb_simulations = 100,
  return_best_model_fit = TRUE,
  return_model_selection_df = TRUE,
  plot_map = FALSE)

# Plot result
plot_densityMaps_overlay(Ponerinae_cat_2lvl_data_old_calib$densityMaps,
                         colors_per_levels = colors_per_states)

# ## Export Ponerinae_cat_2lvl_data_old_calib in deepSTRAPP
# usethis::use_data(Ponerinae_cat_2lvl_data_old_calib, overwrite = TRUE)

# # Load directly output
# data(Ponerinae_cat_2lvl_data_old_calib, package = "deepSTRAPP")

## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
# nb_time_steps <- 5
time_step_duration <- 5
time_range <- c(0, 40)

## Run deepSTRAPP on net diversification rates across time-steps.
Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40 <- run_deepSTRAPP_over_time(
  densityMaps = Ponerinae_cat_2lvl_data_old_calib$densityMaps,
  nb_simulations = 100,
  ace = Ponerinae_cat_2lvl_data_old_calib$ace,
  # tip_data = Ponerinae_data,
  # tip_data = fake_cat_tip_data,
  # tip_data = fake_cat_2lvl_tip_data,
  tip_data = Ponerinae_cat_2lvl_tip_data,
  trait_data_type = "categorical",
  BAMM_object = Ponerinae_BAMM_object_old_calib,
  time_range = time_range,
  time_step_duration = time_step_duration,
  uncertainty_strategy = "paired",
  rate_type = "net_diversification",
  seed = 1234,
  posthoc_pairwise_tests = TRUE,
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40, max.level = 1)

# Display test summaries
# Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
# showing the evolution of the test results across time.
Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40$pvalues_summary_df

# Plot p-values of overall Kruskal-Wallis test across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40,
  alpha = 0.05)

# Plot rates through time
plot_rates_through_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40,
  colors_per_levels = colors_per_states,
  plot_CI = TRUE)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40, overwrite = TRUE)


### 17/ Generate Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 ####

## Load trait evolution data
data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")

## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
# nb_time_steps <- 5
time_step_duration <- 5
time_range <- c(0, 40)

## Run deepSTRAPP on net diversification rates across time-steps.
Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- run_deepSTRAPP_over_time(
   densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
   nb_simulations = 100,
   ace = Ponerinae_cat_3lvl_data_old_calib$ace,
   tip_data = Ponerinae_cat_3lvl_tip_data,
   trait_data_type = "categorical",
   BAMM_object = Ponerinae_BAMM_object_old_calib,
   time_range = time_range,
   time_step_duration = time_step_duration,
   uncertainty_strategy = "paired",
   rate_type = "net_diversification",
   seed = 1234,
   alpha = 0.10, # Select a generous level of significance for the sake of the example
   posthoc_pairwise_tests = TRUE,
   return_perm_data = TRUE,
   extract_trait_data_melted_df = TRUE,
   extract_diversification_data_melted_df = TRUE,
   return_STRAPP_results = TRUE,
   return_updated_Maps = TRUE,
   return_updated_BAMM_objects = TRUE,
   verbose = TRUE,
   verbose_extended = TRUE)

## Explore output
str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, max.level = 1)

# Display test summaries
# Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
# showing the evolution of the test results across time.
Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40$pvalues_summary_df
# Results for posthoc pairwise Dunn's tests over time-steps
Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40$pvalues_summary_df_for_posthoc_pairwise_tests

# Plot p-values of overall Kruskal-Wallis test across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
  alpha = 0.10)

# Plot rates through time
plot_rates_through_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
  colors_per_levels = colors_per_states,
  plot_CI = TRUE)

# Plot p-values of post hoc pairwise Dunn's tests between pairs of tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
  plot_posthoc_tests = TRUE)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, overwrite = TRUE)



### 18/ Generate Ponerinae_deepSTRAPP_biogeo_old_calib_0_40 ####

## deepSTRAPP output for Ponerinae over time from 0 to 40My (steps = 10 My)
# Tests for New World vs. Old World ranges
# Using the old ill-calibrated phylogeny

## Load data

# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
# Load range data
data(Ponerinae_binary_range_table, package = "deepSTRAPP")
# Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")

## Prepare range data for Old World vs. New World

# No overlap in ranges
table(Ponerinae_binary_range_table$Old_World, Ponerinae_binary_range_table$New_World)

Ponerinae_ON_tip_data <- stats::setNames(object = Ponerinae_binary_range_table$Old_World,
                                     nm = Ponerinae_binary_range_table$Taxa)
Ponerinae_ON_tip_data <- as.character(Ponerinae_ON_tip_data)
Ponerinae_ON_tip_data[Ponerinae_ON_tip_data == "TRUE"] <- "O" # O = Old World
Ponerinae_ON_tip_data[Ponerinae_ON_tip_data == "FALSE"] <- "N" # N = New World
names(Ponerinae_ON_tip_data) <- Ponerinae_binary_range_table$Taxa
table(Ponerinae_ON_tip_data)

colors_per_ranges <- c("mediumpurple2", "peachpuff2")
names(colors_per_ranges) <- c("N", "O")

# ## Run evolutionary models for biogeographic inferences
# Ponerinae_biogeo_data_old_calib <- prepare_trait_data(
#   tip_data = Ponerinae_ON_tip_data,
#   trait_data_type = "biogeographic",
#   phylo = Ponerinae_tree_old_calib,
#   evolutionary_models = "DEC+J", # Default = "DEC" for biogeographic
#   prefix_for_files = "Ponerinae",
#   max_range_size = 2,
#   # split_multi_area_ranges = TRUE, # Set to TRUE to display the two outputs
#   nb_simulations = 100, # Reduce to save time (Default = '1000')
#   colors_per_levels = colors_per_ranges,
#   return_model_selection_df = TRUE,
#   verbose = TRUE)

## Load directly range evolution data
data(Ponerinae_biogeo_data_old_calib, package = "deepSTRAPP")

## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
time_range <- c(0, 40)
time_step_duration <- 5

## Run deepSTRAPP on net diversification rates for time-steps = 0 to 40 Mya.
Ponerinae_deepSTRAPP_biogeo_old_calib_0_40 <- run_deepSTRAPP_over_time(
   densityMaps = Ponerinae_biogeo_data_old_calib$densityMaps,
   nb_simulations = 100,
   ace = Ponerinae_biogeo_data_old_calib$ace,
   tip_data = Ponerinae_ON_tip_data,
   trait_data_type = "biogeographic",
   BAMM_object = Ponerinae_BAMM_object_old_calib,
   time_range = time_range,
   time_step_duration = time_step_duration,
   seed = 1234, # Set seed for reproducibility
   alpha = 0.10, # Select a generous level of significance for the sake of the example
   uncertainty_strategy = "paired",
   rate_type = "net_diversification",
   return_perm_data = TRUE,
   extract_trait_data_melted_df = TRUE,
   extract_diversification_data_melted_df = TRUE,
   return_STRAPP_results = TRUE,
   return_updated_Maps = TRUE,
   return_updated_BAMM_objects = TRUE,
   verbose = TRUE,
   verbose_extended = TRUE)

## Explore output
str(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40, max.level = 1)

# Display test summaries
# Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
# showing the evolution of the test results across time.
Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$pvalues_summary_df

# Access bioregeographic range data in a melted data.frame
head(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$trait_data_df_over_time)

# Access the diversification data in a melted data.frame
head(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$diversification_data_df_over_time)

# Access details of deepSTRAPP results
str(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$STRAPP_results, max.level = 2)

# Plot updated densityMaps for time step n°2 = 10My
densityMaps_2 <- Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$updated_Maps_over_time[[2]]
plot_densityMaps_overlay(densityMaps_2$densityMaps)

# Plot diversification rates on updated phylogeny for time step n°2 = 10 My
plot_BAMM_rates(BAMM_object = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$updated_BAMM_objects_over_time[[2]], legend = TRUE, labels = FALSE,
                colorbreaks = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$updated_BAMM_objects_over_time[[2]]$initial_colorbreaks$net_diversification)

# Plot mapped trait and diversification rates on updated phylogeny for time step n°2 = 10 My
BAMM_object_10 <- Ponerinae_deepSTRAPP_biogeo_old_calib_0_40$updated_BAMM_objects_over_time[[2]]
plot_traits_vs_rates_on_phylogeny_for_focal_time(
   deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
   focal_time = 10,
   legend = TRUE, labels = FALSE,
   colorbreaks = BAMM_object_10$initial_colorbreaks$net_diversification)

# Plot p-values of Mann-Whitney-Wilcoxon tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
  alpha = 0.10)

# Plot histogram of Mann-Whitney-Wilcoxon test stats for time step n°2 = 10My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
  focal_time = 10)

# Plot histograms of STRAPP overall test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
  display_plots = TRUE,
  plot_posthoc_tests = FALSE)

# Plot rates through time
plot_rates_through_time(
  deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
  colors_per_levels = colors_per_levels,
  plot_CI = TRUE)

# Export in deepSTRAPP
usethis::use_data(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40, overwrite = TRUE)












