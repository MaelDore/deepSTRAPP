
##### Script to convert package between development versions and CRAN releases #####

### Development versions

# Contain all datasets, including deepSTRAPP and BAMM outputs for all four types of data, and eel_biogeo_data with BioGeoBEARS classes
# Include doc for all datasets
# Produce vignette visual outputs using code evaluation
# Run all examples including loading of deepSTRAPP outputs
# Add BioGeoBEARS/contsimmap as Imports with link through Remotes

### CRAN Releases ###

# Needs to reduce size to comply with CRAN policies
# Do not contain datasets of deepSTRAPP and BAMM output
# Contains a modified version of eel_biogeo_data without BioGeoBEARS classes
# Do not include doc for deepSTRAPP datasets
# Produce vignette visual outputs based on pre-rendered PNG images
# Do not run examples including loading of deepSTRAPP outputs
# Add BioGeoBEARS/contsimmap as Suggests with link through Additional_repositories

### 1/ From development version to CRAN release ####

## 1.1/ Remove light deepSTRAPP datasets from data/
unlink(x = "./data/Ponerinae_BAMM_object.rda", force = T)
unlink(x = "./data/Ponerinae_BAMM_object_old_calib.rda", force = T)
unlink(x = "./data/Ponerinae_BAMM_object_10My.rda", force = T)

## 1.2/ Remove heavy deepSTRAPP datasets from inst/extdata/
unlink(x = "./inst/extdata/Ponerinae_cont_data_old_calib.rds", force = T)
unlink(x = "./inst/extdata/Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", force = T)
unlink(x = "./inst/extdata/Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40.rds", force = T)
unlink(x = "./inst/extdata/Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", force = T)
unlink(x = "./inst/extdata/Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rds", force = T)

## 1.3/ Replace eel_biogeo_data
unlink(x = "./data/eel_biogeo_data.rda", force = T)
file.copy(from = "./for_CRAN/eel_biogeo_data_for_CRAN.rda", to = "./data/eel_biogeo_data.rda", overwrite = T)

## 1.4/ Replace datasets documentation
unlink(x = "./R/datasets_doc.R", force = T)
file.copy(from = "./for_CRAN/datasets_doc_for_CRAN.R", to = "./R/datasets_doc_for_CRAN.R", overwrite = T)

## 1.5/ Add pre-rendered visual outputs to the vignette folder
prerendered_vignette_outputs_path <- list.files(path = "./for_CRAN/Figures_for_vignettes_Low_resolution/", pattern = ".PNG")
file.copy(from = paste0("./for_CRAN/Figures_for_vignettes_Low_resolution/",prerendered_vignette_outputs_path), to = paste0("./vignettes/figures/",prerendered_vignette_outputs_path), overwrite = T)

## 1.6/ Replace DESCRIPTION with BioGeoBEARS as Suggests with link through Additional_repositories
unlink(x = "./DESCRIPTION", force = T)
file.copy(from = "./for_CRAN/DESCRIPTION_for_CRAN", to = "./DESCRIPTION", overwrite = T)

## 1.7/ Increment version number
usethis::use_version(which = "minor") # Increase as minor update  1.X.0 when introducing new features that are compatible with previous usages
usethis::use_version(which = "major") # Increase as major update  X.0.0 when introducing changes that are not backward compatible

## 1.8/ Rerun check once done
# Devtools check
devtools::check()
# Local check as in CRAN check tool
devtools::check(cran = TRUE)


### 2/ From CRAN release to development version ####

## 2.1/ Add light deepSTRAPP datasets to data
file.copy(from = "./for_CRAN/Ponerinae_BAMM_object.rda", to = "./data/Ponerinae_BAMM_object.rda", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_BAMM_object_old_calib.rda", to = "./data/Ponerinae_BAMM_object_old_calib.rda", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_BAMM_object_10My.rda", to = "./data/Ponerinae_BAMM_object_10My.rda", overwrite = T)

## 2.2/ Add heavy deepSTRAPP datasets to inst/extdata
file.copy(from = "./for_CRAN/Ponerinae_cont_data_old_calib.rds", to = "./inst/extdata/Ponerinae_cont_data_old_calib.rds", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", to = "./inst/extdata/Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40.rds", to = "./inst/extdata/Ponerinae_deepSTRAPP_cat_2lvl_old_calib_0_40.rds", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", to = "./inst/extdata/Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rds", to = "./inst/extdata/Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rds", overwrite = T)

## 2.3/ Replace eel_biogeo_data
unlink(x = "./data/eel_biogeo_data.rda", force = T)
file.copy(from = "./for_CRAN/eel_biogeo_data.rda", to = "./data/eel_biogeo_data.rda", overwrite = T)

## 2.4/ Replace datasets documentation
unlink(x = "./R/datasets_doc_for_CRAN.R", force = T)
file.copy(from = "./for_CRAN/datasets_doc.R", to = "./R/datasets_doc.R", overwrite = T)

## 2.5/ Remove pre-rendered visual outputs to the vignette folder
prerendered_vignette_outputs_path <- list.files(path = "./for_CRAN/Figures_for_vignettes_Low_resolution/", pattern = ".PNG")
unlink(x = paste0("./vignettes/figures/",prerendered_vignette_outputs_path), force = T)

## 2.6/ Replace DESCRIPTION with BioGeoBEARS as Imports with link through Remotes
unlink(x = "./DESCRIPTION", force = T)
file.copy(from = "./for_CRAN/DESCRIPTION", to = "./DESCRIPTION", overwrite = T)

## 2.7/ Turn to the next development version
usethis::use_dev_version()

## 2.8/ Rerun check once done
devtools::check()

