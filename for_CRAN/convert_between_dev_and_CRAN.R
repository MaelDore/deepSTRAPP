
##### Script to convert package between development versions and CRAN releases #####

### Development versions

# Contain all datasets, including deepSTRAPP outputs for all four types of data, and eel_biogeo_data with BioGeoBEARS classes
# Include doc for all datasets
# Produce vignette visual outputs using code evaluation
# Run all examples including loading of deepSTRAPP outputs

### CRAN Releases ###

# Needs to reduce size to comply with CRAN policies
# Do not contain datasets of deepSTRAPP outputs (only for binary data as used in the Main tutorial)
# Contains a modified version of eel_biogeo_data without BioGeoBEARS classes
# Do not include doc for deepSTRAPP datasets (only for binary data as used in the Main tutorial)
# Produce vignette visual outputs based on pre-rendered PNG images
# Do not run examples including loading of deepSTRAPP outputs

### 1/ From development version to CRAN release ####

## 1.1/ Remove deepSTRAPP datasets from data
unlink(x = "./data/Ponerinae_deepSTRAPP_cont_old_calib_0_40.rda", force = T)
unlink(x = "./data/Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rda", force = T)
unlink(x = "./data/Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rda", force = T)

## 1.2/ Replace eel_biogeo_data
unlink(x = "./data/eel_biogeo_data.rda", force = T)
file.copy(from = "./for_CRAN/eel_biogeo_data_for_CRAN.rda", to = "./data/eel_biogeo_data.rda", overwrite = T)

## 1.3/ Replace datasets documentation
unlink(x = "./R/datasets_doc.R", force = T)
file.copy(from = "./for_CRAN/datasets_doc_for_CRAN.R", to = "./R/datasets_doc_for_CRAN.R", overwrite = T)

## 1.4/ Add pre-rendered visual outputs to the vignette folder
prerendered_vignette_outputs_path <- list.files(path = "./for_CRAN/", pattern = ".PNG")
file.copy(from = paste0("./for_CRAN/",prerendered_vignette_outputs_path), to = paste0("./vignettes/figures/",prerendered_vignette_outputs_path), overwrite = T)

## 1.5/ Rerun check once done
devtools::check()


### 2/ From CRAN release to development version ####

## 2.1/ Add deepSTRAPP datasets to data
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_cont_old_calib_0_40.rda", to = "./data/Ponerinae_deepSTRAPP_cont_old_calib_0_40.rda", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rda", to = "./data/Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rda", overwrite = T)
file.copy(from = "./for_CRAN/Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rda", to = "./data/Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rda", overwrite = T)

## 2.2/ Replace eel_biogeo_data
unlink(x = "./data/eel_biogeo_data.rda", force = T)
file.copy(from = "./for_CRAN/eel_biogeo_data.rda", to = "./data/eel_biogeo_data.rda", overwrite = T)

## 2.3/ Replace datasets documentation
unlink(x = "./R/datasets_doc_for_CRAN.R", force = T)
file.copy(from = "./for_CRAN/datasets_doc.R", to = "./R/datasets_doc.R", overwrite = T)

## 2.4/ Remove pre-rendered visual outputs to the vignette folder
prerendered_vignette_outputs_path <- list.files(path = "./for_CRAN/", pattern = ".PNG")
unlink(x = paste0("./vignettes/figures/",prerendered_vignette_outputs_path), force = T)

## 2.5/ Rerun check once done
devtools::check()

