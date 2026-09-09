# Compare trait evolutionary model fits with AICc and Akaike's weights

Compare evolutionary models fit with
[`geiger::fitContinuous()`](https://rdrr.io/pkg/geiger/man/fitContinuous.html)
or
[`geiger::fitDiscrete()`](https://rdrr.io/pkg/geiger/man/fitDiscrete.html)
using AICc and Akaike's weights. Generate a data.frame summarizing
information. Identify the best model and extract its results.

## Usage

``` r
select_best_trait_model_from_geiger(list_model_fits)
```

## Arguments

- list_model_fits:

  Named list with the results of a model fit with
  [`geiger::fitContinuous()`](https://rdrr.io/pkg/geiger/man/fitContinuous.html)
  or
  [`geiger::fitDiscrete()`](https://rdrr.io/pkg/geiger/man/fitDiscrete.html)
  in each element.

## Value

The function returns a list with three elements.

- `$model_comparison_df` Data.frame summarizing information to compare
  model fits. It includes the model name (`$model`), the log-likelihood
  (`$logLik`), the number of free-parameters (`$k`), the AIC (`$AIC`),
  the corrected AIC (`$AICc`), the delta to the best/lowest AICc
  (`$delta_AICc`), the Akaike weights (`$Akaike_weights`), and their
  rank based on AICc (`$rank`).

- `$best_model_name` Character string. Name of the best model.

- `$best_model_fit` List containing the output of
  [`geiger::fitContinuous()`](https://rdrr.io/pkg/geiger/man/fitContinuous.html)
  or
  [`geiger::fitDiscrete()`](https://rdrr.io/pkg/geiger/man/fitDiscrete.html)
  for the model with the best fit.

## See also

[`geiger::fitContinuous()`](https://rdrr.io/pkg/geiger/man/fitContinuous.html)
[`geiger::fitDiscrete()`](https://rdrr.io/pkg/geiger/man/fitDiscrete.html)

## Author

Maël Doré

## Examples

``` r

# ----- Example 1: Continuous data ----- #

# Load phylogeny and tip data
library(phytools)
data(eel.tree)
data(eel.data)

# Extract body size
eel_data <- stats::setNames(eel.data$Max_TL_cm,
                            rownames(eel.data))

 # (May take several minutes to run)
# Fit BM model
BM_fit <- geiger::fitContinuous(phy = eel.tree, dat = eel_data, model = "BM")
# Fit EB model
EB_fit <- geiger::fitContinuous(phy = eel.tree, dat = eel_data, model = "EB")
#> Warning: 
#> Parameter estimates appear at bounds:
#>  a
# Fit OU model
OU_fit <- geiger::fitContinuous(phy = eel.tree, dat = eel_data, model = "OU")

# Store models
list_model_fits <- list(BM = BM_fit, EB = EB_fit, OU = OU_fit)

# Compare models
model_comparison_output <- select_best_trait_model_from_geiger(list_model_fits = list_model_fits)

# Explore output
str(model_comparison_output, max.level = 2)
#> List of 3
#>  $ models_comparison_df:'data.frame':    3 obs. of  8 variables:
#>   ..$ model         : chr [1:3] "BM" "EB" "OU"
#>   ..$ logL          : num [1:3] -339 -339 -329
#>   ..$ k             : num [1:3] 2 3 3
#>   ..$ AIC           : num [1:3] 682 684 665
#>   ..$ AICc          : num [1:3] 682 684 665
#>   ..$ delta_AICc    : num [1:3] 17.1 19.4 0
#>   ..$ Akaike_weights: 'aic.w' num [1:3] 0 0 100
#>   ..$ rank          : int [1:3] 2 3 1
#>  $ best_model_name     : chr "OU"
#>  $ best_model_fit      :List of 4
#>   ..$ lik:function (pars, ...)  
#>   .. ..- attr(*, "argn")= chr [1:2] "alpha" "sigsq"
#>   .. ..- attr(*, "cache")=List of 26
#>   .. ..- attr(*, "class")= chr [1:2] "bm" "function"
#>   .. ..- attr(*, "model")= chr "OU"
#>   ..$ bnd:'data.frame':  2 obs. of  2 variables:
#>   ..$ res: num [1:100, 1:4] 4.34e-02 7.12e-218 1.74e-217 4.34e-02 4.34e-02 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ opt:List of 8
#>   ..- attr(*, "class")= chr [1:2] "gfit" "list"

# Print comparison
print(model_comparison_output$models_comparison_df)
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> BM    BM -338.9352 2 681.8704 682.0773   17.14864              0    2
#> EB    EB -338.9355 3 683.8710 684.2921   19.36339              0    3
#> OU    OU -329.2538 3 664.5076 664.9287    0.00000            100    1

# Print best model fit
print(model_comparison_output$best_model_fit) 
#> GEIGER-fitted comparative model of continuous data
#>  fitted ‘OU’ model parameters:
#>  alpha = 0.043393
#>  sigsq = 278.055245
#>  z0 = 92.068580
#> 
#>  model summary:
#>  log-likelihood = -329.253814
#>  AIC = 664.507628
#>  AICc = 664.928681
#>  free parameters = 3
#> 
#> Convergence diagnostics:
#>  optimization iterations = 100
#>  failed iterations = 0
#>  number of iterations with same best fit = 13
#>  frequency of best fit = 0.130
#> 
#>  object summary:
#>  'lik' -- likelihood function
#>  'bnd' -- bounds for likelihood search
#>  'res' -- optimization iteration summary
#>  'opt' -- maximum likelihood parameter estimates

# ----- Example 2: Categorical data ----- #

# Load phylogeny and tip data
library(phytools)
data(eel.tree)
data(eel.data)

# Extract feeding mode data
eel_data <- stats::setNames(eel.data$feed_mode, rownames(eel.data))
table(eel_data)
#> eel_data
#>    bite suction 
#>      31      30 

 # (May take several minutes to run)
# Fit ER model
ER_fit <- geiger::fitDiscrete(phy = eel.tree, dat = eel_data, model = "ER")
# Fit ARD model
ARD_fit <- geiger::fitDiscrete(phy = eel.tree, dat = eel_data, model = "ARD")

# Store models
list_model_fits <- list(ER = ER_fit, ARD = ARD_fit)

# Compare models
model_comparison_output <- select_best_trait_model_from_geiger(list_model_fits = list_model_fits)

# Explore output
str(model_comparison_output, max.level = 2)
#> List of 3
#>  $ models_comparison_df:'data.frame':    2 obs. of  8 variables:
#>   ..$ model         : chr [1:2] "ER" "ARD"
#>   ..$ logL          : num [1:2] -37 -37
#>   ..$ k             : int [1:2] 1 2
#>   ..$ AIC           : num [1:2] 76 78
#>   ..$ AICc          : num [1:2] 76.1 78.2
#>   ..$ delta_AICc    : num [1:2] 0 2.1
#>   ..$ Akaike_weights: 'aic.w' num [1:2] 74.1 25.9
#>   ..$ rank          : int [1:2] 1 2
#>  $ best_model_name     : chr "ER"
#>  $ best_model_fit      :List of 4
#>   ..$ lik:function (pars, ...)  
#>   .. ..- attr(*, "class")= chr [1:4] "constrained" "mkn" "dtlik" "function"
#>   .. ..- attr(*, "argn")= chr "q12"
#>   .. ..- attr(*, "formulae")=List of 1
#>   .. ..- attr(*, "func")=function (pars, ...)  
#>   .. .. ..- attr(*, "class")= chr [1:3] "mkn" "dtlik" "function"
#>   .. .. ..- attr(*, "argn")= chr [1:2] "q12" "q21"
#>   .. .. ..- attr(*, "levels")= chr [1:2] "bite" "suction"
#>   .. ..- attr(*, "levels")= chr [1:2] "bite" "suction"
#>   .. ..- attr(*, "trns")= logi TRUE
#>   .. ..- attr(*, "transform")= chr "none"
#>   ..$ bnd:'data.frame':  2 obs. of  2 variables:
#>   ..$ res: num [1:100, 1:3] 0.0157 0.0157 0.0157 0.0157 0.0157 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ opt:List of 7
#>   ..- attr(*, "class")= chr [1:2] "gfit" "list"

# Print comparison
print(model_comparison_output$models_comparison_df)
#>     model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER     ER -37.02368 1 76.04736 76.11516   0.000000           74.1    1
#> ARD   ARD -37.00289 2 78.00577 78.21267   2.097512           25.9    2

# Print best model fit
print(model_comparison_output$best_model_fit) 
#> GEIGER-fitted comparative model of discrete data
#>  fitted Q matrix:
#>                    bite     suction
#>     bite    -0.01569437  0.01569437
#>     suction  0.01569437 -0.01569437
#> 
#>  model summary:
#>  log-likelihood = -37.023680
#>  AIC = 76.047361
#>  AICc = 76.115158
#>  free parameters = 1
#> 
#> Convergence diagnostics:
#>  optimization iterations = 100
#>  failed iterations = 0
#>  number of iterations with same best fit = 100
#>  frequency of best fit = 1.000
#> 
#>  object summary:
#>  'lik' -- likelihood function
#>  'bnd' -- bounds for likelihood search
#>  'res' -- optimization iteration summary
#>  'opt' -- maximum likelihood parameter estimates
```
