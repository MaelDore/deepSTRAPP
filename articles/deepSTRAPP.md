# deepSTRAPP: All tutorials

## Main tutorial

### Quick-to-run example

A **simple use-case** that shows how deepSTRAPP can be used to **test
for differences in diversification rates between two trait states along
evolutionary times** is available
[here](https://maeldore.github.io/deepSTRAPP/articles/main_tutorial.html)
and within R:
[`vignette("main_tutorial")`](https://maeldore.github.io/deepSTRAPP/articles/main_tutorial.md).

This tutorial presents the main functions in a typical **deepSTRAPP
workflow**.  
For more advanced uses, please refer to the vignettes/tutorials below.

## Advanced uses / tutorials

This vignette points to **tutorials** detailing how to use the
\[deepSTRAPP\] package beyond the **simple use-case** presented in the
README file and also available here:
[`vignette("main_tutorial")`](https://maeldore.github.io/deepSTRAPP/articles/main_tutorial.md).

The following tutorials present more **advanced usages** of deepSTRAPP.
They provide explanations on available arguments and interpretations of
results of deepSTRAPP across multiple types of data.

### **1/ Full deepSTRAPP workflows on different types of data**

    ●   [1.1/ Full deepSTRAPP workflow for **continuous** trait
data](https://maeldore.github.io/deepSTRAPP/articles/deepSTRAPP_continuous_data.html):
[`vignette("deepSTRAPP_continuous_data")`](https://maeldore.github.io/deepSTRAPP/articles/deepSTRAPP_continuous_data.md).

    ●   [1.2/ Full deepSTRAPP workflow for **categorical** trait data
with
3-levels](https://maeldore.github.io/deepSTRAPP/articles/deepSTRAPP_categorical_data.html):
[`vignette("deepSTRAPP_categorical_data")`](https://maeldore.github.io/deepSTRAPP/articles/deepSTRAPP_categorical_data.md).

    ●   [1.3/ Full deepSTRAPP workflow for **biogeographic** range
data](https://maeldore.github.io/deepSTRAPP/articles/deepSTRAPP_biogeographic_data.html):
[`vignette("deepSTRAPP_biogeographic_data")`](https://maeldore.github.io/deepSTRAPP/articles/deepSTRAPP_biogeographic_data.md).

### **2/ Explore options for trait evolution**

    ●   [2.1/ Model evolution of **continuous** trait
data](https://maeldore.github.io/deepSTRAPP/articles/model_continuous_trait_evolution.html):
[`vignette("model_continuous_trait_evolution")`](https://maeldore.github.io/deepSTRAPP/articles/model_continuous_trait_evolution.md).

    ●   [2.2/ Model evolution of **categorical** trait
data](https://maeldore.github.io/deepSTRAPP/articles/model_categorical_trait_evolution.html):
[`vignette("model_categorical_trait_evolution")`](https://maeldore.github.io/deepSTRAPP/articles/model_categorical_trait_evolution.md).

    ●   [2.3/ Model evolution of **biogeographic** range
data](https://maeldore.github.io/deepSTRAPP/articles/model_biogeographic_range_evolution.html):
[`vignette("model_biogeographic_range_evolution")`](https://maeldore.github.io/deepSTRAPP/articles/model_biogeographic_range_evolution.md).

### **3/ Explore options for BAMM**

    ●   [Model **diversification dynamics** with BAMM within
deepSTRAPP](https://maeldore.github.io/deepSTRAPP/articles/model_diversification_dynamics.html):
[`vignette("model_diversification_dynamics")`](https://maeldore.github.io/deepSTRAPP/articles/model_diversification_dynamics.md).

### **4/ Explore the STRAPP test options**

    ●   [Test different
hypotheses](https://maeldore.github.io/deepSTRAPP/articles/explore_STRAPP_test_types.html):
[`vignette("explore_STRAPP_test_types")`](https://maeldore.github.io/deepSTRAPP/articles/explore_STRAPP_test_types.md).

- Type of STRAPP tests: **two-tailed** vs. **one-tailed**.
- Continuous: “negative” or “positive” correlation.
- Binary with hypothesis: (A \> B) vs. (B \> A).
- Multinominal: Hypotheses for all post hoc tests.

### **5/ Plot rates through time (RTT)**

    ●   [Explore options for plotting diversification **rates through
time** in relation to trait
data](https://maeldore.github.io/deepSTRAPP/articles/plot_rates_through_time.html):
[`vignette("plot_rates_through_time")`](https://maeldore.github.io/deepSTRAPP/articles/plot_rates_through_time.md).

### **6/ Handle uncertainty**

    ●   [Handle **uncertainty** in trait and rate
estimates](https://maeldore.github.io/deepSTRAPP/articles/handle_uncertainty.html):
[`vignette("handle_uncertainty")`](https://maeldore.github.io/deepSTRAPP/articles/handle_uncertainty.md).

Explore the three strategies available: - ‘rates_only’: Only accounts
for diversification-rate uncertainty across BAMM posterior samples. -
‘paired’: Accounts for both diversification-rate and ancestral
trait/range reconstruction uncertainty by pairing BAMM posterior samples
with stochastic maps. - ‘full’: Accounts for both diversification-rate
and ancestral reconstruction uncertainty by evaluating every combination
of BAMM posterior samples and stochastic maps.  

### **7/ Import external analyses**

    ●   [Import **external
analyses**](https://maeldore.github.io/deepSTRAPP/articles/import_external_analyses.html):
[`vignette("import_external_analyses")`](https://maeldore.github.io/deepSTRAPP/articles/import_external_analyses.md).

Import and format results of external analyses of trait-evolution
histories and diversification dynamics, and make them ready-to-use as
inputs for a deepSTRAPP run.  

### **8/ Cut phylogenies**

    ●   [Cut different types of **(mapped) phylogenies** for a given
focal-time](https://maeldore.github.io/deepSTRAPP/articles/cut_phylogenies.html):
[`vignette("cut_phylogenies")`](https://maeldore.github.io/deepSTRAPP/articles/cut_phylogenies.md).

- time-calibrated phylogenies.
- contMap(s) for continuous traits.
- densityMap(s) for categorical and biogeographic traits.
- simmaps for categorical and biogeographic traits.
- BAMM_object for diversification dynamics.  
