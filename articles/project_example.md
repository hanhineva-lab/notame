# Project example

In this project example, core functionality of notame, notameViz and
notameStats is demonstrated in unison with the associated protocol
article (Klåvus et al. 2020).

## Project setup

### Set up path and logging

Let’s set up a path and start logging our project. Many functions in the
notame packages will log information when they have finished. This helps
in monitoring the analysis process in case something goes wrong and in
reviewing the results later. The log will also print to console,
although we’ll exclude console output for the rest of the document for
brevity.

``` r

library(notame)
library(notameViz)
library(notameStats)
library(dplyr)
ppath <- tempdir()
init_log(log_file = file.path(ppath, "log.txt"))
```

### Read data

The data is mock data with two balanced groups and two balanced time
points. The time points also correspond to two batches. There are 50
samples with ten regularly interspersed QC samples and 80 features
divided into four analytical modes.

``` r

data <- import_from_excel(
  file = system.file("extdata", "toy_notame_set.xlsx", package = "notame"), 
  sheet = 1, corner_row = 11, corner_column = "H", 
  split_by = c("Column", "Ion_mode"))
  
names(assays(data)) <- "abundances"
```

The function
[`import_from_excel()`](https://hanhineva-lab.github.io/notame/reference/import_from_excel.md)
returns a SummarizedExperiment object with three main parts and
corresponding accessors for streamlined data wrangling.

- `assay`: feature abundances across the samples\
- `colData`: sample information\
- `rowData`: feature information

The example data contains four analytical modes, so we separate them
using
[`fix_object()`](https://hanhineva-lab.github.io/notame/reference/fix_object.md),
which also conveniently checks the object for compatibility with all of
notame.

``` r

modes <- fix_object(data, split_data = TRUE, assay.type = "abundances")
```

Now we have a list of objects, one per mode (LC column x ionization
mode) or however many unique values were specified in
`rowData(data)$Split`.

``` r

names(modes)
## [1] "HILIC_neg" "HILIC_pos" "RP_neg"    "RP_pos"
sapply(modes, class)
##              HILIC_neg              HILIC_pos                 RP_neg 
## "SummarizedExperiment" "SummarizedExperiment" "SummarizedExperiment" 
##                 RP_pos 
## "SummarizedExperiment"
```

## Tabular data preprocessing

### Preprocessing by mode

Tabular data preprocessing is performed to complete the dataset by way
of reducing unwanted variation and data preparation dependent on
downstream methods. Steps performed separately for each mode include
marking missing values as NA, flagging features with low detection rate
(Broadhurst et al. 2018), drift correction (Kirwan et al. 2013), and
flagging of low-quality features (Broadhurst et al. 2018).
Visualizations are drawn to monitor the process and explore the data
globally. The
[`save_QC_plots()`](https://rdrr.io/pkg/notameViz/man/save_QC_plots.html)
wrapper does not include feature-wise plots as it may not be feasible to
inspect them at this stage of the analysis. For example, drift
correction visualizations are best drawn for a subset of interesting
features after feature selection.

``` r

# Initialize empty list for processed objects
processed <- list()
for (i in seq_along(modes)) {
  name <- names(modes)[i]
  mode <- modes[[i]]
  # Set all zero abundances to NA
  mode <- mark_nas(mode, value = 0)
  # Flag features with low detection rate
  mode <- flag_detection(mode, group = "Group")
  # Visualize data before drift correction
  save_QC_plots(mode, prefix = paste0(ppath, "figures/", name, "_ORIG"),
                 perplexity = 5, group = "Group", time = "Time", 
                 id = "Subject_ID", color = "Group")
  # Correct drift
  corrected <- correct_drift(mode)
  # Visualize data after drift correction
  save_QC_plots(corrected, prefix = paste0(ppath, "figures/", name, "_DRIFT"),
                perplexity = 5, group = "Group", time = "Time", 
                id = "Subject_ID", color = "Group")
  # Flag low-quality features
  corrected <- corrected %>% 
    assess_quality() %>% 
    flag_quality()
  # Visualize data after removal of low-quality features
  save_QC_plots(corrected, prefix = paste0(ppath, "figures/", name, "_CLEAN"),
                perplexity = 5, group = "Group", time = "Time",
                id = "Subject_ID", color = "Group")
  # Save result of iteration
  processed[[i]] <- corrected
}
```

Feature-wise flagging information, quality metrics, and brief drift
correction notes are included in the feature information after the above
steps.

``` r

rowData(processed[[1]])$DC_note
##  [1] "Drift_corrected" "Drift_corrected" "Drift_corrected" "Drift_corrected"
##  [5] "Drift_corrected" "Drift_corrected" "Drift_corrected" "Drift_corrected"
##  [9] "Drift_corrected" "Drift_corrected" "Drift_corrected" "Drift_corrected"
## [13] "Drift_corrected" "Drift_corrected" "Drift_corrected" "Drift_corrected"
## [17] "Drift_corrected" "Drift_corrected" "Drift_corrected" "Drift_corrected"
```

### Preprocessing for the complete dataset

Next, it is time to merge the modes together and visualize the complete
dataset.

``` r

merged <- merge_notame_sets(processed)
save_QC_plots(merged, prefix = paste0(ppath, "figures/_FULL"),
              group = "Group", time = "Time",
              id = "Subject_ID", color = "Group")
```

Then, QC samples can (and should) be removed, as they are no longer
needed, and the dataset is visualized anew.

``` r

merged_no_qc <- drop_qcs(merged)
save_QC_plots(merged_no_qc, prefix = paste0(ppath, "figures/FULL_NO_QC"),
              group = "Group", time = "Time",
              id = "Subject_ID", color = "Group")
```

If there are any missing values in the data, they need to be imputed.
Let’s use random forest imputation to impute these values and visualize
the dataset one last time before statistical analysis. Seed number
should be set before random forest imputation to guarantee
reproducibility of results.

``` r

set.seed(2024)
imputed <- impute_rf(merged_no_qc)
save_QC_plots(imputed, prefix = paste0(ppath, "figures/FULL_IMPUTED"),
              group = "Group", time = "Time",
              id = "Subject_ID", color = "Group")
```

By default, the imputation procedure only operates on good quality
features, i.e. those that have not been flagged. To use flagged features
in statistical analysis, they should be imputed as well. This can be
achieved through a second round of imputation, now with all features
included. This two-step imputation makes sure that low-quality features
don’t affect the imputation of quality features. Imputation could also
be performed separately for each mode to reduce execution time,
especially if the amounts of features still allows for good imputation
results.

``` r

base <- impute_rf(imputed, all_features = TRUE)
```

It is a good idea to save the merged and processed data, so
experimenting with different statistical analyses becomes easier.

``` r

save(base, file = paste0(ppath, "full_data.RData"))
```

## Feature selection

### Univariate analysis

First, we’ll use a linear model to evaluate the features in terms of the
probability of obtaining the observed abundances given the null
hypothesis, namely that there is no difference in feature abundance
across study groups.

``` r

assay(base, "log") <- log(assay(base))
lm_results <- perform_lm(base, assay.type = "log", 
  formula_char = "Feature ~ Group")
base <- join_rowData(base, lm_results)
```

Next, univariate and supervised rankings of features could be combined
in a final ranking. Such a final ranking would combine the qualities of
univariate analysis and supervised learning. Next, carefully selected
features would undergo labour-intensive scrutiny for biological meaning,
for example by identification and pathway analysis.

That concludes our project example. The last thing to do is to write the
results to an Excel file and finish the logging. Thanks!

``` r

write_to_excel(base, file = paste0(ppath, "results.xlsx"))
finish_log()
```

## Session information

    ## R Under development (unstable) (2026-02-22 r89452)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] dplyr_1.2.0                 notameStats_1.1.1          
    ##  [3] notameViz_1.1.4             notame_1.1.4               
    ##  [5] SummarizedExperiment_1.41.1 Biobase_2.71.0             
    ##  [7] GenomicRanges_1.63.1        Seqinfo_1.1.0              
    ##  [9] IRanges_2.45.0              S4Vectors_0.49.0           
    ## [11] BiocGenerics_0.57.0         generics_0.1.4             
    ## [13] MatrixGenerics_1.23.0       matrixStats_1.5.0          
    ## [15] ggplot2_4.0.2               BiocStyle_2.39.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1     viridisLite_0.4.3    farver_2.1.2        
    ##  [4] S7_0.2.1             fastmap_1.2.0        digest_0.6.39       
    ##  [7] lifecycle_1.0.5      statmod_1.5.1        magrittr_2.0.4      
    ## [10] compiler_4.6.0       rngtools_1.5.2       rlang_1.1.7         
    ## [13] sass_0.4.10          tools_4.6.0          yaml_2.3.12         
    ## [16] knitr_1.51           lambda.r_1.2.4       doRNG_1.8.6.3       
    ## [19] S4Arrays_1.11.1      labeling_0.4.3       htmlwidgets_1.6.4   
    ## [22] DelayedArray_0.37.0  RColorBrewer_1.1-3   abind_1.4-8         
    ## [25] BiocParallel_1.45.0  Rtsne_0.17           withr_3.0.2         
    ## [28] purrr_1.2.1          itertools_0.1-3      desc_1.4.3          
    ## [31] grid_4.6.0           iterators_1.0.14     scales_1.4.0        
    ## [34] MASS_7.3-65          cli_3.6.5            rmarkdown_2.30      
    ## [37] ragg_1.5.0           otel_0.2.0           cachem_1.1.0        
    ## [40] parallel_4.6.0       formatR_1.14         BiocManager_1.30.27 
    ## [43] XVector_0.51.0       vctrs_0.7.1          Matrix_1.7-4        
    ## [46] jsonlite_2.0.0       bookdown_0.46        systemfonts_1.3.1   
    ## [49] foreach_1.5.2        limma_3.67.0         jquerylib_0.1.4     
    ## [52] tidyr_1.3.2          ggdendro_0.2.0       missForest_1.6.1    
    ## [55] glue_1.8.0           pkgdown_2.2.0        codetools_0.2-20    
    ## [58] cowplot_1.2.0        stringi_1.8.7        gtable_0.3.6        
    ## [61] futile.logger_1.4.9  tibble_3.3.1         pillar_1.11.1       
    ## [64] pcaMethods_2.3.0     htmltools_0.5.9      randomForest_4.7-1.2
    ## [67] R6_2.6.1             Rdpack_2.6.6         textshaping_1.0.4   
    ## [70] evaluate_1.0.5       lattice_0.22-9       backports_1.5.0     
    ## [73] rbibutils_2.4.1      futile.options_1.0.1 openxlsx_4.2.8.1    
    ## [76] broom_1.0.12         bslib_0.10.0         Rcpp_1.1.1          
    ## [79] zip_2.3.3            SparseArray_1.11.10  ranger_0.18.0       
    ## [82] xfun_0.56            fs_1.6.6             pkgconfig_2.0.3

## References

Broadhurst, David, Royston Goodacre, Stacey N Reinke, et al. 2018.
“Guidelines and Considerations for the Use of System Suitability and
Quality Control Samples in Mass Spectrometry Assays Applied in
Untargeted Clinical Metabolomic Studies.” *Metabolomics* 14: 1–17.

Kirwan, JA, DI Broadhurst, RL Davidson, and MR Viant. 2013.
“Characterising and Correcting Batch Variation in an Automated Direct
Infusion Mass Spectrometry (DIMS) Metabolomics Workflow.” *Analytical
and Bioanalytical Chemistry* 405: 5147–57.

Klåvus, Anton, Marietta Kokla, Stefania Noerman, et al. 2020. “‘Notame’:
Workflow for Non-Targeted LC–MS Metabolic Profiling.” *Metabolites* 10
(4): 135.
