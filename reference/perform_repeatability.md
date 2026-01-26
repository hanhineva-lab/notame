# Compute repeatability measures

Computes repeatability for each feature with the following formula:
\$\$\frac{\sigma^2\_{between}}{\sigma^2\_{between} +
\sigma^2\_{within}}\$\$ The repeatability ranges from 0 to 1. Higher
repeatability depicts less variation between batches.

## Usage

``` r
perform_repeatability(object, group, assay.type = NULL)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- group:

  column name of pheno data giving the group labels

- assay.type:

  character, assay to be used in case of multiple assays

## Value

A data frame with one row per feature with the repeatability measure.

## Examples

``` r
data(toy_notame_set)
# Batch correction
replicates <- list(which(toy_notame_set$QC == "QC"))
batch_corrected <- ruvs_qc(toy_notame_set, replicates = replicates)
# Evaluate batch correction
rep_orig <- perform_repeatability(toy_notame_set, group = "Group")
mean(rep_orig$Repeatability, na.rm = TRUE)
#> [1] 0.3658633
rep_corr <- perform_repeatability(batch_corrected, group = "Group")
mean(rep_corr$Repeatability, na.rm = TRUE)
#> [1] 0.4051463
```
