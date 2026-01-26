# Remove Unwanted Variation (RUV) between batches

An interface for `link[RUVSeq]{RUVs}` method

## Usage

``` r
ruvs_qc(object, replicates, k = 3, assay.type = NULL, name = NULL, ...)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- replicates:

  list of numeric vectors, indexes of replicates

- k:

  The number of factors of unwanted variation to be estimated from the
  data.

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay in case of multiple assays

- ...:

  other parameters passed to `link[RUVSeq]{RUVs}`

## Value

A SummarizedExperiment object with the normalized data.

## Examples

``` r
data(toy_notame_set)
# Batch correction
replicates <- list(which(toy_notame_set$QC == "QC"))
batch_corrected <- ruvs_qc(toy_notame_set, replicates = replicates)
# Evaluate batch correction
pca_bhattacharyya_dist(toy_notame_set, batch = "Batch")
#>          [,1]     [,2]
#> [1,]       NA 11.28682
#> [2,] 11.28682       NA
pca_bhattacharyya_dist(batch_corrected, batch = "Batch")
#> Warning: Variance is below eps for 1 variables. Not scaling them.
#>          [,1]     [,2]
#> [1,]       NA 1.516404
#> [2,] 1.516404       NA
```
