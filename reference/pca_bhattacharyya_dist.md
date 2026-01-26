# Bhattacharyya distance between batches in PCA space

Computes Bhattacharyya distance between all pairs of batches using
[`fpc:bhattacharyya.matrix`](https://rdrr.io/pkg/fpc/man/bhattacharyya.matrix.html)
after projecting the samples into PCA space with
[`pca`](https://rdrr.io/pkg/pcaMethods/man/pca.html).

## Usage

``` r
pca_bhattacharyya_dist(
  object,
  batch,
  all_features = FALSE,
  center = TRUE,
  scale = "uv",
  nPcs = 3,
  assay.type = NULL,
  ...
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- batch:

  column name of pheno data giving the batch labels

- all_features:

  logical, should all features be used? If FALSE (the default), flagged
  features are removed before imputation.

- center:

  logical, should the data be centered prior to PCA? (usually yes)

- scale:

  scaling used, as in
  [`prep`](https://rdrr.io/pkg/pcaMethods/man/prep.html). Default is
  "uv" for unit variance

- nPcs:

  the number of principal components to use

- assay.type:

  character, assay to be used in case of multiple assays

- ...:

  other parameters to
  [`pca`](https://rdrr.io/pkg/pcaMethods/man/pca.html)

## Value

A matrix of Bhattacharyya distances between batches.

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
