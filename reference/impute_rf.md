# Impute missing values using random forest

Impute the missing values in the peak table of the object using a random
forest. The estimated error in the imputation is logged. It is
recommended to set the seed number for reproducibility (it is called
random forest for a reason). This a wrapper around
[`missForest`](https://rdrr.io/pkg/missForest/man/missForest.html). Use
parallelize = "variables" to run in parallel for faster testing. NOTE:
running in parallel prevents user from setting a seed number.

## Usage

``` r
impute_rf(object, all_features = FALSE, assay.type = NULL, name = NULL, ...)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- all_features:

  logical, should all features be used? If FALSE (the default), flagged
  features are removed before imputation.

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay in case of multiple assays

- ...:

  passed to
  [`missForest`](https://rdrr.io/pkg/missForest/man/missForest.html)

## Value

An object as the one supplied, with missing values imputed.

## See also

[`missForest`](https://rdrr.io/pkg/missForest/man/missForest.html) for
detail about the algorithm and the parameters

## Examples

``` r
data(toy_notame_set)
missing <- mark_nas(toy_notame_set, 0)
set.seed(38)
imputed <- impute_rf(missing)
#> INFO [2026-02-26 10:48:36] 
#> Starting random forest imputation at 2026-02-26 10:48:36.923678
#> INFO [2026-02-26 10:48:37] Out-of-bag error in random forest imputation: 0.466
#> INFO [2026-02-26 10:48:37] Random forest imputation finished at 2026-02-26 10:48:37.964655 
```
