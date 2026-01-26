# Simple imputation

Impute missing values using a simple imputation strategy. All missing
values of a feature are imputed with the same value. It is possible to
only impute features with a large number of missing values this way.
This can be useful for using this function before random forest
imputation to speed things up. The imputation strategies available are:

- a numeric value: impute all missing values in all features with the
  same value, e.g. 1

- "mean": impute missing values of a feature with the mean of observed
  values of that feature

- "median": impute missing values of a feature with the median of
  observed values of that feature

- "min": impute missing values of a feature with the minimum observed
  value of that feature

- "half_min": impute missing values of a feature with half the minimum
  observed value of that feature

- "small_random": impute missing values of a feature with random numbers
  between 0 and the minimum of that feature (uniform distribution,
  remember to set the seed number!).

## Usage

``` r
impute_simple(object, value, na_limit = 0, assay.type = NULL, name = NULL)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- value:

  the value used for imputation, either a numeric or one of '"min",
  "half_min", "small_random", see above

- na_limit:

  only impute features with the proportion of NAs over this limit. For
  example, if `na_limit = 0.5`, only features with at least half of the
  values missing are imputed.

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay in case of multiple assays

## Value

A SummarizedExperiment object with imputed peak table.

## Examples

``` r
data(toy_notame_set)
missing <- mark_nas(toy_notame_set, 0)
imputed <- impute_simple(missing, value = "min")
```
