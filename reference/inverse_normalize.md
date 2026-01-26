# Inverse-rank normalization

Applies inverse rank normalization to all features to approximate a
normal distribution.

## Usage

``` r
inverse_normalize(object, assay.type = NULL, name = NULL)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay in case of multiple assays

## Value

An object as the one supplied, with normalized features.

## Examples

``` r
data(toy_notame_set)
normalized <- inverse_normalize(toy_notame_set)
```
