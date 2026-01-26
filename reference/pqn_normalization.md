# Probabilistic quotient normalization

Apply probabilistic quotient normalization (PQN) to the peak table of a
SummarizedExperiment object. By default, reference is calculated from
high-quality QC samples and the median of the reference is used for
normalization. Check parameters for more options.

## Usage

``` r
pqn_normalization(
  object,
  ref = c("qc", "all"),
  method = c("median", "mean"),
  all_features = FALSE,
  assay.type = NULL,
  name = NULL
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- ref:

  character, the type of reference samples to use for normalization.

- method:

  character, the method to use for calculating the reference sample.

- all_features:

  logical, should all features be used for calculating the reference
  sample?

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay in case of multiple assays

## Value

A SummarizedExperiment object with altered feature abundances.

## Examples

``` r
data(toy_notame_set)
pqn_set <- pqn_normalization(toy_notame_set)
#> INFO [2026-01-26 11:40:41] Starting PQN normalization
#> INFO [2026-01-26 11:40:41] Using median of qc samples as reference spectrum
```
