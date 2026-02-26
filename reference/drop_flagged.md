# Drop flagged features

Removes all features that have been flagged by quality control
functions. Only features that do not have a flag (Flag == NA) are
retained.

## Usage

``` r
drop_flagged(object, all_features = FALSE)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- all_features:

  logical, should all features be retained? Mainly used by internal
  functions

## Value

A SummarizedExperiment object without the previously flagged features.

## Examples

``` r
data(toy_notame_set)
dim(toy_notame_set)
#> [1] 80 50
flagged <- flag_quality(toy_notame_set)
#> INFO [2026-02-26 10:48:32] 
#> 92% of features flagged for low quality
noflags <- drop_flagged(flagged)
dim(noflags)
#> [1]  6 50
```
