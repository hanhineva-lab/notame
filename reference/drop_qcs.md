# Drop QC samples

Drop QC samples

## Usage

``` r
drop_qcs(object)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

## Value

A SummarizedExperiment object as the one supplied, without QC samples.

## Examples

``` r
data(toy_notame_set)
dim(toy_notame_set)
#> [1] 80 50
noqc <- drop_qcs(toy_notame_set)
dim(noqc)
#> [1] 80 40
```
