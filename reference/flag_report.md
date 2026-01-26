# A report of flagged features

Computes the number of features at each stage of flagging for each mode.

## Usage

``` r
flag_report(object)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

## Value

A data frame with the number of features at each stage of flagging.

## Examples

``` r
data(toy_notame_set)
flagged <- toy_notame_set |>
  mark_nas(0) |>
  flag_detection(group = "Group") |>
  flag_quality()
#> INFO [2026-01-26 11:40:28] 
#> 1% of features flagged for low detection rate
#> INFO [2026-01-26 11:40:28] 
#> 91% of features flagged for low quality
flag_report(flagged)
#>       Split Kept Low_quality Total Flagged Low_qc_detection
#> 1 HILIC_neg    3          17    20      17               NA
#> 2 HILIC_pos    1          18    20      19                1
#> 3    RP_neg    2          18    20      18               NA
#> 4    RP_pos    0          20    20      20               NA
```
