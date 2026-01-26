# Join new columns to pheno data

Join a new data frame of information to pheno data of a
SummarizedExperiment object.

## Usage

``` r
join_colData(object, dframe)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- dframe:

  a data frame with the new information

## Value

A SummarizedExperiment object with the new information added to
colData(object).

## Examples

``` r
data(toy_notame_set)
new_info <- data.frame(
  Sample_ID = colnames(toy_notame_set),
  BMI = stats::runif(ncol(toy_notame_set), 22, 26)
)
with_new_info <- join_colData(toy_notame_set, new_info)
colnames(colData(with_new_info))
#>  [1] "Sample_ID"          "Injection_order"    "Subject_ID"        
#>  [4] "Group"              "QC"                 "Time"              
#>  [7] "Batch"              "HILIC_neg_Datafile" "HILIC_pos_Datafile"
#> [10] "RP_neg_Datafile"    "RP_pos_Datafile"    "BMI"               
```
