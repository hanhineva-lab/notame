# Join new columns to feature data

Join a new data frame of information to feature data of a
SummarizedExperiment object. The data frame needs to have a column
"Feature_ID". This function is usually used internally by some of the
functions in the package, but can be useful.

## Usage

``` r
join_rowData(object, dframe)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- dframe:

  a data frame with the new information

## Value

A SummarizedExperiment object with the new information added to
rowData(object).

## Examples

``` r
data(toy_notame_set)
new_info <- data.frame(
  Feature_ID = rownames(toy_notame_set),
  Feature_number = seq_len(nrow(toy_notame_set))
)
with_new_info <- join_rowData(toy_notame_set, new_info)
colnames(rowData(with_new_info))
#> [1] "Feature_ID"     "Split"          "Alignment"      "Average_Mz"    
#> [5] "Average_Rt_min" "Column"         "Ion_mode"       "Flag"          
#> [9] "Feature_number"
```
