# Get and set the values in the flag column

Get and set the values in the flag column

## Usage

``` r
flag(object)

flag(object) <- value
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- value:

  character vector, values for flag column

## Value

Character vector of feature flags.

For the endomorphism, an object with a modified flag column.

## Examples

``` r
data(toy_notame_set)
# Get values in flag column of rowData
flag(toy_notame_set)
#>  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [76] NA NA NA NA NA

data(toy_notame_set)
# Flag a suspicious feature manually
flag(toy_notame_set)[1] <- "Contaminant, known from experience"
```
