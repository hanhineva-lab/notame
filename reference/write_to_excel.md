# Write results to Excel file

Writes all the data in a SummarizedExperiment object to an Excel
spreadsheet. The format is similar to the one used to read data in,
except for the fact that EVERYTHING NEEDS TO BE WRITTEN AS TEXT. To fix
numeric values in Excel, choose any cell with a number, press Ctrl + A,
then go to the dropdown menu in upper left corner and choose "Convert to
Number". This will fix the file, but can take quite a while.

## Usage

``` r
write_to_excel(object, file, ...)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- file:

  path to the file to write

- ...:

  Additional parameters passed to
  [`write.xlsx`](https://rdrr.io/pkg/openxlsx/man/write.xlsx.html)

## Value

None, the function is invoked for its side effect.

## Examples

``` r
data(toy_notame_set)
write_to_excel(toy_notame_set, file = "toy_notame_set.xlsx")
#> INFO [2026-01-26 11:40:43] Moved RP_pos_Datafile column to last to get meaningful column names for abundances
```
