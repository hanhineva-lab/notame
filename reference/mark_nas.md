# Mark specified values as missing

Replaces all values in the peak table that equal the specified value
with NA. For example, vendor software might use 0 or 1 to signal a
missing value, which is not understood by R.

## Usage

``` r
mark_nas(object, value, assay.type = NULL, name = NULL)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- value:

  the value to be converted to NA

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay in case of multiple assays

## Value

SummarizedExperiment object as the one supplied, with missing values
correctly set to NA.

## Examples

``` r
data(toy_notame_set)
nas_marked <- mark_nas(toy_notame_set, value = 0)
```
