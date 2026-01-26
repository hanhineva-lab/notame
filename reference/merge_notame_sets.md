# Merge SummarizedExperiment objects together

Merges two or more SummarizedExperiment objects together. Can be used to
merge analytical modes or batches.

## Usage

``` r
merge_notame_sets(..., merge = c("features", "samples"), assay.type = NULL)
```

## Arguments

- ...:

  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  objects or a list of objects

- merge:

  what to merge? features is used for combining analytical modes,
  samples is used for batches

- assay.type:

  character, assay to be used in case of multiple assays. The same assay
  needs to be present in all objects to be merged, and the resultant
  object contains this single assay.

## Value

A merged SummarizedExperiment object.

## Details

When merging samples, sample IDs that begin with "QC" or "Ref" are
combined so that they have running numbers on them. This means that if
both batches have samples called "QC_1", this will not result in an
error, but the sample IDs will be adjusted so that they are unique.

## Examples

``` r
# Merge analytical modes
data(hilic_neg_sample, hilic_pos_sample, rp_neg_sample, rp_pos_sample)
merged <- merge_notame_sets(
  hilic_neg_sample, hilic_pos_sample,
  rp_neg_sample, rp_pos_sample
)
# Merge batches
batch1 <- toy_notame_set[, toy_notame_set$Batch == 1]
batch2 <- toy_notame_set[, toy_notame_set$Batch == 2]
merged <- merge_notame_sets(batch1, batch2, merge = "samples")
```
