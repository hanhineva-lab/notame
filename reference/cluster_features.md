# Cluster correlated features originating from the same metabolite

Clusters features potentially originating from the same compound.
Features with high Pearson correlation coefficient and small retention
time difference are linked together. Then clusters are formed by setting
a threshold for the relative degree that each node in a cluster needs to
fulfil. Each cluster is named after the feature with the highest median
peak area (median abundance). This is a wrapper around numerous
functions that are based on the MATLAB code by David Broadhurst.

## Usage

``` r
cluster_features(
  object,
  mz_col = NULL,
  rt_col = NULL,
  all_features = FALSE,
  rt_window = 1/60,
  corr_thresh = 0.9,
  d_thresh = 0.8,
  assay.type = NULL
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- mz_col:

  the column name in feature data that holds mass-to-charge ratios

- rt_col:

  the column name in feature data that holds retention times

- all_features:

  logical, should all features be included in the clustering? If FALSE,
  as the default, flagged features are not included in clustering

- rt_window:

  the retention time window for potential links NOTE: use the same unit
  as the retention time

- corr_thresh:

  the correlation threshold required for potential links between
  features

- d_thresh:

  the threshold for the relative degree required by each node

- assay.type:

  character, assay to be used in case of multiple assays

## Value

a SummarizedExperiment object, with median peak area (MPA), the cluster
ID, the features in the cluster, and cluster size added to feature data.

## Examples

``` r
data(toy_notame_set)
# The parameters are really weird because example data is imaginary
clustered <- cluster_features(toy_notame_set, rt_window = 1,
  corr_thresh = 0.5, d_thresh = 0.6)
#> INFO [2026-01-26 11:40:22] Identified m/z column Average_Mz and retention time column Average_Rt_min
#> INFO [2026-01-26 11:40:22] Identified m/z column Average_Mz and retention time column Average_Rt_min
#> INFO [2026-01-26 11:40:22] Identified m/z column Average_Mz and retention time column Average_Rt_min
#> INFO [2026-01-26 11:40:22] 
#> Starting feature clustering at 2026-01-26 11:40:22.423585
#> INFO [2026-01-26 11:40:22] Finding connections between features in HILIC_neg
#> INFO [2026-01-26 11:40:22] Found 1 connections in HILIC_neg
#> INFO [2026-01-26 11:40:22] Finding connections between features in HILIC_pos
#> INFO [2026-01-26 11:40:22] Found 4 connections in HILIC_pos
#> INFO [2026-01-26 11:40:22] Finding connections between features in RP_neg
#> INFO [2026-01-26 11:40:22] Found 1 connections in RP_neg
#> INFO [2026-01-26 11:40:22] Finding connections between features in RP_pos
#> INFO [2026-01-26 11:40:22] Found 2 connections in RP_pos
#> INFO [2026-01-26 11:40:22] Found 8 connections
#> 5 components found
#> 1 components found
#> INFO [2026-01-26 11:40:22] Found 5 clusters of 2 or more features, clustering finished at 2026-01-26 11:40:22.446714
```
