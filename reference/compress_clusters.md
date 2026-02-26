# Compress clusters of features to a single feature

This function compresses clusters found by cluster_features, keeping
only the feature with the highest median peak area. The features that
were discarded are recorded in feature data, under Cluster_features.

## Usage

``` r
compress_clusters(object)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

## Value

A SummarizedExperiment object with only one feature per cluster.

## See also

[`cluster_features`](https://hanhineva-lab.github.io/notame/reference/cluster_features.md)

## Examples

``` r
data(toy_notame_set)
clustered <- cluster_features(toy_notame_set, 
  rt_window = 1, corr_thresh = 0.5, d_thresh = 0.6)
#> INFO [2026-02-26 10:48:29] Identified m/z column Average_Mz and retention time column Average_Rt_min
#> INFO [2026-02-26 10:48:29] Identified m/z column Average_Mz and retention time column Average_Rt_min
#> INFO [2026-02-26 10:48:29] Identified m/z column Average_Mz and retention time column Average_Rt_min
#> INFO [2026-02-26 10:48:29] 
#> Starting feature clustering at 2026-02-26 10:48:29.391987
#> INFO [2026-02-26 10:48:29] Finding connections between features in HILIC_neg
#> INFO [2026-02-26 10:48:29] Found 1 connections in HILIC_neg
#> INFO [2026-02-26 10:48:29] Finding connections between features in HILIC_pos
#> INFO [2026-02-26 10:48:29] Found 4 connections in HILIC_pos
#> INFO [2026-02-26 10:48:29] Finding connections between features in RP_neg
#> INFO [2026-02-26 10:48:29] Found 1 connections in RP_neg
#> INFO [2026-02-26 10:48:29] Finding connections between features in RP_pos
#> INFO [2026-02-26 10:48:29] Found 2 connections in RP_pos
#> INFO [2026-02-26 10:48:29] Found 8 connections
#> 5 components found
#> 1 components found
#> INFO [2026-02-26 10:48:29] Found 5 clusters of 2 or more features, clustering finished at 2026-02-26 10:48:29.412303
compressed <- compress_clusters(clustered)
#> INFO [2026-02-26 10:48:29] Clusters compressed, left with 73 features
```
