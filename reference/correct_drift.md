# Correct drift using cubic spline

A wrapper function for applying cubic spline drift correction and saving
before and after plots.

## Usage

``` r
correct_drift(
  object,
  log_transform = TRUE,
  spar = NULL,
  spar_lower = 0.5,
  spar_upper = 1.5,
  check_quality = FALSE,
  condition = "RSD_r < 0 & D_ratio_r < 0",
  file = NULL,
  width = 16,
  height = 8,
  color = "QC",
  shape = color,
  color_scale = getOption("notame.color_scale_dis"),
  shape_scale = scale_shape_manual(values = c(15, 16)),
  assay.type = NULL,
  name = NULL
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- log_transform:

  logical, should drift correction be done on log-transformed values?
  See Details

- spar:

  smoothing parameter as in
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html)

- spar_lower, spar_upper:

  lower and upper limits for the smoothing parameter

- check_quality:

  logical, whether quality should be monitored.

- condition:

  a character specifying the condition used to decide whether drift
  correction works adequately, see Details

- file:

  path to the PDF file where the plots should be saved

- width, height:

  width and height of the plots in inches

- color:

  character, name of the column used for coloring the points

- shape:

  character, name of the column used for shape

- color_scale, shape_scale:

  the color and shape scales as returned by a ggplot function

- assay.type:

  character, assay to be used in case of multiple assays

- name:

  character, name of the resultant assay

## Value

A SummarizedExperiment object as the one supplied, with drift corrected
features.

## Details

If `log_transform = TRUE`, the correction will be done on
log-transformed values. The correction formula depends on whether the
correction is run on original values or log-transformed values. In
log-space: \\corrected = original + mean of QCs - prediction by cubic
spline\\. In original space: \\corrected = original \* prediction for
first QC / prediction for current point\\. We recommend doing the
correction in the log-space since the log-transfomred data better
follows the assumptions of cubic spline regression. The drift correction
in the original space also sometimes results in negative values, and
results in rejection of the drift corrrection procedure. If
`check_quality = TRUE`, the `condition` parameter should be a character
giving a condition compatible with
[`filter`](https://dplyr.tidyverse.org/reference/filter.html). The
condition is applied on the **changes** in the quality metrics RSD,
RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r \< 0
and D_ratio_r \< 0", meaning that both RSD_r and D_ratio_r need to
decrease in the drift correction, otherwise the drift corrected feature
is discarded and the original is retained. By default, the column used
for color is also used for shape.

## See also

[`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html) for
details about the regression

## Examples

``` r
data(toy_notame_set)
corrected <- correct_drift(mark_nas(toy_notame_set[1:5, ], value = 0))
#> INFO [2026-02-26 10:48:29] Starting drift correction
#> INFO [2026-02-26 10:48:31] Recomputing quality metrics for drift corrected data
#> INFO [2026-02-26 10:48:31] Drift correction performed
#> INFO [2026-02-26 10:48:31] Inspecting drift correction results
#> INFO [2026-02-26 10:48:31] Original quality metrics missing, recomputing
#> INFO [2026-02-26 10:48:31] Drift correction results inspected: Drift_corrected: 100%
```
