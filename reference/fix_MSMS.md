# Transform the MS/MS output to publication ready

Change the MS/MS output from MS-DIAL format to publication-ready format.
Original spectra is sorted according to abundance percentage and
clarified. See the example below.

## Usage

``` r
fix_MSMS(
  object,
  ms_ms_spectrum_col = "MS_MS_spectrum",
  peak_num = 10,
  min_abund = 5,
  deci_num = 3
)
```

## Arguments

- object:

  a
  ` `[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object

- ms_ms_spectrum_col:

  name of column with original MS/MS spectra

- peak_num:

  maximum number of peak that is kept (Recommended: 4-10)

- min_abund:

  minimum relative abundance to be kept (Recommended: 1-5)

- deci_num:

  maximum number of decimals to m/z value (Recommended: \>2)

## Value

A SummarizedExperiment object as the one supplied, with
publication-ready MS/MS peak information.

## Details

Original MS/MS spectra from MS-DIAL: m/z:Raw Abundance

23.193:254 26.13899:5 27.50986:25 55.01603:82 70.1914:16 73.03017:941
73.07685:13 73.13951:120

Spectra after transformation: m/z (Abundance)

73.03 (100), 23.193 (27), 73.14 (12.8), 55.016 (8.7)

## Examples

``` r
data(toy_notame_set)
# Spectra before fixing
ex_set <- toy_notame_set
rowData(ex_set)$MS_MS_spectrum <- NA
rowData(ex_set)[1, ]$MS_MS_spectrum <- 
  "28.769:53 44.933:42 52.106:89 69.518:140"
rowData(ex_set)$MS_MS_spectrum[
  !is.na(rowData(ex_set)$MS_MS_spectrum)]
#> [1] "28.769:53 44.933:42 52.106:89 69.518:140"
# Fixing spectra with default settings
fixed_MSMS_peaks <- fix_MSMS(ex_set)
#> INFO [2026-01-26 11:40:26] Saving fixed MS/MS spectra to column 'MS_MS_Spectrum_clean' in rowData
# Spectra after fixing
rowData(fixed_MSMS_peaks)$MS_MS_Spectrum_clean[
  !is.na(rowData(fixed_MSMS_peaks)$MS_MS_Spectrum_clean)]
#> [1] "69.518 (100), 52.106 (63.6), 28.769 (37.9), 44.933 (30)"
```
