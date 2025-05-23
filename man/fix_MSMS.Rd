% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{fix_MSMS}
\alias{fix_MSMS}
\title{Transform the MS/MS output to publication ready}
\usage{
fix_MSMS(
  object,
  ms_ms_spectrum_col = "MS_MS_spectrum",
  peak_num = 10,
  min_abund = 5,
  deci_num = 3
)
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{ms_ms_spectrum_col}{name of column with original MS/MS spectra}

\item{peak_num}{maximum number of peak that is kept (Recommended: 4-10)}

\item{min_abund}{minimum relative abundance to be kept (Recommended: 1-5)}

\item{deci_num}{maximum number of decimals to m/z value (Recommended: >2)}
}
\value{
A SummarizedExperiment or MetaboSet object as the one supplied, with 
publication-ready MS/MS peak information.
}
\description{
Change the MS/MS output from MS-DIAL format to publication-ready format.
Original spectra is sorted according to abundance percentage and clarified. 
See the example below.
}
\details{
Original MS/MS spectra from MS-DIAL:
m/z:Raw Abundance

23.193:254 26.13899:5 27.50986:25 55.01603:82 70.1914:16 73.03017:941 
73.07685:13 73.13951:120

Spectra after transformation:
m/z  (Abundance)

73.03 (100), 23.193 (27), 73.14 (12.8), 55.016 (8.7)
}
\examples{
data(example_set)
# Spectra before fixing
ex_set <- example_set
rowData(ex_set)$MS_MS_spectrum <- NA
rowData(ex_set)[1, ]$MS_MS_spectrum <- 
  "28.769:53 44.933:42 52.106:89 69.518:140"
rowData(ex_set)$MS_MS_spectrum[
  !is.na(rowData(ex_set)$MS_MS_spectrum)]
# Fixing spectra with default settings
fixed_MSMS_peaks <- fix_MSMS(ex_set)
# Spectra after fixing
rowData(fixed_MSMS_peaks)$MS_MS_Spectrum_clean[
  !is.na(rowData(fixed_MSMS_peaks)$MS_MS_Spectrum_clean)]

}
