% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drift_correction.R
\name{correct_drift}
\alias{correct_drift}
\title{Correct drift using cubic spline}
\usage{
correct_drift(
  object,
  log_transform = TRUE,
  spar = NULL,
  spar_lower = 0.5,
  spar_upper = 1.5,
  check_quality = FALSE,
  condition = "RSD_r < 0 & D_ratio_r < 0",
  plotting = FALSE,
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
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{log_transform}{logical, should drift correction be done on 
log-transformed values? See Details}

\item{spar}{smoothing parameter as in 
\code{\link[stats]{smooth.spline}}}

\item{spar_lower, spar_upper}{lower and upper limits for the smoothing 
parameter}

\item{check_quality}{logical, whether quality should be monitored.}

\item{condition}{a character specifying the condition used to decide whether 
drift correction
works adequately, see Details}

\item{plotting}{logical, whether plots should be drawn}

\item{file}{path to the PDF file where the plots should be saved}

\item{width, height}{width and height of the plots in inches}

\item{color}{character, name of the column used for coloring the points}

\item{shape}{character, name of the column used for shape}

\item{color_scale, shape_scale}{the color and shape scales as returned by a 
ggplot function}

\item{assay.type}{character, assay to be used in case of multiple assays}

\item{name}{character, name of the resultant assay in case of multiple assays}
}
\value{
A SummarizedExperiment or MetaboSet object as the one supplied, with 
drift corrected features.
}
\description{
A wrapper function for applying cubic spline drift correction and saving
before and after plots.
}
\details{
If \code{log_transform = TRUE}, the correction will be done on 
log-transformed values.
The correction formula depends on whether the correction is run on original 
values or log-transformed values.
In log-space: \eqn{corrected = original + mean of QCs - prediction by cubic 
spline}.
In original space: \eqn{corrected = original * prediction for first QC / 
prediction for current point}.
We recommend doing the correction in the log-space since the log-transfomred 
data better follows the assumptions of cubic spline regression. The drift 
correction in the original space also sometimes results
in negative values, and results in rejection of the drift corrrection 
procedure.
If \code{spar} is set to \code{NULL} (the default), the smoothing parameter 
will be separately chosen for each feature from the range
[\code{spar_lower, spar_upper}] using cross validation. 
If  \code{check_quality = TRUE}, the \code{condition} parameter should be a 
character giving a condition compatible with \code{\link[dplyr]{filter}}. 
The condition is applied on the \strong{changes} in the quality metrics
RSD, RSD_r, D_ratio and D_ratio_r. For example, the default is "RSD_r < 0 
and D_ratio_r < 0", meaning that both RSD_r and D_ratio_r need to decrease 
in the drift correction, otherwise the drift corrected feature is discarded 
and the original is retained.
By default, the column used for color is also used for shape.
}
\examples{
data(example_set)
corrected <- correct_drift(example_set)

}
\seealso{
\code{\link[stats]{smooth.spline}} for details about the regression
}
