% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_metrics.R
\name{assess_quality}
\alias{assess_quality}
\alias{assess_quality,MetaboSet-method}
\alias{assess_quality,SummarizedExperiment-method}
\title{Assess quality information of features}
\usage{
assess_quality(object, ...)

\S4method{assess_quality}{MetaboSet}(object)

\S4method{assess_quality}{SummarizedExperiment}(object, assay.type = NULL)
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{...}{additional arguments passed to methods}

\item{assay.type}{character, assay to be used in case of multiple assays}
}
\value{
A SummarizedExperiment or MetaboSet object with quality metrics in 
feature data.
}
\description{
Assess features using the quality metrics defined in (Broadhurst 
2018). The quality metrics are described in Details section of 
\code{\link{flag_quality}}
}
\examples{
data(example_set)
ex_set <- assess_quality(example_set)
rowData(ex_set)

}
