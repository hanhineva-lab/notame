% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{flag}
\alias{flag}
\alias{flag<-}
\alias{flag,SummarizedExperiment-method}
\alias{flag<-,SummarizedExperiment-method}
\title{Get and set the values in the flag column}
\usage{
flag(object)

flag(object) <- value

\S4method{flag}{SummarizedExperiment}(object)

\S4method{flag}{SummarizedExperiment}(object) <- value
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{value}{character vector, values for flag column}
}
\value{
Character vector of feature flags.

For the endomorphism, an object with a modified flag column.
}
\description{
Get and set the values in the flag column
}
\examples{
data(example_set)
# Get values in flag column of rowData
flag(example_set)

data(example_set)
# Flag a suspicious feature manually
flag(example_set)[1] <- "Contaminant, known from experience"
}
