% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{subject_col}
\alias{subject_col}
\alias{subject_col<-}
\title{Get and set the name of special column for subject identifiers}
\usage{
subject_col(object)

subject_col(object) <- value
}
\arguments{
\item{object}{a \code{\link{MetaboSet}} object}

\item{value}{string, name of column to be designated for holding subject 
identifiers}
}
\value{
Character, the name of the subject variable.

For the endomorphism, an object with the subject variable set to the 
specified variable.
}
\description{
Get and set the name of special column for subject identifiers
}
\examples{
data(example_set)
ex_set <- as(example_set, "MetaboSet")
# Get name of subject variable
subject_col(ex_set)
data(example_set)
ex_set <- as(example_set, "MetaboSet")
# Set subject variable
subject_col(ex_set) <- "Subject_ID"
}
