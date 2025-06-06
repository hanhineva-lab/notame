% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{write_to_excel}
\alias{write_to_excel}
\title{Write results to Excel file}
\usage{
write_to_excel(object, file, ...)
}
\arguments{
\item{object}{a \code{
\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
or \code{\link{MetaboSet}} object}

\item{file}{path to the file to write}

\item{...}{Additional parameters passed to
\code{\link[openxlsx]{write.xlsx}}}
}
\value{
None, the function is invoked for its side effect.
}
\description{
Writes all the data in a SummarizedExperiment or MetaboSet object to an 
Excel spreadsheet.
The format is similar to the one used to read data in, except for the fact 
that EVERYTHING NEEDS TO BE WRITTEN AS TEXT. To fix numeric values in Excel,
choose any cell with a number, press Ctrl + A, then go to the dropdown menu
in upper left corner and choose "Convert to Number". This will fix the file,
but can take quite a while.
}
\examples{
\dontshow{.old_wd <- setwd(tempdir())}
data(example_set)
write_to_excel(example_set, file = "example_set.xlsx")
\dontshow{setwd(.old_wd)}

}
