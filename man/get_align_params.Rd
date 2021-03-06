% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conservationFuncs.R
\name{get_align_params}
\alias{get_align_params}
\title{Get alignment dimensions}
\usage{
get_align_params(alignment)
}
\arguments{
\item{alignment}{data loaded with \code{\link[seqinr]{read.alignment}}}
}
\value{
\item{row_no }{number of sequences}
\item{col_no }{length of aligned sequences}
}
\description{
This function returns size of alignment, which facilitates the convenient performing upcoming steps of analysis.
}
\details{
Function returns list of two elements row_no(number of rows, sequences) and col_no(number of columns,length of aligned sequences)
}
\examples{
data("alignment")
parameters=get_align_params(alignment);
}
