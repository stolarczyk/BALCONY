\name{align_params}
\alias{align_params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get alignment dimensions
}
\description{
This function returns size of alignment, which facilitates the convenient performing upcoming steps of analysis.
}
\usage{
align_params(alignment)
}

\arguments{
  \item{alignment}{ alignment loaded with \code{\link[seqinr]{read.alignment}}}

}
\details{
Function returns list of two elements row_no(number of rows, sequences) and col_no(number of columns,length of aligned sequences)
}
\value{

  \item{row_no }{number of sequences}
  \item{col_no }{length of aligned sequences}

}
\references{
seqinr
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}


\seealso{
\code{\link[seqinr]{read.alignment}}
}
\examples{
data("alignment")
parameters=align_params(alignment);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ alignment }
\keyword{ dimension }% __ONLY ONE__ keyword per line
