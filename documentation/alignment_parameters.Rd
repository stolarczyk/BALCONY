\name{alignment_parameters}
\alias{alignment_parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get alignment dimensions
}
\description{
This function returns size of alignment, which facilitates upcoming steps of analysis with BALCONY package.
}
\usage{
alignment_parameters(alignment_sequences)
}

\arguments{
  \item{alignment_sequences}{ alignment[[3]] loaded with read.alignment()}

}
}
\details{
Function returns list of two elements row_no(number of rows, sequences) and col_no(number of columns,length of aligned sequences)
}
\value{

  \item{row_no }{number of sequences}
  \item{col_no }{length of aligned sequences}

}
\references{

}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{}

\seealso{
read.alignment {seqinr}
}
\examples{
aligned_sequences=file[[3]]
parametry=alignment_parameters(aligned_sequences);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ alignment }
\keyword{ dimension }% __ONLY ONE__ keyword per line
