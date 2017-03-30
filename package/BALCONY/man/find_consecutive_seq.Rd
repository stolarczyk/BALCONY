\name{find_consecutive_seq}
\alias{find_consecutive_seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find sequences of numbers in a numeric vector
}
\description{
This function finds sequences of consecutive numbers in numeric vectors
}
\usage{
find_consecutive_seq(vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vector}{A numeric vector to be analyzed}
}
\details{
 Out of the following vector: 1,2,3,4,5,6,7,20,21,140,141 the function will find values starting the sequences: 1,20,140 and their lengths 7,2,2 respectively}
\value{
\item{values}{A vector of values starting the consecutive sequences}
\item{lengths}{A vector of lengths of identified sequences}

}
\author{
Michal Stolarczyk & Alicja Pluciennik
}
\examples{
find_consecutive_seq(c(1,2,3,4,5,6,7,20,21,140,141,300,301,302))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{find}% use one of  RShowDoc("KEYWORDS")
\keyword{sequences}% __ONLY ONE__ keyword per line
\keyword{lengths}
\keyword{consecutive}
