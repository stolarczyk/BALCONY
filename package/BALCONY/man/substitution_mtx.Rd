\name{substitution_mtx}
\alias{substitution_mtx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read a substitution matrix
}
\description{
This function facilitates reading of substitution matrices for further use
}
\usage{
substitution_mtx(matrix_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix_name}{
A string with path to the substitution matrix in a text file to be read
}
}

\value{
\item{names}{A vector of characters with amino acid names included in the matrix}
\item{matrix}{A numeric matrix with values}
}

\author{
Michal Stolarczyk & Alicja Pluciennik
}

\examples{
path = system.file("extdata", "GONNET.txt", package = "BALCONY")
sub_mat = substitution_mtx(path)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ substitution matrix }% use one of  RShowDoc("KEYWORDS")
\keyword{ read }% __ONLY ONE__ keyword per line
