\name{schneider_conservativity}
\alias{schneider_conservativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Schneider conservation metric
}
\description{
This function facilitates the calculation of Schneider conservation metric.
}
\usage{
schneider_conservativity(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{Alignment data read with \code{\link[seqinr]{read.alignment}} function}
}

\value{
\item{conservation_score}{A vector of length equal to the length of aligned sequences}
}

\author{
Alicja Plucennik & Michal Stolarczyk
}
\note{
Please note that the Schneider matric formula can be found in the paper listed in "See Also" section below.
Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
\deqn{conservation = 1 - entropy}
}


\seealso{
\url{http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract}
}
\examples{
data("alignment")
conservation_score = schneider_conservativity(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conservation}% use one of  RShowDoc("KEYWORDS")
\keyword{mertic}
\keyword{alignment}% __ONLY ONE__ keyword per line
