\name{RealValET_conservativity}
\alias{RealValET_conservativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate real-value Evolutionary Trace (ET)
}
\description{
This function allows to calculate real-valued ET for MSA.
}
\usage{
realval_ET = RealValET_conservativity(alignment)
}

\arguments{
  \item{alignment}{Alignment data read with read.alignment() function}
}
\details{
Here, the real-valued ET is calculated using an evolutionary tree calculated for given alignment. The tree is calculated using UPGMA method.
Real-valued ET score can be used as complimentary analysis of evolutionary entropy measures.
}
\value{\
item{realValuedET}{A vector of real valued ET score corresponding to each MSA column }
}
\references{ Mihalek, Res, Lichtarge, 2004}
\author{
Alicja Plucennik & Michal Stolarczyk
}
\examples{
data("alignment")
conservation_score = schneider_conservativity(alignment)
realval_ET = RealValET_conservativity(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conservation}% use one of  RShowDoc("KEYWORDS")
\keyword{Evolutionary Trace}% __ONLY ONE__ keyword per line
