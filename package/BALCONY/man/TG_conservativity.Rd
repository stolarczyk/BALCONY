\name{TG_conservativity}
\alias{TG_conservativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Tunneling Group (TG) conservation metric
}
\description{
This function facilitates the calculation of TG conservation metric.
}
\usage{
TG_conservativity(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{Alignment data read with read.alignment() function}
}
\details{
The conservativity score is calculated according to the following formula:
\deqn{P(i) = max(p(i))/n(i)}
\deqn{Pnorm(i) = P(i)/max(P)}
\deqn{score = -ln(P_norm(i))/max(-ln(P_norm))}
\cr where:\cr
\eqn{p(i)} - amino acids frequency on i-th position where gaps are included \cr
\eqn{n(i)} - amino acids count on i-th position where gaps are excluded
}
\value{
\item{conservation_score}{A vector of length equal to the length of aligned sequences}
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{
Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
\deqn{conservation = 1 - entropy}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract}
}
\examples{
data("alignment")
conservation_score = TG_conservativity(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conservation}% use one of  RShowDoc("KEYWORDS")
\keyword{mertic}
\keyword{alignment}% __ONLY ONE__ keyword per line
