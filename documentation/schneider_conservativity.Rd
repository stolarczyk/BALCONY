\name{schneider_conservation}
\alias{schneider_conservation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Schneider conservation metric
}
\description{
This function facilitates the calculation of Schneider conservation metric.
}
\usage{
conservation_score = schneider_conservation(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{Alignment data read with read.alignment() function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{conservation_score}{A vector of length equal to the length of aligned sequences}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Plucennik & Michal Stolarczyk
}
\note{
Please note that the Schneider matric formula can be found in the paper listed in "See Also" section below.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract
}
\examples{
load(alignment)
conservation_score = schneider_conservation(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conservation}% use one of  RShowDoc("KEYWORDS")
\keyword{mertic}
\keyword{alignment}% __ONLY ONE__ keyword per line
