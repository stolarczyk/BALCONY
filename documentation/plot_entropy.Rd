\name{plot_entropy}
\alias{plot_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot entropies for protein
}
\description{
This function allows to plot entropies for protein. Plots might be superimposed or not.
}
\usage{
plot_entropy(ProtConservation,description, colors, impose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ProtConservation}{list or a vector of protein conservation, entropies}
  \item{description}{vector of names for entropy/conservation scores}
  \item{colors}{(optional) vector of colors for each plot, default-rainbow}
  \item{impose}{(optional)boolean, if True/T plots are superimposed, if False/F plots are printed separately}


}
\details{
Function makes plots for given values, on X axis are further amino acids, on Y axis are values of entropy/conservation. Legend contains score names form description values. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Płuciennik & Michał Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
FIX ME
NAMES=c("Shannon","Schneider", "Kabat", "Landgraf", "TG")
ProtConservation=get_prot_entropy(structure_index$proteinIndices,SCORE_LIST, NAMES)
colors=c("red", "green", "cyan", "hotpink", "hotpink2")
plot_entropy(ProtConservation,NAMES, colors, T)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
