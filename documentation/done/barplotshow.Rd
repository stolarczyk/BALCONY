\name{barplotshow}
\alias{barplotshow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shows barplot with amino acid variation on the specified position
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
barplotshow(position, AA_variation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{A cumber of column of alignment to be visualized}
  \item{AA_variation}{A percentage frequency of amino acids in the alignment, calculated with display_AA_variation() function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
Alicja Pluciennik & Michal Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
display_AA_variation(), calculate_AA_variation()
}
\examples{
load(alignment)
position = 13
var_aa = alignment
barplotshow(position, var_aa)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{barplot}
\keyword{amino acids variation}% __ONLY ONE__ keyword per line
