\name{barplotshow}
\alias{barplotshow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shows barplot with amino acid variation on the specified position
}
\description{
This function facilitates a visual inspection of multiple sequence alignment (MSA) position variablity.
}
\usage{
barplotshow(position, AA_variation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{A number of column of alignment to be visualized}
  \item{AA_variation}{A percentage frequency of amino acids in the alignment, calculated with \code{\link{calculate_AA_variation}} function}
}
\value{
This function produces a barchart
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\seealso{
\code{\link{calculate_AA_variation}}
}
\examples{
data("small_alignment")
position = 100
threshold = 0.01
var_aa = calculate_AA_variation(small_alignment,threshold)
barplotshow(position, var_aa)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{barplot}
\keyword{amino acids variation}% __ONLY ONE__ keyword per line
