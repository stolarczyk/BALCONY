\name{calculate_AA_variation}
\alias{calculate_AA_variation}

\title{
Calculate AA variations on each position of the multiple sequence alignment
}
\description{
 This function calculates AA variations on each position of the alignment which may be further used for the conservativity study of the set of sequences in quiestion
}
\usage{
calculate_AA_variation(alignment, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{The data loaded with read.alignment() function }
  \item{threshold}{A number in range 0-100. A of minimal frequency of occurences of amino acids at each position}
}
\details{
The output consists of amino acids and their fractions on each position of alignment. Amino acids with occurence frequencies lower than the threshold of user's choice are excluded.
}
\value{
Returns list of two matrices with tabelarized symbols of the most common AA in alignment column and percentage values for contributed AA
 \item{var_aa$AA }{A matrix of AA on all alignment positions with decreasing frequencies in columns}
 \item{var_aa$per}{The percentage of AA frequencies corresponding to the $AA}
 \item{var_aa$matrix}{A combination of this two. The best suited element for visual inspection of the variability at each position}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Pluciennik & Michal Stoalrczyk
}
\note{
}

\seealso{
alignment_parameters(), distplay_AA_variation()
}
\examples{
load(alignment)
threshold=10
var_aa=calculate_AA_variation(alignment, threshold)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
