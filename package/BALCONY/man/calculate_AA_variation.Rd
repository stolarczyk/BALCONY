\name{calculate_AA_variation}
\alias{calculate_AA_variation}

\title{
Calculate AA variations on each position of the multiple sequence alignment
}
\description{
 This function calculates AA variations on each position of the alignment which may be further used for the conservativity study of the set of sequences in quiestion
}
\usage{
calculate_AA_variation(alignment, threshold, grouped, grouping_method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{The data loaded with \code{\link[seqinr]{read.alignment}} function }
  \item{threshold}{(optional) A number in range 0-100. A of minimal frequency of occurences of amino acids at each position. Default: all the residues are visualized.}
  \item{grouped}{(optional) A logical indicating if the grouping of amino acids should be applied. Default: FALSE}
  \item{grouping_method}{(optional) A string which specifies the grouping method to be used. One of following: 'substitution_matrix', 'polarity', 'size', 'aromaticity'. Default: 'substitution_matrix'}
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
\author{
Alicja Pluciennik & Michal Stoalrczyk
}
\seealso{
\code{\link{alignment_parameters}}
}
\examples{
data("alignment")
alignment = delete_isoforms(alignment)
threshold=10
grouped = FALSE
var_aa=calculate_AA_variation(alignment, threshold, grouped)
}
\keyword{ amino acids }
\keyword{ variation }% __ONLY ONE__ keyword per line
