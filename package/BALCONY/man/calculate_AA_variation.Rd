\name{calculate_AA_variation}
\alias{calculate_AA_variation}

\title{
Calculate AA variations on each position of the multiple sequence alignment
}
\description{
 This function calculates AA variations on each position of the alignment which may be further used for the conservativity study of the set of sequences in quiestion
}
\usage{
calculate_AA_variation(alignment, threshold, grouped,
  grouping_method, weights, pseudo_counts=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{The data loaded with \code{\link[seqinr]{read.alignment}} function }
  \item{threshold}{(optional) A number in range 0-1. A of minimal frequency of occurences of amino acids at each position. Default: all the residues are visualized.}
  \item{grouped}{(optional) A logical indicating if the grouping of amino acids should be applied. Default: FALSE}
  \item{grouping_method}{(optional) A string which specifies the grouping method to be used. One of following: 'substitution_matrix', 'polarity', 'size', 'aromaticity'. Default: 'substitution_matrix'. Default: 'substitution_matrix' if grouped is TRUE.}
  \item{weights}{(optional) A vector of length equal number of sequences in the alignment object with weights to overcome the taxonomic bias in the conservation analysis.}
    \item{pseudo_counts}{(optional) A logical indicating if pseudo-counts should be added to the MSA. Pseudo-counts can be used only in non-group mode and without weights. Using these options with pseudo-counts will be suppressed. Default: FALSE}
}
\details{
The output consists of amino acids and their fractions on each position of alignment. Amino acids with occurence frequencies lower than the threshold of user's choice are excluded.
}
\value{
Returns list of thre matrices with tabelarized symbols of the most common AA in alignment column, percentage values for contributed AA and combined one.
 \item{var_aa$AA }{A matrix of AA on all alignment positions with decreasing frequencies in columns}
 \item{var_aa$per}{The percentage of AA frequencies corresponding to the $AA}
 \item{var_aa$matrix}{A combination of this two. The best suited element for visual inspection of the variability at each position}
}
\author{
Michal Stolarczyk & Alicja Pluciennik
}
\seealso{
\code{\link{align_params}},\code{\link{calculate_pseudo_counts}}
}
\examples{
data("small_alignment")
alignment = delete_isoforms(small_alignment)
threshold=10
grouped = FALSE
var_aa=calculate_AA_variation(small_alignment, threshold, grouped)
}
\keyword{ amino acids }
\keyword{ variation }% __ONLY ONE__ keyword per line
