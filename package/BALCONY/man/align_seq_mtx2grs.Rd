\name{align_seq_mtx2grs}
\alias{align_seq_mtx2grs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert amino acid symbols to groups according to their properties of user's choice.
}
\description{
This function performs a conversion of amino acid symbols to group symbols according to their properties. Implemented grouping methods are: substitution_matrix (majority of properties taken into account), polarity, size and aromaticity. "GX", where X stands for group number, are group symbols.
}
\usage{
align_seq_mtx2grs(aligned_sequences_matrix,grouping_method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aligned_sequences_matrix}{A matrix that contains aligned sequences. It is an output of \code{\link{alignment2matrix}} function}
  \item{grouping_method}{A string which specifies the grouping method to be used. One of following: 'substitution_matrix', 'polarity', 'size', 'aromaticity' }
}
\value{
\item{grouped_aligned_sequences_matrix}{A matrix of size of the input matrix but with group symbols instead of amino acid symbols}
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}

\seealso{
\code{\link{alignment2matrix}}, \code{\link[seqinr]{read.alignment}}
}
\examples{
data(alignment)
alignment = delete_isoforms(alignment)
grouping_method = "general"
aligned_sequences_matrix = alignment2matrix(alignment)
grouped = align_seq_mtx2grs(aligned_sequences_matrix,grouping_method)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ groups }% use one of  RShowDoc("KEYWORDS")
\keyword{ matrix }% __ONLY ONE__ keyword per line
