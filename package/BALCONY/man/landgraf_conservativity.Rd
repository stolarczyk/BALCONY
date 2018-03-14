\name{landgraf_conservativity}
\alias{landgraf_conservativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Landgraf conservation score
}
\description{
This function calculates Landgraf conservarion score
}
\usage{
landgraf_conservativity(matrix_name = NULL, alignment, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix_name}{A string with path to the file woith substitution matrix to be used to calculate the Landgraf conservation score. Optional parameter, if not provided the Gonnet substitution matrix is used (according to author's suggestion)}
  \item{alignment}{An alignment object read with \code{\link[seqinr]{read.alignment}} function}
  \item{weights}{A vector with weight for each sequence in the alignment  to be used to calculate the Landgraf conservation score e.g. each sequence similarity to the consensus sequence from the alignment - output from \code{\link{cons2seqs_ident}} fuction}
}

\value{
\item{score}{A vector of length equal to the length of aligned sequences}
}

\author{
Michal Stolarczyk & Alicja Pluciennik
}
\note{
Please note that the Shannon matric formula can be found in the paper listed in "See Also" section below.
Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
\deqn{conservation = 1 - entropy}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{consensus}}, \code{\link{cons2seqs_ident}}, \code{\link[seqinr]{read.alignment}}

\url{http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract}
}
\examples{
data("small_alignment")
alignment = small_alignment
threshold_consensus = 30
consensus_seq=consensus(alignment, threshold_consensus);
consensus_sequences_identity=cons2seqs_ident(alignment, consensus_seq)
score = landgraf_conservativity(alignment = alignment, weights = consensus_sequences_identity)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Landgraf }% use one of  RShowDoc("KEYWORDS")
\keyword{ conservation }% __ONLY ONE__ keyword per line
