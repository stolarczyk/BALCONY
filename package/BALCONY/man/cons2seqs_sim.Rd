\name{cons2seqs_sim}
\alias{cons2seqs_sim}

\title{
Group consensus to each sequence in the alignment similarity
}
\description{
The function calculates similarity of group consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment. Grouping amino acids allows to check similiarity between sequences by amino acids properties of user's choice.
}
\usage{
cons2seqs_sim(alignment, consensus_seq, grouping_method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{The output of read.alignment() function}
  \item{consensus_seq}{The output of consensus() function}
  \item{grouping_method}{One of following: 'hydrophobicity','size','aromaticity','general'}
}
\details{
AA in  consensus sequences and aligned sequences are converted into groups symbols according to method of user's choice. Returned values are percentage of similar amino acids considering the properties in consensus sequence and aligned sequence.

}
\value{
\item{percentage}{numeric vector of identity score (percentage); positions in the numeric vector correspond to sequences in alignment, respectively}
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}


\seealso{
read.alignment(), consensus(), alignment_parameters()
}
\examples{
data("alignment")
alignment = delete_isoforms(alignment)
threshold_consensus = 30
grouping_method = "general"
consensus_seq = consensus(alignment, threshold_consensus)
consensus_to_seqs_similarity = cons2seqs_sim(alignment, consensus_seq, grouping_method)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ amino acid }
\keyword{ properties }% __ONLY ONE__ keyword per line
