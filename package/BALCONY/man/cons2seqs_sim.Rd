\name{cons2seqs_sim}
\alias{cons2seqs_sim}

\title{
Group consensus to each sequence in the alignment similarity
}
\description{
The function calculates similarity of group consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment. Grouping amino acids allows to check similiarity between sequences by amino acids properties of user's choice.
}
\usage{
cons2seqs_sim(grouped_alignment, grouped_consensus_seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grouped_alignment}{The output of \code{\link[seqinr]{read.alignment}} function}
  \item{grouped_consensus_seq}{A string of amino acids, the output of  \code{\link{consensus}} function}
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
\code{\link[seqinr]{read.alignment}}, \code{\link{consensus}}, \code{\link{align_params}}
}
\examples{
data("small_alignment")
alignment = delete_isoforms(small_alignment)
threshold_consensus = 30
grouping_method = "substitution_matrix"
alignment_grouped = align_seq_mtx2grs(alignment2matrix(alignment),grouping_method)
consensus_seq_grouped = consensus(alignment_grouped, threshold_consensus)
consensus_to_seqs_similarity = cons2seqs_sim(alignment_grouped, consensus_seq_grouped)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ amino acid }
\keyword{ properties }% __ONLY ONE__ keyword per line
