\name{noteworthy_seqs}
\alias{noteworthy_seqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find noteworthy sequences in the dataset (aligned sequences)
}
\description{
This function detects noteworthy sequences (most common, closest to the consensus and most different from the consesus) to facilitate convenient detection of outlying sequences that might be excluded from the further analysis.
}
\usage{
noteworthy_seqs(percentage, alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{percentage}{The identity of each sequence in the alignment to the consensus sequence. Output of the \code{\link{cons2seqs_ident}} function}
  \item{alignment}{Alignment loaded with \code{\link[seqinr]{read.alignment}} function}

}
\value{
\item{best_consensus}{Sequence closest to the consensus}
\item{worst_consensus}{Sequence most different to the consensus}
\item{most_common}{Most common sequence in the alignment}
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
data("alignment")
consensus_seq = consensus(alignment, 30)
consensus_to_sequences_identity=cons2seqs_ident(alignment,consensus_seq)
noteworthy_seqs(consensus_to_sequences_identity, alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Consensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ Sequences }% __ONLY ONE__ keyword per line
