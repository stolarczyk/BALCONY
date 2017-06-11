\name{consensus}
\alias{consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consensus sequence determination
}
\description{
Function calculates consensus sequence for given alignment with a threshold of user's choice.
}
\usage{
consensus(alignment, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{output of of \code{\link[seqinr]{read.alignment}} function or grouped alignment created with: \code{\link{aligned_sequences_matrix2groups}} and \code{\link{alignment2matrix}}}
  \item{threshold}{minimal fraction of amino acids on the certain position in all sequences of the alignment to be taken for consensus letter on this position; number in range 0-100.  }
  }
%%     ~~Describe \code{x} here~~
\details{
%%  ~~ If necessary, more details than the description above ~~
If maximum fraction of any amino acid on the certain position is lower than a threshold then "*" is printed instead.
}
\value{
\item{consensus_sequence}{A character vector of length of the aligned sequence containing consesus sequence based on the input alignment}}

\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{
Please note that this function masks the seqinr package function  \code{\link[seqinr]{consensus}}
}

\seealso{
\code{\link[seqinr]{s2c}}
}
\examples{
data("alignment")
alignment = delete_isoforms(alignment)
threshold=80 # Set the consensus threshold
consensus_sequence=consensus(alignment, threshold)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ consensus_sequence }
\keyword{ alignment }
