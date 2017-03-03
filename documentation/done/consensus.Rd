\name{consensus}
\alias{consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
consensus
}
\description{
Function calculates consensus sequence for given alignment with a threshold of user's choice.
}
\usage{
consensus_sequence = consensus(alignment, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{alignment, threshold
  \item{alignment}{output of of read.alignment() function}
  \item{threshold}{minimal fraction of amino acids on the certain position in all sequences of the alignment to be taken for consensus letter on this position; number in range 0-100.  }
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If maximum fraction of any amino acid on the certain position is lower than a threshold then "*" is printed instead. 
}
\value{consensus_sequence}{A character vector of length of the aligned sequence containing consesus sequence based on the input alignment}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{
Please note that this function masks the seqinr package function "consensus"
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
seqinr
}
\examples{
load(alignment) #Load the alignment data provided with the package
threshold=80 # Set the consensus threshold
consensus_sequence=consensus(alignment, threshold)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ consensus_sequence }
\keyword{ alignment }
