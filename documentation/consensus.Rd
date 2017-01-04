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
consensus(alignment, threshold)
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
\value{
character vector of length of the aligned sequence containing consesus sequence based on the input alignment}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
seqinr
}
\examples{
alignment=read.alignment(file="alignment.mfa", format="fasta", forceToLower=F) 
threshold=80
cons=consensus(alignment, threshold)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ consensus_sequence }
