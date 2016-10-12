\name{consensus}
\alias{consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculates consensus based on read file
}
\description{
Function calclates consensus sequence for given alignment  with choosen threshold, required: {seqinr}.
}
\usage{
consensus(alignment, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{alignment, threshold
  \item{alignment}{from output of read.alignment()}
  \item{threshold}{minimal fraction of amino acids on the certain position in all sequences of the alignment to be taken for consensus letter on this position; number in range 0-100.  }
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If maximum fraction of any amino acid on the certain position is lower than a threshold then "*" is printed instead. 
}
\value{

  \item{consensus }{character of length of the aligned sequence}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja P³uciennik & Micha³ Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
