\name{cons2seqs_ident}
\alias{cons2seqs_ident}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identity of each sequence in the alignment to the consensus sequence.
}
\description{
The function calculates identity of consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment. Also, it can be used to weight conservativity metrics results in further steps of analysis with BALCONY package.
}
\usage{
cons2seqs_ident(alignment, consensus)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{Data loaded with read.alignment() function }
  \item{consensus}{Consensus sequence (output of consensus() function)}
%%     ~~Describe \code{x} here~~
}
}
\details{
Returned values are percentage of identical symbols (AA and "-") in consensus sequence and aligned sequence.
}
\value{

 \item{percentage }{Numeric vector of identity score (percentage); positions in the numeric vector correspond to sequences in alignment, respectively}

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
 \code{\link{consensus()}}
 \code{\link{cons2seq_sim()}}
}
\examples{
load(alignment)
threshold=60
consensus=consensus(alignment, threshold)
true_consensus=cons2seqs_ident(alignment, consensus)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ identity }% __ONLY ONE__ keyword per line
