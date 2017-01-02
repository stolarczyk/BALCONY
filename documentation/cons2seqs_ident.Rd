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
cons2seqs_ident(alignment_sequences, number_of_seq, consensus)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment_sequences}{the third element of the list loaded with read.alignment() function }
  \item{number_of_seq}{number of sequences in the alignment. The first element of list yielded by alignment_parameters() function}
  \item{consensus}{consensus sequence (output of consensus() function)}
%%     ~~Describe \code{x} here~~
}
}
\details{
Returned values are percentage of identical symbols (AA and "-") in consensus sequence and aligned sequence.
}
\value{

 \item{percentage }{numeric vector of identity score (percentage); positions in the numeric vector correspond to sequences in alignment, respectively}

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
alignment=read.alignment(file="alignment.mfa", format="fasta", forceToLower=F) 
threshold=60
cons=consensus(alignment, threshold)
aligned_sequences=alignment[[3]]
number=alignment[[1]]
true_consensus=cons2seqs_ident(aligned_sequences,number, consensus)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ identity }% __ONLY ONE__ keyword per line
