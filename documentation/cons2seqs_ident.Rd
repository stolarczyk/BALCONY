\name{cons2seqs_ident}
\alias{cons2seqs_ident}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consensus to each sequence in the alignment identity
}
\description{
The function calculates identity of consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment.
}
\usage{
cons2seqs_ident(alignment_sequences, number_of_seq, consensus)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ , 
  \item{alignment_sequences}{ alignment[[3]] loaded with read.alignment()}
  \item{number_of_seq}{number of sequences in the alignment}
  \item{consensus}{calculated consensus sequences (output of consensus())}
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
Alicja P³uciennik & Micha³ Stolarczyk
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

hist(true_consensus,breaks=15, col="firebrick3",cex.main=1.8,cex.lab=1.4,main="",xlim=c(0.3,1),ylim=c(0,250),xlab="Sequence to consensus identity score") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ identity }% __ONLY ONE__ keyword per line
