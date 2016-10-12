\name{cons2seqs_sim}
\alias{cons2seqs_sim}

\title{
Group consensus to each sequence in the alignment similarity
}
\description{
The function calculates similarity of group consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment. Grouping amino acids allows to check similiarity between sequences by amino acids properties of user's choice.
}
\usage{
cons2seq_sim(prmt, aligned_sequences_matrix, consensus, grouping_method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prmt}{list of alignment dimensions(number of sequences and length of sequences in alignment)}
  \item{alignment_sequences_matrix}{ output of alignment2matrix()}
  \item{consensus}{output of consensus()}
  \item{grouping_method}{one of following: 'hydrophobicity','size','aromaticity','general'}

}
}
\details{
AA in  consensus sequences and aligned sequences are converted into groups symbols according to method of user's choice. Returned values are percentage of similar amino acids considering the properties in consensus sequence and aligned sequence.

}
\value{
\item{percentage}{numeric vector of identity score (percentage); positions in the numeric vector correspond to sequences in alignment, respectively}
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
read.alignment()[seqinr] consensus() alignment_parameters() cons2seqs_ident()
}
\examples{
alignment=read.alignment(file="alignment.mfa", format="fasta", forceToLower=F) 
threshold=60
cons=consensus(alignment, threshold)
aligned_sequences=alignment[[3]]
number=alignment[[1]]
parameters=alignment_parameters(aligned_sequences)
true_consensus=cons2seqs_sim(parameters, aligned_sequences_matrix, consensus,'hydrophobicity')

hist(true_consensus,breaks=15, col="firebrick3",cex.main=1.8,cex.lab=1.4,main="",xlim=c(0.3,1),ylim=c(0,250),xlab="Sequence to consensus identity score") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }
\keyword{ amino acid }
\keyword{ properties }% __ONLY ONE__ keyword per line
