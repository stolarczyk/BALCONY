\name{outlying_sequences}
\alias{outlying_sequences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outlying sequences detection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
outlying_sequences(percentage, alignment_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{percentage}{The identity of each sequence in the alignment to the consensus sequence}
  \item{alignment_file}{Alignment loaded with read.alignment() function}
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{best_consensus}{Sequence closest to the consensus}
\item{worst_consensus}{Sequence most different to the consensus}
\item{most_common}{Most common sequence in the alignment ???}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
file=read.alignment(file="alignment.fasta", format="fasta", forceToLower=F)
aligned_sequences=file[[3]];
seq_count=file[[1]];
threshold_consensus=30;
consensus_seq=consensus(file, threshold_consensus);
consensus_sequences_identity=cons2seqs_ident(aligned_sequences, seq_count, consensus_seq)
outliers=outlying_sequences(consensus_sequences_identity, file);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Consensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ Outliers }% __ONLY ONE__ keyword per line
