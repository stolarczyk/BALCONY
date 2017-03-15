\name{landgraf_conservativity}
\alias{landgraf_conservativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Landgraf conservation score
}
\description{
This function calculates Landgraf conservarion score
}
\usage{
score = landgraf_conservativity(matrix_name = NULL, alignment_file, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix_name}{A string with path to the file woith substitution matrix to be used to calculate the Landgraf conservation score. Optional parameter, if not provided the Gonnet substitution matrix is used (according to author's suggestion)}
  \item{alignment_file}{An alignment object read with read.alignment() function}
  \item{weights}{A vector with weight for each sequence in the alignment  to be used to calculate the Landgraf conservation score e.g. each sequence similarity to the consensus sequence from the alignment - output from cons2seqs_ident() fuction}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michal Stolarczyk & Alicja Pluciennik
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
consensus(), cons2seqs_ident(), read.alignment()
}
\examples{
load(alignment)
threshold_consensus = 30
consensus_seq=consensus(file, threshold_consensus);
consensus_sequences_identity=cons2seqs_ident(alignment, consensus_seq)
score = landgraf_conservativity(alignment_file = alignment, weights = consensus_sequences_identity)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Landgraf }% use one of  RShowDoc("KEYWORDS")
\keyword{ conservation }% __ONLY ONE__ keyword per line
