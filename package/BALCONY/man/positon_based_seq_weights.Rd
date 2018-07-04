\name{get_pos_based_seq_weights}
\alias{get_pos_based_seq_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get position based weights of sequences in alignment
}
\description{
This function calculates position based weights of sequences based on Heinkoff & Heinkoff (1994) for given MSA. The score is calculated as sum of scores for each sequence position c. Score for position c is equal 1/r if there is r different residues at column c in MSA but 1/rs if r symbol is repeated in s sequences.
 
}
\usage{
get_pos_based_seq_weights(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{alignment loaded with \code{\link[seqinr]{read.alignment}}
  \item{gap}{logical parameter, if TRUE(default) the gaps in MSA are included}
  \item{normalized}{logical parameter, if TRUE (default) weights for all sequences are divided by number of columns in alignment (when gap = TRUE weights sum up to 1)}
}
}
\details{
  The weights might be calculated only for amino acids symbols or for all symbols (including gaps). Also weights can be normalized by number of columns in MSA, then the sum of weights for all sequences is 1.
}
\value{
 \item{weights}{a vector of position based weights for each sequence in given alignment}
}
\references{
Steven Henikoff and Jorja G. Henikoff,Position-based Sequence Weights,J: Mol. Biol. (1994) 243, 574-578
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
data("alignment")
pos_based_weights <- positon_based_seq_weights(alignment)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
