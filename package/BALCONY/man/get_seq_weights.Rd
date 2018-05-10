\name{get_seq_weights}
\alias{get_seq_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get sequences weights
}
\description{
This function returns weights of the sequneces in the alignment object
}
\usage{
get_seq_weights(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{alignment loaded with \code{\link[seqinr]{read.alignment}}}
}
\details{
The weights are calculated as shown in: \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/1097-0134\%2820010101\%2942\%3A1\%3C108\%3A\%3AAID-PROT110\%3E3.0.CO\%3B2-O}{Valdar and Thronton (2001)}

\strong{According to the following formulas:}

\deqn{W_{j} = \frac{\sum_{k\neq j}^{N}Dist(s_{j},s_{k}))}{N-1}}{Wj = (\sum Dist(sj,sk))/N-1}\cr
where:\cr
\eqn{W_{j}}{Wj} is the weight of sequence \eqn{s_{j}}{sj}, and is defined as the average evolutionary
distance between \eqn{s_{j}}{sj} and all other sequences in the alignment\cr
\eqn{N} is the number of sequences in the alignment.

\deqn{Dist(s_{j},s_{k})) = 1 = \frac{\sum_{i\epsilon Aligned_{jk}}Mut(s_{j},s_{k}))}{n(Aligned_{jk}))}}{Dist(sj,sk)) = 1 = {\sum Mut(sj,sk)/n(Alignedjk))}}
where:\cr
\eqn{Dist(s_{j},s_{k})}{Dist(sj,sk)}, the evolutionary distance between sequences \eqn{s_{j}}{sj} and \eqn{s_{k}}{sk}\cr
\eqn{Aligned_{jk}}{Alignedjk} is the set of all non-gap positions in \eqn{s_{j}}{sj} or \eqn{s_{k}}{sk}, \eqn{n(Aligned_{jk})}{n(Alignedjk)} is the number of such positions.

\deqn{Mut(a,b) = \frac{m(a,b) - min(m)}{max(m) - min(m)}}{Mut(a,b) = (m(a,b) - min(m))/(max(m) - min(m))}
where:\cr
\eqn{Mut(a,b)} measures the similarity between amino acids \eqn{a} and \eqn{b} as derived from \eqn{a} mutation data matrix \eqn{m}
}
\value{
A vector with weights of length equal to the number of sequences in the alignment
}
\references{
Valdar, W. S. J. & Thornton, J. M. Protein–protein interfaces: Analysis of amino acid conservation in homodimers. Proteins: Structure, Function, and Bioinformatics 42, 108–124 (2001).

}
\author{
Michal Stolarczyk & Alicja Pluciennik
}

\examples{
data("small_alignment")
alignment = small_alignment
weights = get_seq_weights(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
