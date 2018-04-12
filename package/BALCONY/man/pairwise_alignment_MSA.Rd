\name{pairwise_alignment_MSA}
\alias{pairwise_alignment_MSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate pairwise alignment for whole MSA
}
\description{
For given alignment calculate pariwise alignments and returns alignment score.
}
\usage{
pairwise_alignment_MSA(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{An alignment object read with \code{\link[seqinr]{read.alignment}} function}
}

\value{
 \item{score_mtx }{Matrix of alignment scores}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Michal Stolarczyk & Alicja Pluciennik
}

\examples{
data("small_alignment")
\donttest{pairwiseAlignemnt_scores=pairwise_alignment_MSA(small_alignment)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pairwiseAlignment }% use one of  RShowDoc("KEYWORDS")

