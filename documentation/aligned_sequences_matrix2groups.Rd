\name{aligned_sequences_matrix2groups}
\alias{aligned_sequences_matrix2groups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert amino acid symbols to groups according to their properties of user's choice.
}
\description{
This function performs a conversion of amino acid symbols to group symbols according to their properties. Implemented grouping methods are: general (majority of properties taken into account), hydrophobicity, size and aromaticity. "GX", where X stands for group number, are group symbols.  
}
\usage{
grouped_aligned_sequences_matrix = aligned_sequences_matrix2groups(aligned_sequences_matrix,grouping_method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aligned_sequneces_matrix}{A matrix that contains aligned sequences. It is an output of alignment2matrix() function}
  \item{grouping_method}{A string which specifies the grouping method to be used. One of following: 'general', 'hydrophobicity', 'size', 'aromaticity' }
}
\details{
}
\value{grouped_aligned_sequences_matrix}{A matrix of size of the input matrix but with group symbols instead of amino acid symbols}
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
alignment2matrix()
}
\examples{
FIXME!!!!!!!!!!!!!!!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
