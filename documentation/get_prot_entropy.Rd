\name{get_prot_entropy}
\alias{get_prot_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get MSA-based calculated entropy for choosen protein.
}
\description{
This function allows to obtain vector of entropies for one complete protein sequence from MSA (gaps introduced in alignment are omitted)
}
\usage{
ProtConservation=get_prot_entropy(whole_prot,entropyScores_list,score_descript)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{whole_prot}{indices of given protein aminoacids in aligned sequence}
  \item{entropyScores_list}{list of entropy scores calculated for MSA }
  \item{score_descript}{list of score identifiers (names) for each score from entropyScores_list}
}
\details{
This function can be used on list of entropies or on vector for one entropy score.
}
\value{

  \item{entropy }{vector or a list of entropy values for further amino acids in given protein sequence }
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
FIX ME!!!!!!!!!!!!!!!!!!!!!
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
