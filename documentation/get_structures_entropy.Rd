\name{get_structures_entropy}
\alias{get_structures_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get entropy/conservation for previously choosen amino acids (which for structure/region of interest) in given protein.
}
\description{
This function allows to get values of entropy/conservation for amino acids dispersed in sequence of given protein. It works well with list of dispersed amino acids in one protein.
}
\usage{
EntropyStructures=get_structures_entropy(structure_index,Score_list,names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_index}{A is a list of indices in alignment of protein and structures}
  \item{score_list}{A list of entropies for whole alignment}
}
\details{
Function allows to obtain entropy (calculated on MSA) for dispersed amino acids in protein e.g. surface, binding site, tunnels, etc. The input might be a list of few structure indices in given protein sequence. Function calculates position of those in aligned sequence and returns vector/matrix or a list of matrices with entropy values. 
}
\value{

 \item{entropies for structures}{list of matrices(one entropy score for each row)}

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
create_structure_seq(), read_structure()
\examples{
FIXME!!!!!!!!!!!!!!!
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
