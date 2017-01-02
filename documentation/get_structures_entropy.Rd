\name{get_structures_entropy}
\alias{get_structures_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get entropy/conservation for previously choosen amino acids in given protein.
}
\description{
This function allows to get vlaues of entropy/conservation for amino acids dispersed in sequence of given protein. It works well with list of dispersed amino acids in one protein.
}
\usage{
EntropyStructures=get_structures_entropy(structure_index,Score_list,names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_index}{ is a list of indexes in alignment of protein and structures}
  \item{Score_list}{ list of entropies for whole alignment}
  \item{names}{names of entropy scores}
}
\details{
Function allows to obtain entropy (calculated on MSA) for dispersed amino acids in protein (like surface, binding site, tunnels, etc). The input might be a list of few structure indexes in given protein sequence. Function calculates position of those residues in aligned sequence and returns vector/matrix or a list of matrices with entropy values. 
}
\value{

 \item{entropies for structures}{list of matrices(one entropy score for each row)}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Płuciennik & Michał Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
create_structure_seq(), read_structure()
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
