\name{get_structure_idx}
\alias{get_structure_idx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get IDs of structure elements from aligned sequences (MSA)
}
\description{
This function allows to obtain positions in aligned sequences for analyzed structure(eg. functionally connected amino acids dispersed in sequence) based on sequence corresponding to crystal structure.
}
\usage{
structure_loop=read_structure("1CQZ_wt_petla.txt", uniprot,file,3)
structure_index=get_structure_idx(structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{
structure - the aligned quasi-sequence of protein (contains gaps "-",and marked aminoacids which are part of crystal structure "T" or not "N")
}
}
\details{
It facilitate to menage and operate on the entropy values calculated for given MSA
}
\value{
Output is a list of two elements:
  \item{$proteinIndices}{ sorted vector of structure amino acids in MSA}
  \item{$strucureIndices}{list of sorted vectors of amino acids in aligned sequence for each structure }
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
FIXME
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
