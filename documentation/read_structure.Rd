\name{read_structure}
\alias{read_structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read structure from text file
}
\description{
Using this function you can read text file and create  aligned structure quasi-sequence.Text file schould contain 2 lines. In the first line are placed IDs (separated with single spaces). The scecond line contains one letter amino acid codes (separated with single space), respectively.
\usage{

structure_loop=read_structure(structure_file, sequence_id, alignment_file, shift)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_file}{filename eg. "1CQZ_loop.txt"}
  \item{sequence_id}{id of reference protein-structure sequence used in MSA}
  \item{alignment_file}{read in alignment file}
  \item{shift}{ between reference sequence and sequence known from crystal structure}
}
\details{
The shift between reference structure and crystal structure is easy to check using sequence alignment for the pair of sequences.This must be checked before analysis. 
}
\value{

  \item{strucure_sequence }{an aligned quasi sequence of structure with an extra first element (description) what provide compatibility with other formants in this package}

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
FIX ME!
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
\keyword{ CREATE }
\keyword{ structure }% __ONLY ONE__ keyword per line
