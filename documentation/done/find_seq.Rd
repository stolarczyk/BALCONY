\name{find_seq}
\alias{find_seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find sequence by id in alignment.
}
\description{
This function allows to search for a sequence with its id. Useful for browsing a larg multiple sequence alignment data or for automatization purposes.
}
\usage{
find_seq(sequence_id, alignment_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence_id}{identifier of desired sequence from alignment}
  \item{alignment_file}{alignment file loaded with read.alignment() }
}
\details{

}
\value{
 \item{sequence}{desired aligned sequence form alignment}

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
load(alignment)
#creating library uniprot - PDB
lib=list(c("Q6U6J0","4QLA"), 
         c("Q9UR30","3G0I","3G02","1QO7"))
sequence_id=find_seqid("4QLA",lib)
sequence=find_seq(sequence_id, alignment)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sequence}
\keyword{ alignment}% __ONLY ONE__ keyword per line
