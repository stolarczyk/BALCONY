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
find_seq(sequence_id, alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence_id}{identifier of desired sequence from alignment}
  \item{alignment}{alignment file loaded with \code{\link[seqinr]{read.alignment}} }
}

\value{
 \item{sequence}{A string, the desired aligned sequence form alignment}

}
\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
data("alignment")
#creating library uniprot - PDB
lib=list(  c("Q84HB8","4I19","4QA9"),
  c("P34913","4JNC"),
  c("P34914","1EK2","1CR6","1EK1","1CQZ"))
sequence_id=find_seqid("1CQZ",lib)
sequence=find_seq(sequence_id, alignment)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sequence}
\keyword{ alignment}% __ONLY ONE__ keyword per line
