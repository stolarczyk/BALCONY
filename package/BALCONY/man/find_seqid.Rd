\name{find_seqid}
\alias{find_seqid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find sequence identifier by other sequence identifier in given alignment within a specified library
}
\description{
This function allows to find sequence id from alignment file corresponding to the given sequence id. Function requires library of equivalent sequences id defined by user and it is useful to find sequences from other databases in alignment for examined sequence from other database (like PDB sequence for structure and UniProt sequences in alignment).
%
}
\usage{
find_seqid(sequence_id, library)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence_id}{A string. An ID of e.g. PDB structure identifier}
  \item{library}{A list of vectors which contain a defined by user library e.g. of UniProt ids <-> PDB ids. See examples}
}

\value{

\item{seqid}{A string. The equivalent ID to the one provided as the input.}

}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
#creating library uniprot - PDB
lib=list(  c("Q84HB8","4I19","4QA9"),
  c("P34913","4JNC"),
  c("P34914","1EK2","1CR6","1EK1","1CQZ"))
PDB_name = "1CQZ"
find_seqid(PDB_name,lib)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ find }
\keyword{ seqid }% __ONLY ONE__ keyword per line
