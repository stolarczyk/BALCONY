\name{find_seqid}
\alias{find_seqid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find sequence identifier by other sequence identifier in given alignment within a specified library
}
\description{
This function allows to find sequence id form alignment file corresponding to the given sequence id. Function requires library of equivalent sequences id defined by user and it is useful to find sequences from other databases in alignment for examined sequence from other database (like PDB sequence for structure and UniProt sequences in alignment).
%
}
\usage{
find_seqid(sequence_id, library)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence_id}{id of e.g. PDB sequence}
  \item{library}{defined by user library e.g. of UniProt ids and PDB ids}

}

\value{

\item{seqid}{function gives sequence name (character)}

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
