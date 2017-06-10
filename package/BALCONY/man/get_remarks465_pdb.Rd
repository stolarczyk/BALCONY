\name{get_remarks465_pdb}
\alias{get_remarks465_pdb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get "REMARK 465" data from PDB file
}
\description{
This function extracts the data concerning missing amino acids in PDB protein structure from the PDB file
}
\usage{
get_remarks465_pdb(pdb_path, chain_identifier)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pdb_path}{A string specifying the path tp the PDB file}
  \item{chain_identifier}{A character specifying the chain to be considered}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{aa_numbers}{A numeric vector of indices of missing amino acids}
\item{chain}{A character specifying the chain which was considered in remark 465 data extraction}
%% ...
}

\author{
Michal Stolarczyk & Alicja Pluciennik
}
\seealso{
\code{\link[Rpdb]{read.pdb}}
}
\examples{
require(Rpdb)
chain_identifier = "A"
pdb_path = system.file("extdata", "4jnc.pdb", package = "BALCONY")
print(pdb_path)
#pdb_file_path = "path_to_file"
remark465_data = get_remarks465_pdb(pdb_path,chain_identifier)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{REMARK 465}% use one of  RShowDoc("KEYWORDS")
\keyword{missing amino acids}% __ONLY ONE__ keyword per line
