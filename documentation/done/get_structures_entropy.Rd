\name{get_structures_entropy}
\alias{get_structures_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get entropy of amino acids (for region of interest) in given protein
}
\description{
This function allows to get values of entropy/conservation for amino acids dispersed in sequence of given protein. It works well with a list of dispersed amino acids in one protein.
}
\usage{
structure_entropy=get_structures_entropy(structure_index,score_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_index}{A is a list of indices in alignment of protein and structures. Output output of get_structure_idx() function}
  \item{score_list}{A list of entropies for whole alignment}
}
\details{
This function allows to obtain entropy (calculated on MSA) for dispersed amino acids in protein e.g. surface, binding site, tunnels etc. The input is a list of few structure indices in given protein sequence. Function calculates position of those in aligned sequence and returns a vector/matrix or a list of matrices with entropy values. 
}
\value{
 \item{structure_entropies}{A list of matrices. Rows are entropy scores, columns are}
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
create_structure_seq(), read_structure()}
\examples{
load(structure)

#creating library uniprot - PDB
pdb_path = "1CQZ" #A string with path to PDB file
uniprot="P34914"
chain_identifier = "B"
tunnel=create_structure_seq(structure,uniprot,alignment,pdb_path,chain_identifier)
structure_index=get_structures_idx(tunnel)
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),TG_entropy = TG_conservativity(alignment))

structure_entropy=get_structures_entropy(structure_index, entropy_score_list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ structure }
\keyword{ entropy }% __ONLY ONE__ keyword per line
