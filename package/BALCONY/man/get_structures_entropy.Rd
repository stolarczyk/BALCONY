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
get_structures_entropy(structure_index,score_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_index}{A is a list of indices in alignment of protein and structures. Output output of \code{\link{get_structures_idx}} function}
  \item{score_list}{A list of entropies for whole alignment}
}
\details{
This function allows to obtain entropy (calculated on MSA) for dispersed amino acids in protein e.g. surface, binding site, tunnels etc. The input is a list of few structure indices in given protein sequence. Function calculates position of those in aligned sequence and returns a vector/matrix or a list of matrices with entropy values.
}
\value{
 \item{structure_entropies}{A list of matrices. Rows are entropy scores, columns are}
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\seealso{
\code{\link{create_structure_seq}}, \code{\link{read_structure}}}
\examples{
data("structure")
data("alignment")

#creating library uniprot - PDB
uniprot="P34914"
tunnel=create_structure_seq(structure,uniprot,alignment)
indices=get_structures_idx(structure)
protein_index = indices$proteinIndices
structure_index = indices$structureIndices
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
                        Escore_entropy = Escore_conservativity(alignment))
structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ structure }
\keyword{ entropy }% __ONLY ONE__ keyword per line
