\name{get_structures_idx}
\alias{get_structures_idx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get IDs of structure(s) elements from aligned sequences (MSA)
}
\description{
This function allows to obtain positions in aligned sequences for analyzed structure (e.g. functionally related amino acids dispersed in sequence) based on sequence corresponding to the crystal structure.
}
\usage{
get_structures_idx(structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{
  The output of create_structure_seq() function
}
}
\details{
It facilitates the management and oparation on the entropy values calculated for given MSA.
}
\value{
Output is a list of two elements:
  \item{proteinIndices}{A sorted vector of amino acids of analyzed sequence in MSA}
  \item{strucureIndices}{A list of sorted vectors of amino acids indices in aligned sequence for each structure}
%% ...
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\examples{
data("structure")

#creating library uniprot - PDB
lib=list(c("Q84HB8","4I19","4QA9"),
        c("P34913","4JNC"),
        c("P34914","1EK2","1CR6","1EK1","1CQZ"))
pdb_name = "1CQZ" #A string with path to PDB file
uniprot=find_seqid(pdb_name,lib)
tunnel=create_structure_seq(structure,uniprot,alignment)
structure_index=get_structures_idx(tunnel)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ structure }
\keyword{ indices }% __ONLY ONE__ keyword per line
