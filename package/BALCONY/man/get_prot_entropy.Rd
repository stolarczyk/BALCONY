\name{get_prot_entropy}
\alias{get_prot_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get MSA-based calculated entropy for chosen protein.
}
\description{
This function allows to obtain vector of entropies for one complete protein sequence from MSA (gaps introduced in alignment are omitted)
}
\usage{
get_prot_entropy(protein_index, score_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_index}{Indices of given protein aminoacids in aligned sequence}
  \item{score_list}{A list of entropy scores calculated for MSA }
}
\details{
This function can be used on list of entropies or list with one element for one entropy score.
}
\value{

  \item{entropy}{A list where each element is a vector of entropy values provided in entropy_scores_list}
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
data("structure")
data("alignment")
pdb_name = "1CQZ" #A string with path to PDB file
uniprot="P43914"
chain_identifier = "B"
structure_index=get_structures_idx(structure)
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
                        TG_entropy = TG_conservativity(alignment))
prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)

# In case of one entropy score
entropy_scores_list = list()
entropy_scores_list[[1]] = Schneider_entropy = schneider_conservativity(alignment)
prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ get }
\keyword{ entropy }% __ONLY ONE__ keyword per line
