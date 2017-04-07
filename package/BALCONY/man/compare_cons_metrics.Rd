\name{compare_cons_metrics}
\alias{compare_cons_metrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compare_cons_metrics
}
\description{
This function is designed to compare the conservation metrics used in the analysis. This way the user can notice the significant correlation or differences between these to evaluate their performance in a specific case.
}
\usage{
compare_cons_metrics(protein_entropy, structure_profile, pdb_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_entropy}{List of entropy scores values for a whole protein (output of \code{\link{get_prot_entropy}}).}
  \item{structure_profile}{Each element is a list of entropy values (matrix of entropy scores) and indices of residues building structure in protein of interest  (output of \code{\link{prepare_structure_profile}}).}
  \item{pdb_name}{The name of the analyzed protein.}
%%     ~~Describe \code{x} here~~
}
\details{
This function allows to show the scatterplots of an entropy scores. The protein is marked as gray points, the structures are marked with symbols. It is useful to visualise differences between entropy scores, and choose the best one for further analysis.}
\value{
This function produces a set of scatter plots facilitating the visual inspection of entropy metrics dependancies.
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
data("alignment")
alignment = delete_isoforms(alignment)
data("structure")
uniprot="P34913"
#pdb_path = system.file("extdata", "1cqz.pdb", package = "BALCONY") #A string with path to PDB file
#chain_identifier = "B"
#tunnel=create_structure_seq(structure,uniprot,alignment,pdb_path,chain_identifier)
indices=get_structures_idx(structure)
protein_index = indices$proteinIndices
structure_index = indices$structureIndices
structure=create_structure_seq(structure,uniprot,alignment)
entropy_scores_list=list(
                        Schneider_entropy = schneider_conservativity(alignment),
                        TG_entropy = TG_conservativity(alignment)
                        )
structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)

structure_profile = prepare_structure_profile(structure, structure_entropy)
protein_entropy=get_prot_entropy(protein_index, entropy_scores_list)
compare_cons_metrics(protein_entropy, structure_profile, "1CQZ")

}

