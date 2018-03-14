\name{plot_structure_on_protein}
\alias{plot_structure_on_protein}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot structure entropy on protein background
}
\description{
This function enables to visually asses the stucture(s) entropy in comparison with protein's entropy
}
\usage{
plot_structure_on_protein(protein_entropy, structure_profiles,
                          pdb_name, colors, structure_names, legend_pos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_entropy}{A list of entropy values for protein of interest. Output of \code{\link{get_prot_entropy}} function}
  \item{structure_profiles}{Output of \code{\link{prepare_structure_profile}} function}
  \item{pdb_name}{(optional) A string with protein's name e.g. PDB ID, default: none}
  \item{colors}{(optional) A vector of strings with colors to be used to plot the stucture markers of length equal to number of structures in structure profiles, default: rainbow()}
  \item{structure_names}{(optional) A vector of strings to be displayed as names in the legend, default: "stru <no>"}
  \item{legend_pos}{(optional) A string witch legend position - one of following: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default: "bottomleft" }
}
\details{
For each entropy score from structure_profiles (these must correspond to prot_entropy) this function plots separate plots. Each plot presents entropy score for whole protein each structure is marked as one of 21 symbols available in generic \code{\link[graphics]{plot}} function.
}
\value{
This function produces plot
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}



\examples{
data("alignment")
data("structure")
indices=get_structures_idx(structure)
protein_index = indices$proteinIndices
structure_index = indices$structureIndices
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
                        Escore_entropy = Escore_conservativity(alignment))
structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)
structure_profile = prepare_structure_profile(structure, structure_entropy)
prot_entropy=get_prot_entropy(protein_index, entropy_scores_list)

plot_structure_on_protein(prot_entropy, structure_profile)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }% use one of  RShowDoc("KEYWORDS")
\keyword{ profile }% __ONLY ONE__ keyword per line
