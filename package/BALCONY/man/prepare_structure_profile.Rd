\name{prepare_structure_profile}
\alias{prepare_structure_profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function combines the entropy data for structure building amino acids with their indices
}
\description{
This function combines the entropy data for structure building amino acids with its indices. It prepares the data for convenient visualization or processing.
}
\usage{
structure_profile = prepare_structure_profile(structure, structure_entropy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{A structure data, the output of create_structure_seq() function}
  \item{structure_entropy}{The entropy values for the structure building residues, the output of get_structures_entropy() function}
}

\value{
\item{List of structures}{Each element is a list of entropy values (matrix of entropy scores) and indices of residues building structure in protein of interest.}
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
load(alignment)
load(structure)
uniprot="P34914"
pdb_path = "1CQZ" #A string with path to PDB file
chain_identifier = "B"
tunnel=create_structure_seq(structure,uniprot,alignment,pdb_path,chain_identifier)
structure_index=get_structures_idx(tunnel)
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),TG_entropy = TG_conservativity(alignment))
structure_entropy=get_structures_entropy(structure_index, entropy_score_list)

structure_profile = prepare_structure_profile(tunnel, structure_entropy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ profile }% use one of  RShowDoc("KEYWORDS")
\keyword{ structure }% __ONLY ONE__ keyword per line
