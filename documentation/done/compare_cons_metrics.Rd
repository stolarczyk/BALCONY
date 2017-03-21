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
compare_cons_metrics(prot_entropy, structure_profile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_entropy}{List of entropy scores values for a whole protein (output of get_prot_entropy()).}
  \item{structure_profile}{Each element is a list of entropy values (matrix of entropy scores) and indices of residues building structure in protein of interest  (output of prepare_structure_profile()).}
  \item{pdb_name}{The name of the analyzed protein.}
%%     ~~Describe \code{x} here~~
}
}
\details{
This function allows to show the scatterplots of an entropy scores. The protein is marked as gray points, the structures are marked with symbols. It is useful to visualise differences between entropy scores, and choose the best one for further analysis.}
\value{

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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),TG_entropy = TG_conservativity(alignment))
prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)
compare_cons_metrics(prot_entropy, structure_profile, "1CQZ")

}

