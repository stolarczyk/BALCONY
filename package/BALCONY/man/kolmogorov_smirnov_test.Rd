\name{kolmogorov_smirnov_test}
\alias{kolmogorov_smirnov_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Kolmogorov-Smirnov test for structural data
}
\description{
This function facilitates the comparison of conservativity of structure of interest with the rest of the protein. For example comparison of tunnel conservativity with overall protein conservativity.
}
\usage{
kolmogorov_smirnov_test(protein_entropy, structure_entropy, alternative,
                        pdb_name = "Reference", range = NULL, make_plot = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_entropy}{A list of calculated entropy scores (vectors of numeric values). The output of \code{\link{get_prot_entropy}} function}
  \item{structure_entropy}{A list where each element is a list of structure indices in the protein and matrix with corresponding entropy values (each row is a separate score metric)}
  \item{alternative}{A numeric value indicating the character of alternative hypothesis of the test to be performed: 1 - two sided test, 2 - less, 3 - greater, following the generic \code{\link[stats]{ks.test}} function.}
  \item{pdb_name}{(optional) A string with name of the reference protein, default: "Reference"}
  \item{range}{(optional) A numeric vector with region of reference protein to be excluded from the data set. Useful when protein consists of additional chains with outstandingly low/high entropy values which may distort result of the test, default: NULL}
  \item{make_plot}{(optional) A logical indicating if cumulative distribution functions should be displayed, default: TRUE}
  }
\value{A matrix of p-values for each entropy metric (rows) and structure (columns)}
\author{
Michal Stolarczyk & Alicja Pluciennik
}
\seealso{
\code{\link[stats]{ecdf}}, \code{\link[stats]{ks.test}}
}

\examples{
data("alignment")
data("structure")
entropy_data=list(Schneider.entropy=schneider_conservativity(alignment),
              Escore.entropy = Escore_conservativity(alignment),
              Kabat.entropy =  kabat_conservativity(alignment))
indices=get_structures_idx(structure)
protein_index = indices$proteinIndices
structure_index = indices$structureIndices
prot_cons=get_prot_entropy(protein_index,entropy_data)
stru_entropy=get_structures_entropy(structure_index,entropy_data)
profiles_for_structure=prepare_structure_profile(structure, stru_entropy)
EQUAL=kolmogorov_smirnov_test(protein_entropy = prot_cons,
                              structure_entropy = profiles_for_structure,
                              alternative = 1,range = c(1:233),make_plot = TRUE)
}

\keyword{ test }% use one of  RShowDoc("KEYWORDS")
\keyword{ structure }% __ONLY ONE__ keyword per line
