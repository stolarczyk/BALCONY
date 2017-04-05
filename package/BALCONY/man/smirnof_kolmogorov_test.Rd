\name{smirnof_kolmogorov_test}
\alias{smirnof_kolmogorov_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Smirnof-Kolmogorov test for structural data
}
\description{
This function facilitates the comparion of conservativity of structure of interest with the rest of the protein. For example comparison of tunnel conservativity with overall protein conservativity.
}
\usage{
smirnof_kolmogorov_test(protein_cons, structure_cons, alternative,
                        pdb_name = "Reference", range = NULL, make_plot = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_cons}{A list of calculated entropy scores (vectors of numeric values). The output of \code{\link{get_prot_entropy}} function}
  \item{structure_cons}{A list where each element is a list of structure indices in the protein and matrix with corresponding entropy values (each row is a separate score metric)}
  \item{alternative}{A numeric value indicating the character of alternative hypothesis of the test to be performed: 1 - two sided test, 2 - less, 3 - greater, following the generic \code{\link[stats]{ks.test}} function.}
  \item{pdb_name}{(optional) A string with name of the reference protein, default: "Reference"}
  \item{range}{(optional) A numeric vector with region of reference protein to be excluded from the dataset. Useful when protein consists of additional chains with outstandingly low/high entropy values which may distort result of the test, default: NULL}
  \item{make_plot}{(optional) A logical indicating if cumulative distribution functions should be displayed, default: TRUE}
  }
\value{A matrix of p-values for each entropy metric (rows) and structure (columns)}
\author{
Michal Stolarczyk & Alicja Pluciennik
}
\note{
\code{\link[stats]{ecdf}}, \code{\link[stats]{ks.test}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("alignment")
data("structure")
entropy_data=list(Schneider.entropy=schneider_conservativity(alignment),
              TG.entropy = TG_conservativity(alignment),
              Kabat.entropy =  kabat_conservativity(alignment))
tunnel_index=get_structures_idx(structure)
prot_cons=get_prot_entropy(tunnel_index[[1]],entropy_data)
stru_entropy=get_structures_entropy(tunnel_index,entropy_data)
profiles_for_structure=prepare_structure_profile(structure, stru_entropy)
EQUAL=smirnof_kolmogorov_test(protein_cons = prot_cons,
                              structure_cons = profiles_for_structure,
                              alternative = 1,range = c(1:233),make_plot = TRUE)
}

\keyword{ test }% use one of  RShowDoc("KEYWORDS")
\keyword{ structure }% __ONLY ONE__ keyword per line
