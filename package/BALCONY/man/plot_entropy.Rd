\name{plot_entropy}
\alias{plot_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot entropies for protein
}
\description{
This function plots entropies of protein. Plots might be superimposed or not.
}
\usage{
plot_entropy(protein_conservation, colors, impose = NULL,
             prot_name = NULL, legend_pos = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_conservation}{A list or a vectors of protein conservation/entropies. The output of \code{\link{get_prot_entropy}} function}
  \item{colors}{(optional) A vector of colors for each plot, default: rainbow}
  \item{impose}{(optional) A boolean, if True/T plots are superimposed, if False/F plots are printed separately, default: T}
  \item{prot_name}{(optional) A string with structure name (to be used in the tile of the plot), default: none}
  \item{legend_pos}{(optional) A string witch legend position - one of following: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default: "bottomleft" }
}
\details{
This function produces plots for given values, on X axis are amino acids, on Y axis are values of entropy/conservation. Legend contains score names for description values.
}
\value{
This function produces plots
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
data("alignment")
data("structure")
uniprot="P34914"
structure_index=get_structures_idx(structure)
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
                        TG_entropy = TG_conservativity(alignment))
prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)

plot_entropy(prot_entropy, colors = c("red","green","blue"),
            impose = TRUE, prot_name = "Murine Epoxide Hydrolase",
            legend_pos = "bottomright")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ entropy }% __ONLY ONE__ keyword per line
