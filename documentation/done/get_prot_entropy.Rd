\name{get_prot_entropy}
\alias{get_prot_entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get MSA-based entropy for chosen protein.
}
\description{
This function allows to obtain vector of entropies for one complete protein sequence from MSA (gaps introduced in alignment are omitted)
}
\usage{
prot_entropy=get_prot_entropy(whole_prot, entropy_scores_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{whole_prot}{Indices of given protein amino acids in aligned sequence}
  \item{entropy_scores_list}{A list of entropy scores calculated for MSA }
}
\details{
This function can be used on list of entropies or list with one element for one entropy score.
}
\value{

  \item{entropy}{A list where each element is a vector of entropy values provided in entropy_scores_list}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
load(structure)
load(alignment)
pdb_name = "1CQZ" #A string with path to PDB file
uniprot="P43914"
chain_identifier = "B"
tunnel=create_structure_seq(structure,uniprot,alignment,pdb_path,chain_identifier)
structure_index=get_structures_idx(tunnel)
entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),TG_entropy = TG_conservativity(alignment))
prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)

# In case of one entropy score
entropy_scores_list = list()
entropy_scores_list[[1]] = Schneider_entropy = schneider_conservativity(alignment)
prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
