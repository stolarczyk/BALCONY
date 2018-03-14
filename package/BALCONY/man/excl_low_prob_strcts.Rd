\name{excl_low_prob_strcts}
\alias{excl_low_prob_strcts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exclude low probability structural data
}
\description{
This function facilitates the exclusion of low probability structural data from the downstream conservativity analysis, which helps to reduce the effect of non-consistent structural amino acids on the conservativity analysis of the structure of interest
}
\usage{
excl_low_prob_strcts(structure, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{A structure object generated with \code{\link{create_structure_seq}} function}
  \item{threshold}{The threshold for the structural data exclusion}
}

\value{
\item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
\item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
\item{structure_probabilities}{A matrix of numeric values: probabilities of corresponding to the structural information from first element of the output}
}
\author{
Michal Stolarczyk & Alicja Pluciennik
}

\seealso{
\code{\link{create_structure_seq}}
}
\examples{
data("alignment")
structure_files = c(system.file("extdata", "T1_4JNC.structure", package = "BALCONY"),
                    system.file("extdata", "T2_4JNC.structure", package = "BALCONY"),
                    system.file("extdata", "T3_4JNC.structure", package = "BALCONY")
                    )
structure_list = read_structure(structure_files)
#creating library uniprot - PDB
lib=list(c("Q84HB8","4I19","4QA9"),
        c("P34913","4JNC"),
        c("P34914","1EK2","1CR6","1EK1","1CQZ"))
pdb_name = "4JNC"
uniprot=find_seqid(pdb_name,lib)
tunnel=create_structure_seq(structure_list,uniprot,alignment)
tunnel_excluded = excl_low_prob_strcts(tunnel, 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ structure }% use one of  RShowDoc("KEYWORDS")
