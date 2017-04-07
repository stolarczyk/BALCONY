\name{create_structure_seq}
\alias{create_structure_seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Superimpose structural data of interest on sequence after the alignmment
}
\description{
Create sequence of a protein structure model based on numbers of amino acids given in a text file (list of IDs and numbers in protein)
}
\usage{
create_structure_seq(structure_list, sequence_id, alignment,
                     pdb_path = NULL, chain_identifier = NULL, shift = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_list}{A list of structure data used for further evolutionary analysis. It can be text file(s) read by the \code{\link{read_structure}} function (text file with 2 columns: numbers of amino acids and 3-letters codes of AA; First row needs to contain markers)}
  \item{sequence_id}{The id/name of the target sequence in alignment which will be a base of structure sequence}
  \item{alignment}{An alignment object read with \code{\link[seqinr]{read.alignment}} function, must contain the target sequence}
  \item{pdb_path}{A string specifying the path to the PDB file with structural information. Optional parameter, required if the structure is incomplete e.g. fragments such as loops are missing}
  \item{chain_identifier}{A character specifying the chain of interest e.g. "A" or "B"}
  \item{shift}{A numeric value. In case there is a need to adjust the amino acids numeration due to missing amino acids at the beginning of the structure (that are not considered in the PDB file REMARK465 section)}
}
\details{
This function is useful to create sequence covered with structural data provided in a .txt file. This sequence can be compared with alignment to check the conservation for interesting amino acid(s). Additionally, if path to the PDB file is provided the function corrects the output accordingly to the information in REMARK465 on missing amino acids.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
\item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
\item{structure_probabilities}{A matrix of numeric values: probabilities of corresponding to the structural information from first element of the output, which helps to reduce the effect of non-consistent structural amino acids on the conservativity analysis of the structure of interest}
}

\author{
Alicja Pluciennik & Michal Stolarczyk}

\seealso{
\code{\link{get_remarks465_pdb}}, \code{\link{find_consecutive_seq}}, \code{\link{read_structure}}, \code{\link[seqinr]{read.alignment}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{protein structure}
\keyword{incomplete structure correction}
\keyword{PDB file}

