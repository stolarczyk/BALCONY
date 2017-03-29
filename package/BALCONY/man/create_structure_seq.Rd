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
create_structure_seq(structure_list, sequence_id, alignment, pdb_path = NULL, chain_identifier = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_list}{A list of strucure data used for further evolutionary analysis. It can be text file(s) read by the read.table() function (text file with 2 columns: numbers of amino acids and 3-letters codes of AA; First row needs to contain markers)}
  \item{sequence_id}{The id/name of the target sequence in alignment which will be a base of structure sequence}
  \item{alignment}{ alignment read with read.alignment() function , must contain the target sequence}
  \item{pdb_path}{A string specifying the path to the PDB file with structural information. Optional parameter, required if the structure is incomplete e.g. fragments such as loops are missing}
  \item{chain_identifier}{A character specyfiying the chain of interest e.g. "A" or "B"}
}
\details{
This function is useful to create sequence covered with structural data provided in a .txt file. This sequence can be compared with aligment to check the conservaion for interesting amino acid(s). Addionally, if path to the PDB file is provided the function corrects the output accrdingly to the information in REMARK465 on missing amino acids.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
\item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
}

\author{
Alicja Pluciennik & Michal Stolarczyk}

\seealso{
\code{\link{get_remarks465_pdb}}, \code{\link{find_consecutive_seq}}
}
\examples{
data("alignment")
#creating library uniprot - PDB
lib=list(c("Q84HB8","4I19","4QA9"),
        c("P34913","4JNC"),
        c("P34914","1EK2","1CR6","1EK1","1CQZ"))
pdb_name = "1CQZ"
uniprot=find_seqid(pdb_name,lib)
chain_identifier = "B"
tunnel=create_structure_seq(structure,uniprot,alignment)
tunnel_numbers = tunnel$structure_numbers
tunnel_matrix = tunnel$structure_matrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{protein structure}
\keyword{incomplete structure correction}
\keyword{PDB file}

