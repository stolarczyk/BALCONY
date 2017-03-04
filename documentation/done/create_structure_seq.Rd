\name{create_structure_seq}
\alias{create_structure_seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create aligned sequence of examined structure in protein.
}
\description{
Create sequence of a protein structure model based on numbers of amino acids given in a text file (list of IDs and numbers in protein).
}
\usage{
create_structure_seq(structure_list, sequence_id, alignment, shift)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure_list}{A list of strucure data used for furhter evolutionary analysis. It can be text file(s) read by the read.table() function (text file with 2 columns: numbers of amino acids and 3-letters codes of AA; First row needs to contain markers)}
  \item{sequence_id}{The id/name of the target sequence in alignment which will be a base of structure sequence}
  \item{alignment}{ alignment read with read.alignment() function , must contain the target sequence}
  \item{shift}{Sequence shift. In case of lack of amino acid(s) at the beginning of the crystal structure sequence (PDB) compared with one used in sequences alignment the shift must be provided.}
}
\details{
This function is useful to create sequence covered with structure element. This sequence can be compared with aligment to check the conservaion for interesting amino acid(s).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{structure$structure_matrix}{A matrix of characters T and N marking on sequence the structural element; T-amino acid which forms analyzed structure, N-amino acid which not build a structure. Number of rows of the matrix corresponds to the number of structures analyzed.}
\item{structure$structure_numbers}{A list of vector containing the numbers of just the amino acids in the aligned sequences (no gaps)}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Pluciennik & Michal Stolarczyk}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(alignment)
load(structure_data)
#creating library uniprot - PDB
lib=list(c("Q84HB8","4I19","4QA9"),
        c("P34913","4JNC"),
        c("P34914","1EK2","1CR6","1EK1","1CQZ"))
uniprot=find_seqid(pdb_name,lib)
shift=3
tunnel=create_structure_seq(structure_data,uniprot,alignment,shift)
tunnel_numbers = tunnel$structure_numbers
tunnel_matrix = tunnel$structure_matrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{protein structure }

