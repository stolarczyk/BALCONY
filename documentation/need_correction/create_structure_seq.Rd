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
create_structure_seq(tunnel_file, sequence_id, alignment_file, shift)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tunnel_file}{ text file with 2 columns: numbers of amino acids and 3-letters codes of AA; First row needs to contain markers}
  \item{sequence_id}{ id/name of sequence in alignment which will be a base of tunnel sequence}
  \item{alignment_file}{ alignment, must contain target sequence}
  \item{shift}{INFO NEEDED!}
}
\details{
This function is useful to create sequence covered with structure element. This sequence can be compared with aligment to check conservaion rate for interesting amino acid.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{tunel }{a vector of characters T and N to mark on sequence structural element; T-amino acid which forms analyzed structure, N-amino acid which not build a structure}
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
alignment=read.alignment(file="test.mfa", format="fasta", forceToLower=F)
#creating library uniprot - PDB
lib=list(c("Q6U6J0","4QLA"), 
         c("Q9UR30","3G0I","3G02","1QO7"))
sequence_id=find_UPid(3G0I,lib)
my_seq=find_uniprot_seq(sequence_id, alignment_file)
tunnel_file=read.table(tunel1.txt) #??
tunnel=create_tunel_seq(tunnel_file,uniprot,file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ structure }

