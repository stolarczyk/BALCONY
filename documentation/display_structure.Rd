\name{display_structure}
\alias{display_structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display structure on alignment matrix
}
\description{
This function allows to display desired structures in matrix as an aligned sequences. It is useful for visualisation of amino acids conservation which create structures of interest and compare the values with the whole alignment background.
}
\usage{
display_structure(structure, structure_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{created structure sequences, the output of create_structure_seq() function}
  \item{structure_file}{list of structures, where in the first column of each element of the list are amino acid numbers in the protein and in the second column their symbols}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

{information about structures on aligned sequence, "-" means alignment gap, "N" means amino acid which do not create structure, "T" means amino acid forming the structure}

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
#reading file(s) containing information on tunnels
structure_file1=read.table("T1_1cqz.txt");
structure_file2=read.table("T2_1cqz.txt");
structure_file3=read.table("T3_1cqz.txt");
structure_file=list(structure_file1,structure_file2,structure_file3)
seqid="P34913"
structure=create_structure_seq(tunnel_file,seqid,file)
display_structure_matrix=display_structure(structure,structure_file)
}

