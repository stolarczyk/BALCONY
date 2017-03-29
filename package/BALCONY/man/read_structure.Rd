\name{read_structure}
\alias{read_structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read structure data from a text file
}
\description{
By using this function you can read text file and create an stucture list which can be used in furter evolutionary analysis with BALCONY package. Text file should comprise 2 columns: first one should contain indices (positions) of amino acids in the protein, the scecond one should contain amino acid symbols on specified positions.}
\usage{
structure = read_structure(file_names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_names}{A vector of strings with stucture file(s) name(s)}
}
\details{
The files should be formatted as follows:\cr

\code{
res AA\cr
2  ASP\cr
6   TYR\cr
11  PHE\cr
6   TYR}
}
\value{
  \item{strucure_list}{A list with read structure data. Number of elements of this list equals to the number of files specified.}

}

\author{
Alicja Pluciennik & Michal Stolarczyk
}

\examples{
#Generating exemplary input files for the function
fileConn<-file("exemplary_input1.txt")
writeLines(c("2 TYR","3 LEU", "7 VAL", "10 PHE", "20 SER"), fileConn)
close(fileConn)
fileConn<-file("exemplary_input2.txt")
writeLines(c("5 ALA","6 ILE", "18 GLY", "40 PHE"), fileConn)
close(fileConn)

structure_list = read_structure(file_names = c("exemplary_input1.txt", "exemplary_input2.txt"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read }
\keyword{ structure file }% __ONLY ONE__ keyword per line
