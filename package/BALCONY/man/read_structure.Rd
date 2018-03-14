\name{read_structure}
\alias{read_structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read structure data from a text file
}
\description{
By using this function you can read text file and create an structure list which can be used in further evolutionary analysis with BALCONY package. Text file should comprise 2 or 3 columns: first one should contain indices (positions) of amino acids in the protein, the second one should contain amino acid symbols on specified positions and the third one (optionally) the numeric property of given residue. }
\usage{
read_structure(file_names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_names}{A vector of strings with structure file(s) name(s)}
}
\details{
The files should be formatted as follows:\cr

2  ASP 100\cr
6   TYR 80\cr
11  PHE 30\cr
6   TYR 30
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
writeLines(c("2 TYR 100","3 LEU 100", "7 VAL 50", "10 PHE 30", "20 SER 20"), fileConn)
close(fileConn)
fileConn<-file("exemplary_input2.txt")
writeLines(c("5 ALA 100","6 ILE 100", "18 GLY 100", "40 PHE 100"), fileConn)
close(fileConn)

structure_list = read_structure(file_names = c("exemplary_input1.txt", "exemplary_input2.txt"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read }
\keyword{ structure file }% __ONLY ONE__ keyword per line
