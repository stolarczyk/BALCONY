\name{structure}
\alias{structure}
\docType{data}
\title{
A sample structure data
}
\description{
This sample structure data consists of the amino acids names forming tunnels and their numbers is analyzed protein. The data is a result of CAVER which is a software tool for analysis and visualization of tunnels and channels in protein structures.
}
\format{
  A structure object with three elements:
  \describe{
\item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
\item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
\item{structure_probabilities}{A matrix of numeric values: probabilities of corresponding to the structural information from first element of the output, which helps to reduce the effect of non-consistent structural amino acids on the conservativity analysis of the structure of interest}
  }
}
\details{
The tunnel analysis with CAVER was performed on human epoxide hydrolase structure (PDB ID: 4JNC) 50ns MD simulation.
}
\seealso{
CAVER: \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002708}\cr
4JNC: \url{http://www.sciencedirect.com/science/article/pii/S0960894X13004885}
}
\examples{
data("structure")
structure
}
\keyword{datasets}
