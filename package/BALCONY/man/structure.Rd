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
  A list of three lists with two elements each.
  \describe{
    \item{T1_mysz_1cqz}{A list of two elements: numbers of amino acids and names of amino acids. Tunnel T1}
    \item{T2_mysz_1cqz}{A list of two elements: numbers of amino acids and names of amino acids. Tunnel T2}
    \item{T3_mysz_1cqz}{A list of two elements: numbers of amino acids and names of amino acids. Tunnel T3}
  }
}
\details{
The tunnel analysis with CAVER was performed on murine epoxide hydrolase structure (PDB ID: 1CQZ)
}
\seealso{
\url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002708}
}
\examples{
load(structure)
}
\keyword{datasets}
