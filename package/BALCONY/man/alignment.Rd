\name{alignment}
\alias{alignment}
\docType{data}
\title{
Sample alignment of soluable epoxide hydrolase family
}
\description{
This alignment consists of 52 proteins which belong to the soluable epoxide hydrolase family. The amino acid sequences were aligned using MUSCLE algorithm with default settings.
}

\format{
  An alignment object read with read.alignment() function from seqinr package.
  \describe{
    \item{alignment$nb}{A numeric: number of sequences}
    \item{alignment$nam}{A vector of characters: names of the sequences}
    \item{alignment$seq}{A vector of characters: amino acid sequences}
  }
}
\examples{
data("alignment")
alignment
}
\keyword{datasets}
