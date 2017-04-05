\name{small_alignment}
\alias{small_alignment}
\docType{data}
\title{
Sample small alignment of soluable epoxide hydrolase family
}
\description{
This alignment consists of 10 proteins which belong to the soluable epoxide hydrolase family. The amino acid sequences were aligned using MUSCLE algorithm with default settings.
}

\format{
  An alignment object read with \code{\link[seqinr]{read.alignment}} function from seqinr package.
  \describe{
    \item{alignment$nb}{A numeric: number of sequences}
    \item{alignment$nam}{A vector of characters: names of the sequences}
    \item{alignment$seq}{A vector of characters: amino acid sequences}
  }
}
\details{
This is a smaller version of sample alignment which facilitates faster presentation of the functions capabilities.
}
\examples{
data("small_alignment")
small_alignment
}
\keyword{datasets}
