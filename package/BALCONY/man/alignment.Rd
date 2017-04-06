\name{alignment}
\alias{alignment}
\docType{data}
\title{
Sample alignment of soluable epoxide hydrolase family
}
\description{
was performed on a dataset comprising 311 soluble epoxide hydrolase peptide orthologous sequences acquired from UniProtKB. The alignment was performed and edited with MUSCLE algorithm in JALVIEW, respectively.
}

\format{
  An alignment object read with \code{\link[seqinr]{read.alignment}} function from seqinr package.
  \describe{
    \item{alignment$nb}{A numeric: number of sequences}
    \item{alignment$nam}{A vector of characters: names of the sequences}
    \item{alignment$seq}{A vector of characters: amino acid sequences}
  }
}
\references{
MUSCLE: \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113}\cr
JALVIEW: \url{https://academic.oup.com/bioinformatics/article/25/9/1189/203460/Jalview-Version-2-a-multiple-sequence-alignment}
}
\examples{
data("alignment")
alignment
}
\keyword{datasets}

