\name{calculate_pseudo_counts}
\alias{calculate_pseudo_counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate pseudo counts for alignment }
\description{
This function calculates pseudo-counts (as shown in \href{https://doi.org/10.1093/bioinformatics/12.2.135}{Henikoff et al. (1996)}) for an alignment with the use of substitution matrices. It is recommended to estimate amino acid frequencies for alignments with small number of sequences (in order to calculate reliable entropy scores)
}
\usage{
calculate_pseudo_counts(alignment, substitution_mtx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{alignment loaded with \code{\link[seqinr]{read.alignment}}}
  \item{substitution_mtx}{Matrix with amino acids substitution frequencies. Default: GONNET}
}

\value{
 \item{pseudoCounts}{Matrix with pseudo counts of size 21x number of alignment columns}
 }
\references{
Henikoff et al.(1996) Using substitution probabilities to improve position-specific scoring matrices, Bioinformatics, 12, 135–143\cr
Claverie (1994) Some useful statistical properties of position-weight matrices.
Comput. Chem., 18, 287-293
}
\author{
Alicja Płuciennik & Michał Stolarczyk
}
\note{
Please note that when using other scoring matrix user needs to make sure that all alignment symbols are present there. Missing symbol will issue an error.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("alignment")
PC <- calculate_pseudo_counts(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pseudo counts}% use one of  RShowDoc("KEYWORDS")
\keyword{ substitution matrix}% __ONLY ONE__ keyword per line
