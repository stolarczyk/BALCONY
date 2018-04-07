\name{pseudo_counts}
\alias{pseudo_counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate pseudo counts for alignment 
}
\description{
This function calculates pseudo counts for alignment using substitution matrices. It is recommended to calculate amino acid frequencies for MSA with small number of sequences. This allows to calculate reliable entropy scores.
}
\usage{
pseudo_counts(alignment, substitution_mtx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment }{Alignment data read with read.alignment() function}
  \item{substitution_mtx }{Matrix with amino acids substitution frequencies. Default: GONNET.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{pseudoCounts }{Matrix with pseudo counts of size 21x number of alignment columns}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Henikoff et al.(1996) Using substitution probabilities to improve position-specific scoring matrices, Bioinformatics, 12, 135–143
[https://doi.org/10.1093/bioinformatics/12.2.135]
Claverie (1994) Some useful statistical properties of position-weight matrices.
Comput. Chem., 18, 287-293
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{
Using custom scoring matrix make sure that it takes all alignment symbols into consideration! Missing symbol will issue an error.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(alignment)
PC <- pseudo_counts_sub(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pseudo counts}% use one of  RShowDoc("KEYWORDS")
\keyword{ substitution matrix }% __ONLY ONE__ keyword per line
