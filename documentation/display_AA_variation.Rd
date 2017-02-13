\name{display_AA_variation}
\alias{display_AA_variation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates table of AA and respective precentage of occurences.
}
\description{
This function concatenates matrices from the input list (calculate_AA_variation()).
}
\usage{
display_AA_variation(AA_variation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AA_variaton}{list of matrices with AA and their frequencies}

}
\details{
}
\value{

\item{table}{Matrix of AA and their frequencies}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja Pluciennik & Michal Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
calculate_AA_variation()
}
\examples{
alignment=read.alignment(file="alignment.mfa", format="fasta", forceToLower=F)
aligned_sequences=alignment[[3]]
prmt=alignment_parameters(aligned_sequences)
threshold=10

var_aa=calculate_AA_variation(prmt,aligned_sequences,threshold)
table = display_AA_variation(var_aa)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ variation}

