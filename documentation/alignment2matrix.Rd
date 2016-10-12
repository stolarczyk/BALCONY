\name{alignment2matrix}
\alias{alignment2matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load alignment into matrix}
\description{The function loads alignment into matrix to facilitate more convenient operations
}
\usage{
alignment2matrix(prmt,alignment_squences)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{alignment_sequences}{ alignment[[3]] loaded with read.alignment()}
 \item{prmt}{output of alignment_parameter()}
}
\value{aligned sequences matrix where number of rows equals to number of aligned sequences and number of columns equals to the length length of aligned sequence}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja P?uciennik & Micha? Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
alignment_parameters() read.alignment() 
}
\examples{
alignment=read.alignment(file="alignment.mfa", format="fasta", forceToLower=F)
aligned_sequences=alignment[[3]]
prmt=alignment_parameters(aligned_sequences)
matrix=alignment2matrix(prmt,aligned_sequences)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ alignment }
\keyword{ conversion }
