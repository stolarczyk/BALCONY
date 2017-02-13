\name{alignment2matrix}
\alias{alignment2matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load alignment into matrix}
\description{The function loads alignment into matrix to facilitate a convenient with BALCONY package.
}
\usage{
alignment2matrix(prmt,alignment_squences)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{alignment_sequences}{the third element of the list loaded with read.alignment() function }
 \item{prmt}{output of alignment_parameters() function. Number of sequences in the alignment and length of aligned sequences}
}
\value{Aligned sequences matrix where number of rows equals to number of aligned sequences and number of columns equals to the length length of aligned sequences}
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
