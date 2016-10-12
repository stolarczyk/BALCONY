\name{get_seq_names}
\alias{get_seq_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get names of sequences from alignment
}
\description{
This simple function allows to get sequence names/identifies form alignment, useful to forward analysis. Alignment is obtained with read.alignment() {seqinr}.
}
\usage{
get_seq_names(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{alignment structure (list)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{names} list of names from alignment

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja P³uciennik & Micha³ Stolarczyk
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
alignment=read.alignment(file="test.mfa", format="fasta", forceToLower=F)
list_names=get_seq_names(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
