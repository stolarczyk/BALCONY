\name{get_seq_names}
\alias{get_seq_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get names of sequences from alignment
}
\description{
This function allows to get sequence names/identifires from alignment.
}
\usage{
get_seq_names(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{The alignment object read with read.alignment() function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{names}{A vector of characters with names of each sequence from the alignment}

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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(alignment)
sequences_names=get_seq_names(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ alignment }
\keyword{ sequences names }% __ONLY ONE__ keyword per line
