\name{worst_cons_for}
\alias{worst_cons_for}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consensus the worst fitted to sequence.
}
\description{
Find the sequence to which consensus is the worst fitted.
}
\usage{
worst_cons_for(percentage, alignment_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{percentage}{the percentage of similarity/identity to sequence}
  \item{alignment_file}{loaded alignment file (read.alignment()[seqinr])}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{name}{name of the worst fitted sequence}
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
percentage=cons2seqs_sim(aligned_sequences,number, consensus,'hydrophobicity')
best_consensus = worst_cons_for(percentage, file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
