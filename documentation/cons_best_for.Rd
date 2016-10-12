\name{cons_best_for}
\alias{cons_best_for}

\title{
Consensus the best fitted to sequence
}
\description{
This function allows to find seqences which are the best represented with calculated consensus to better explore dataset.
}
\usage{
cons_best_for(percentage, alignment_file)
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

 \item{name}{name of the best fitted sequence}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
worst_cons_for()
}
\examples{
percentage=cons2seqs_sim(aligned_sequences,number, consensus,'hydrophobicity')
best_consensus = cons_best_for(percentage, file)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
