\name{most_common}
\alias{most_common}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List sequences with the most frequent percentage of fitting consensus.
}
\description{
This function lists sequences which has the same, the most common level of similarity with consensus sequence. This allows to find dominating, group of probable closely related sequences in dataset. (There is no waranty that consensus match to the same positions with sequences)
}
\usage{

most_common(percentage, alignment_file)
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
  \item{list_of_sequences}{list of sequences with the most common similarity/identity percentage to the consensus}
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
list_homologs=most_common(percentage, consensus_file)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ homologs }

