\name{barplotshow}
\alias{barplotshow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shows barplot with amino acid variation on the specified position.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
barplotshow(position, AA_variation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position}{choosen column of alignment}
  \item{AA_variation}{percentage frequency of amino acids on alignment, calculated with display_AA_variation()}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
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
display_AA_variation(), calculate_AA_variation()
}
\examples{
seq=read.alignment("aligment.mfa")
aligned_sequences=seq[[3]]
prmt=alignnment_parameters(aligned_sequences)
theshold=4 # exemplary threshold, percentage, 4\%, below this value AA are not included into displayed matrix
var_aa=calculate_AA_variation(prmt,aligned_sequences,4)
AA_variation=display_AA_variation(var_aa)
barplotshow(AA_variation)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
