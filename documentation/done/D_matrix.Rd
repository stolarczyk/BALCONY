\name{D_matrix}
\alias{D_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate substitution rate matrix between two amino acids
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to calculate Landgraf conservation metric. D_matrix contains substitution rates between two aminoacids in the alignment, according to the following formula:\cr \cr
\deqn{D(a,b)= (d(a,a)-d(a,b))/d(a,a)}\cr
where:\cr
\eqn{d(a,a)} is a probability of AA substitution by itself\cr
\eqn{d(a,b)} is a probability of substitution of amino acid a with other amino acid.
}
\usage{

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frequencies}{A matrix with probablity of substitutions, e.g.  Gonnet substitution matrix
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

 \item{distance}{A matrix of substitution probablities for all amino acids }
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
load(gonnet)
distance=D_matrix(gonnet)
}