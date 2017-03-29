\name{is_upper}
\alias{is_upper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check if the letter is uppercase.
}
\description{
This function facilitates the detection of uppercase strings/characters.
}
\usage{
is_upper(string)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{string}{
A string or character
}
}
\details{
All letters of a string must be uppercase for the string to be identified as an uppercase one
}
\value{
A logical value indicating if the string/character is an uppercase one
}

\author{
Michal Stolarczyk & Alicja Pluciennik
}

\examples{
string = "ABCD"
is_upper(string)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ uppercase }% use one of  RShowDoc("KEYWORDS")
\keyword{ detect }% __ONLY ONE__ keyword per line
