\name{convert_AA_symbol}
\alias{convert_AA_symbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Amino acids symbols conversion
}
\description{
This function facilitates the conversion of three letter amino acids' codes to one letter equivalents.
}
\usage{
convert_AA_symbol(amino_acids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amino_acids}{
A character or vector of characters with amino acid(s) three letter code(s)
}
}
\details{
In case a vector of amino acid three letter codes is provided the function returns a vector of their one letter equivalents.
}
\value{
A chracter or vector of characters with amino acids one letter code(s)
}

\author{
Michal Stolarczyk & Alicja Pluciennik
}

\examples{
three_letter_codes = c("LEU", "VAL", "ALA")
convert_AA_symbol(three_letter_codes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ amino acids }% use one of  RShowDoc("KEYWORDS")
\keyword{ convert }% __ONLY ONE__ keyword per line
\keyword{ symbols }
