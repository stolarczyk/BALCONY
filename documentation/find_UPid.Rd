nd\name{find_UPid}
\alias{find_UPid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find UniProt sequence for PDB sequence in given alignment
}
\description{
This function allows to find UniProt sequence form alignment file corresponding to the given PDBid. Function requires library of equivalent sequences id defined by user.
%
}
\usage{
find_UPid(sequence="", library)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{id of PDB sequence}
  \item{library}{defined by user library of UniProt ids and PDB ids}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{uniprot}{function gives uniprot name (character)}

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
lib=list(c("Q6U6J0","4QLA"), #creating library uniprot - PDB
         c("Q84HB8","4I19","4QA9"),
         c("P34913","4JNC"),
         c("Q41413","3CXU","2CJP"))

uniprot=read_file(filename="4JNC");
UPid=find_UPid(sequence="2CJP", library);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ find }
\keyword{ UniProtID }% __ONLY ONE__ keyword per line
