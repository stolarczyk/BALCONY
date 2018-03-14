\name{alignment2matrix}
\alias{alignment2matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load alignment into matrix}
\description{The function loads alignment into matrix to facilitate a convenient data manipulation
}
\usage{
alignment2matrix(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{alignment}{data loaded with \code{\link[seqinr]{read.alignment}} function }
 }
\value{
\item{matrix}{Aligned sequences matrix where number of rows equals to number of aligned sequences and number of columns equals to the length length of aligned sequences}
}

\author{
Alicja Pluciennik & Michal Stolarczyk
}


\seealso{
\code{\link{align_params}}, \code{\link[seqinr]{read.alignment}}
}
\examples{
data("alignment")
alignment = delete_isoforms(alignment)
matrix=alignment2matrix(alignment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matrix }
\keyword{ alignment }
\keyword{ conversion }
