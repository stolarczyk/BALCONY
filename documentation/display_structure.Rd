\name{display_structure}
\alias{display_structure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display structure on alignment matrix
}
\description{
This function allows to display desired structures in matrix as an aligned sequences. It is usefull for visualisation of amino acids conservation which create structures on a whole alignment background. For realiable analysis alignment should contains only close homologs of examined protein.
}
\usage{
display_structure(structure, structure_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{created structure sequences (chceck function: create_structure_seq())}
  \item{structure_file}{structure or list of structures, where in column there are amino acid numbers and amino acid symbols (one letter code) }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{output }{information about structures on aligned sequence, "-" means alignment gap, "N" means amino acid which do not create structure, "T" means amino acid forming structure}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja P³uciennik & Micha³ Stolarczyk
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#reading file(s) containing information on tunnels
tunnel_file1=read.table("T1_1cqz.txt");
tunnel_file2=read.table("T2_1cqz.txt");
tunnel_file3=read.table("T3_1cqz.txt");
tunnel_file=list(tunnel_file1,tunnel_file2,tunnel_file3)
seqid="P34913"
tunnel=create_tunel_seq(tunnel_file,seqid,file)
display_tunnels_matrix=display_tunnels(tunnel,tunnel_file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
