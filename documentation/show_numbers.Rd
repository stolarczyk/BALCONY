\name{show_numbers}
\alias{show_numbers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show number of amino acid in examined sequence 
}
\description{
This function allows to show on alignment table numbers of amino acids in given structure to identify and compare it with examined sequence.
}
\usage{
show_numbers(structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{structure}{structure sequence 
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{numbers}{vector numbers of amino acids in aligned, examined structure}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
alignment=read.alignment("aligment.mfa")
structure_file1=read.table("T1_mysz_1cqz.txt");
structure_file2=read.table("T2_mysz_1cqz.txt");
structure_file3=read.table("T3_mysz_1cqz.txt");
structure_file=list(structure_file1,structure_file2,structure_file3)
aligned_sequences=seq[[3]]
seqid="P34914"
structure=create_structure_seq(structure_file,seqid,alignment)
display_tunnels_matrix=display_structure(structure,structure_file)
show_numbers=show_numbers(structure);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
