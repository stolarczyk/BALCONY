\name{calculate_AA_variation}
\alias{calculate_AA_variation}

\title{
Calculate AA variations.
}
\description{
 This function visualises AA diversity on each position of the alignment.
}
\usage{
calculate_AA_variation(prmt, sequence_alignment, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prmt}{output of alignment_parameters())}
  \item{sequence_alignment}{ alignment[[3]] loaded with read.alignment()}
  \item{threshold}{number in range 0-100 of minimal frequency of occurences}
  

}
\details{
The output consists of amino acids and their fractions on each position of alignment. Amino acids with occurence frequencies lower than the threshold of user's choice are excluded.
}
\value{
Returns list of matrices with tabelarised symbols of the most common AA in alignment column and percentage values for contributed AA
 \item{var_aa$AA }{Matrix of AA on all alignment positions with decreasing frequencies in columns}
 \item{var_aa$per}{percentage of AA frequencies corresponding to the $AA}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alicja P³uciennik & Micha³ Stoalrczyk
}
\note{

}


\seealso{
alignment_parameters(), distplay_AA_variation()
}
\examples{
alignment=read.alignment(file="alignment.mfa", format="fasta", forceToLower=F)
aligned_sequences=alignment[[3]]
prmt=alignment_parameters(aligned_sequences)
threshold=10

var_aa=calculate_AA_variation(prmt,aligned_sequences,threshold)
display_AA_variation(var_aa)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
