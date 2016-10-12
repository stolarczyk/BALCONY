\name{create_final_csv}
\alias{create_final_csv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create CSV file to save results
}
\description{
Create_final_csv() saves results as table into csv file. Combination of given variation allows to compare protein structure with evolutionary data content from alignment. Each position on alignment has its own column in csv file. WARNING: for long aligned sequences last colums can exceed csv file dimension!
}
\usage{
create_final_csv(variations_matrix,structure_matrix,structure_numbers)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variations_matrix}{matrix which contains alignment and frequencies of occurences each amino acids on each position of alignment}
  \item{structure_matrix}{matrix of aligned, examined protein sequence covered structure markers (T/N) }
  \item{structure_numbers}{vector of unaligned amino acids numbers(contains gaps) }
}
\details{

}
\value{

  \item{csv_file }{ file has tittled "final_csv.csv"}

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
show_numbers(), display_structure(), display_AA_variation()
}
\examples{
file=read.alignment(file="test.mfa", format="fasta", forceToLower=F)
aligned_sequences=file[[3]]
parameters=alignment_parameters(aligned_sequences);
var_aa=calculate_AA_variation(parameters,aligned_sequences,threshold=1)
variations_matrix=display_AA_variation(var_aa)
tunnel_file=read.table("structure1_4jnc.txt");
structure=create_structure_seq(tunnel_file,uniprot="P34914",file)
structure_matrix=display_structure(structure,tunnel_file)
structure_numbers=show_numbers(structure);

final_csv=create_final_csv(variations_matrix,structure_matrix,structure_numbers)
}
\keyword{ save }
\keyword{ output }% __ONLY ONE__ keyword per line
