\name{create_final_CSV}
\alias{create_final_CSV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create CSV file to save results
}
\description{
Create_final_CSV() saves results as table into csv file. Combination of given variation allows to compare protein structure with evolutionary data content from alignment. Each position on alignment has its own column in csv file. If the length of the alignmnet exceeds 1000 characters, the output is divided into separate files with suffixes corresponing to the number of file produced by this function.

}
\usage{
create_final_CSV(filename,variations_matrix,structure,sequence_id,alignment,score_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{name of the output file produced by the function}
  \item{variations_matrix}{An object which contains alignment and frequencies of occurences each amino acids on each position of alignment. Output of \code{\link{calculate_AA_variation}}}
  \item{structure}{An strcture object - matrix of aligned, examined protein sequence covered by structure markers (S/N). Output of \code{\link{create_structure_seq}}.}
  \item{sequence_id}{the Uniprot code of the sequence of interest}
  \item{alignment}{the output of \code{\link[seqinr]{read.alignment}} function. A variable containing alignment data. One of the sequences must be the sequence of interest}
  \item{score_list}{list of calculated entropy/conservation scores. Optional parameter. If not provided, this rows are not present in the output file}
}
\value{

  \item{csv_file}{A comma separated variable file containing information provided to this function. It is also written in the current directory.}

}

\author{
Alicja Pluciennik & Michal Stolarczyk
}


\seealso{
\code{\link{create_structure_seq}}, \code{\link{schneider_conservativity}}, \code{\link{TG_conservativity}}, \code{\link{landgraf_conservativity}}, \code{\link[seqinr]{read.alignment}}
}
\examples{
data("alignment")
data("structure")
uniprot="P34914"
alignment = delete_isoforms(alignment)
threshold = 1
var_aa=calculate_AA_variation(alignment,threshold)
entropy_data=list(Schneider.entropy=schneider_conservativity(alignment),
                  TG.entropy = TG_conservativity(alignment))
create_final_CSV("my_filename",var_aa,structure,uniprot,alignment,entropy_data)
}
\keyword{ save }
\keyword{ output }% __ONLY ONE__ keyword per line
