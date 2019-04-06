#' Sample alignment of soluable epoxide hydrolase family
#'
#' was performed on a dataset comprising 311 soluble epoxide hydrolase peptide orthologous sequences acquired from UniProtKB. The alignment was performed and edited with MUSCLE algorithm in JALVIEW, respectively.
#'
#' @name alignment
#' @docType data
#'
#' @format An alignment object read with \code{\link[seqinr]{read.alignment}} function from seqinr package.
#' \describe{
#'   \item{nb}{A numeric: number of sequences}
#'   \item{nam}{A vector of characters: names of the sequences}
#'   \item{seq}{A vector of characters: amino acid sequences}
#' }
#'
#' @return A character vector of length of the aligned sequence containing consesus sequence based on the input alignment
#'
#' @keywords data
#'
#' @references
#' MUSCLE: \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113}\cr
#' JALVIEW: \url{https://academic.oup.com/bioinformatics/article/25/9/1189/203460/Jalview-Version-2-a-multiple-sequence-alignment}
#'
#' @examples
#' data("alignment")
#' alignment
NULL

#' Sample small alignment of soluable epoxide hydrolase family
#'
#' This alignment consists of 10 proteins which belong to the soluable epoxide hydrolase family. The amino acid sequences were aligned using MUSCLE algorithm with default settings.
#'
#' This is a smaller version of sample alignment which facilitates faster presentation of the functions capabilities.
#'
#' @name small_alignment
#' @docType data
#'
#' @format An alignment object read with \code{\link[seqinr]{read.alignment}} function from seqinr package.
#' \describe{
#'   \item{nam}{A vector of characters: names of the sequences}
#'   \item{nb}{A numeric: number of sequences}
#'   \item{seq}{A vector of characters: amino acid sequences}
#' }
#'
#' @return A character vector of length of the aligned sequence containing consesus sequence based on the input alignment
#'
#' @keywords data
#'
#' @examples
#' data("small_alignment")
#' small_alignment
NULL

#' Gonnet substitution matrix
#'
#' This dataset comprises the Gonnet substitution matrix which facilitates e.g. the calculation of Landgraf conservation score
#'
#' @name gonnet
#' @docType data
#'
#' @format An alignment object read with \code{\link[seqinr]{read.alignment}} function from seqinr package.
#' A data frame with 0 observations on the following 2 variables.
#' \describe{
#'   \item{AA names}{Names of amino acids included in the matrix}
#'   \item{matrix}{The substitution matrix itself}
#' }
#'
#' @return A character vector of length of the aligned sequence containing consesus sequence based on the input alignment
#'
#' @note Please note that this function masks the seqinr package function  \code{\link[seqinr]{consensus}}
#'
#' @keywords data
#'
#' @source
#' \url{http://imed.med.ucm.es/Tools/sias_help.html}
#'
#' @examples
#' data("gonnet")
NULL

#' Sample structure data
#'
#' This sample structure data consists of the amino acids names forming tunnels and their numbers is analyzed protein. The data is a result of CAVER which is a software tool for analysis and visualization of tunnels and channels in protein structures.
#'
#' The tunnel analysis with CAVER was performed on human epoxide hydrolase structure (PDB ID: 4JNC) 50ns MD simulation.
#'
#' @name structure
#' @docType data
#'
#' @format A structure object with three elements:
#' \describe{
#'   \item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
#'   \item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
#'   \item{structure_probabilities}{A matrix of numeric values: probabilities of corresponding to the structural information from first element of the output, which helps to reduce the effect of non-consistent structural amino acids on the conservativity analysis of the structure of interest}
#' }
#'
#' @return A character vector of length of the aligned sequence containing consesus sequence based on the input alignment
#'
#' @note Please note that this function masks the seqinr package function  \code{\link[seqinr]{consensus}}
#'
#' @keywords data
#'
#' @references
#' CAVER: \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002708}\cr
#' 4JNC: \url{http://www.sciencedirect.com/science/article/pii/S0960894X13004885}
#'
#' @examples
#' data("structure")
#' structure
NULL
