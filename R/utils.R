.check_alignment_characters <- function(aligned_sequences_matrix) {
  total_char = length(unique(c(aligned_sequences_matrix)))
  if(total_char >  AA_COUNT +1)
    warning("There are ", total_char," character (AA) types in your alignment, which is more than prueba", AA_COUNT,": (", AA_COUNT-1 ,"AA + '-')")
  return(total_char)
}


#' Remove gaps and perform MSA
#'
#' @param seqj sequence to align
#' @param seqi sequence to align
#'
#' @return matrix of MSA scores
#' 
#' @export
#' 
#' @importFrom Biostrings pairwiseAlignment
#' 
#' @examples 
#' data("small_alignment")
#' seq = small_alignment$seq[1]
#' seq2 = small_alignment$seq[2]
#' .convertAlign(seq, seq2)
.convertAlign <- function(seqj, seqi) {
    return(pairwiseAlignment(.removeGaps(seqi),
                             .removeGaps(seqj), scoreOnly = T))
}

#' Remove gaps from string
#'
#' Can be used to remove gaps from an alignment sequence
#' 
#' @param string string with gaps decoded as "-"
#'
#' @return string without gaps
#' 
#' @export
#' 
#' @importFrom seqinr c2s s2c
#' 
#' @examples
#' data("small_alignment")
#' seq = small_alignment$seq[1]
#' .removeGaps(seq)
.removeGaps <- function(string) {
    char = s2c(string)
    return(c2s(char[char != "-"]))
}

#' Find sequences of numbers in a numeric vector
#'
#' This function finds sequences of consecutive numbers in numeric vectors
#'
#' Out of the following vector: 1,2,3,4,5,6,7,20,21,140,141 the function will find values starting the sequences: 1,20,140 and their lengths 7,2,2 respectively}
#'
#' @param vector A numeric vector to be analyzed
#'
#' @return \item{values}{A vector of values starting the consecutive sequences}
#' \item{lengths}{A vector of lengths of identified sequences}
#'
#' @export
#'
#' @examples
#' find_consecutive_seq(c(1,2,3,4,5,6,7,20,21,140,141,300,301,302))
find_consecutive_seq <- function(vector) {
  vector = append(vector, vector[length(vector)] + 2)
  diffs = diff(vector)
  ind = append(vector[1], vector[which(diffs != 1) + 1])
  cnt = 1
  lengths = c()
  for (i in seq(2, length(ind))) {
    lengths[cnt] = which(vector == ind[i]) - which(vector == ind[i - 1])
    cnt = cnt + 1
  }
  ind = ind[-length(ind)] # usuniecie sztucznej wartosci
  return(list(values = ind, lengths = lengths))
}

#' apply the pseudocounts to the position
#'
#' @param alignment The data loaded with \code{\link[seqinr]{read.alignment}} function
#' @param df_table a data.frame object
#' @param position position id
#'
#' @return a matrix
#'
#' @importFrom dplyr %>%
#'
#' @examples
#' #no example
.apply_pseudocounts_position <- function(alignment, df_table, position) {
  pc = calculate_pseudo_counts(alignment)
  names(df_table) = c("AAs", "freq")
  pseudo_counts = pc[, position]
  df_pc = as.data.frame(pseudo_counts)
  df_pc = dplyr::mutate(df_pc, AAs = rownames(df_pc))
  df_merged_table = suppressMessages(df_pc %>% dplyr::full_join(df_table))
  df_merged_table[is.na(df_merged_table)] = 0
  df_merged_table =
    df_merged_table %>%
    dplyr::mutate(Frequency = pseudo_counts + freq) %>%
    dplyr::select(AAs, Frequency)
  return(df_merged_table)
}

.merge_matrices <- function(m1, m2) {
  row_count = NROW(m1)
  col_count = NCOL(m1)
  output = matrix(NA, row_count * 2, col_count)
  counter = 1
  for (row in seq(1, row_count * 2, 2)) {
    output[row, ] = m1[counter, ]
    output[row + 1, ] = m2[counter, ]
    counter = counter + 1
  }
  return(output)
}

.is_upper <- function(string) {
  return(all(grepl("[[:upper:]]", strsplit(string, "")[[1]])))
}

#' Get sequence identifier by other sequence identifier in given alignment within a specified library
#'
#' This function allows to find sequence id from alignment file corresponding to the given sequence id. Function requires library of equivalent sequences id defined by user and it is useful to find sequences from other databases in alignment for examined sequence from other database (like PDB sequence for structure and UniProt sequences in alignment).
#'
#' @param sequence_id A string. An ID of e.g. PDB structure identifier
#' @param library A list of vectors which contain a defined by user library e.g. of UniProt ids <-> PDB ids. See examples
#'
#' @return A string. The equivalent ID to the one provided as the input.
#' @export
#'
#' @examples
#' #creating library uniprot - PDB
#' lib=list(  c("Q84HB8","4I19","4QA9"),
#'            c("P34913","4JNC"),
#'            c("P34914","1EK2","1CR6","1EK1","1CQZ"))
#' PDB_name = "1CQZ"
#' find_seqid(PDB_name,lib)
find_seqid <- function(sequence_id, library) {
  found = c()
  sequence_id = toupper(sequence_id)
  for (i in seq(1, length(library))) {
    found[i] = length(which((library[[i]] == sequence_id) == TRUE))
    found1 = which(found == 1)
  }
  seqid = library[[found1]][1]
  return(seqid)
}

#' Get sequence by id in alignment.
#'
#' This function allows to search for a sequence with its id. Useful for browsing a larg multiple sequence alignment data or for automatization purposes.
#'
#' @param sequence_id identifier of desired sequence from alignment
#' @param alignment data loaded with \code{\link[seqinr]{read.alignment}}
#'
#' @return A string, the desired aligned sequence form alignment
#' @export
#'
#' @examples
#' data("alignment")
#' #creating library uniprot - PDB
#' lib=list(  c("Q84HB8","4I19","4QA9"),
#'            c("P34913","4JNC"),
#'            c("P34914","1EK2","1CR6","1EK1","1CQZ"))
#' sequence_id=find_seqid("1CQZ",lib)
#' sequence=find_seq(sequence_id, alignment)
find_seq <- function(sequence_id, alignment) {
  mat = as.matrix(as.character(alignment[[3]]))
  nazwy_mat = .get_seq_names(alignment)
  which_uniprot = c()
  for (i in seq(1, length(nazwy_mat))) {
    #finding indices of uniprot ID in seq
    which_uniprot[i] = length(
      grep(
        sequence_id,
        nazwy_mat[i],
        ignore.case = FALSE,
        perl = FALSE,
        value = FALSE,
        fixed = FALSE,
        useBytes = FALSE,
        invert = FALSE
      )
    )
  }
  seqs = mat[which(which_uniprot == 1)]
  #seq of protein uniprot ID
  seqs_character = length(which(seqinr::s2c(seqs) != "-") == T)
  seq = list(sequence = seqs, len = seqs_character)
  return(seq)
}

excl_low_prob_strcts <-
  function(structure, threshold) {
    if (length(structure) == 3) {
      for (i in seq(1, dim(structure[[3]])[1], by = 1)) {
        to_exclude = which(structure[[3]][i, ] < threshold)
        structure[[3]][i, to_exclude] = NaN
        structure[[1]][i, to_exclude] = "N"
      }
    } else{
      warning("No probability information provided for the structure.")
    }
    return(structure)
  }

.get_seq_names <- function(alignment) {
  names = alignment[[2]]
  return(names)
}

.preprocess_hmm_output <- function(hmm_out) {
  size = dim(hmm_out)
  hmm_out = hmm_out[seq(from = 1, to = size[1] - 1, by = 3),-1]
  sites = as.numeric(unlist(hmm_out[, 21]))
  probs = hmm_out[, 1:20]
  #calculate probabilities; prob = -log(x)
  probs = exp(-probs)
  return(list(probabilities = probs, alignment_positions = sites))
}
