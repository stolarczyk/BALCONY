.check_alignment_characters <- function(aligned_sequences_matrix) {
  total_char = length(unique(c(aligned_sequences_matrix)))
  if(total_char != AA_COUNT-1)
    warning("There are ", total_char," character (AA) types in your alignment, which is more than ", AA_COUNT,": (", AA_COUNT-1 ,"AA + '-')")
  return(total_char)
}

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

.find_seqid <- function(sequence_id, library) {
  found = c()
  sequence_id = toupper(sequence_id)
  for (i in seq(1, length(library))) {
    found[i] = length(which((library[[i]] == sequence_id) == TRUE))
    found1 = which(found == 1)
  }
  seqid = library[[found1]][1]
  return(seqid)
}

.get_seq <- function(sequence_id, alignment) {
  mat = as.matrix(as.character(alignment[[3]]))
  nazwy_mat = get_seq_names(alignment)
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

.excl_low_prob_strcts <-
  function(structure, threshold) {
    if (length(structure) == 3) {
      for (i in seq(1, dim(structure[[3]])[1], by = 1)) {
        to_exclude = which(structure[[3]][i, ] < threshold)
        structure[[3]][i, to_exclude] = NaN
        structure[[1]][i, to_exclude] = "N"
      }
    } else{
      print("No probability information provided for the structure.")
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
