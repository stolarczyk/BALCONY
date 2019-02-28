








if (getRversion() >= "2.15.1")
  utils::globalVariables(c("sequence", "gonnet"))

# Conservation analysis ---------------------------------------------------
read_structure <- function(file_names) {
  structure_names = c()

  for (i in seq(1, length(file_names), by = 1)) {
    structure_names = append(structure_names, strsplit(file_names, "[.]")[[i]][1])

  }
  structure_names_list = as.list(structure_names)
  structure_list = list()

  i = 1

  for (structure in file_names) {
    temp = read.table(structure)
    structure_list[[i]] = temp

    i = i + 1

  }
  names(structure_list) = structure_names
  return(structure_list)
}

is_upper <- function(string) {
  return(all(grepl("[[:upper:]]", strsplit(string, "")[[1]])))
}

delete_isoforms <- function(alignment) {
  # This function searches for isoforms in the alignment (entries with "-digit|" in the name) and deletes them
  # As input it takes the alignment file
  # Output: alignment without isoforms
  lines_to_delete = c()

  pattern = "(-\\d+\\|)"

  for (i in seq(1, length(alignment$nam))) {
    test = grep(perl = T,
                pattern = pattern,
                x = alignment$nam[i])

    if (length(test) == 1) {
      lines_to_delete = append(lines_to_delete, i)
    }
  }
  new = list()
  if (length(lines_to_delete >= 1)) {
    alignment$nb = alignment$nb - (length(lines_to_delete))
    alignment$nam = alignment$nam[-lines_to_delete]
    alignment$seq = alignment$seq[-lines_to_delete]
    output = seqinr::as.alignment(
      nb = alignment$nb,
      nam = alignment$nam,
      seq = alignment$seq,
      com = NA
    )
    no_deleted = length(lines_to_delete)
    warning(paste(no_deleted, "isoforms were deleted"))
  }
  else{
    output = alignment
    warning("No isoforms detected")
  }
  return(output)
}
consensus <-  function(alignment, threshold) {
  #Function which calculates consensus
  #alignment-output of read.alignment() or alignment2matrix() function
  #threshold-given threshold of conservation (%)
  if (!is.matrix(alignment)) {
    alignment_matrix = alignment2matrix(alignment = alignment)
  } else{
    alignment_matrix = alignment
  }
  #alignment_matrix = as.matrix(as.character(alignment[[3]]))

  count_cols = dim(alignment_matrix)[1]
  #count_cols = length(alignment_matrix)
  vec = seq(1, dim(alignment_matrix)[2], 1)
  #vec = seq(1, length(seqinr::s2c(alignment_matrix[1,])))
  len_of_vers = length(vec)
  mat = matrix("-", count_cols, len_of_vers)

  for (i in seq(1, count_cols)) {
    #temp = seqinr::s2c(alignment_matrix[i])
    temp = alignment_matrix[i,]
    for (j in vec) {
      mat[i, j] = temp[j]
    }
  }

  consens = c(rep("*", each = len_of_vers))
  for (j in vec) {
    alignment_col = mat[, j]
    m = names(sort(table(alignment_col), decreasing = T))
    num = sort(table(alignment_col), decreasing = T)[1]
    value = num / length(alignment_col) * 100
    if (value >= threshold) {
      consens[j] = m[1]
    }
  }
  return(consens)
}
cons2seqs_ident <-  function(alignment, consensus_seq) {
  #alignment_sequence- file[[3]]
  # number_of_seq- file[[1]]
  #consensus_seq - calculated consensus (output of consensusus())
  true_percentage = c()

  for (i in seq(1, align_params(alignment)$row_no)) {
    true_percentage[i] = length(which((
      consensus_seq == seqinr::s2c(alignment$seq[i])
    ) == TRUE)) / length(consensus_seq)
    # true percentage calculation (number of the same AAs in consensusus and in each sequence/number of all AAs)
  }
  return(true_percentage)
}
align_params <- function(alignment) {
  #alignment data
  #return list of alignment size [row_numbers, col_numbers]
  aligned_sequences = alignment$seq
  row_num = length(aligned_sequences)
  col_num = length(seqinr::s2c(aligned_sequences[1]))
  param = list(row_no = row_num, col_no = col_num)
  return(param)
}

align_seq_mtx2grs <-
  function(aligned_sequences_matrix,
           grouping_method) {
    rows = dim(aligned_sequences_matrix)[1]
    cols = dim(aligned_sequences_matrix)[2]
    aligned_sequences_matrixG = aligned_sequences_matrix
    if (grouping_method == 'substitution_matrix') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (BLOSUM62 distance clustering based similarity)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "P" ||
              aligned_sequences_matrix[i, j] == "S" ||
              aligned_sequences_matrix[i, j] == "T")
            aligned_sequences_matrixG[i, j] = "G1"
          if (aligned_sequences_matrix[i, j] == "R" ||
              aligned_sequences_matrix[i, j] == "Q" ||
              aligned_sequences_matrix[i, j] == "E" ||
              aligned_sequences_matrix[i, j] == "K")
            aligned_sequences_matrixG[i, j] = "G2"
          if (aligned_sequences_matrix[i, j] == "N" ||
              aligned_sequences_matrix[i, j] == "D" ||
              aligned_sequences_matrix[i, j] == "G" ||
              aligned_sequences_matrix[i, j] == "H")
            aligned_sequences_matrixG[i, j] = "G3"
          if (aligned_sequences_matrix[i, j] == "I" ||
              aligned_sequences_matrix[i, j] == "L" ||
              aligned_sequences_matrix[i, j] == "M" ||
              aligned_sequences_matrix[i, j] == "V")
            aligned_sequences_matrixG[i, j] = "G4"
          if (aligned_sequences_matrix[i, j] == "F" ||
              aligned_sequences_matrix[i, j] == "W" ||
              aligned_sequences_matrix[i, j] == "Y")
            aligned_sequences_matrixG[i, j] = "G5"
          if (aligned_sequences_matrix[i, j] == "C")
            aligned_sequences_matrixG[i, j] = "G6"
        }
      }
    }

    if (grouping_method == 'polarity') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to the polarity and charges of the residues)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "W" ||
              aligned_sequences_matrix[i, j] == "S" ||
              aligned_sequences_matrix[i, j] == "T" ||
              aligned_sequences_matrix[i, j] == "Y" ||
              aligned_sequences_matrix[i, j] == "N" ||
              aligned_sequences_matrix[i, j] == "Q")
            #POLAR
            aligned_sequences_matrixG[i, j] = "G1"
          if (aligned_sequences_matrix[i, j] == "H" ||
              aligned_sequences_matrix[i, j] == "R" ||
              aligned_sequences_matrix[i, j] == "K")
            #POLAR AND CHARGED POSITIVELY
            aligned_sequences_matrixG[i, j] = "G1A"
          if (aligned_sequences_matrix[i, j] == "E" ||
              aligned_sequences_matrix[i, j] == "D")
            #POLAR AND CHARGED NEGATIVELY
            aligned_sequences_matrixG[i, j] = "G1B"
          if (aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "P" ||
              aligned_sequences_matrix[i, j] == "I" ||
              aligned_sequences_matrix[i, j] == "L" ||
              aligned_sequences_matrix[i, j] == "M" ||
              aligned_sequences_matrix[i, j] == "F" ||
              aligned_sequences_matrix[i, j] == "V" ||
              aligned_sequences_matrix[i, j] == "G" ||
              aligned_sequences_matrix[i, j] == "K")
            #OTHERS
            aligned_sequences_matrixG[i, j] = "G4"
        }
      }
    }
    if (grouping_method == 'size') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to size)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "G" ||
              aligned_sequences_matrix[i, j] == "S")
            #TINY
            aligned_sequences_matrixG[i, j] = "G1"
          if (aligned_sequences_matrix[i, j] == "V" ||
              aligned_sequences_matrix[i, j] == "C" ||
              aligned_sequences_matrix[i, j] == "D" ||
              aligned_sequences_matrix[i, j] == "N" ||
              aligned_sequences_matrix[i, j] == "T" ||
              aligned_sequences_matrix[i, j] == "P")
            #SMALL
            aligned_sequences_matrixG[i, j] = "G2"
          if (aligned_sequences_matrix[i, j] == "R" ||
              aligned_sequences_matrix[i, j] == "E" ||
              aligned_sequences_matrix[i, j] == "H" ||
              aligned_sequences_matrix[i, j] == "I" ||
              aligned_sequences_matrix[i, j] == "L" ||
              aligned_sequences_matrix[i, j] == "K" ||
              aligned_sequences_matrix[i, j] == "M" ||
              aligned_sequences_matrix[i, j] == "F" ||
              aligned_sequences_matrix[i, j] == "W" ||
              aligned_sequences_matrix[i, j] == "Y" ||
              aligned_sequences_matrix[i, j] == "Q")
            #OTHERS
            aligned_sequences_matrixG[i, j] = "G3"
        }
      }
    }
    if (grouping_method == 'aromaticity') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to the aromaricity)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "F" ||
              aligned_sequences_matrix[i, j] == "W" ||
              aligned_sequences_matrix[i, j] == "H" ||
              aligned_sequences_matrix[i, j] == "Y")
            #AROMATIC
            aligned_sequences_matrixG[i, j] = "G1"
          #NON-AROMATIC
          if (aligned_sequences_matrix[i, j] == "R" ||
              aligned_sequences_matrix[i, j] == "E" ||
              aligned_sequences_matrix[i, j] == "I" ||
              aligned_sequences_matrix[i, j] == "L" ||
              aligned_sequences_matrix[i, j] == "K" ||
              aligned_sequences_matrix[i, j] == "M" ||
              aligned_sequences_matrix[i, j] == "Q" ||
              aligned_sequences_matrix[i, j] == "V" ||
              aligned_sequences_matrix[i, j] == "C" ||
              aligned_sequences_matrix[i, j] == "D" ||
              aligned_sequences_matrix[i, j] == "N" ||
              aligned_sequences_matrix[i, j] == "T" ||
              aligned_sequences_matrix[i, j] == "P" ||
              aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "G" ||
              aligned_sequences_matrix[i, j] == "S")
          aligned_sequences_matrixG[i, j] = "G2"
        }
      }
    }
    return(aligned_sequences_matrixG)
  }

cons2seqs_sim <-
  function(grouped_alignment,
           grouped_consensus_seq) {
    row_no = dim(grouped_alignment)[1]
    col_no = dim(grouped_alignment)[2]
    aligned_sequences_matrixG = grouped_alignment
    consensusG = grouped_consensus_seq

    true_percentage_G = c()

    for (i in seq(from = 1, to = row_no, by = 1)) {
      true_percentage_G[i] = length(which((
        consensusG == (aligned_sequences_matrixG[i,])
      ) == TRUE)) / col_no
      # true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
    }
    return(true_percentage_G)
  }

alignment2matrix <- function(alignment) {
  #alignment data
  #returns alignment as a matrix
  prmt = align_params(alignment = alignment)
  aligned_sequences_matrix = matrix("-", prmt$row_no, prmt$col_no)

  for (i in seq(1, prmt$row_no)) {
    #Putting aligned seqs into matrix
    temp = toupper(seqinr::s2c(alignment$seq[i]))
    for (j in seq(1, prmt$col_no)) {
      aligned_sequences_matrix[i, j] = temp[j]

    }
  }
  return(aligned_sequences_matrix)
}

calculate_AA_variation <-
  function(alignment,
           threshold = NULL,
           grouped = F,
           grouping_method = "substitution_matrix",
           weights = NULL,
           pseudo_counts = F) {
    #prmt- size of alignment (output of get_parameter())
    #sequence_alignment-file[[3]]
    #threshold-threshold for detecting key amino acids (the percentage of all at the given position)
    #returns list of matrices with tabelarised symbols of the most common AA in alignment column and percentage values for contributed AA
    freq = c()
    prmt = align_params(alignment = alignment)
    if (pseudo_counts) {
      grouped = F
      weights = NULL
    }
    if (is.null(threshold)) {
      # keyaas_treshold = 1 / prmt$row_no
      keyaas_treshold = 0.0000001
    } else{
      if (threshold <= 0) {
        stop("The threshold must be > 0! Do not specify any to inspect all the residues in MSA.")
      } else{
        keyaas_treshold = prmt$row_no * (threshold / 100)
      }
    }

    aligned_sequences_matrix = alignment2matrix(alignment = alignment)
    if (grouped == T) {
      aligned_sequences_matrix = align_seq_mtx2grs(aligned_sequences_matrix, grouping_method = grouping_method)
    }
    keyaas = matrix("n", dim(aligned_sequences_matrix)[2], 21 * 2)
    keyaas_per = matrix("n", dim(aligned_sequences_matrix)[2], 21 * 2)
    for (i in seq(1, dim(aligned_sequences_matrix)[2])) {
      table = (sort(table(aligned_sequences_matrix[, i]), decreasing = T))
      if (pseudo_counts) {
        pc = calculate_pseudo_counts(alignment)
        df_table = as.data.frame(table)
        names(df_table) = c("AAs", "freq")
        pseudo_counts = pc[, i]
        df_pc = as.data.frame(pseudo_counts)
        df_pc = dplyr::mutate(df_pc, AAs = rownames(df_pc))
        df_merged_table = df_pc %>% dplyr::full_join(df_table)
        df_merged_table[is.na(df_merged_table)] = 0
        df_merged_table =
          df_merged_table %>%
          dplyr::mutate(Frequency = pseudo_counts + freq) %>%
          dplyr::select(AAs, Frequency)
      } else{
        if(length(table) == 1){
          df_merged_table = data.frame(col1=names(table), col2=table[1],row.names = 1)
        }else{
          df_merged_table = as.data.frame(table)
        }
        names(df_merged_table) = c("AAs", "Frequency")
      }
      #AA types and frequencies
      aas = names(table)
      over_threshold = which(df_merged_table$Frequency > keyaas_treshold)
      AAs = df_merged_table$AAs[over_threshold]
      Frequency = df_merged_table$Frequency[over_threshold]
      order_decresasing = order(Frequency, decreasing = T)
      Frequency = Frequency[order_decresasing]
      AAs = AAs[order_decresasing]

      keyaas[i, 1:length(matrix(AAs, 1, length(AAs)))] = matrix(AAs, 1, length(AAs))
      #if there are any AA (=always) then these are introduced to the keyaas matrix
      keyaas_per[i, 1:length(matrix(Frequency, 1, length(Frequency)))] = matrix(round(Frequency / sum(Frequency), 6) * 100, 1, length(AAs))
      #similarly in the percentages case
    }
    i = 1
    while (length(which(keyaas[, i] != "n")) > 0) {
      i = i + 1
    }
    i = i - 1
    keyaas = t(keyaas[, 1:i])
    #transpose matrix
    keyaas_per = t(keyaas_per[, 1:i]) #transpose matrix
    #merging key AAs symbols table with key AAs percentages table

    size = dim(keyaas)
    output = matrix("-", size[1] * 2, size[2])
    if (!is.null(weights)) {
      if (length(weights) != prmt$row_no)
        stop("The length of weights vector must equal the number of sequences in the alignment!")
      total_char = length(unique(c(aligned_sequences_matrix)))
      if(total_char != 21)
        warning("There are more than 20 AA types in your alignment.")
      weight = matrix("n",
                      ncol =  dim(aligned_sequences_matrix)[2],
                      nrow =  total_char)
      weights = weights / mean(weights)
      for (i in seq_len(ncol(aligned_sequences_matrix))) {
        representatives = unique(keyaas[, i])
        representatives = representatives[!representatives == "n"]
        for (j in seq_len(length(representatives))) {
          which_representative = which(aligned_sequences_matrix[, i] == representatives[j])
          weight[j, i] = mean(weights[which_representative])
        }
        a = weight[,i][which(weight[,i]!="n")]
        b = keyaas_per[,i][which(keyaas_per[,i]!="n")]
        keyaas_per[1:length(a), i] = as.numeric(b) * as.numeric(a)
        }
      }
    j = 1
    for (i in seq(1, size[1] * 2, 2)) {
      output[i, ] = keyaas[j, ]
      output[i + 1, ] = keyaas_per[j, ]
      j = j + 1
    }
    return(list(
      AA = keyaas,
      percentage = keyaas_per,
      matrix = output
    ))
  }

noteworthy_seqs <- function(percentage, alignment) {
  max = which.max(percentage)
  namelist = alignment[[2]]
  out.max = list(c(namelist[max], max)) #output is a name of sequence and position in alignment
  min = which.min(percentage) #percentage is an output of calculate_group_consensus or calculate_tru_consensus
  out.min = list(c(namelist[min], min)) #output is a name of sequence and position in alignment
  value = which.max(table(percentage))
  vector = which(percentage %in% names(value))
  names = match(vector, alignment[[2]])
  out.common = namelist[vector]
  out = list(
    best_consensus = out.max,
    worst_consensus = out.min,
    most_common = out.common
  )

  return(out)
}
convert_AA_symbol <- function(amino_acids) {
  map <-
    c(
      "H",
      "H",
      "H",
      "A",
      "R",
      "N",
      "D",
      "C",
      "E",
      "Q",
      "G",
      "H",
      "I",
      "L",
      "K",
      "M",
      "F",
      "P",
      "S",
      "T",
      "W",
      "Y",
      "V"
    )

  names(map) <-
    c(
      "HIP",
      "HID",
      "HIE",
      "ALA",
      "ARG",
      "ASN",
      "ASP",
      "CYS",
      "GLU",
      "GLN",
      "GLY",
      "HIS",
      "ILE",
      "LEU",
      "LYS",
      "MET",
      "PHE",
      "PRO",
      "SER",
      "THR",
      "TRP",
      "TYR",
      "VAL"
    )

  sapply(strsplit(amino_acids, "(?<=[A-Z]{3})", perl = TRUE),
         function(x)
           paste(map[x], collapse = ""))
}
find_seqid <- function(sequence_id, library) {
  ktory = c()
  sequence_id = toupper(sequence_id)
  for (i in seq(1, length(library))) {
    ktory[i] = length(which((library[[i]] == sequence_id) == TRUE))

    ktory1 = which(ktory == 1)

  }
  seqid = library[[ktory1]][1]

  return(seqid)
}
get_seq_names <- function(alignment) {
  names = alignment[[2]]

  return(names)
}
find_seq <- function(sequence_id, alignment) {
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


excl_low_prob_strcts <-
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

barplotshow <- function (position, AA_variation)
{
  row = max(which(AA_variation$percentage[, position] != "n"))
  x_names = AA_variation$AA[seq(1, row), position]
  x_val = round(as.numeric(AA_variation$percentage[seq(1,
                                                       row), position]), digits = 1)
  x_val = format(x_val, nsmall = 1)

  x_val_raw = as.numeric(AA_variation$percentage[seq(1,
                                                     row), position])

  ind = order(-x_val_raw)

  x_names = x_names[ind]
  x_val = x_val[ind]
  x_val_raw = x_val_raw[ind]



  plot = barplot(
    x_val_raw,
    names.arg = x_names,
    cex.main = 1.8,
    cex.lab = 1.0,
    cex.names = 1.2,
    cex.axis = 1.2,
    width = rep(0.95,
                row),
    xlim = c(0, 1 + dim(AA_variation[[1]])[1]),
    ylim = c(0,
             100),
    ylab = "Percentage",
    main = bquote(Amino ~ acid ~ variation ~ on ~ the ~ .(position) ~ position),
    col = "lightgrey"
  )

  z <- function(s) {
    if (s == "0.0") {
      "<0.1"
    } else {
      s
    }
  }

  x_val_raw <- x_val_raw + 9

  text(
    x = plot,
    y = x_val_raw,
    cex = 0.9,
    labels = paste(sapply(format(x_val, nsmall = 1), z), "%", sep = ""),
    xpd = T,
    srt = 90
  )
}



get_seq_weights <- function(alignment) {
  weights = c()
  gonnet_mtx = BALCONY::gonnet
  maxm = max(as.numeric(as.matrix(gonnet_mtx[[2]])))
  minm = min(as.numeric(as.matrix(gonnet_mtx[[2]])))
  sumMut = 0
  dist = matrix(
    data = NaN,
    nrow = length(alignment$seq),
    ncol = length(alignment$seq)
  )
  pb <- progress_bar$new(
    format = paste("Distance calculation", " [:bar] :percent eta: :eta"),
    total = length(alignment$seq),
    clear = T,
    width = 80
  )
  for (i in seq_len(length(alignment$seq))) {
    pb$tick()
    for (j in seq_len(length(alignment$seq))) {
      if (i != j) {
        seqi = s2c(alignment$seq[i])
        seqj = s2c(alignment$seq[j])
        non_gaps = intersect(which(seqi != "-"), which(seqj != "-"))
        for (pos in non_gaps) {
          aai = which(gonnet_mtx[[1]] == seqi[pos])
          aaj = which(gonnet_mtx[[1]] == seqj[pos])
          m = as.numeric(as.character(gonnet_mtx[[2]][aai, aaj]))
          Mut = (m - minm) / (maxm - minm)
          sumMut = sumMut + Mut
        }
        dist[i, j] = 1 - (sumMut / length(non_gaps))
        sumMut = 0
      } else{
        dist[i, j] = 0
      }
    }
  }
  for (sq in seq_len(length(alignment$seq))) {
    weights[sq] = sum(dist[sq,]) / (length(alignment$seq) - 1)
  }
  return(weights)
}

get_pos_based_seq_weights <- function(alignment, gap = TRUE, normalized = TRUE){
  align_param = align_params(alignment)
  weights_mtx = matrix(NA,nrow = align_param$row_no, ncol = align_param$col_no)
  alignment_mtx = alignment2matrix(alignment)

  if (gap == TRUE){
    # treat gap as other residues
    for(i in seq(1, ncol(alignment_mtx))){
      temp_table <- table(alignment_mtx[,i])
      for(j in 1:align_param$row_no){
        weights_mtx[j,i] = 1/(length(temp_table)*temp_table[alignment_mtx[j,i]])
      }
    }

    # chcecking if weights should be normalized
    if(normalized ==TRUE){
      weights = apply(weights_mtx,MARGIN = 1, sum)/align_param$col_no

    }else{
      weights = apply(weights_mtx,MARGIN = 1, sum)
    }
    return(weights)


  }else{

    # skip gaps
    for(i in seq(1, ncol(alignment_mtx))){
      temp_table <- table(alignment_mtx[,i])
      if("-" %in% names(temp_table)){
        temp_table <- temp_table[-which(names(temp_table)=="-")]
      }
      for(j in 1:align_param$row_no){
        if(alignment_mtx[j,i]== "-"){
          weights_mtx[j,i] = 0
        } else{
          weights_mtx[j,i] = 1/(length(temp_table)*temp_table[alignment_mtx[j,i]])

        }
      }
    }
    # chcecking if weights should be normalized
    if(normalized ==TRUE){
      weights = apply(weights_mtx,MARGIN = 1, sum)/align_param$col_no

    }else{
      weights = apply(weights_mtx,MARGIN = 1, sum)
    }
    return(weights)

  }
}


create_final_CSV <-
  function(filename,
           variations_matrix,
           structure,
           sequence_id,
           alignment,
           score_list = NULL) {
    sequence = seqinr::s2c(find_seq(sequence_id, alignment)[[1]])
    structure_output = rbind(structure$structure_matrix, structure$structure_numbers)
    structure_output_names = append(rownames(structure$structure_matrix), "Structure numbers")
    rownames(structure_output) = structure_output_names
    rownames(variations_matrix$matrix) = rep(c("AA name", "Percentage"), (dim(variations_matrix$matrix)[1] /
                                                                            2))
    if (is.null(score_list)) {
      alignment_position = seq(1, dim(variations_matrix$matrix)[2], by = 1)
      final_output = rbind(alignment_position,
                           variations_matrix$matrix,
                           sequence,
                           structure_output)
    }
    else{
      #Dodawanie wynikow konserwatywnosci tylko takich, jakie podal user
      scores_mtx = matrix(NA,
                          nrow = length(score_list),
                          ncol = length(score_list[[1]]))
      scores_mtx_names = c()

      for (i in seq(1, length(score_list))) {
        scores_mtx[i,] = score_list[[i]]
        scores_mtx_names[i] = names(score_list)[i]

      }
      rownames(scores_mtx) = scores_mtx_names
      alignment_position = seq(1, dim(variations_matrix$matrix)[2], by = 1)
      final_output = rbind(
        alignment_position,
        variations_matrix$matrix,
        sequence,
        structure_output,
        scores_mtx
      )

    }
    files_no = ceiling(dim(final_output)[2] / 1000)
    for (i in seq(1, files_no, 1)) {
      if (!i == files_no) {
        output = final_output[, ((i - 1) * 1000 + 1):(i * 1000)]
      }
      else{
        output = final_output[, ((i - 1) * 1000 + 1):dim(final_output)[2]]
      }
      write.table(
        output,
        file = paste(filename, "_", i, ".csv", sep = ""),
        row.names = T,
        col.names = F,
        sep = ","
      )
    }
    if (Sys.info()[[1]] == "Linux") {
      warning(paste("Output written to: ", getwd(), "/", filename, ".csv", sep = ""))
    }
    return(final_output)
  }

Escore_conservativity <-
  function(alignment,
           grouping_method = NULL,
           weights = NULL,
           pseudo_counts = F) {
    if (is.null(grouping_method)) {
      var_aa = calculate_AA_variation(
        alignment,
        threshold = 0.01,
        weights = weights,
        pseudo_counts = pseudo_counts
      )
    } else {
      var_aa = calculate_AA_variation(
        alignment,
        threshold = 0.01,
        grouped = T,
        grouping_method = grouping_method,
        weights = weights,
        pseudo_counts = F
      )
    }
    max_cons = c()

    for (i in seq(1, length(var_aa$matrix[1, ]), 1)) {
      if (is.na(as.numeric(var_aa$matrix[2, i])) == FALSE) {
        max_cons[i] = as.numeric(var_aa$matrix[2, i])
      }
    }
    AA = which(max_cons != 0)

    ile_var = c()

    for (i in seq(1, dim(var_aa$AA)[2], 1)) {
      ile_var[i] = length(which(var_aa$AA[, i] != "n" &
                                  var_aa$AA[, i] != "-"))
    }
    pre_conservativity = max_cons / ile_var

    pre_conservativity[which(is.na(pre_conservativity))] = 0
    part_con = pre_conservativity
    part_con[which(is.infinite(part_con))] = 0
    part_conserv = part_con / max(part_con)
    E = -(log(part_conserv))
    E[which(is.infinite(E))] = 0 # change Infs to 0
    E_score = (E / max(E))
    return_data = E_score

    return(return_data)
  }

kabat_conservativity <-
  function(alignment,
           weights = NULL,
           pseudo_counts = F) {
    if (!is.matrix(alignment)) {
      aligned_sequences_matrix = alignment2matrix(alignment = alignment)
    } else{
      aligned_sequences_matrix = alignment
    }
    Kabat = rep(NaN, dim(aligned_sequences_matrix)[2])

    var_aa = suppressWarnings(
      calculate_AA_variation(alignment, weights = weights, pseudo_counts = pseudo_counts)
    )
    N = length(aligned_sequences_matrix[, 1])
    for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
      column = aligned_sequences_matrix[, rep]
      #KABAT
      K = c()
      p = c()
      n = c()
      tab = c()
      #AA occurrences on the alignment position
      occureneces = sort(table(column), decreasing = T)
      #no of classes
      K = length(as.numeric(table(column)))
      AAs = var_aa$AA[, rep]
      PER = suppressWarnings(as.numeric(var_aa$percentage[, rep]))
      per = PER[!is.na(PER)]
      aas = AAs[!is.na(PER)]
      names(per) = aas
      #no of times the most frequent AA appeared
      n1 = occureneces[1]
      Kabat[rep] = (K / n1) * N
    }
    Kabat_entropy_normalized = Kabat / max(Kabat)

    return(Kabat_entropy_normalized)
  }

schneider_conservativity <-
  function(alignment,
           weights = NULL,
           pseudo_counts = F) {
    if (!is.matrix(alignment)) {
      aligned_sequences_matrix = alignment2matrix(alignment = alignment)
    } else{
      aligned_sequences_matrix = alignment
    }
    symbols = 21
    sum_schneider = rep(NaN, dim(aligned_sequences_matrix)[2])
    var_aa = suppressWarnings(
      calculate_AA_variation(alignment, weights = weights, pseudo_counts = pseudo_counts)
    )
    for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
      column = aligned_sequences_matrix[, rep]
      AAs = var_aa$AA[, rep]
      PER = suppressWarnings(as.numeric(var_aa$percentage[, rep]))
      per = PER[!is.na(PER)]
      aas = AAs[!is.na(PER)]
      names(per) = aas
      tab = as.numeric(per)
      K = c()
      p = c()
      n = c()
      N = length(column)
      K = length(as.numeric(table(column)))
      sum_internal_schneider = 0
      for (i in seq(1, K, by = 1)) {
        n[i] = tab[i]
        p[i] = n[i] / 100
        res = p[i] * log2(p[i])
        res_schneider = p[i] * log(p[i]) * (1 / log(symbols))
        sum_internal_schneider = sum_internal_schneider + res_schneider
      }
      sum_schneider[rep] = -sum_internal_schneider
    }
    return(sum_schneider)
  }

shannon_conservativity <-
  function(alignment,
           weights = NULL,
           pseudo_counts = F) {
    if (!is.matrix(alignment)) {
      aligned_sequences_matrix = alignment2matrix(alignment = alignment)
    }
    else{
      aligned_sequences_matrix = alignment
    }
    sum = rep(NaN, dim(aligned_sequences_matrix)[2])
    var_aa = suppressWarnings(
      calculate_AA_variation(
        alignment,
        threshold = 0.01,
        weights = weights,
        pseudo_counts = pseudo_counts
      )
    )

    for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
      column = aligned_sequences_matrix[, rep]
      AAs = var_aa$AA[, rep]
      PER = suppressWarnings(as.numeric(var_aa$percentage[, rep]))
      per = PER[!is.na(PER)]
      aas = AAs[!is.na(PER)]
      names(per) = aas
      tab = as.numeric(per)
      K = c()
      p = c()
      n = c()
      N = length(column)
      K = length(as.numeric(table(column)))
      sum_internal = 0
      for (i in seq(1, K, by = 1)) {
        n[i] = tab[i]
        p[i] = n[i] / 100
        res = p[i] * log2(p[i])
        sum_internal = sum_internal + res
      }
      sum[rep] = -sum_internal
    }
    return(sum)
  }

pairwise_alignment_MSA <- function(alignment) {
  #define a function to use in sapply
  convert_then_align <- function(seqj, seqi) {
    seqi = seqinr::s2c(seqi)
    seqi = seqinr::c2s(seqi[seqi != "-"])
    seqj = seqinr::s2c(seqj)
    seqj = seqinr::c2s(seqj[seqj != "-"])
    return(Biostrings::pairwiseAlignment(seqi, seqj, scoreOnly = T))
  }

  score_mtx = matrix(NA, nrow = alignment$nb, ncol = alignment$nb)
  for (i in seq_len(alignment$nb)) {
    seqi = alignment$seq[i]
    score_mtx[i,] = sapply(alignment$seq, convert_then_align, seqi)
    print(i)
  }
  return(score_mtx)
}

preprocess_hmm_output <- function(hmm_out) {
  size = dim(hmm_out)
  hmm_out = hmm_out[seq(from = 1, to = size[1] - 1, by = 3),-1]
  sites = as.numeric(unlist(hmm_out[, 21]))
  probs = hmm_out[, 1:20]
  #calculate probabilities; prob = -log(x)
  probs = exp(-probs)
  return(list(probabilities = probs, alignment_positions = sites))
}

CRE_conservativity <-
  function(alignment,
           hmmbuild_path = NULL,
           pairwiseAlignemnt_scores = NULL) {
    if (is.null(pairwiseAlignemnt_scores)) {
      pairwiseAlignemnt_scores = pairwise_alignment_MSA(alignment)
    }
    #Make sure the hmmbuild is there
    if (.Platform$OS.type == "unix") {
      if (is.null(hmmbuild_path)) {
        hmmbuild_path = system("which hmmbuild", intern = T)
      }
      if (length(hmmbuild_path) == 0) {
        stop("You need to specify the hmmbuild path as it cannot be located in your system")
      }
    } else {
      if (is.null(hmmbuild_path)) {
        stop("You need to specify the hmmbuild path as it cannot be located in your system")
      }
    }

    #perform hierarchical clustering on the distance matrix obtained from pairwise alignemnt matrix
    dendro = stats::hclust(stats::dist(pairwiseAlignemnt_scores), method = "average")
    #cut the tree and get clusters
    clusters = stats::cutree(dendro, h = 8000)
    #get the number of clusters
    no_clusters = max(clusters)
    RE = matrix(data = NA,
                nrow = no_clusters,
                ncol = length(seqinr::s2c(alignment$seq[1])))
    CRE = c()
    Z = c()
    #iterate over each clusters
    for (i in seq_len(no_clusters)) {
      #create copies of the alignemnt object to work on
      sub_alignment = alignment
      leftover_alignment = alignment
      #get the indices of the sequences in the cluster of interest
      to_delete = which(clusters == clusters[i])
      if (length(to_delete) > 1) {
        #delete appropriate sequences and their names in the objects and adjust the sequeces counts
        sub_alignment$seq = sub_alignment$seq[to_delete]
        leftover_alignment$seq = leftover_alignment$seq[-to_delete]
        sub_converted_seq = lapply(sub_alignment$seq, seqinr::s2c)
        leftover_converted_seq = lapply(leftover_alignment$seq, seqinr::s2c)
        sub_alignment$nam = sub_alignment$nam[to_delete]
        leftover_alignment$nam = leftover_alignment$nam[-to_delete]
        leftover_alignment$nb = leftover_alignment$nb - length(to_delete)
        sub_alignment$nb = length(to_delete)
        #write updated objects to files
        leftover_fasta_name = "leftover_MSA.fasta"
        sub_fasta_name = "sub_MSA.fasta"
        seqinr::write.fasta(sequences = leftover_converted_seq,
                            names = leftover_alignment$nam,
                            file.out = leftover_fasta_name)
        seqinr::write.fasta(sequences = sub_converted_seq,
                            names = sub_alignment$nam,
                            file.out = sub_fasta_name)
        #prepare commands to run HMMER
        sub_hmm_command = paste(hmmbuild_path, "sub_hmm.out", sub_fasta_name)
        leftover_hmm_command = paste(hmmbuild_path, "leftover_hmm.out", leftover_fasta_name)
        system(command = sub_hmm_command, wait = T)
        system(command = leftover_hmm_command, wait = T)
        #read the HMMER outputs in
        leftover_hmm <-
          readr::read_table("leftover_hmm.out",
                            col_names = FALSE,
                            skip = 21)
        sub_hmm <-
          readr::read_table("sub_hmm.out", col_names = FALSE, skip = 21)
        #preprocess the data
        leftover_prob = preprocess_hmm_output(leftover_hmm)$probabilities
        leftover_pos = preprocess_hmm_output(leftover_hmm)$alignment_positions
        sub_prob = preprocess_hmm_output(sub_hmm)$probabilities
        sub_pos = preprocess_hmm_output(sub_hmm)$alignment_positions
        #get the indices of alignemnt positions avaialble in both groups
        intersection_pos = intersect(sub_pos, leftover_pos)
        for (pos in seq_len(length(seqinr::s2c(alignment$seq[1])))) {
          which_pos = which(intersection_pos == pos)
          if (length(which_pos) > 0) {
            pre_RE = c()
            for (aa in seq_len(dim(sub_prob)[2])) {
              sub_prob_pos = which(sub_pos == pos)
              leftover_prob_pos = which(leftover_pos == pos)
              pre_RE[aa] = sub_prob[sub_prob_pos, aa] * log(sub_prob[sub_prob_pos, aa] /
                                                              leftover_prob[leftover_prob_pos, aa])
            }
            RE[i, pos] = sum(pre_RE, na.rm = T)
          } else{
            RE[i, pos] = 0
          }
        }
      } else{
        RE[i,] = 0
      }
    }
    for (pos in seq_len(length(seqinr::s2c(alignment$seq[1])))) {
      CRE[pos] = sum(RE[, pos])
    }
    return(CRE / max(CRE))
  }


substitution_mtx <- function (matrix_name) {
  # function can read .txt file substitution matrix and convert it
  # into a list of alphabet (the range of letters in matrix) and
  # mtx(valuse into substitution matrix)
  matrix = read.table(matrix_name)
  mtx = matrix[-1,-1]
  alphabet = as.character(unlist(matrix[1,-1]))
  sub_mtx = list(names = alphabet, matrix = mtx)
  return(sub_mtx)
}
D_matrix <- function(substitution_matrix) {
  values = as.numeric(as.matrix(substitution_matrix[[2]]))
  dim(values) <- dim(substitution_matrix[[2]])
  k = dim(substitution_matrix[[2]])[1]
  distance = matrix(0, k, k)
  #CREATE DISTANCE MATRIX
  for (i in 1:k) {
    index = i

    for (j in 1:k) {
      if (i != j) {
        distance[i, j] = (values[i, i] - values[i, j] / values[i, i])
      }
      if (i == j) {
        distance[i, j] <- 0
      }

    }
  }
  output = list(substitution_matrix[[1]], distance)
  return(output)
}
landgraf_conservativity <-
  function(matrix_name = NULL, alignment, weights) {
    if (is.null(matrix_name)) {
      gonnet = BALCONY::gonnet
      pre_dissim_mtx = gonnet
    } else{
      pre_dissim_mtx = substitution_mtx(matrix_name)
    }
    aligned_sequences_matrix = alignment2matrix(alignment = alignment)
    dissim_mtx = D_matrix(pre_dissim_mtx)
    conservation = rep(NaN, dim(aligned_sequences_matrix)[2])
    status = 0

    pb <- progress_bar$new(
      format = paste("Landgraf calculation", " [:bar] :percent eta: :eta"),
      total = dim(aligned_sequences_matrix)[2],
      clear = T,
      width = 80
    )
    for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
      pb$tick()
      column = aligned_sequences_matrix[, rep]

      values = as.numeric(as.matrix(dissim_mtx[[2]]))
      alpha = dissim_mtx[[1]]
      dim(values) <- dim(dissim_mtx[[2]])
      iterator = seq(1:length(column))
      sum_of_dist = 0

      global_sum = 0

      for (i in iterator) {
        iterator2 = c((i + 1):length(column))
        posa = match(column[i], alpha)
        if (!i == length(column)) {
          for (j in iterator2) {
            posb = match(column[j], alpha)
            tempi = weights[i] * values[posa, posb]
            tempj = weights[j] * values[posb, posa]
            sum_of_dist = tempi + tempj
            global_sum = global_sum + sum_of_dist
          }
        }
      }
      conservation[rep] = global_sum / length(column)
    }
    Landgraf_normalized_entropy = conservation / max(conservation)
    return(Landgraf_normalized_entropy)
  }

calculate_pseudo_counts <-
  function(alignment, substitution_mtx = NULL) {
    if (is.null(substitution_mtx)) {
      gonnet = BALCONY::gonnet
      substitution_mtx <-
        apply(
          X = gonnet[[2]],
          MARGIN = 2,
          FUN = function(X)
            as.numeric(as.character(X))
        )
      substitution_mtx = exp(substitution_mtx)
      colnames(substitution_mtx) <- gonnet[[1]]
      rownames(substitution_mtx) <- gonnet[[1]]
    }
    #warning - other substitution matrices may have different symbols (like '*', or other letters)
    mtx_alignment = alignment2matrix(alignment)
    pseudoCounts = matrix(NA,
                          nrow = length(append(Biostrings::AA_STANDARD, "-")),
                          ncol = dim(mtx_alignment)[2])
    B = apply(
      X = mtx_alignment,
      MARGIN = 2,
      FUN = function(X)
        (5 * length(unique(X)))
    )
    calc_ba <- function(column, B, substitution_mtx) {
      pseudocounts = matrix(NA, nrow = length(append(Biostrings::AA_STANDARD, "-")), ncol = 1)
      names(pseudocounts) <- append(Biostrings::AA_STANDARD, "-")
      N = sum((table(column)))
      AA = table(column)
      for (i in names(pseudocounts)) {
        to_sum <- rep(NA, length(AA))
        cnt = 1
        for (j in names(AA)) {
          to_sum[cnt] <- as.numeric(AA[j] / N) * substitution_mtx[j, i]
          cnt = cnt + 1
        }
        pseudocounts[i] <- sum(to_sum)
      }
      output = pseudocounts * (B / sum(pseudocounts))
      return(output)
    }
    for (column in seq_len(ncol(mtx_alignment))) {
      pseudoCounts[, column] <-
        calc_ba(mtx_alignment[, column], B[column], substitution_mtx)
    }
    rownames(pseudoCounts) = append(Biostrings::AA_STANDARD, "-")
    return(pseudoCounts)
  }



# Structure analysis ------------------------------------------------------


create_structure_seq <-
  function(structure_list,
           sequence_id,
           alignment,
           pdb_path = NULL,
           chain_identifier = NULL,
           shift = NULL) {
    #structure_list-> list of structures in protein
    #sequence_id -> uniprot id  which has been found by read_file(filename="PDBid") with PDB indentifier ;
    #alignment-> file with alignment (alignment.fasta)
    #shift-> if there is missing domain in structure
    seqs = list()
    probs = list()
    structure = list()
    probability = list()
    if (!is.null(pdb_path)) {
      if (is.null(chain_identifier)) {
        print("Chain identifier required!")
        stop()
      }
      missing = find_consecutive_seq(
        get_remarks465_pdb(pdb_path = pdb_path, chain_identifier = chain_identifier)$aa_numbers
      )
    } else{
      missing = NULL
    }
    if (length(missing$values) == 0) {
      missing = NULL
    }

    #finds an uniprot name from alignent- extract uniprot seq from alignment
    base_seq = find_seq(sequence_id, alignment)
    structure_probabilities = list()
    for (i in seq(1, length(structure_list), by = 1)) {
      structures_indices = as.vector(structure_list[[i]][[1]])
      structures_names = as.vector(structure_list[[i]][[2]])
      if (length(structure_list[[i]]) == 3) {
        prob_data = T
        structures_probability = as.vector(structure_list[[i]][[3]])
        structures_probability = structures_probability[-1]
        structures_probability = structures_probability / max(structures_probability)# normalization

      } else{
        prob_data = F
        structures_probability = NULL
      }
      structures_indices = structures_indices[-1]
      structure_idx = as.numeric(structures_indices)
      structures_names = structures_names[-1]

      if (!is.null(shift)) {
        structure_idx = structure_idx + rep(shift, length(structure_idx))
      }

      if (length(missing) != 0) {
        for (z in seq(1, length(missing$values))) {
          for (l in seq(1, length(structure_idx))) {
            if (structure_idx[l] >= missing$values[z]) {
              structure_idx[l] = structure_idx[l] + missing$lengths[z]
            }
          }
        }
      }
      seqs[[i]] = rep("N", each = base_seq$len)
      seqs[[i]][structure_idx] = "S"
      aa_positions = which(seqinr::s2c(base_seq$sequence) != "-")
      just_align = alignment[[3]]
      length_alignment = align_params(alignment)$col_no
      structure[[i]] = rep("-", each = length_alignment)
      probability[[i]] = rep(NaN, each = length_alignment)
      if (prob_data == T) {
        probs[[i]] = rep(NaN, each = base_seq$len)
        probs[[i]][structure_idx] = structures_probability
      }
      j = 1
      for (a in aa_positions) {
        #aligning the structures information with the alignment sequence
        structure[[i]][a] = seqs[[i]][j]
        if (prob_data == T) {
          probability[[i]][a] = probs[[i]][[j]]
        }
        j = j + 1
      }
    }
    struc_length = length(structure[[1]])
    count = length(structure_list)
    struct_output = matrix("-", count, struc_length)
    probability_output = matrix(NaN, count, struc_length)
    v = seq(1, count, by = 1)
    for (i in v) {
      struct_output[i, ] = structure[[i]]
      if (prob_data == T) {
        probability_output[i, ] = probability[[i]]
      }
    }
    nr_stru = rep("-", length(structure[[1]]))
    j = 1
    for (i in seq(1, length(structure[[1]]), 1)) {
      if (structure[[1]][i] != "-") {
        nr_stru[i] = j
        j = j + 1
      }
    }
    rownames(struct_output) = names(structure_list)
    if (prob_data == T) {
      rownames(probability_output) = names(structure_list)
    }
    if (prob_data == T) {
      return(
        list(
          structure_matrix = struct_output,
          structure_numbers = nr_stru,
          structure_probabilities = probability_output
        )
      )
    } else{
      return(list(
        structure_matrix = struct_output,
        structure_numbers = nr_stru
      ))
    }
  }

get_remarks465_pdb <- function(pdb_path, chain_identifier) {
  #pdb_file_path: a path to the pdb file
  #chain_identifier: a character spcifying the chain to analyze
  #returns: list of 1) numbers of aminocids which are missing 2) chain identifier
  pdb_file = Rpdb::read.pdb(
    file = pdb_path,
    REMARK = T,
    ATOM = T,
    CRYST1 = F,
    TITLE = T,
    MODEL = F,
    HETATM = F,
    CONECT = F
  )
  pdb_file = pdb_file$remark
  pattern = "REMARK 465"
  indices = grep(perl = T,
                 pattern = pattern,
                 x = pdb_file)
  remarks = pdb_file[indices[8:length(indices)]]
  aa_number = c()
  chain = c()
  for (i in seq(1, length(remarks), by = 1)) {
    remark_line = strsplit(remarks[i], split = " ")
    x = remark_line[[1]]
    aa_number[i] = x[which(x != "")][length(x[which(x != "")])]
    chain[i] = x[which(x != "")][length(x[which(x != "")]) - 1]
  }
  aa_number = as.numeric(aa_number)
  chain_indices = which(chain == chain_identifier)
  aa_number = aa_number[chain_indices]
  chain = chain[chain_indices]
  return(list(aa_numbers = aa_number, chain = unique(chain)))
}

find_consecutive_seq <- function(vector) {
  vector = append(vector, vector[length(vector)] + 2) # wydluzenie o sztuczna wartosc wieksza od ostatniej rzeczywistej o 2
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

get_structures_idx <- function(structure) {
  #documentation get_structure_idx.Rd
  #get idx of structure in alignemnt
  #return sorted list of indices in MSA, one list element is one structure
  #the first element is indices of whole protein in alignment
  structure = structure$structure_matrix
  stru_index = list()
  for (i in seq(1:length(structure[, 1]))) {
    stru_index[[i]] = which(structure[i, ] == "S")

  }
  whole_prot = which(structure[1, ] != "-")

  out = list(proteinIndices = whole_prot, structureIndices = stru_index)
  names(out[[2]]) = rownames(structure)
  return(out)
}
get_prot_entropy <- function(protein_index, score_list) {
  #documentation get_prot_entropy.Rd
  #allows to get idx of whole protein in alignment
  #returns list of entropy for protein
  if (is.list(score_list))
  {
    prot_cons = list()
    for (i in seq(1, length(score_list))) {
      prot_cons[[i]] = score_list[[i]][protein_index]
    }
    names(prot_cons) <- names(score_list)
  }
  else{
    print("score_list is not a list!")
    stop()
  }
  return(prot_cons)
}

plot_entropy <-
  function(protein_conservation,
           colors,
           impose = NULL,
           prot_name = NULL,
           legend_pos = NULL) {
    #plots scores on one plot, if colors are not specified plot as rainbow
    #automagically uses name of pdb FIXME
    # recive list of entropy scores and list of Names in this list
    if (missing(colors)) {
      colors <- rainbow(length(protein_conservation))
    }
    if (is.null(impose)) {
      impose <- T
    }
    if (!is.null(prot_name)) {
      title_str = paste("for ", prot_name)
    }
    else
    {
      title_str = ""
    }
    if (is.null(legend_pos)) {
      legend_pos = "bottomleft"
    }
    plot(
      protein_conservation[[1]],
      ylim = c(0, 1),
      col = colors[1],
      xlab = "amino acid",
      ylab = "entropy score",
      main = paste("Entropy score for entire protein", title_str) ,
      type = "l"
    )
    for (i in seq(2, length(protein_conservation))) {
      par(new = impose)
      plot(
        protein_conservation[[i]],
        ylim = c(0, 1),
        col = colors[i],
        xlab = "",
        ylab = "",
        main = "",
        type = "l"
      )
    }
    legend(legend_pos,
           names(protein_conservation),
           col = colors,
           lty = c(1))
  }

get_structures_entropy <- function(structure_index, score_list) {
  #structure_index output of get_structures_idx  is a list of indexes in alignment of protein and structures
  #SCORE_LIST list of entropies for whole alignment
  #output is a list of matrixes where each row contains values of entropy for AA in structure
  t_index = structure_index
  Entropy = list()
  lengths = list()
  for (i in seq(1:length(t_index))) {
    lengths[[i]] = length(t_index[[i]])
    output = matrix(NA, nrow = length(score_list), ncol = lengths[[i]])
    for (j in seq(1:length(score_list))) {
      output[j, ] = score_list[[j]][t_index[[i]]]
    }
    rownames(output) <- names(score_list)
    Entropy[[i]] = output
  }
  names(Entropy) = names(structure_index[[2]])
  return(Entropy)
}

prepare_structure_profile <-
  function(structure, structure_entropy) {
    #structure-> output of create_structure_seq
    #structure_entropy -> list of entropy scores for alignment
    megalist = list()
    score_count = dim(structure_entropy[[1]])[1] # ????
    names <- rownames(structure[[1]])
    for (i in seq(1, length(structure[[1]][, 1]))) {
      StruEnt = list(entropy = c(), idx = c())
      StruEnt$idx = as.numeric(structure[[2]][which(structure[[1]][i, ] ==
                                                      "S")])
      entropy_mtx = matrix(NaN, nrow = score_count, ncol = length(StruEnt$idx))
      for (j in seq(1, score_count)) {
        entropy_mtx[j, ] = structure_entropy[[i]][j,]
      }
      rownames(entropy_mtx) <- rownames(structure_entropy[[1]])
      StruEnt$entropy = entropy_mtx
      megalist[[i]] = StruEnt
    }
    names(megalist) <- names
    return(megalist)
  }
plot_structure_on_protein <-
  function(protein_entropy,
           structure_profiles,
           pdb_name = NULL,
           colors = NULL,
           structure_names = NULL,
           legend_pos = NULL) {
    #protein entropy- list of an entropy values (with different scores)
    #structure_profils - list of entropy scores for structures (as list) where [[1]] entropy value
    # and [[2]] index in protein


    if (length(protein_entropy) == length(structure_profiles[[1]][[1]][, 1])) {
      score_names = names(protein_entropy)
      prot_lenght = length(protein_entropy[[1]])
      StruLen = length(structure_profiles)
      if (is.null(pdb_name)) {
        pdb_name = " "
      }
      if (is.null(colors)) {
        colors = rainbow(StruLen)
      }
      if (is.null(structure_names)) {
        structure_names = c()
        for (i in seq(1, StruLen)) {
          #structure_names[i] = paste("stru", i)
          structure_names[i] = names(structure_profiles)[i]

        }
      }
      if (is.null(legend_pos)) {
        legend_pos = "bottomleft"
      }
      for (i in seq(1, length(score_names))) {
        if (isRStudio <- Sys.getenv("RSTUDIO") == "1") {

        } else{
          #X11()
          dev.new()
        }
        plot(
          protein_entropy[[i]],
          col = "black",
          type = "l",
          main = paste(score_names[i], " score for ", pdb_name),
          xlim = c(0, prot_lenght),
          ylim = c(0.0, 1.0),
          xlab = 'Amino Acid',
          ylab = 'Entropy'
        )
        for (j in seq(1, StruLen)) {
          par(new = T)
          plot(
            structure_profiles[[j]][[2]],
            structure_profiles[[j]][[1]][i, ],
            col = colors[[j]],
            main = "",
            pch = j,
            xlim = c(0, prot_lenght),
            ylim = c(0.0, 1.0),
            xlab = '',
            ylab = ''
          )
          legend(
            legend_pos,
            c(pdb_name, structure_names),
            lty = c(1, rep(0, j)),
            pch = c(-1, seq(1, j)),
            lwd = c(2.5, 2.5),
            col = c("black", colors)
          )
        }
      }
    } else{
      print("The lists contain different number of conservation/entropy scores!")
    }
  }
compare_cons_metrics <-
  function(protein_entropy,
           structure_profile,
           pdb_name) {
    metrics_count = length(protein_entropy)
    structures_count = length(structure_profile)
    colors = rainbow(structures_count)
    structure_names = c()
    for (i in seq(1, structures_count)) {
      #structure_names[i] = paste("stru", i)
      structure_names[i] = names(structure_profile)[i]
    }
    for (i in seq(1, metrics_count)) {
      for (j in seq(1, metrics_count)) {
        if (i != j) {
          if (isRStudio <- Sys.getenv("RSTUDIO") == "1") {
            plot(
              protein_entropy[[i]],
              protein_entropy[[j]],
              main = paste(
                'Scatterplot of',
                names(protein_entropy)[i],
                "vs. ",
                names(protein_entropy)[j]
              ),
              xlim = c(0, 1),
              ylim = c(0, 1),
              xlab = names(protein_entropy)[i],
              ylab = names(protein_entropy)[j],
              pch = 20,
              col = alpha("slategray", 0.5)
            )
            for (k in seq(1, structures_count)) {
              par(new = T)
              plot(
                structure_profile[[k]][[1]][i, ],
                structure_profile[[k]][[1]][j, ],
                col = alpha(colors[k], 0.7),
                pch = k,
                main = "",
                xlim = c(0, 1),
                ylim = c(0, 1),
                xlab = "",
                ylab = ""
              )
            }
            legend(
              'bottomright',
              c(pdb_name, structure_names),
              pch = c(20, seq(1, k)),
              col = scales::alpha(c("slategray", colors), 0.5)
            )
          } else{
            #X11()
            dev.new()
            plot(
              protein_entropy[[i]],
              protein_entropy[[j]],
              main = paste(
                'Scatterplot of',
                names(protein_entropy)[i],
                "vs. ",
                names(protein_entropy)[j]
              ),
              xlim = c(0, 1),
              ylim = c(0, 1),
              xlab = names(protein_entropy)[i],
              ylab = names(protein_entropy)[j],
              pch = 20,
              col = alpha("slategray", 0.5)
            )
            for (k in seq(1, structures_count)) {
              par(new = T)
              plot(
                structure_profile[[k]][[1]][i, ],
                structure_profile[[k]][[1]][j, ],
                col = alpha(colors[k], 0.7),
                pch = k,
                main = "",
                xlim = c(0, 1),
                ylim = c(0, 1),
                xlab = "",
                ylab = ""
              )
            }
            legend(
              'bottomright',
              c(pdb_name, structure_names),
              pch = c(20, seq(1, k)),
              col = scales::alpha(c("slategray", colors), 0.5)
            )
          }
        }
      }
    }

  }

kolmogorov_smirnov_test <-
  function(protein_entropy,
           structure_entropy,
           alternative,
           pdb_name = "Reference",
           range = NULL,
           make_plot = NULL) {
    if (is.null(make_plot)) {
      make_plot <- T
    }
    metrics_count = length(protein_entropy)
    structures_count = length(structure_entropy)
    structure_names = names(structure_entropy)

    alt_hip = c("two.sided", "less", "greater")[alternative]

    pval_mtx = matrix(NA, nrow = metrics_count, ncol = structures_count)
    rownames(pval_mtx) <- names(protein_entropy)
    colnames(pval_mtx) <- structure_names
    for (i in seq(1, metrics_count)) {
      for (j in seq(1, structures_count)) {
        reference = protein_entropy[[i]][-c(range, structure_entropy[[j]][[2]])]
        temp = structure_entropy[[j]][[1]][i,]
        pval_mtx[i, j] = stats::ks.test(reference, temp, alternative = alt_hip)$p.value
      }
    }

    if (make_plot == T) {
      colors = rainbow(structures_count)
      for (i in seq(1, metrics_count)) {
        par(new = F)
        for (j in seq(1, structures_count)) {
          reference = protein_entropy[[i]][-c(range, structure_entropy[[j]][[2]])]
          cumulative_distribution = stats::ecdf(reference)
          temp = ecdf(structure_entropy[[j]][[1]][i,])
          plot(
            cumulative_distribution,
            xlim = c(0, 1),
            ylim = c(0, 1),
            xlab = "Entropy",
            ylab = "Cumulative probability",
            col = scales::alpha("slategray", 0.7),
            main = paste(
              "CDF of",
              names(protein_entropy)[i],
              " for",
              pdb_name,
              "structures"
            )
          )
          par(new = T)
          plot(
            temp,
            xlim = c(0, 1),
            ylim = c(0, 1),
            col = scales::alpha(colors[j], 0.7),
            ylab = "",
            xlab = "",
            main = ""
          )
          legend(
            'bottomright',
            c(pdb_name, structure_names[j]),
            pch = c(16, 16),
            col = scales::alpha(c("slategray", colors[j]), 0.8)
          )
        }
      }
    }
    return(pval_mtx)
  }
RealValET_conservativity <- function(alignment) {
  ET_groups = seq(1, length(alignment$nam))
  dist <- seqinr::dist.alignment(alignment)
  clust <- stats::hclust(dist, method = "average") #upgma
  groups <- list() # groups for all nodes
  for (i in ET_groups) {
    groups[[i]] <- stats::cutree(clust, i)
  }
  aligned_matrix <- alignment2matrix(alignment = alignment)
  shannon_group <- function(group) {
    counts = table(group)
    len = length(group)
    freq = counts / len
    Shannon = -sum(freq * log(freq))
  }
  realValET <- function(column, group_tree) {
    nodes = length(group_tree) - 1
    node_group = rep(0, nodes)
    for (i in 1:nodes) {
      temp = split(column, group_tree[[i]])
      group = rep(0, length(temp))
      for (j in 1:length(temp)) {
        group[j] = shannon_group(temp[[j]])
      }
      node_group[i] = (1 / i) * sum(group)
    }
    RO = 1 + sum(node_group)
    return(RO)
  }
  x = rep(0, dim(aligned_matrix)[2])
  pb <- progress_bar$new(
    format = paste("Real valued ET score calculation", " [:bar] :percent eta: :eta"),
    total = dim(aligned_matrix)[2],
    clear = T,
    width = 80
  )
  for (i in 1:dim(aligned_matrix)[2]) {
    pb$tick()
    x[i] = realValET(aligned_matrix[, i], groups)
  }
  return(x / max(x))
}
