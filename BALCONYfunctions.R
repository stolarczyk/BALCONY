



# Conservation analysis ---------------------------------------------------
isUpper <- function(s) {
  return (all(grepl("[[:upper:]]", strsplit(s, "")[[1]])))
}

delete_isoforms <- function(alignment) {
  # This functions searches for isoforms in the alignment (entries with "-digit|" in the name) and deletes them
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
    new$nb = alignment$nb - (length(lines_to_delete))
    new$nam = alignment$nam[-lines_to_delete]
    new$seq = alignment$seq[-lines_to_delete]
    output = new
  }
  else{
    output = alignment
    warning("No isoforms detected")
  }
  return(output)
}
consensus <-  function(alignment, thresh) {
  #Function which calculates consensus
  #alignment-output of read.alignment()
  #thresh-given threshold of conservation (%)
  alignment_matrix = as.matrix(as.character(alignment[[3]]))
  
  count_cols = length(alignment_matrix)
  vec = seq(1, length(s2c(alignment_matrix[1, ])))
  len_of_vers = length(vec)
  mat = matrix("-", count_cols, len_of_vers)
  
  for (i in seq(1, count_cols)) {
    temp = s2c(alignment_matrix[i])
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
    if (value >= thresh) {
      consens[j] = m[1]
    }
  }
  return(consens)
}
cons2seqs_ident <-  function(alignment, consensus_seq) {
  #alignment_sequence- file[[3]]
  # number_of_seq- file[[1]]
  #consensus- calculated consensusus (output of my_consensusus())
  true_percentage = c()
  
  for (i in seq(1, alignment_parameters(alignment)$row_no)) {
    true_percentage[i] = length(which((consensus_seq == s2c(alignment$seq[i])) == TRUE)) / length(consensus_seq)
    # true percentage calculation (number of the same AAs in consensusus and in each sequence/number of all AAs)
  }
  return(true_percentage)
}
alignment_parameters <- function(alignment) {
  #alignment data
  #return list of alignment size [row_numbers, col_numbers]
  aligned_sequences = alignment$seq
  row_num = length(aligned_sequences)
  col_num = length(s2c(aligned_sequences[1]))
  param = list(row_no = row_num, col_no = col_num)
  return(param)
}

aligned_sequences_matrix2groups <-
  function(aligned_sequences_matrix,
           grouping_method) {
    rows = dim(aligned_sequences_matrix)[1]
    cols = dim(aligned_sequences_matrix)[2]
    aligned_sequences_matrixG = aligned_sequences_matrix
    if (grouping_method == 'general') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (general similarity)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "I" ||
              aligned_sequences_matrix[i, j] == "L" ||
              aligned_sequences_matrix[i, j] == "V")
            aligned_sequences_matrixG[i, j] = "G1"
          if (aligned_sequences_matrix[i, j] == "F" ||
              aligned_sequences_matrix[i, j] == "Y" ||
              aligned_sequences_matrix[i, j] == "W")
            aligned_sequences_matrixG[i, j] = "G2"
          if (aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "G")
            aligned_sequences_matrixG[i, j] = "G3"
          if (aligned_sequences_matrix[i, j] == "P")
            aligned_sequences_matrixG[i, j] = "G4"
          if (aligned_sequences_matrix[i, j] == "N" ||
              aligned_sequences_matrix[i, j] == "Q" ||
              aligned_sequences_matrix[i, j] == "S" ||
              aligned_sequences_matrix[i, j] == "T")
            aligned_sequences_matrixG[i, j] = "G5"
          if (aligned_sequences_matrix[i, j] == "R" ||
              aligned_sequences_matrix[i, j] == "H" ||
              aligned_sequences_matrix[i, j] == "K")
            aligned_sequences_matrixG[i, j] = "G6"
          if (aligned_sequences_matrix[i, j] == "D" ||
              aligned_sequences_matrix[i, j] == "E")
            aligned_sequences_matrixG[i, j] = "G7"
          if (aligned_sequences_matrix[i, j] == "M" ||
              aligned_sequences_matrix[i, j] == "C")
            aligned_sequences_matrixG[i, j] = "G8"
        }
      }
    }
    
    if (grouping_method == 'hydrophobicity') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to the hydrophobicity scale)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "I" ||
              aligned_sequences_matrix[i, j] == "L" ||
              aligned_sequences_matrix[i, j] == "V")
            aligned_sequences_matrixG[i, j] = "G1"
          if (aligned_sequences_matrix[i, j] == "F" ||
              aligned_sequences_matrix[i, j] == "C" ||
              aligned_sequences_matrix[i, j] == "M" ||
              aligned_sequences_matrix[i, j] == "A")
            aligned_sequences_matrixG[i, j] = "G2"
          if (aligned_sequences_matrix[i, j] == "G" ||
              aligned_sequences_matrix[i, j] == "S" ||
              aligned_sequences_matrix[i, j] == "T" ||
              aligned_sequences_matrix[i, j] == "W" ||
              aligned_sequences_matrix[i, j] == "P" ||
              aligned_sequences_matrix[i, j] == "Y")
            aligned_sequences_matrixG[i, j] = "G3"
          if (aligned_sequences_matrix[i, j] == "H" ||
              aligned_sequences_matrix[i, j] == "N" ||
              aligned_sequences_matrix[i, j] == "E" ||
              aligned_sequences_matrix[i, j] == "Q" ||
              aligned_sequences_matrix[i, j] == "D" ||
              aligned_sequences_matrix[i, j] == "K")
            aligned_sequences_matrixG[i, j] = "G4"
        }
      }
    }
    if (grouping_method == 'size') {
      for (i in seq(1, rows)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to size)
        for (j in seq(1, cols)) {
          if (aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "C" ||
              aligned_sequences_matrix[i, j] == "G" ||
              aligned_sequences_matrix[i, j] == "S")
            #TINY
            aligned_sequences_matrixG[i, j] = "G1"
          if (aligned_sequences_matrix[i, j] == "P" ||
              aligned_sequences_matrix[i, j] == "V" ||
              aligned_sequences_matrix[i, j] == "T" ||
              aligned_sequences_matrix[i, j] == "D" ||
              aligned_sequences_matrix[i, j] == "N")
            #SMALL
            aligned_sequences_matrixG[i, j] = "G2"
          if (aligned_sequences_matrix[i, j] == "A" ||
              aligned_sequences_matrix[i, j] == "G")
            aligned_sequences_matrixG[i, j] = "G3"
          else
            aligned_sequences_matrixG[i, j] = "G4"
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
            aligned_sequences_matrixG[i, j] = "G1"
          else
            aligned_sequences_matrixG[i, j] = "G2"
        }
      }
    }
    return(aligned_sequences_matrixG)
  }

cons2seqs_sim <-
  function(alignment,
           consensus_seq,
           grouping_method) {
    prmt = alignment_parameters(alignment = alignment)
    aligned_sequences_matrix = alignment2matrix(alignment)
    aligned_sequences_matrixG = aligned_sequences_matrix2groups(aligned_sequences_matrix, grouping_method)
    consensusG = rep(x = "-", length(consensus_seq))
    
    true_percentage_G = c()
    
    
    if (grouping_method == 'general') {
      for (i in seq(1, length(consensus_seq))) {
        if (consensus_seq[i] == "I" ||
            consensus_seq[i] == "L" || consensus_seq[i] == "V")
          consensusG[i] = "G1"
        if (consensus_seq[i] == "F" ||
            consensus_seq[i] == "Y" || consensus_seq[i] == "W")
          consensusG[i] = "G2"
        if (consensus_seq[i] == "A" || consensus_seq[i] == "G")
          consensusG[i] = "G3"
        if (consensus_seq[i] == "P")
          consensusG[i] = "G4"
        if (consensus_seq[i] == "N" || consensus_seq[i] == "Q")
          consensusG[i] = "G5"
        if (consensus_seq[i] == "R" ||
            consensus_seq[i] == "H" || consensus_seq[i] == "K")
          consensusG[i] = "G6"
        if (consensus_seq[i] == "D" || consensus_seq[i] == "E")
          consensusG[i] = "G7"
        if (consensus_seq[i] == "M" || consensus_seq[i] == "C")
          consensusG[i] = "G8"
      }
      for (i in seq(1, prmt$row_no)) {
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i, ])
        ) == TRUE)) / prmt$col_no
        # true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    
    if (grouping_method == 'hydrophobicity') {
      for (i in seq(1, length(consensus_seq))) {
        if (consensus_seq[i] == "I" ||
            consensus_seq[i] == "L" || consensus_seq[i] == "V")
          consensusG[i] = "G1"
        if (consensus_seq[i] == "F" ||
            consensus_seq[i] == "C" ||
            consensus_seq[i] == "M" || consensus_seq[i] == "A")
          consensusG[i] = "G2"
        if (consensus_seq[i] == "G" ||
            consensus_seq[i] == "S" ||
            consensus_seq[i] == "T" ||
            consensus_seq[i] == "W" ||
            consensus_seq[i] == "P" || consensus_seq[i] == "Y")
          consensusG[i] = "G3"
        if (consensus_seq[i] == "H" ||
            consensus_seq[i] == "N" ||
            consensus_seq[i] == "E" ||
            consensus_seq[i] == "Q" ||
            consensus_seq[i] == "D" || consensus_seq[i] == "K")
          consensusG[i] = "G4"
      }
      for (i in seq(1, prmt$row_no)) {
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i, ])
        ) == TRUE)) / prmt$col_no
        # true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    
    if (grouping_method == 'size') {
      for (i in seq(1, length(consensus_seq))) {
        if (consensus_seq[i] == "A" ||
            consensus_seq[i] == "C" ||
            consensus_seq[i] == "G" || consensus_seq[i] == "S")
          #TINY
          consensusG[i] = "G1"
        if (consensus_seq[i] == "P" ||
            consensus_seq[i] == "V" ||
            consensus_seq[i] == "T" ||
            consensus_seq[i] == "D" || consensus_seq[i] == "N")
          #SMALL
          consensusG[i] = "G2"
        if (consensus_seq[i] == "A" || consensus_seq[i] == "G")
          consensusG[i] = "G3"
        else
          consensusG[i] = "G4"
      }
      for (i in seq(1, prmt$row_no)) {
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i, ])
        ) == TRUE)) / prmt$col_no
        # true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    
    if (grouping_method == 'aromaticity') {
      for (i in seq(1, length(consensus))) {
        if (consensus_seq[i] == "F" ||
            consensus_seq[i] == "W" ||
            consensus_seq[i] == "H" || consensus_seq[i] == "Y")
          consensusG[i] = "G1"
        else
          consensusG[i] = "G2"
      }
      for (i in seq(1, prmt$row_no)) {
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i, ])
        ) == TRUE)) / prmt$col_no
        # true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    return(true_percentage_G)
  }

alignment2matrix <- function(alignment) {
  #alignment data
  #returns alignment as a matrix
  prmt = alignment_parameters(alignment = alignment)
  aligned_sequences_matrix = matrix("-", prmt$row_no, prmt$col_no)
  
  for (i in seq(1, prmt$row_no)) {
    #Putting aligned seqs into matrix
    temp = s2c(alignment$seq[i])
    for (j in seq(1, prmt$col_no)) {
      aligned_sequences_matrix[i, j] = temp[j]
      
    }
  }
  return(aligned_sequences_matrix)
}
calculate_AA_variation <-
  function(alignment,
           threshold,
           grouped = F,
           grouping_method = "general") {
    #prmt- size of alignment (output of get_parameter())
    #sequence_alignment-file[[3]]
    #threshold-threshold for detecting key amino acids (the percentage of all at the given position)
    #returns list of matrices with tabelarised symbols of the most common AA in alignment column and percentage values for contributed AA
    prmt = alignment_parameters(alignment = alignment)
    keyaas_treshold = prmt$row_no * (threshold / 100)
    
    aligned_sequences_matrix = alignment2matrix(alignment = alignment)
    
    if (grouped == T) {
      aligned_sequences_matrix = aligned_sequences_matrix2groups(aligned_sequences_matrix, grouping_method = grouping_method)
    }
    keyaas = matrix("n", dim(aligned_sequences_matrix)[2], 20 * 2)
    
    keyaas_per = matrix("n", dim(aligned_sequences_matrix)[2], 20 * 2)
    
    
    for (i in seq(1, dim(aligned_sequences_matrix)[2])) {
      table = (sort(table(aligned_sequences_matrix[, i]), decreasing = T))
      #Czestosci i rodzaje wystepowania AA
      aas = names(table)
      
      keyaas[i, 1:length(matrix(aas[which(table >= keyaas_treshold)], 1, length(aas[which(table >= keyaas_treshold)])))] = matrix(aas[which(table >= keyaas_treshold)], 1, length(aas[which(table >= keyaas_treshold)]))
      #jeżeli są jakieś AA (=zawsze) to wpisywane są one do macierzy keyaas (odpowiendia ilość za sprawą sprawdzenia)
      keyaas_per[i, 1:length(matrix((table)[which(table >= keyaas_treshold)], 1, length((table)[which(table >= keyaas_treshold)])))] = matrix(round((table)[which(table >= keyaas_treshold)] / prmt$row_no, 3) * 100, 1, length(aas[which(table >= keyaas_treshold)]))
      #podobne wpisanie w odpowiedi wiersz macierzy dla ilości AA
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
noteworthy_sequences <- function(percentage, alignment) {
  max = which.max(percentage)
  namelist = alignment[[2]]
  out.max = list(c(namelist[max], max)) #output is a name of sequence and position in alignment
  min = which.min(percentage) #percentage is an output of calculate_group_consensus or calculate_tru_consensus
  out.min = list(c(namelist[min], min)) #output is a name of sequence and position in alignment
  value = which.max(table(percentage))
  vector = which(percentage %in% names(value))
  names = match(vector, file[[2]])
  out.common = namelist[vector]
  out = list(
    best_consensus = out.max,
    worst_consensus = out.min,
    most_common = out.common
  )
  
  return(out)
}
convert <- function(amino_acids) {
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
find_seqid <- function(sequence, library) {
  ktory = c()
  
  for (i in seq(1, length(library))) {
    ktory[i] = length(which((library[[i]] == sequence) == TRUE))
    
    ktory1 = which(ktory == 1)
    
  }
  seqid = library[[ktory1]][1]
  
  return(seqid)
}
get_seq_names <- function(alignment) {
  names = alignment[[2]]
  
  return(names)
}
find_seq <- function(sequence_id, alignment_file) {
  mat = as.matrix(as.character(alignment_file[[3]]))
  
  nazwy_mat = get_seq_names(alignment_file)
  
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
  seqs_character = length(which(s2c(seqs) != "-") == T)
  
  seq = list(sequence = seqs, len = seqs_character)
  return(seq)
}
create_structure_seq <-
  function(structure_file,
           sequence_id,
           alignment,
           shift) {
    #structure_file-> list of structures in protein
    #sequence_id -> uniprot id  which has been found by read_file(filename="PDBid") with PDB indentifier ;
    #alignment-> file with alignment (alignment.fasta)
    seq = list()
    
    structure = list()
    
    #finds an uniprot name from alignent- extract uniprot seq from alignment
    base_seq = find_seq(sequence_id, alignment)
    for (i in seq(1, length(structure_file), by = 1)) {
      structures_indices = as.vector(structure_file[[i]][[1]])
      
      structures_names = as.vector(structure_file[[i]][[2]])
      
      structures_indices = structures_indices[-1]
      
      structure_idx = as.numeric(structures_indices)
      structures_names = structures_names[-1]
      
      if (shift != 0) {
        structure_idx = structure_idx + rep(shift, length(structures_indices))
      }
      seq[[i]] = rep("N", each = base_seq$len)
      
      seq[[i]][structure_idx] = "T"
      aa_positions = which(s2c(base_seq$sequence) != "-")
      just_align = alignment[[3]]
      length_alignment = alignment_parameters(alignment)$col_no
      structure[[i]] = rep("-", each = length_alignment)
      
      j = 1
      
      for (a in aa_positions) {
        #aligning the structures information with the alignment sequence
        structure[[i]][a] = seq[[i]][j]
        
        j = j + 1
        
      }
    }
    
    struc_length = length(structure[[1]])
    count = length(structure_file)
    
    output = matrix("-", count, struc_length)
    
    v = seq(1, count, by = 1)
    for (i in v) {
      output[i, ] = structure[[i]]
    }
    nr_stru = rep("-", length(structure[[1]]))
    
    j = 1
    
    for (i in seq(1, length(structure[[1]]), 1)) {
      if (structure[[1]][i] != "-") {
        nr_stru[i] = j
        j = j + 1
      }
    }
    rownames(output) = names(structure_file)
    return(list(structure_matrix = output, structure_numbers = nr_stru))
  }

barplotshow <- function(position, AA_variation) {
  #position=position-1;
  row = max(which(AA_variation$percentage[, position] != "n"))
  #checking how many different AA is on the position
  x_names = AA_variation$AA[seq(1, row), position]
  x_val = round(as.numeric(AA_variation$percentage[seq(1, row), position]))
  plot = barplot(
    x_val,
    names.arg = x_names,
    cex.main = 1.8,
    cex.lab = 1.3,
    cex.names =
      1.5,
    cex.axis = 1.5,
    width = rep(1.1, row),
    xlim = c(0, dim(AA_variation[[1]])[1]),
    ylim = c(0, max(x_val) +
               10),
    ylab = "Percentage",
    main = bquote(Amino ~
                    acid ~ variation ~ on ~ the ~ .(position) ~ position)
  )
  
  text(
    x = plot,
    y = x_val / 2,
    cex = 1.0,
    labels = paste(as.character(x_val), "%", sep = ""),
    xpd =
      T
  )
}
create_final_CSV <-
  function(FILENAME,
           AA_variations,
           structure,
           uniprot,
           alignment,
           list_of_scores = NULL) {
    sequence = s2c(find_seq(uniprot, alignment)$sequence)
    structure_output = rbind(structure$structure_matrix, structure$structure_numbers)
    structure_output_names = append(rownames(structure$structure_matrix), "Structure numbers")
    rownames(structure_output) = structure_output_names
    rownames(AA_variations$matrix) = rep(c("AA name", "Percentage"), dim(variations_matrix)[1] /
                                           2)
    if (is.null(list_of_scores)) {
      alignment_position = seq(1,dim(AA_variations$matrix)[2],by = 1)
      final_output = rbind(alignment_position, AA_variations$matrix, sequnece, structure_output)
    }
    else{
      #Dodawanie wynikow konserwatywnosci tylko takich, jakie podal user
      scores_mtx = matrix(NA,
                          nrow = length(list_of_scores),
                          ncol = length(list_of_scores[[1]]))
      scores_mtx_names = c()
      
      for (i in seq(1, length(list_of_scores))) {
        scores_mtx[i, ] = list_of_scores[[i]]
        scores_mtx_names[i] = names(list_of_scores)[i]
        
      }
      rownames(scores_mtx) = scores_mtx_names
      alignment_position = seq(1,dim(AA_variations$matrix)[2],by = 1)
      final_output = rbind(alignment_position,
                           AA_variations$matrix,
                           sequence,
                           structure_output,
                           scores_mtx)
      
    }
    files_no = ceiling(dim(final_output)[2] / 1000)
    for (i in seq(1, files_no, 1)) {
      if (!i == files_no) {
        output = final_output[, ((i - 1) * 1000 + 1):(i * 1000)]
      }
      else{
        output = final_output[, ((i - 1) * 1000 + 1):dim(final_output)[2]]
      }
      write.table(output,file = paste(FILENAME, "_", i, ".csv", sep = ""),
                row.names = T,col.names = F,sep = ",")
    }
    if (Sys.info()[[1]] == "Linux") {
      print(paste("Output written to: ", getwd(), "/", FILENAME, ".csv", sep = ""))
    }
    return(final_output)
  }

TG_conservativity <- function(alignment) {
  var_aa = calculate_AA_variation(alignment, threshold = 0.01)
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
  # change NaNs to 0
  part_con = pre_conservativity
  
  part_conserv = part_con / max(part_con)
  TG = -(log(part_conserv))
  TG[which(is.infinite(TG[1]))] = 0 # change Infs to 0
  TG_score = (TG / max(TG))
  return_data = TG_score
  
  return(return_data)
}
conservativity <- function(aligned_sequences_matrix) {
  suma = rep(NaN, dim(aligned_sequences_matrix)[2])
  suma_schneider = rep(NaN, dim(aligned_sequences_matrix)[2])
  Kabat = rep(NaN, dim(aligned_sequences_matrix)[2])
  for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
    column = aligned_sequences_matrix[, rep]
    #KABAT
    K = c()
    sum_wew = c()
    p = c()
    n = c()
    tab = c()
    sum_wew_schneider = c()
    
    N = length(column)
    #no of objects
    K = length(as.numeric(table(column)))
    #no of classes
    n1 = as.numeric(sort(table(column))[length(table(column))])
    
    Kabat[rep] = (K / n1) * N
    
    #SCHNEIDER, SHANNON
    tab = as.numeric(table(column))
    
    suma_wew = 0
    sum_wew_schneider = 0
    
    for (i in seq(1, K, by = 1)) {
      n[i] = tab[i]
      
      p[i] = n[i] / N
      
      res = p[i] * log2(p[i])
      
      res_schneider = p[i] * log(p[i]) * (1 / log(21))
      
      sum_wew = suma_wew + res
      
      sum_wew_schneider = sum_wew_schneider + res_schneider
      
      
    }
    suma[rep] = -sum_wew
    suma_schneider[rep] = -sum_wew_schneider
    
  }
  Kabat_entropy_normalized = Kabat / max(Kabat)
  
  ret = list(Shannon = suma,
             Schneider = suma_schneider,
             Kabat = Kabat_entropy_normalized)
  
  return(ret)
}

kabat_conservativity <- function(alignment_file) {
  aligned_sequences_matrix = alignment2matrix(alignment = alignment_file)
  Kabat = rep(NaN, dim(aligned_sequences_matrix)[2])
  for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
    column = aligned_sequences_matrix[, rep]
    #KABAT
    K = c()
    p = c()
    n = c()
    tab = c()
    N = length(column)
    #no of objects
    K = length(as.numeric(table(column)))
    #no of classes
    n1 = as.numeric(sort(table(column))[length(table(column))])
    
    Kabat[rep] = (K / n1) * N
    
  }
  Kabat_entropy_normalized = Kabat / max(Kabat)
  
  return(Kabat_entropy_normalized)
}

schneider_conservativity <- function(alignment_file) {
  aligned_sequences_matrix = alignment2matrix(alignment = alignment_file)
  sum_schneider = rep(NaN, dim(aligned_sequences_matrix)[2])
  for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
    column = aligned_sequences_matrix[, rep]
    #SCHNEIDER, SHANNON
    tab = as.numeric(table(column))
    K = c()
    p = c()
    n = c()
    N = length(column)
    K = length(as.numeric(table(column)))
    sum_internal_schneider = 0
    for (i in seq(1, K, by = 1)) {
      n[i] = tab[i]
      
      p[i] = n[i] / N
      
      res = p[i] * log2(p[i])
      
      res_schneider = p[i] * log(p[i]) * (1 / log(21))
      
      sum_internal_schneider = sum_internal_schneider + res_schneider
      
    }
    sum_schneider[rep] = -sum_internal_schneider
    
  }
  return(sum_schneider)
}

shannon_conservativity <- function(alignment_file) {
  aligned_sequences_matrix = alignment2matrix(alignment_file)
  sum = rep(NaN, dim(aligned_sequences_matrix)[2])
  for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
    column = aligned_sequences_matrix[, rep]
    #SCHNEIDER, SHANNON
    tab = as.numeric(table(column))
    K = c()
    p = c()
    n = c()
    N = length(column)
    K = length(as.numeric(table(column)))
    sum_internal = 0
    for (i in seq(1, K, by = 1)) {
      n[i] = tab[i]
      
      p[i] = n[i] / N
      
      res = p[i] * log2(p[i])
      
      sum_internal = sum_internal + res
      
    }
    sum[rep] = -sum_internal
  }
  return(sum)
}

weights <- function(file,
                    aligned_sequences_matrix,
                    threshold) {
  consensus_seq = consensus(file, threshold)
  
  consensus_seq = consensus_seq[-which(consensus_seq == "\r")]
  
  length_cons_seq = length(consensus_seq)
  
  for (i in seq(1, dim(aligned_sequences_matrix)[1], by = 1)) {
    eq = length(which(consensus_seq == aligned_sequences_matrix[i, ]))
    
    score[i] = eq / length_cons_seq
    
  }
  return(score)
}
substitution_mtx <- function (matrix_name) {
  # function can read .txt file substitution matrix and convert it
  # into a list of alphabet (the range of letters in matrix) and
  # mtx(valuse into substitution matrix)
  matrix = read.table(matrix_name)
  mtx = matrix[-1, -1]
  alphabet = as.character(unlist(matrix[1, -1]))
  sub_mtx = list(alphabet, mtx)
  return(sub_mtx)
}
D_matrix <- function(sub_mtx) {
  values = as.numeric(as.matrix(sub_mtx[[2]]))
  dim(values) <- dim(sub_mtx[[2]])
  k = dim(sub_mtx[[2]])[1]
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
  output = list(sub_mtx[[1]], distance)
  return(output)
}
landgraf_conservativity <-
  function(matrix_name = NULL, alignment_file, weights) {
    if (is.null(matrix_name)) {
      load("Gonnet_mtx.rd") # Proper path?
      pre_dissim_mtx = gonnet_mtx
    }
    else{
      pre_dissim_mtx = substitution_mtx(matrix_name)
    }
    aligned_sequences_matrix = alignment2matrix(alignment = alignment_file)
    dissim_mtx = D_matrix(pre_dissim_mtx)
    conservation = rep(NaN, dim(aligned_sequences_matrix)[2])
    status = 0
    
    for (rep in seq(1, dim(aligned_sequences_matrix)[2], 1)) {
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
      
      if (round((rep / dim(aligned_sequences_matrix)[2]) * 100) != status) {
        print(paste("Position: ", rep, ", ", round((
          rep / dim(aligned_sequences_matrix)[2]
        ) * 100), "% DONE", sep = ""))
        
        status = round((rep / dim(aligned_sequences_matrix)[2]) * 100)
      }
    }
    Landgraf_normalized_entropy = conservation / max(conservation)
    return(Landgraf_normalized_entropy)
  }
sequence_stats <-
  function(alignment_file,
           uniprot,
           landgraf,
           schneider,
           TG) {
    sequence = s2c(find_seq(uniprot, alignment_file)$sequence)
    
    alignment_positions = which(sequence != "-")
    sequence = sequence[alignment_positions]
    
    protein_positions = seq(1, length(sequence))
    
    landgraf_metric = landgraf[alignment_positions]
    
    schneider_metric = schneider[alignment_positions]
    
    TG_metric = TG[alignment_positions]
    
    return_table = matrix("", 6, length(sequence))
    
    return_table[1, ] = sequence
    
    return_table[2, ] = protein_positions
    
    return_table[3, ] = alignment_positions
    
    return_table[4, ] = landgraf_metric
    
    return_table[5, ] = schneider_metric
    
    return_table[6, ] = TG_metric
    
    rownames(return_table) = c(
      "sequence",
      "protein positions",
      "alignment positions",
      "landraf metric entropy score",
      "schneider metric entropy score",
      "TG metric entropy score"
    )
    
    return_list = list()
    
    return_list[[1]] = sequence
    
    return_list[[2]] = protein_positions
    
    return_list[[3]] = alignment_positions
    
    return_list[[4]] = landgraf_metric
    
    return_list[[5]] = schneider_metric
    
    return_list[[6]] = TG_metric
    
    names(return_list) = c(
      "sequence",
      "protein positions",
      "alignment positions",
      "landraf metric entropy score",
      "schneider metric entropy score",
      "TG metric entropy score"
    )
    
    return_set = list(return_list, return_table)
    
    names(return_set) = c("list", "matrix")
    
    write.csv(x = return_table, file = "sequence_stats.csv")
    print(paste("Table saved in: ", getwd(), "sequence_stats.csv", sep = ""))
    return(return_set)
    
  }
entropy_profile <-
  function(tunnel_file,
           sequence_id,
           alignment_file,
           shift,
           prot_entropy,
           index) {
    #tunnel_file-> list of tunnels in protein
    #sequence_id -> uniprot id  which has been found by read_file(filename="PDBid") with PDB indentifier ;
    #alignment_file-> file wiht alignment (alignment.fst)
    #Shift-> shift of AA from caver to reference UniProt sequence
    #index of tunne;
    seq = list()
    
    
    index
    #finds an uniprot name from alignent- extract uniprot seq from alignment
    base_seq = find_seq(sequence_id, alignment_file)
    tunnels_indices = as.vector(tunnel_file[[index]][[1]])
    
    tunnels_names = as.vector(tunnel_file[[index]][[2]])
    
    tunnels_indices = tunnels_indices[-1]
    
    tunnel_idx = as.numeric(tunnels_indices)
    tunnels_names = tunnels_names[-1]
    
    if (shift != 0) {
      tunnel_idx = tunnel_idx + rep(shift, length(tunnels_indices))
    }
    #profile=tunnel_idx
    #for(i in seq(1:length(profile))){
    # id=tunnel_idx[i]
    # profile[i]=prot_entropy[id]
    profile = prot_entropy[tunnel_idx]
    output = list(profile, tunnel_idx)
    
    return(output)
  }


# Structure analysis ------------------------------------------------------
get_remarks465_pdb <- function(pdb_path,chain_identifier){
  #pdb_path: a path to the pdb file
  #chain_identifier: a character spcifying the chain to analyze
  #returns: list of 1) numbers of aminocids which are missing 2) chain identifier
  require(Rpdb)
  pdb_file = read.pdb(file = pdb_path, REMARK = T,ATOM = T, CRYST1 = F,TITLE = T,MODEL = F,HETATM = F,CONECT = F)
  pdb_file = pdb_file$remark
  pattern = "REMARK 465"
  indices = grep(perl = T,
              pattern = pattern,
              x = pdb_file)
  remarks = pdb_file[indices[8:length(indices)]]
  aa_number=c()
  chain=c()
  for(i in seq(1,length(remarks),by = 1)){
    remark_line = strsplit(remarks[i],split = " ")
    x = remark_line[[1]]
    aa_number[i] = x[which(x!="")][length(x[which(x!="")])]
    chain[i] = x[which(x!="")][length(x[which(x!="")])-1]
  }
  aa_number=as.numeric(aa_number)
  chain_indices = which(chain==chain_identifier)
  aa_number=aa_number[chain_indices]
  chain=chain[chain_indices]
  return(list(aa_numbers = aa_number,chain=chain))
}

correct_structure_seq_missing <- function(structure,pdb_path,chain_identifier){
  
}
missing_aa2list <- function(remark465_data){
  #cdn
}
get_structure_idx <- function(structure) {
  #documentation get_structure_idx.Rd
  #get idx of structure in alignemnt
  #return sorted list of indices in MSA, one list element is one structure
  #the first element is indices of whole protein in alignment
  stru_index = list()
  for (i in seq(1:length(structure))) {
    stru_index[[i]] = which(structure[[i]] == "T")
  }
  whole_prot = which(structure[[1]] != "-")
  
  out = list(proteinIndices = whole_prot, structureIndices = stru_index)
  names(out[[2]]) = names(structure)
  return(out)
}

get_prot_entropy <- function(whole_prot, score_list) {
  #documentation get_prot_entropy.Rd
  #allows to get idx of whole protein in alignment
  #returns list of entropy for protein
  
  prot_cons = list()
  for (i in seq(1, length(score_list))) {
    prot_cons[[i]] = score_list[[i]][whole_prot]
  }
  names(prot_cons) <- names(score_list)
  return(prot_cons)
}
plot_entropy <- function(prot_cons, colors, impose = NULL) {
  #plots scores on one plot, if colors are not specified plot as rainbow
  #automagically uses name of pdb FIXME
  # recive list of entropy scores and list of Names in this list
  if (missing(colors)) {
    colors <- rainbow(length(prot_cons))
  }
  if (is.null(impose)) {
    impose <- T
  }
  
  plot(
    prot_cons[[1]],
    ylim = c(0, 1),
    col = colors[1],
    xlab = "amino acid",
    ylab = "entropy score",
    main = paste("entropy score for ", pdb_name),
    type = "l"
  )
  for (i in seq(2, length(prot_cons))) {
    par(new = impose)
    plot(
      prot_cons[[i]],
      ylim = c(0, 1),
      col = colors[i],
      xlab = "",
      ylab = "",
      main = "",
      type = "l"
    )
  }
  legend("topleft",
         names(prot_cons),
         col = colors,
         lty = c(1))
}
get_structures_entropy <- function(structure_index, score_list) {
  #structure_index is a list of indexes in alignment of protein and structures
  #score_list list of entropies for whole alignment
  #output is a list of matrixes where each row contains values of entropy for AA in structure
  t_index = structure_index$structureIndices
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
  return(Entropy)
}
entropy_for_all <-
  function(tunnel_file,
           uniprot,
           file,
           shift,
           prot_cons) {
    profilet = list()
    stru_numb = length(tunnel_file)
    names <- paste(names(prot_cons))
    #k- possible entropy score iterator
    megalist = list()
    for (k in seq(1, length(prot_cons))) {
      for (i in seq(1, stru_numb)) {
        profilet[[i]] = entropy_profile(tunnel_file, uniprot, file, shift, prot_cons[[k]], i)
        #profilet[[i]]<- setNames(mapply(paste("stru",i))
      }
      megalist[[names[k]]] <- profilet
      
    }
    return(megalist)
  }

plot_structure_on_protein <-
  function(protein_entropy,
           structure_profils,
           pdb_name,
           colors,
           structure_names = NULL) {
    #protein entropy- list of an entropy values (with different scores)
    #structure_profils - list of entropy scores for structures (as list) where [[1]] entropy value
    # and [[2]] index in protein
    protein_entropy = prot_cons
    structure_profils = profils_for_structure
    
    if (length(protein_entropy) == length(structure_profils)) {
      score_names = names(protein_entropy)
      prot_lenght = length(protein_entropy[[1]])
      StruLen = length(structure_profils[[1]])
      if (missing(colors)) {
        colors = rainbow(StruLen)
      }
      if (is.null(structure_names)) {
        structure_names = c()
        for (i in seq(1, StruLen)) {
          structure_names[i] = paste("stru", i)
        }
      }
      for (i in seq(1, length(protein_entropy))) {
        plot(
          protein_entropy[[i]],
          col = "black",
          type = "l",
          main = paste(score_names[i], " score for ", pdb_name),
          pch = 20,
          xlim = c(0, prot_lenght),
          ylim = c(0.0, 1.0),
          xlab = 'Amino Acid',
          ylab = 'Entropy'
        )
        for (j in seq(1, StruLen)) {
          par(new = T)
          plot(
            structure_profils[[i]][[j]][[2]],
            structure_profils[[i]][[j]][[1]],
            col = colors[[j]],
            main = "",
            pch = j,
            xlim = c(0, prot_lenght),
            ylim = c(0.0, 1.0),
            xlab = '',
            ylab = ''
          )
        }
        legend(
          'topleft',
          c(pdb_name, structure_names),
          lty = c(1, rep(0, j)),
          pch = c(-1, seq(1, j)),
          lwd = c(2.5, 2.5),
          col = c("black", colors)
        )
        
        
      }
    }
    else
      print("The lists contain different number of conservation/entropy scores!")
  }

compare_cons_metrics <-
  function(protein, structure_cons, pdb_name) {
    #protein- list of conservation scores for proteins
    #structure cons- list of conservation scores for [[list]] of strustures with [[list]] scores and indices
    #pdb_name - protein name
    metrics_count = length(protein)
    structures_count = length(structure_cons[[1]])
    colors = rainbow(structures_count)
    structure_names = c()
    for (i in seq(1, structures_count)) {
      structure_names[i] = paste("stru", i)
    }
    for (i in seq(1, metrics_count)) {
      for (j in seq(1, metrics_count)) {
        if (i != j) {
          plot(
            protein[[i]],
            protein[[j]],
            main = paste(
              'Scatterplot of',
              names(protein)[i],
              "vs. ",
              names(protein)[j]
            ),
            xlim = c(0, 1),
            ylim = c(0, 1),
            xlab = names(protein)[i],
            ylab = names(protein)[j],
            pch = 20,
            col = alpha("slategray", 0.5)
          )
          for (k in seq(1, structures_count)) {
            par(new = T)
            plot(
              structure_cons[[i]][[k]][[1]],
              structure_cons[[j]][[k]][[1]],
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
            col = alpha(c("slategray", colors), 0.5)
          )
        }
      }
    }
  }
