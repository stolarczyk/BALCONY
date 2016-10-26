consensus <-
  function(alignment, thresh) {
    #Function which calculates consensus
    #alignment-output of read.alignment()
    #thresh-given threshold of conservation (%)
    alignment_matrix = as.matrix(as.character(alignment[[3]]));
    count_cols = length(alignment_matrix)
    vec = seq(1,length(s2c(alignment_matrix[1,])))
    len_of_vers = length(vec)
    mat = matrix("-",count_cols, len_of_vers)
    
    for (i in seq(1,count_cols)) {
      temp = s2c(alignment_matrix[i])
      for (j in vec) {
        mat[i,j] = temp[j]
      }
    }
    
    consens = c(rep("*", each = len_of_vers))
    for (j in vec) {
      alignment_col = mat[,j]
      m = names(sort(table(alignment_col), decreasing = T))
      num = sort(table(alignment_col),decreasing = T)[1]
      value = num / length(alignment_col) * 100
      if (value >= thresh) {
        consens[j] = m[1]
      }
    }
    return(consens)
  }
cons2seqs_ident <-
  function(alignment_sequence, number_of_seq, consensusus_seq) {
    #alignment_sequence- file[[3]]
    # number_of_seq- file[[1]]
    #consensus- calculated consensusus (output of my_consensusus())
    true_percentage = c();
    for (i in seq(1,number_of_seq)) {
      true_percentage[i] = length(which((
        consensus_seq == s2c(alignment_sequence[i])
      ) == TRUE)) / length(consensus_seq);# true percentage calculation (number of the same AAs in consensusus and in each sequence/number of all AAs)
    }
    return(true_percentage)
  }
alignment_parameters <- function(alignment_sequences) {
  #alignment_sequences= file[[3]]
  #return list of alignment size [row_numbers, col_numbers]
  row_num = length(alignment_sequences)
  col_num = length(s2c(alignment_sequences[1]))
  param = list(row_no = row_num, col_no = col_num)
  return(param)
}
cons2seqs_sim <-
  function(prmt, aligned_sequences_matrix, consensus_seq, grouping_method) {
    aligned_sequences_matrixG = aligned_sequences_matrix
    consensusG = rep(x = "-",length(consensus_seq));
    true_percentage_G = c();
    
    if (grouping_method == 'general') {
      for (i in seq(1,length(consensus_seq))) {
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
      for (i in seq(1,prmt$row_no)) {
        #Dividing AAs into groups of similar AAs in the matrices (general similarity)
        for (j in seq(1,prmt$col_no)) {
          if (aligned_sequences_matrix[i,j] == "I" ||
              aligned_sequences_matrix[i,j] == "L" ||
              aligned_sequences_matrix[i,j] == "V")
            aligned_sequences_matrixG[i,j] = "G1"
          if (aligned_sequences_matrix[i,j] == "F" ||
              aligned_sequences_matrix[i,j] == "Y" ||
              aligned_sequences_matrix[i,j] == "W")
            aligned_sequences_matrixG[i,j] = "G2"
          if (aligned_sequences_matrix[i,j] == "A" ||
              aligned_sequences_matrix[i,j] == "G")
            aligned_sequences_matrixG[i,j] = "G3"
          if (aligned_sequences_matrix[i,j] == "P")
            aligned_sequences_matrixG[i,j] = "G4"
          if (aligned_sequences_matrix[i,j] == "N" ||
              aligned_sequences_matrix[i,j] == "Q")
            aligned_sequences_matrixG[i,j] = "G5"
          if (aligned_sequences_matrix[i,j] == "R" ||
              aligned_sequences_matrix[i,j] == "H" ||
              aligned_sequences_matrix[i,j] == "K")
            aligned_sequences_matrixG[i,j] = "G6"
          if (aligned_sequences_matrix[i,j] == "D" ||
              aligned_sequences_matrix[i,j] == "E")
            aligned_sequences_matrixG[i,j] = "G7"
          if (aligned_sequences_matrix[i,j] == "M" ||
              aligned_sequences_matrix[i,j] == "C")
            aligned_sequences_matrixG[i,j] = "G8"
        }
        
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i,])
        ) == TRUE)) / prmt$col_no;# true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    
    if (grouping_method == 'hydrophobicity') {
      for (i in seq(1,length(consensus_seq))) {
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
      for (i in seq(1,prmt$row_no)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to the hydrophobicity scale)
        for (j in seq(1,prmt$col_no)) {
          if (aligned_sequences_matrix[i,j] == "I" ||
              aligned_sequences_matrix[i,j] == "L" ||
              aligned_sequences_matrix[i,j] == "V")
            aligned_sequences_matrixG[i,j] = "G1"
          if (aligned_sequences_matrix[i,j] == "F" ||
              aligned_sequences_matrix[i,j] == "C" ||
              aligned_sequences_matrix[i,j] == "M" ||
              aligned_sequences_matrix[i,j] == "A")
            aligned_sequences_matrixG[i,j] = "G2"
          if (aligned_sequences_matrix[i,j] == "G" ||
              aligned_sequences_matrix[i,j] == "S" ||
              aligned_sequences_matrix[i,j] == "T" ||
              aligned_sequences_matrix[i,j] == "W" ||
              aligned_sequences_matrix[i,j] == "P" ||
              aligned_sequences_matrix[i,j] == "Y")
            aligned_sequences_matrixG[i,j] = "G3"
          if (aligned_sequences_matrix[i,j] == "H" ||
              aligned_sequences_matrix[i,j] == "N" ||
              aligned_sequences_matrix[i,j] == "E" ||
              aligned_sequences_matrix[i,j] == "Q" ||
              aligned_sequences_matrix[i,j] == "D" ||
              aligned_sequences_matrix[i,j] == "K")
            aligned_sequences_matrixG[i,j] = "G4"
        }
        
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i,])
        ) == TRUE)) / prmt$col_no;# true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    
    if (grouping_method == 'size') {
      for (i in seq(1,length(consensus_seq))) {
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
      
      for (i in seq(1,prmt$row_no)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to size)
        for (j in seq(1,prmt$col_no)) {
          if (aligned_sequences_matrix[i,j] == "A" ||
              aligned_sequences_matrix[i,j] == "C" ||
              aligned_sequences_matrix[i,j] == "G" ||
              aligned_sequences_matrix[i,j] == "S")
            #TINY
            aligned_sequences_matrixG[i,j] = "G1"
          if (aligned_sequences_matrix[i,j] == "P" ||
              aligned_sequences_matrix[i,j] == "V" ||
              aligned_sequences_matrix[i,j] == "T" ||
              aligned_sequences_matrix[i,j] == "D" ||
              aligned_sequences_matrix[i,j] == "N")
            #SMALL
            aligned_sequences_matrixG[i,j] = "G2"
          if (aligned_sequences_matrix[i,j] == "A" ||
              aligned_sequences_matrix[i,j] == "G")
            aligned_sequences_matrixG[i,j] = "G3"
          else
            aligned_sequences_matrixG[i,j] = "G4"
        }
        
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i,])
        ) == TRUE)) / prmt$col_no;# true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    
    if (grouping_method == 'aromaticity') {
      for (i in seq(1,length(consensus))) {
        if (consensus_seq[i] == "F" ||
            consensus_seq[i] == "W" ||
            consensus_seq[i] == "H" || consensus_seq[i] == "Y")
          consensusG[i] = "G1"
        else
          consensusG[i] = "G2"
      }
      for (i in seq(1,prmt$row_no)) {
        #Dividing AAs into groups of similar AAs in the matrices (according to the aromaricity)
        for (j in seq(1,prmt$col_no)) {
          if (aligned_sequences_matrix[i,j] == "F" ||
              aligned_sequences_matrix[i,j] == "W" ||
              aligned_sequences_matrix[i,j] == "H" ||
              aligned_sequences_matrix[i,j] == "Y")
            aligned_sequences_matrixG[i,j] = "G1"
          else
            aligned_sequences_matrixG[i,j] = "G2"
        }
        
        true_percentage_G[i] = length(which((
          consensusG == (aligned_sequences_matrixG[i,])
        ) == TRUE)) / prmt$col_no;# true percentage calculation (number of the same AAs in consensus (grouped AAs) and in each sequence/number of all AAs)
      }
    }
    return(true_percentage_G)
  }
alignment2matrix <- function(prmt, aligned_sequences) {
  #prmt= output of get_parameters
  #alignment_sequences=file[[3]]
  #returns alignment as a matrix
  aligned_sequences_matrix = matrix("-", prmt$row_no, prmt$col_no);
  for (i in seq(1,prmt$row_no)) {
    #Putting aligned seqs into matrix
    temp = s2c(aligned_sequences[i])
    for (j in seq(1,prmt$col_no)) {
      aligned_sequences_matrix[i,j] = temp[j];
    }
  }
  return(aligned_sequences_matrix)
}
calculate_AA_variation <-
  function(prmt, sequence_alignment, threshold) {
    #prmt- size of alignment (output of get_parameter())
    #sequence_alignment-file[[3]]
    #threshold-threshold for detecting key amino acids (the percentage of all at the given position)
    #returns list of matrices with tabelarised symbols of the most common AA in alignment column and percentage values for contributed AA
    keyaas_treshold = prmt$row_no * (threshold / 100);
    aligned_sequences_matrix = alignment2matrix(prmt,sequence_alignment);
    keyaas = matrix("n",dim(aligned_sequences_matrix)[2],20 * 2);
    keyaas_per = matrix("n",dim(aligned_sequences_matrix)[2],20 * 2);
    
    for (i in seq(1,dim(aligned_sequences_matrix)[2])) {
      table = (sort(table(aligned_sequences_matrix[,i]),decreasing = T)); #Czestosci i rodzaje wystepowania AA
      aas = names(table);
      keyaas[i,1:length(matrix(aas[which(table >= keyaas_treshold)],1,length(aas[which(table >= keyaas_treshold)])))] = matrix(aas[which(table >= keyaas_treshold)],1,length(aas[which(table >= keyaas_treshold)])); #jeżeli są jakieś AA (=zawsze) to wpisywane są one do macierzy keyaas (odpowiendia ilość za sprawą sprawdzenia)
      keyaas_per[i,1:length(matrix((table)[which(table >= keyaas_treshold)],1,length((table)[which(table >= keyaas_treshold)])))] = matrix(round((table)[which(table >= keyaas_treshold)] / prmt$row_no,3) * 100,1,length(aas[which(table >= keyaas_treshold)])); #podobne wpisanie w odpowiedi wiersz macierzy dla ilości AA
    }
    i = 1;
    while (length(which(keyaas[,i] != "n")) > 0) {
      i = i + 1;
    }
    i = i - 1;
    
    keyaas = t(keyaas[,1:i]); #transpose matrix
    keyaas_per = t(keyaas_per[,1:i]) #transpose matrix
    return(list(AA = keyaas,per = keyaas_per))
  }
calculate_GROUP_variation <-
  function(prmt, sequence_alignment, threshold) {
    #prmt- size of alignment (output of get_parameter())
    #sequence_alignment-file[[3]]
    #threshold-threshold for detecting key amino acids (the percentage of all at the given position)
    #returns list of matrices with tabelarised symbols of the most common AA in alignment column and percentage values for contributed AA
    keyaas_treshold = prmt$row_no * (threshold / 100);
    #returns list of matrices with tabled symbols of the most common AA in alignment column and percentage values for contributed AA
    aligned_sequences_matrix = alignment2matrix(prmt,sequence_alignment);
    keyaas_gr = matrix("n",dim(aligned_sequences_matrix)[2],20 * 2);
    keyaas_per_gr = matrix("n",dim(aligned_sequences_matrix)[2],20 * 2);
    
    grupy = matrix("",dim(aligned_sequences_matrix)[1],dim(aligned_sequences_matrix)[2]);
    for (i in seq(1,(dim(aligned_sequences_matrix)[1]),by = 1)) {
      for (j in seq(1,(dim(aligned_sequences_matrix)[2]),by = 1)) {
        if (aligned_sequences_matrix[i,j] == "I" ||
            aligned_sequences_matrix[i,j] == "L" ||
            aligned_sequences_matrix[i,j] == "V")
          grupy[i,j] = "G1"
        if (aligned_sequences_matrix[i,j] == "F" ||
            aligned_sequences_matrix[i,j] == "Y" ||
            aligned_sequences_matrix[i,j] == "W")
          grupy[i,j] = "G2"
        if (aligned_sequences_matrix[i,j] == "A" ||
            aligned_sequences_matrix[i,j] == "G" ||
            aligned_sequences_matrix[i,j] == "S")
          grupy[i,j] = "G3"
        if (aligned_sequences_matrix[i,j] == "P")
          grupy[i,j] = "G4"
        if (aligned_sequences_matrix[i,j] == "N" ||
            aligned_sequences_matrix[i,j] == "Q")
          grupy[i,j] = "G5"
        if (aligned_sequences_matrix[i,j] == "R" ||
            aligned_sequences_matrix[i,j] == "H" ||
            aligned_sequences_matrix[i,j] == "K")
          grupy[i,j] = "G6"
        if (aligned_sequences_matrix[i,j] == "D" ||
            aligned_sequences_matrix[i,j] == "E" ||
            aligned_sequences_matrix[i,j] == "T")
          grupy[i,j] = "G7"
        if (aligned_sequences_matrix[i,j] == "M" ||
            aligned_sequences_matrix[i,j] == "C")
          grupy[i,j] = "G8"
        if (aligned_sequences_matrix[i,j] == "-")
          grupy[i,j] = "-"
      }
    }
    
    for (i in seq(1,dim(grupy)[2])) {
      table = (sort(table(grupy[,i]),decreasing = T)); #Czestosci i ridzaje wystepowania AA
      aas = names(table);
      keyaas_gr[i,1:length(matrix(aas[which(table >= keyaas_treshold)],1,length(aas[which(table >=
                                                                                            keyaas_treshold)])))] = matrix(aas[which(table >= keyaas_treshold)],1,length(aas[which(table >=
                                                                                                                                                                                     keyaas_treshold)])); #jeżeli są jakieś AA (=zawsze) to wpisywane są one do macierzy keyaas (odpowiendia ilość za sprawą sprawdzenia)
      keyaas_per_gr[i,1:length(matrix((table)[which(table >= keyaas_treshold)],1,length((table)[which(table >=
                                                                                                        keyaas_treshold)])))] = matrix(round((table)[which(table >= keyaas_treshold)] /
                                                                                                                                               prmt$row_no,3) * 100,1,length(aas[which(table >= keyaas_treshold)])); #podobne wpisanie w odpowiedi wiersz macierzy dla ilości AA
    }
    
    i = 1;
    while (length(which(keyaas_gr[,i] != "n")) > 0) {
      i = i + 1;
    }
    i = i - 1;
    
    keyaas_gr = t(keyaas_gr[,1:i]); #transpose matrix
    keyaas_per_gr = t(keyaas_per_gr[,1:i]) #transpose matrix
    return(list(AA = keyaas_gr,per = keyaas_per_gr))
  }
display_AA_variation <-
  function(AA_variation) {
    #merging key AAs symbols table with key AAs percentages table
    keyaas = AA_variation$AA; keyaas_per = AA_variation$per;
    size = dim(keyaas);
    output = matrix("-",size[1] * 2,size[2] + 1);
    j = 1;
    for (i in seq(1,size[1] * 2,2)) {
      output[i,-1] = keyaas[j,];
      output[i + 1,-1] = keyaas_per[j,];
      j = j + 1;
    }
    return(output)
  }
cons_best_for <- function(percentage, alignment_file) {
  a = which.max(percentage)
  namelist = alignment_file[[2]]
  out = list(c(namelist[a], a)) #output is a name of sequence and position in alignment
  return(out)
}
worst_cons_for <- function(percentage, alignment_file) {
  a = which.min(percentage) #percentage is an output of calculate_group_consensus or calculate_tru_consensus
  namelist = alignment_file[[2]]
  out = list(c(namelist[a], a)) #output is a name of sequence and position in alignment
  return(out)
}
most_common <- function(percentage, alignment_file) {
  value = which.max(table(percentage))
  vector = which(percentage %in% names(value))
  names = match(vector,file[[2]])
  list_name = alignment_file[[2]]
  out = list_name[vector]
  return(out)
}
convert <- function(amino_acids) {
  map <-
    c(
      "H","H", "H", "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )
  
  names(map) <-
    c(
      "HIP", "HID", "HIE", "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
      "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
      "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    )
  
  sapply(strsplit(amino_acids, "(?<=[A-Z]{3})", perl = TRUE),
         function(x)
           paste(map[x], collapse = ""))
}
find_seqid <- function(sequence,library) {
  ktory = c();
  for (i in seq(1,length(library))) {
    ktory[i] = length(which((library[[i]] == sequence) == TRUE));
    ktory1 = which(ktory == 1);
  }
  seqid = library[[ktory1]][1];
  return(seqid)
}
get_seq_names <- function(alignment) {
  names = file[[2]];
  return(names)
}
find_seq <- function(sequence_id, alignment_file, isoform) {
  mat = as.matrix(as.character(alignment_file[[3]]));
  nazwy_mat = get_seq_names(alignment_file)
  
  which_uniprot = c();
  for (i in seq(1,length(nazwy_mat))) {
    #finding indices of uniprot ID in seq
    which_uniprot[i] = length(
      grep(
        sequence_id, nazwy_mat[i], ignore.case = FALSE, perl = FALSE, value = FALSE,
        fixed = FALSE, useBytes = FALSE, invert = FALSE
      )
    )
  }
  if (sum(which_uniprot) < isoform) {
    isoform = 1;
  }
  seqs = mat[which(which_uniprot == 1)][isoform]; #seq of protein uniprot ID
  seqs_character = length(which(s2c(seqs[isoform]) != "-") == T);
  seq = list(sequence = seqs, len = seqs_character)
  return(seq)
}
create_structure_seq <-
  function(tunnel_file, sequence_id, alignment_file, shift) {
    #tunnel_file-> list of tunnels in protein
    #sequence_id -> uniprot id  which has been found by read_file(filename="PDBid") with PDB indentifier ;
    #alignment_file-> file wiht alignment (alignment.fst)
    
    seq = list();
    tunnel = list();
    
    #finds an uniprot name from alignent- extract uniprot seq from alignment
    base_seq = find_seq(sequence_id, alignment_file,1)
    for (i in seq(1,length(tunnel_file))) {
      tunnels_indices = as.vector(tunnel_file[[i]][[1]]);
      tunnels_names = as.vector(tunnel_file[[i]][[2]]);
      tunnels_indices = tunnels_indices[-1];
      tunnel_idx=as.numeric(tunnels_indices)
      tunnels_names = tunnels_names[-1];
      if(shift!=0){
       tunnel_idx=tunnel_idx+rep(shift,length(tunnels_indices))
      }
      seq[[i]] = rep("N",each = base_seq$len);
      
      #for (z in tunnel_idx) {
        seq[[i]][tunnel_idx] = "T"
      #}
      
      aa_positions = which(s2c(base_seq$sequence) != "-")
      just_align = alignment_file[[3]]
      paramet = alignment_parameters(just_align)
      length_alignment = dim(alignment2matrix(paramet, just_align))[2]
      tunnel[[i]] = rep("-",each = length_alignment);
      j = 1;
      for (a in aa_positions) {
        #aligning the tunnels information with the alignment sequence
        tunnel[[i]][a] = seq[[i]][j];
        j = j + 1;
      }
    }
    return(tunnel)
  }
display_structure <- function(structure,structure_file) {
  struc_length = length(structure[[1]])
  count = length(structure_file);
  output = matrix("-",count,struc_length + 1);
  v = seq(1,count,1)
  for (i in v) {
    output[i,-1] = structure[[i]]
    output[i,1] = as.character(structure_file[[i]][1,2]);
  }
  return((output))
}
isUpper <- function(s) {
  return (all(grepl("[[:upper:]]", strsplit(s, "")[[1]])))
}
barplotshow <- function(position,AA_variation) {
  row = max(which(is.na(as.numeric(
    AA_variation$per[,position]
  )) == FALSE));#checking how many different AA is on the position
  x_names = AA_variation$AA[seq(1,row),position]
  x_val = round(as.numeric(AA_variation$per[seq(1,row),position]))
  plot = barplot(
    x_val,names.arg = x_names,cex.main = 1.8,cex.lab = 1.3,cex.names =
      1.5,cex.axis = 1.5, width = rep(1.1,row),xlim = c(0,dim(AA_variation[[1]])[1]),ylim = c(0,max(x_val)+10),ylab = "Percentage",main = bquote(Amino ~
                                                                                     acid ~ variation ~ on ~ the ~ .(position) ~ position)
  );
  text(
    x = plot,y = x_val / 2,cex = 1.0,labels = paste(as.character(x_val),"%",sep = ""),xpd =
      T
  )
}
show_numbers <- function(structure) {
  nr_stru = rep("-",length(structure[[1]]));
  j = 1;
  for (i in seq(1,length(structure[[1]]),1)) {
    if (structure[[1]][i] != "-") {
      nr_stru[i + 1] = j
      j = j + 1;
    }
  }
  return(nr_stru)
}

create_final_CSV <-
  #do poprawy: powinna wczytywa? list? dodatkowych argument?w
  #dlaczego jak chce stworzy? csv z wybranymi kolumnami to ca?y czas zwraca to samo?
  function(FILENAME,variations_matrix,structure_matrix,structure_numbers,uniprot,alignment_file,list_of_scores=NULL) {
    
    sequence = s2c(find_seq(uniprot,alignment_file,1)$sequence);
    if (is.null(list_of_scores)){
      final_output = rbind(
        variations_matrix,structure_matrix,structure_numbers);
    }
    else{
      landgraf=list_of_scores[[1]]
      schneider=list_of_scores[[2]]
      TG_score= list_of_scores[[3]]
      kabat= list_of_scores[[4]]
      final_output = rbind(
        variations_matrix,structure_matrix,structure_numbers,append("landgraf metric",landgraf),append("schneider metric",schneider),append("TG metric",TG_score),append("Kabat metric",kabat),append("sequence",sequence));
    }
    files_no = ceiling(dim(final_output)[2]/1000);
    for (i in seq(1,files_no,1)){
      if (i == files_no){
        write.csv(final_output[,((i-1)*1000+1):dim(final_output)[2]],file = paste(FILENAME,"_",i,".csv",sep = ""), row.names = F)
      }
      else{
        write.csv(final_output[,((i-1)*1000+1):(i*1000)],file = paste(FILENAME,"_",i,".csv",sep = ""), row.names = F)
      }
    }
    return(final_output)
  }

TG_conservativity <- function(final_output,var_aa,method) {
  max_cons = c();
  for (i in seq(1,length(final_output[1,]),1)) {
    if (final_output[1,i] != "-") {
      #ktore nie są gapami
      max_cons[i] = as.numeric(final_output[2,i]);
    }
    else if (is.na(as.numeric(final_output[4,i])) == FALSE) {
      #czy poza gapem sa jakies inne
      max_cons[i] = as.numeric(final_output[4,i])
    }
    else
      #nie ma innych
      max_cons[i] = 0;
  }
  max_cons = max_cons[-1];#maksymalne wartosci procentowe na kazdej pozycji gdzie nie ma gapa lub poza gapem jest jeszcze inny AA |||| zniwelowanie przesuniecia o 1 kolumne w final_output wzgledem var_aa$AA
  AA = which(max_cons != 0);
  ile_var = c();
  for (i in seq(1,dim(var_aa$AA)[2],1)) {
    ile_var[i] = length(which(var_aa$AA[,i] != "n" &
                                var_aa$AA[,i] != "-"));
  }
  m_i = max_cons / ile_var;#wzgledna konserwatywnosc
  m_i[which(is.nan(m_i) == TRUE)] = 0;#zmiana NaN na 0
  m_i_m = mean(m_i); #srednia wzgledna konserwatywnosc
  m_i_s = sd(m_i);#odchylenie od sredniej wzglednej konserwatywnosci
  position = seq(1,length(m_i),1);
  plot(position,max_cons,main = "Max AA precentage on each alignment position, gaps = 0",ylab =
         "Percentage",pch = "*")
  plot(position,m_i,main = "Relative conservativity on each position",ylab =
         "Realtive conservativity",pch = "*")
  hist(
    m_i,col = "mediumpurple2",cex.main = 1.8,cex.axis = 1.4,cex.lab = 1.4,ylim = c(0,100), main =
      "Relative conservation histogram",xlab = "Relative conservation"
  )
  ret = list(
    relative_conservativity = m_i,relative_conservativity_mean = m_i_m,relative_conservativity_standard_deviation =
      m_i_s
  );
  return(ret)
}# GAPS EXCLUDED - wrong
TG_conservativity1 <- function(final_output,var_aa) {
  max_cons = c();
  for (i in seq(1,length(final_output[1,]),1)) {
    #     if (final_output[1,i] != "-") {
    #       #ktore nie są gapami
    #       max_cons[i] = as.numeric(final_output[2,i]);
    #     }
    if (is.na(as.numeric(final_output[2,i])) == FALSE) {
      #czy poza gapem sa jakies inne
      max_cons[i] = as.numeric(final_output[2,i])
    }}
  #     else
  #       #nie ma innych
  #       max_cons[i] = 0;
  #   }
  max_cons = max_cons[-1];#maksymalne wartosci procentowe na kazdej pozycji gdzie nie ma gapa lub poza gapem jest jeszcze inny AA |||| zniwelowanie przesuniecia o 1 kolumne w final_output wzgledem var_aa$AA
  AA = which(max_cons != 0);
  ile_var = c();
  for (i in seq(1,dim(var_aa$AA)[2],1)) {
    ile_var[i] = length(which(var_aa$AA[,i] != "n" & var_aa$AA[,i] != "-"));
  }
  m_i = max_cons / ile_var;#wzgledna konserwatywnosc
  m_i[which(is.nan(m_i) == TRUE)] = 0;#zmiana NaN na 0
  m_i_m = mean(m_i); #srednia wzgledna konserwatywnosc
  m_i_s = sd(m_i);#odchylenie od sredniej wzglednej konserwatywnosci
  position = seq(1,length(m_i),1);
  plot(position,max_cons,main = "Max AA precentage on each alignment position, gaps = 0",ylab =
         "Percentage",pch = "*")
  plot(position,m_i,main = "Relative conservativity on each position",ylab =
         "Realtive conservativity",pch = "*")
  hist(
    m_i,col = "mediumpurple2",cex.main = 1.8,cex.axis = 1.4,cex.lab = 1.4,ylim = c(0,100), main =
      "Relative conservation histogram",xlab = "Relative conservation"
  )
  ret = list(
    relative_conservativity = m_i,relative_conservativity_mean = m_i_m,relative_conservativity_standard_deviation =
      m_i_s
  );
  return(ret)
}# GAPS INCLUDED 
conservativity <- function(column) {
  #Methods:Shannon, Schneider, Kabat
  K = c();sum_wew = c();p = c();n = c();tab = c();sum_wew_schneider = c();suma =
    c();suma_schneider = c();
    N = length(column); #no of objects
    K = length(as.numeric(table(column))); #no of classes
    n1 = as.numeric(sort(table(column))[length(table(column))]);
    Kabat = (K/n1)*N;
    tab = as.numeric(table(column));
    
    suma_wew = 0;sum_wew_schneider = 0;
    for (i in seq(1,K,by = 1)) {
      n[i] = tab[i];
      p[i] = n[i] / N;
      res = p[i] * log2(p[i]);
      res_schneider = p[i] * log(p[i]) * (1 / log(21));
      sum_wew = suma_wew + res;
      sum_wew_schneider = sum_wew_schneider + res_schneider;
    }
    suma = -sum_wew;
    
    suma_schneider = -sum_wew_schneider;
    ret = list(Shannon = suma,Schneider = suma_schneider,Kabat = Kabat);
    #plot(suma,main="Shannon Entropy of each alignmnet position",xlab="Alignment position",ylab="Shannon entropy",pch="*");
    return(ret)
}
weights <- function(file,aligned_sequences_matrix,threshold) {
  consensus_seq = consensus(file, threshold);
  consensus_seq = consensus_seq[-which(consensus_seq == "\r")];
  length_cons_seq = length(consensus_seq);
  for (i in seq(1,dim(aligned_sequences_matrix)[1],by = 1)) {
    eq = length(which(consensus_seq == aligned_sequences_matrix[i,]));
    score[i] = eq / length_cons_seq;
  }
  return(score)
}

substitution_mtx <- function (matrix_name) {
  # function can read .txt file substitution matrix and convert it
  # into a list of alphabet (the range of letters in matrix) and
  # mtx(valuse into substitution matrix)
  matrix = read.table(matrix_name)
  mtx = matrix[-1,-1]
  alphabet = as.character(unlist(matrix[1,-1]))
  sub_mtx = list(alphabet, mtx)
  return(sub_mtx)
}

D_matrix <- function(sub_mtx) {
  values = as.numeric(as.matrix(sub_mtx[[2]]))
  dim(values) <- dim(sub_mtx[[2]])
  k = dim(sub_mtx[[2]])[1]
  distance = matrix(0,k,k)
  #CREATE DISTANCE MATRIX
  for (i in 1:k) {
    index = i;
    for (j in 1:k) {
      if (i != j) {
        distance[i,j] = (values[i,i] - values[i,j] / values[i,i])
      }
      if (i == j) {
        distance[i,j] <- 0
      }
      
    }
  }
  output = list(sub_mtx[[1]],distance)
  return(output)
}

Landgraf_conservation <-
  function(matrix_name, alignment_column_number, weights) {
    dissim_mtx = substitution_mtx(matrix_name)
    column = aligned_sequences_matrix[,alignment_column_number];
    values = as.numeric(as.matrix(dissim_mtx[[2]]))# values
    alpha = dissim_mtx[[1]]
    dim(values) <- dim(dissim_mtx[[2]])
    iterator = seq(1:length(column))
    sum_of_dist = 0;
    global_sum = 0;
    for (i in iterator) {
      iterator2 = c((i + 1):length(column))
      posa = match(column[i],alpha)
      if (!i == length(column)) {
        for (j in iterator2) {
          posb = match(column[j], alpha)
          tempi = weights[i] * values[posa,posb]
          tempj = weights[j] * values[posb,posa]
          sum_of_dist = tempi + tempj
          global_sum = global_sum + sum_of_dist
        }
        
      }
      
    }
    conservation = global_sum / length(column)
    return(conservation)
  }

sequence_stats <-
  function(alignment_file,uniprot,landgraf,schneider,TG) {
    sequence = s2c(find_seq(uniprot,alignment_file,1)$sequence);
    alignment_positions = which(sequence != "-")
    sequence = sequence[alignment_positions];
    protein_positions = seq(1,length(sequence));
    landgraf_metric = landgraf[alignment_positions];
    schneider_metric = schneider[alignment_positions];
    TG_metric = TG[alignment_positions];
    return_table = matrix("",6,length(sequence));
    return_table[1,] = sequence;
    return_table[2,] = protein_positions;
    return_table[3,] = alignment_positions;
    return_table[4,] = landgraf_metric;
    return_table[5,] = schneider_metric;
    return_table[6,] = TG_metric;
    rownames(return_table) = c(
      "sequence","protein positions","alignment positions","landraf metric entropy score","schneider metric entropy score","TG metric entropy score"
    );
    return_list = list();
    return_list[[1]] = sequence;
    return_list[[2]] = protein_positions;
    return_list[[3]] = alignment_positions;
    return_list[[4]] = landgraf_metric;
    return_list[[5]] = schneider_metric;
    return_list[[6]] = TG_metric;
    names(return_list) = c(
      "sequence","protein positions","alignment positions","landraf metric entropy score","schneider metric entropy score","TG metric entropy score"
    );
    return_set = list(return_list,return_table);
    names(return_set) = c("list","matrix");
    write.csv(x = return_table,file = "sequence_stats.csv");print(paste("Table saved in: ",getwd(),"sequence_stats.csv",sep = ""))
    return(return_set);
  }
entropy_profile<-
  function(tunnel_file, sequence_id, alignment_file, shift,prot_entropy, index) {
    #tunnel_file-> list of tunnels in protein
    #sequence_id -> uniprot id  which has been found by read_file(filename="PDBid") with PDB indentifier ;
    #alignment_file-> file wiht alignment (alignment.fst)
    #Shift-> shift of AA from caver to reference UniProt sequence
    #index of tunne;
    seq = list();
    
    index
    #finds an uniprot name from alignent- extract uniprot seq from alignment
    base_seq = find_seq(sequence_id, alignment_file,1)
    tunnels_indices = as.vector(tunnel_file[[index]][[1]]);
    tunnels_names = as.vector(tunnel_file[[index]][[2]]);
    tunnels_indices = tunnels_indices[-1];
    tunnel_idx=as.numeric(tunnels_indices)
    tunnels_names = tunnels_names[-1];
    if(shift!=0){
      tunnel_idx=tunnel_idx+rep(shift,length(tunnels_indices))
    }
    #profile=tunnel_idx
    #for(i in seq(1:length(profile))){
    # id=tunnel_idx[i]
    # profile[i]=prot_entropy[id]
    profile=prot_entropy[tunnel_idx]
    output=list(profile, tunnel_idx)
    
    return(output)
  }
