#' Read structure data from a text file
#'
#' By using this function you can read text file and create an structure list which can be used in further evolutionary analysis with BALCONY package. Text file should comprise 2 or 3 columns: first one should contain indices (positions) of amino acids in the protein, the second one should contain amino acid symbols on specified positions and the third one (optionally) the numeric property of given residue.
#'
#' The files should be formatted as follows:\cr
#' 2  ASP 100\cr
#' 6   TYR 80\cr
#' 11  PHE 30\cr
#' 6   TYR 30
#'
#' @param file_names A vector of strings with structure file(s) name(s)
#'
#' @return A list with read structure data. Number of elements of this list equals to the number of files specified
#'
#' @keywords structure
#'
#' @export
#'
#' @importFrom utils read.table write.table
#'
#' @examples
#' fileConn<-file("exemplary_input1.txt")
#' writeLines(c("2 TYR 100","3 LEU 100", "7 VAL 50", "10 PHE 30", "20 SER 20"), fileConn)
#' close(fileConn)
#' fileConn<-file("exemplary_input2.txt")
#' writeLines(c("5 ALA 100","6 ILE 100", "18 GLY 100", "40 PHE 100"), fileConn)
#' close(fileConn)
#' structure_list = read_structure(file_names = c("exemplary_input1.txt", "exemplary_input2.txt"))
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

#' Delete protein isoforms from alignment object
#'
#' This function searches for isoforms in the alignment object (entries with "-digit|" in the name) and deletes them
#'
#' The isoforms are detected as entries with \code{"-digit|"} in the sequence name. If no isoforms are detected this function prints a "No isiforms detected" notification instead
#'
#' @param alignment An object class alignment read with \code{\link[seqinr]{read.alignment}} function
#'
#' @return Alignment without isoforms
#'
#' @export
#'
#' @keywords alignment
#'
#' @examples
#' data("alignment")
#' delete_isoforms(alignment)
delete_isoforms <- function(alignment) {
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

#' Determine consensus sequence
#'
#' Function calculates consensus sequence for given alignment with a threshold of user's choice.
#'
#' If maximum fraction of any amino acid on the certain position is lower than a threshold then "*" is printed instead.
#'
#' @param alignment output of of \code{\link[seqinr]{read.alignment}} function or grouped alignment created with: \code{\link{align_seq_mtx2grs}} and \code{\link{alignment2matrix}}
#' @param threshold minimal fraction of amino acids on the certain position in all sequences of the alignment to be taken for consensus letter on this position; number in range 0-100.
#'
#' @return A character vector of length of the aligned sequence containing consesus sequence based on the input alignment
#'
#' @note Please note that this function masks the seqinr package function  \code{\link[seqinr]{consensus}}
#'
#' @keywords consensus
#'
#' @export
#'
#' @examples
#' data("alignment")
#' alignment = delete_isoforms(alignment)
#' threshold=80 # Set the consensus threshold
#' consensus_sequence=consensus(alignment, threshold)
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

#' Calculate identity of each sequence in the alignment to the consensus sequence.
#'
#' The function calculates identity of consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment. Also, it can be used to weight conservativity metrics results in further steps of analysis with BALCONY package.
#'
#' Returned values are percentage of identical symbols (AA and "-") in consensus sequence and aligned sequence.
#'
#' @param alignment Data loaded with \code{\link[seqinr]{read.alignment}} function
#' @param consensus_seq Consensus sequence (output of  \code{\link{consensus}} function)
#'
#' @return Numeric vector of identity score (percentage); positions in the numeric vector correspond to sequences in alignment, respectively
#'
#' @export
#'
#' @keywords consensus
#'
#' @examples
#' data("alignment")
#' alignment = delete_isoforms(alignment)
#' threshold=60
#' consensus=consensus(alignment, threshold)
#' true_consensus=cons2seqs_ident(alignment, consensus)
cons2seqs_ident <-  function(alignment, consensus_seq) {
    true_percentage = c()
    
    for (i in seq(1, get_align_params(alignment)$row_no)) {
        true_percentage[i] = length(which((
            consensus_seq == seqinr::s2c(alignment$seq[i])
        ) == TRUE)) / length(consensus_seq)
        # true percentage calculation (number of the same AAs in consensusus and in each sequence/number of all AAs)
    }
    return(true_percentage)
}

#' Get alignment dimensions
#'
#' This function returns size of alignment, which facilitates the convenient performing upcoming steps of analysis.
#'
#' Function returns list of two elements row_no(number of rows, sequences) and col_no(number of columns,length of aligned sequences)
#'
#' @param alignment data loaded with \code{\link[seqinr]{read.alignment}}
#'
#' @return \item{row_no }{number of sequences}
#' \item{col_no }{length of aligned sequences}
#' @export
#'
#' @examples
#' data("alignment")
#' parameters=get_align_params(alignment);
get_align_params <- function(alignment) {
    #alignment data
    #return list of alignment size [row_numbers, col_numbers]
    aligned_sequences = alignment$seq
    row_num = length(aligned_sequences)
    col_num = length(seqinr::s2c(aligned_sequences[[1]]))
    param = list(row_no = row_num, col_no = col_num)
    return(param)
}

#' Convert amino acid symbols to groups according to their properties of user's choice
#'
#' This function performs a conversion of amino acid symbols to group symbols according to their properties. Implemented grouping methods are: substitution_matrix (majority of properties taken into account), polarity, size and aromaticity. "GX", where X stands for group number, are group symbols.
#'
#' @param aligned_sequences_matrix A matrix that contains aligned sequences. It is an output of \code{\link{alignment2matrix}} function
#' @param grouping_method A string which specifies the grouping method to be used. One of following: 'substitution_matrix', 'polarity', 'size', 'aromaticity'
#'
#' @return A matrix of size of the input matrix but with group symbols instead of amino acid symbols
#' @export
#' @keywords groups matrix
#'
#' @examples
#' data(alignment)
#' alignment = delete_isoforms(alignment)
#' grouping_method = "general"
#' aligned_sequences_matrix = alignment2matrix(alignment)
#' grouped = align_seq_mtx2grs(aligned_sequences_matrix,grouping_method)
align_seq_mtx2grs <-
    function(aligned_sequences_matrix,
             grouping_method="substitution_matrix") {
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

#' Group consensus to each sequence in the alignment similarity
#'
#' The function calculates similarity of group consensus to each sequence in the alignment. It facilitates an assessment of consensus accuracy and identification of outlying sequences in the alignment. Grouping amino acids allows to check similiarity between sequences by amino acids properties of user's choice.
#'
#' AA in  consensus sequences and aligned sequences are converted into groups symbols according to method of user's choice. Returned values are percentage of similar amino acids considering the properties in consensus sequence and aligned sequence.
#'
#' @param grouped_alignment The output of \code{\link[seqinr]{read.alignment}} function
#' @param grouped_consensus_seq A string of amino acids, the output of  \code{\link{consensus}} function
#'
#' @return numeric vector of identity score (percentage); positions in the numeric vector correspond to sequences in alignment, respectively
#'
#' @export
#'
#' @keywords consensus
#'
#' @examples
#' data("small_alignment")
#' alignment = delete_isoforms(small_alignment)
#' threshold_consensus = 30
#' grouping_method = "substitution_matrix"
#' alignment_grouped = align_seq_mtx2grs(alignment2matrix(alignment),grouping_method)
#' consensus_seq_grouped = consensus(alignment_grouped, threshold_consensus)
#' consensus_to_seqs_similarity = cons2seqs_sim(alignment_grouped, consensus_seq_grouped)
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

#' Convert alignment object to a matrix
#'
#' The function loads alignment into matrix to facilitate a convenient data manipulation
#'
#' @param alignment data loaded with \code{\link[seqinr]{read.alignment}} function
#'
#' @return Aligned sequences matrix where number of rows equals to number of aligned sequences and number of columns equals to the length length of aligned sequences
#' @export
#' @keywords alignment conversion matrix
#' @examples
#' data("alignment")
#' alignment = delete_isoforms(alignment)
#' matrix=alignment2matrix(alignment)
alignment2matrix <- function(alignment) {
    #alignment data
    #returns alignment as a matrix
    prmt = get_align_params(alignment = alignment)
    aligned_sequences_matrix = matrix("-", prmt$row_no, prmt$col_no)
    
    for (i in seq(1, prmt$row_no)) {
        #Putting aligned seqs into matrix
        temp = toupper(seqinr::s2c(alignment$seq[[i]]))
        for (j in seq(1, prmt$col_no)) {
            aligned_sequences_matrix[i, j] = temp[j]
        }
    }
    return(aligned_sequences_matrix)
}

#' Calculate AA variations on each position of the multiple sequence alignment
#'
#' This function calculates AA variations on each position of the alignment which may be further used for the conservativity study of the set of sequences in quiestio
#'
#' The output consists of amino acids and their fractions on each position of alignment. Amino acids with occurence frequencies lower than the threshold of user's choice are excluded.
#'
#' @param alignment The data loaded with \code{\link[seqinr]{read.alignment}} function
#' @param threshold (optional) A number in range 0-1. A of minimal frequency of occurences of amino acids at each position. Default: all the residues are visualized.
#' @param grouped (optional) A logical indicating if the grouping of amino acids should be applied. Default: FALSE
#' @param grouping_method (optional) A string which specifies the grouping method to be used. One of following: 'substitution_matrix', 'polarity', 'size', 'aromaticity'. Default: 'substitution_matrix'. Default: 'substitution_matrix' if grouped is TRUE.
#' @param weights (optional) A vector of length equal number of sequences in the alignment object with weights to overcome the taxonomic bias in the conservation analysis.
#' @param pseudo_counts (optional) A logical indicating if pseudo-counts should be added to the MSA. Pseudo-counts can be used only in non-group mode and without weights. Using these options with pseudo-counts will be suppressed. Default: FALSE
#'
#' @return Returns list of three matrices with tabelarized symbols of the most common AA in alignment column, percentage values for contributed AA and combined one.
#' \item{AA }{A matrix of AA on all alignment positions with decreasing frequencies in columns}
#' \item{per}{The percentage of AA frequencies corresponding to the $AA}
#' \item{matrix}{A combination of this two. The best suited element for visual inspection of the variability at each position}
#'
#' @export
#'
#' @examples
#' data("small_alignment")
#' alignment = delete_isoforms(small_alignment)
#' threshold=10
#' grouped = FALSE
#' var_aa=calculate_AA_variation(small_alignment, threshold, grouped)
calculate_AA_variation <-
    function(alignment,
             threshold = NULL,
             grouped = F,
             grouping_method = "substitution_matrix",
             weights = NULL,
             pseudo_counts = F) {
        freq = c()
        prmt = get_align_params(alignment = alignment)
        if (pseudo_counts) {
            # pseudocounts calculation takes a long time, initialize a progress bar
            pb <- progress::progress_bar$new(
                format = paste("Pseudocounts calculation", " [:bar] :percent eta: :eta"),
                total = prmt$col_no,
                clear = T,
                width = 80
            )
            if (grouped)
                message("Using pseudo counts, grouping set to False")
            grouped = F
            if (!is.null(weights))
                message("Using pseudo counts, weights set to NULL")
            weights = NULL
        }
        if (is.null(threshold)) {
            keyaas_treshold = 1e-07
        } else if (threshold > 0){
            keyaas_treshold = prmt$row_no * (threshold / 100)
        } else{
            stop(paste0("The threshold must be > 0. Got: ", threshold, ". Do not specify any to inspect all the residues in MSA."))
        }
        
        aligned_sequences_matrix = alignment2matrix(alignment)
        total_char = .check_alignment_characters(aligned_sequences_matrix)
        if (grouped == T)
            aligned_sequences_matrix = align_seq_mtx2grs(aligned_sequences_matrix, grouping_method)
        keyaas = keyaas_per = matrix(NA, total_char, prmt$col_no)
        for (col in seq_len(prmt$col_no)) {
            aa_counts = sort(table(aligned_sequences_matrix[, col]), decreasing = T)
            if(length(aa_counts) == 1){
                # in case there's only one character type in the column
                df_table = data.frame(col1 = names(aa_counts), col2 = aa_counts[1], row.names = 1)
            } else{
                df_table = as.data.frame(aa_counts)
            }
            if (pseudo_counts) {
                pb$tick()
                df_table = .apply_pseudocounts_position(alignment, df_table, col)
            }
            names(df_table) = c("AAs", "Frequency")
            over_threshold = which(df_table$Frequency > keyaas_treshold)
            aas = df_table$AAs[over_threshold]
            frequency = df_table$Frequency[over_threshold]
            order_decresasing = order(frequency, decreasing = T)
            frequency = frequency[order_decresasing]
            aas = as.character(aas[order_decresasing])
            col_aa_cnt = length(frequency)
            keyaas[1:col_aa_cnt, col] = matrix(aas, col_aa_cnt, 1)
            keyaas_per[1:col_aa_cnt, col] = matrix(round(frequency / sum(frequency), 6) * 100, col_aa_cnt, 1)
        }
        if (!is.null(weights)) {
            if (length(weights) != prmt$row_no)
                stop("The length of weights vector (got: ", length(weights),") must equal the number of sequences in the alignment (", prmt$row_no,")")
            weight = matrix(NA, ncol=prmt$col_no, nrow=total_char)
            weights = weights / mean(weights)
            for (col in seq_len(prmt$col_no)) {
                representatives = unique(keyaas[, col])
                representatives = representatives[!is.na(representatives)]
                for (i in seq_len(length(representatives))) {
                    which_representative = which(aligned_sequences_matrix[, col] == representatives[i])
                    weight[i, col] = mean(weights[which_representative])
                }
                a = weight[,col][which(!is.na(weight[,col]))]
                b = keyaas_per[,col][which(!is.na(keyaas_per[,col]))]
                keyaas_per[1:length(a), col] = round(as.numeric(b) * as.numeric(a), 6)
            }
        }
        return(list(
            AA = keyaas,
            percentage = keyaas_per,
            matrix = .merge_matrices(keyaas, keyaas_per)
        ))
    }

#' Find noteworthy sequences in the dataset (aligned sequences)
#'
#' This function detects noteworthy sequences (most common, closest to the consensus and most different from the consesus) to facilitate convenient detection of outlying sequences that might be excluded from the further analysis.
#'
#' @param percentage The identity of each sequence in the alignment to the consensus sequence. Output of the \code{\link{cons2seqs_ident}} function
#' @param alignment Alignment loaded with \code{\link[seqinr]{read.alignment}} function
#'
#' @return
#' \item{best_consensus}{Sequence closest to the consensus}
#' \item{worst_consensus}{Sequence most different to the consensus}
#' \item{most_common}{Most common sequence in the alignment}
#' @export
#'
#' @examples
#' data("alignment")
#' consensus_seq = consensus(alignment, 30)
#' consensus_to_sequences_identity=cons2seqs_ident(alignment,consensus_seq)
#' noteworthy_seqs(consensus_to_sequences_identity, alignment)
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

#' Amino acids symbols conversion
#'
#' This function facilitates the conversion of three letter amino acids' codes to one letter equivalents.
#'
#' In case a vector of amino acid three letter codes is provided the function returns a vector of their one letter equivalents.
#'
#' @param amino_acids A character or vector of characters with amino acid(s) three letter code(s)
#'
#' @return A chracter or vector of characters with amino acids one letter code(s)
#'
#' @export
#'
#' @keywords amino_acids
#'
#' @examples
#' three_letter_codes = c("LEU", "VAL", "ALA")
#' convert_AA_symbol(three_letter_codes)
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

#' Creates barplot with amino acid variation on the specified position
#'
#' This function facilitates a visual inspection of multiple sequence alignment (MSA) position variablity.
#'
#' @param position A number of column of alignment to be visualized
#' @param AA_variation A percentage frequency of amino acids in the alignment, calculated with \code{\link{calculate_AA_variation}} function
#'
#' @return This function produces a barchart
#' @export
#' @keywords plot
#' @importFrom graphics barplot legend par plot text
#' @examples
#' data("small_alignment")
#' position = 100
#' threshold = 0.01
#' var_aa = calculate_AA_variation(small_alignment,threshold)
#' barplotshow(position, var_aa)
barplotshow <- function (position, AA_variation) {
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

#' Get sequences weights
#'
#' This function returns weights of the sequneces in the alignment object
#'
#' The weights are calculated as shown in: \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/1097-0134\%2820010101\%2942\%3A1\%3C108\%3A\%3AAID-PROT110\%3E3.0.CO\%3B2-O}{Valdar and Thronton (2001)
#'
#' \strong{According to the following formulas:}
#'
#' \deqn{W_{j} = \frac{\sum_{k\neq j}^{N}Dist(s_{j},s_{k}))}{N-1}}{Wj = (\sum Dist(sj,sk))/N-1}\cr
#' where:\cr
#' \eqn{W_{j}}{Wj} is the weight of sequence \eqn{s_{j}}{sj}, and is defined as the average evolutionary
#' distance between \eqn{s_{j}}{sj} and all other sequences in the alignment\cr
#' \eqn{N} is the number of sequences in the alignment.
#'
#' \deqn{Dist(s_{j},s_{k})) = 1 = \frac{\sum_{i\epsilon Aligned_{jk}}Mut(s_{j},s_{k}))}{n(Aligned_{jk}))}}{Dist(sj,sk)) = 1 = {\sum Mut(sj,sk)/n(Alignedjk))}}
#' where:\cr
#' \eqn{Dist(s_{j},s_{k})}{Dist(sj,sk)}, the evolutionary distance between sequences \eqn{s_{j}}{sj} and \eqn{s_{k}}{sk}\cr
#' \eqn{Aligned_{jk}}{Alignedjk} is the set of all non-gap positions in \eqn{s_{j}}{sj} or \eqn{s_{k}}{sk}, \eqn{n(Aligned_{jk})}{n(Alignedjk)} is the number of such positions.
#'
#' \deqn{Mut(a,b) = \frac{m(a,b) - min(m)}{max(m) - min(m)}}{Mut(a,b) = (m(a,b) - min(m))/(max(m) - min(m))}
#' where:\cr
#' \eqn{Mut(a,b)} measures the similarity between amino acids \eqn{a} and \eqn{b} as derived from \eqn{a} mutation data matrix \eqn{m}
#'
#' @param alignment data loaded with \code{\link[seqinr]{read.alignment}}
#'
#' @return A vector with weights of length equal to the number of sequences in the alignment
#'
#' @references Valdar, W. S. J. & Thornton, J. M. Protein–protein interfaces: Analysis of amino acid conservation in homodimers. Proteins: Structure, Function, and Bioinformatics 42, 108–124 (2001).
#'
#' @export
#'
#' @examples
#' data("small_alignment")
#' alignment = small_alignment
#' weights = get_seq_weights(alignment)
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
    pb <- progress::progress_bar$new(
        format = paste("Distance calculation", " [:bar] :percent eta: :eta"),
        total = length(alignment$seq),
        clear = T,
        width = 80
    )
    for (i in seq_len(length(alignment$seq))) {
        pb$tick()
        for (j in seq_len(length(alignment$seq))) {
            if (i != j) {
                seqi = seqinr::s2c(alignment$seq[i])
                seqj = seqinr::s2c(alignment$seq[j])
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

#' Get position based weights of sequences in alignment
#'
#' This function calculates position based weights of sequences based on Heinkoff & Heinkoff (1994) for given MSA. The score is calculated as sum of scores for each sequence position c. Score for position c is equal 1/r if there is r different residues at column c in MSA but 1/rs if r symbol is repeated in s sequences.
#'
#'   The weights might be calculated only for amino acids symbols or for all symbols (including gaps). Also weights can be normalized by number of columns in MSA, then the sum of weights for all sequences is 1.
#'
#' @param alignment alignment loaded with \code{\link[seqinr]{read.alignment}}
#' @param gap (optional) a logical parameter, if TRUE(default) the gaps in MSA are included
#' @param normalized optional) logical parameter, if TRUE (default) weights for all sequences are divided by number of columns in alignment (when gap = TRUE weights sum up to 1)
#'
#' @return a vector of position based weights for each sequence in given alignment
#'
#' @references Henikoff, S. & Henikoff, J. G. Position-based sequence weights. Journal of Molecular Biology 243, 574–578 (1994).
#'
#' @export
#'
#' @examples
#' data("small_alignment")
#' pos_based_weights <- get_pos_based_seq_weights(small_alignment)
get_pos_based_seq_weights <- function(alignment, gap = TRUE, normalized = TRUE){
    align_param = get_align_params(alignment)
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
        } else{
            weights = apply(weights_mtx,MARGIN = 1, sum)
        }
        return(weights)
    } else{
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


#' Create CSV file to save results
#'
#' This function saves results as table into csv file. Combination of given variation allows to compare protein structure with evolutionary data content from alignment. Each position on alignment has its own column in csv file. If the length of the alignmnet exceeds 1000 characters, the output is divided into separate files with suffixes corresponing to the number of file produced by this function.
#'
#' @param filename name of the output file produced by the function
#' @param variations_matrix An object which contains alignment and frequencies of occurences each amino acids on each position of alignment. Output of \code{\link{calculate_AA_variation}}
#' @param structure A strcture object - matrix of aligned, examined protein sequence covered by structure markers (S/N). Output of \code{\link{create_structure_seq}}
#' @param sequence_id the Uniprot code of the sequence of interest
#' @param alignment the output of \code{\link[seqinr]{read.alignment}} function. A variable containing alignment data. One of the sequences must be the sequence of interest
#' @param score_list list of calculated entropy/conservation scores. Optional parameter. If not provided, this rows are not present in the output file
#'
#' @return A comma separated variable file containing information provided to this function. It is also written in the current directory
#' @export
#'
#' @name create_final_CSV
#' @examples
#' data("alignment")
#' data("structure")
#' uniprot="P34914"
#' alignment = delete_isoforms(alignment)
#' threshold = 1
#' \donttest{var_aa=calculate_AA_variation(alignment,threshold)
#' entropy_data=list(Schneider.entropy=schneider_conservativity(alignment),
#'                   Escore.entropy = Escore_conservativity(alignment))
#' create_final_CSV("my_filename",var_aa,structure,uniprot,alignment,entropy_data)}
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
        rownames(variations_matrix$matrix) = rep(c("AA name", "Percentage"), (dim(variations_matrix$matrix)[1] / 2))
        if (is.null(score_list)) {
            alignment_position = seq(1, dim(variations_matrix$matrix)[2], by = 1)
            final_output = rbind(alignment_position,
                                 variations_matrix$matrix,
                                 sequence,
                                 structure_output)
        } else {
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

#' Exclude low probability structural data
#'
#' This function facilitates the exclusion of low probability structural data from the downstream conservativity analysis, which helps to reduce the effect of non-consistent structural amino acids on the conservativity analysis of the structure of interest
#'
#' @param structure A structure object generated with \code{\link{create_structure_seq}} function
#' @param threshold The threshold for the structural data exclusion
#'
#' @return \item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
#' \item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
#' \item{structure_probabilities}{A matrix of numeric values: probabilities of corresponding to the structural information from first element of the output}
#'
#' @export
#'
#' @keywords structure
#'
#' @examples
#' data("alignment")
#' structure_files = c(system.file("extdata", "T1_4JNC.structure", package = "BALCONY"),
#'                     system.file("extdata", "T2_4JNC.structure", package = "BALCONY"),
#'                     system.file("extdata", "T3_4JNC.structure", package = "BALCONY")
#' )
#' structure_list = read_structure(structure_files)
#' #creating library uniprot - PDB
#' lib=list(c("Q84HB8","4I19","4QA9"),
#'          c("P34913","4JNC"),
#'          c("P34914","1EK2","1CR6","1EK1","1CQZ"))
#' pdb_name = "4JNC"
#' uniprot=find_seqid(pdb_name,lib)
#' tunnel=create_structure_seq(structure_list,uniprot,alignment)
#' tunnel_excluded = excl_low_prob_strcts(tunnel, 0.5)
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

#' Calculate the Escore conservation metric
#'
#' This function facilitates the calculation of Escore conservation metric (in amino acid or group mode)
#'
#' The conservativity score is calculated according to the following formula:
#' \deqn{P(i) = max(p(i))/n(i)}
#' \deqn{Pnorm(i) = P(i)/max(P)}
#' \deqn{score = -ln(P_norm(i))/max(-ln(P_norm))}
#' \cr where:\cr
#' \eqn{p(i)} - amino acids frequency on i-th position where gaps are included \cr
#' \eqn{n(i)} - amino acids count on i-th position where gaps are excluded
#'
#' @param alignment data read with \code{\link[seqinr]{read.alignment}} function
#' @param grouping_method (optional) A string which specifies the grouping method to be used. One of following: 'substitution_matrix', 'polarity', 'size', 'aromaticity', default: NULL
#' @param weights (optional) A vector of length equal number of sequences in the alignment object with weights to overcome the taxonomic bias in the conservation analysis.
#' @param pseudo_counts (optional) A logical indicating if pseudo-counts should be added to the MSA. Pseudo-counts can be used only in non-group mode and without weights. Using these options with pseudo-counts will be suppressed. Default: FALSE
#'
#' @return A vector of length equal to the length of aligned sequences
#'
#' @export
#'
#' @keywords conservation_metrics
#'
#' @note Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula: \deqn{conservation = 1 - entropy}
#'
#' @examples
#' data("small_alignment")
#' \donttest{conservation_score = Escore_conservativity(alignment)}
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

#' Calculate Kabat conservation metric
#'
#' This function facilitates the calculation of Kabat conservation metric.
#'
#' @param alignment Alignment data read with \code{\link[seqinr]{read.alignment}} function
#' @param weights (optional) A vector of length equal number of sequences in the alignment object with weights to overcome the taxonomic bias in the conservation analysis.
#' @param pseudo_counts (optional) A logical indicating if pseudo-counts should be added to the MSA. Pseudo-counts can be used only without weights. Using this option with pseudo-counts will be suppressed. Default: FALSE
#'
#'
#' @return A vector of length equal to the length of aligned sequences
#'
#' @note Please note that the Kabat matric formula can be found in the paper listed in "See Also" section below.
#' Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
#'   \deqn{conservation = 1 - entropy}
#'
#' @references http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract
#'
#' @export
#'
#' @keywords conservation_metrics
#'
#' @examples
#' data("small_alignment")
#' \donttest{conservation_score = kabat_conservativity(alignment)}
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
        var_aa = calculate_AA_variation(alignment, weights = weights, pseudo_counts = pseudo_counts)
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

#' Calculate Schneider conservation metric
#'
#' This function facilitates the calculation of Schneider conservation metric.
#'
#' @param alignment Alignment data read with \code{\link[seqinr]{read.alignment}} function
#' @param weights (optional) A vector of length equal number of sequences in the alignment object with weights to overcome the taxonomic bias in the conservation analysis.
#' @param pseudo_counts (optional) A logical indicating if pseudo-counts should be added to the MSA. Pseudo-counts can be used only without weights. Using this option with pseudo-counts will be suppressed. Default: FALSE
#'
#' @return A vector of length equal to the length of aligned sequences
#'
#' @note Please note that the Schneider matric formula can be found in the paper listed in "See Also" section below.
#' Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
#' \deqn{conservation = 1 - entropy}
#'
#' @references http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract
#'
#' @keywords conservation_metrics
#' @export
#'
#' @examples
#' data("small_alignment")
#' \donttest{conservation_score = schneider_conservativity(alignment)}
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
        var_aa = calculate_AA_variation(alignment, weights = weights, pseudo_counts = pseudo_counts)
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

#' Calculate Shannon conservation metric
#'
#' This function facilitates the calculation of Shannon conservation metric.
#'
#' @param alignment Alignment data read with \code{\link[seqinr]{read.alignment}} function
#' @param weights (optional) A vector of length equal number of sequences in the alignment object with weights to overcome the taxonomic bias in the conservation analysis.
#' @param pseudo_counts (optional) A logical indicating if pseudo-counts should be added to the MSA. Pseudo-counts can be used only without weights. Using this option with pseudo-counts will be suppressed. Default: FALSE
#'
#' @return A vector of length equal to the length of aligned sequences
#'
#' @note Please note that the Schneider matric formula can be found in the paper listed in "See Also" section below.
#' Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
#' \deqn{conservation = 1 - entropy}
#'
#' @references http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract
#'
#' @keywords conservation_metrics
#' @export
#'
#' @examples
#' data("small_alignment")
#' \donttest{conservation_score = shannon_conservativity(alignment)}
shannon_conservativity <-
    function(alignment,
             weights = NULL,
             pseudo_counts = F) {
        if (!is.matrix(alignment)) {
            aligned_sequences_matrix = alignment2matrix(alignment = alignment)
        } else{
            aligned_sequences_matrix = alignment
        }
        sum = rep(NaN, dim(aligned_sequences_matrix)[2])
        var_aa = calculate_AA_variation(
            alignment,
            threshold = 0.01,
            weights = weights,
            pseudo_counts = pseudo_counts
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

#' Calculate pairwise alignment for whole MSA
#'
#' For given alignment calculate pariwise alignments and returns alignment score.
#'
#' @param alignment An alignment object read with \code{\link[seqinr]{read.alignment}} function
#' 
#' @return Matrix of alignment scores
#' 
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom progress progress_bar
#'
#' @export
#'
#' @examples
#' data("small_alignment")
#' \donttest{pairwiseAlignemnt_scores=pairwise_alignment_MSA(small_alignment)}
pairwise_alignment_MSA <- function(alignment) {
    pb <- progress_bar$new(
        format = paste("MSA", " [:bar] :percent eta: :eta"),
        total = alignment$nb,
        clear = T,
        width = 80
    )
    score_mtx = matrix(NA, nrow = alignment$nb, ncol = alignment$nb)
    for (i in seq_len(alignment$nb)) {
        pb$tick()
        seqi = alignment$seq[i]
        score_mtx[i,] = sapply(alignment$seq, .convertAlign, seqi)
    }
    return(score_mtx)
}

#' Calculate cumulative relative entropy score
#'
#' This function calculates cumulative relative entropy score according to: \href{https://www.sciencedirect.com/science/article/pii/S0022283600940361?via\%3Dihub}{Hannenhalli and Russell (2000)}.
#'
#' \strong{PSEUDO-ALGORITHM} (According to \href{https://www.sciencedirect.com/science/article/pii/S0022283600940361?via\%3Dihub}{Hannenhalli and Russell (2000)}):
#'   \enumerate{
#'     \item (If score matrix is not provided) Run pairwise alignments for all available sequences in the input MSA and save scores to a matrix
#'     \item (If score matrix is not provided) Calculate a distance matrix based off of the alignment scores one
#'     \item Perform hierarchical clustering on the distance matrix (UPGMA method)
#'     \item Get the sequence clusters
#'     \item Divide the alignment into \code{sub_groups} which are the clusters
#'     \item Run hmmbuild for \code{whole_alignment} without \code{sub-group} and \code{sub_group}
#'     \item Calculate relative entropy using these two as indicated in the Reference and repeat for each \code{sub_group}
#'     \item Calculate the cumulative relative entropy
#'   }
#' \strong{hmmbuild program}:
#' This function uses hmmbuild program of \href{http://www.hmmer.org/}{HMMER} suite for HMM profile generation for MSA.
#' We recommend downloading and installing HMMER by following the instructions and steps in the \href{http://hmmer.org/download.html}{ HMMER installation website }.
#'
#' @references Hannenhalli, S. S. & Russell, R. B. Analysis and prediction of functional sub-types from protein sequence alignments11Edited by J. Thornton. Journal of Molecular Biology 303, 61–76 (2000).
#' @param alignment An alignment object read with \code{\link[seqinr]{read.alignment}} function
#' @param hmmbuild_path (optional if running under UNIX) The aboslute path to the hmmbuild binary
#' @param pairwiseAlignemnt_scores (optional) A matrix with pairwise alignment scores. For example created by \code{\link[Biostrings]{pairwiseAlignment}}. If the matrix is not provideded by the user it is calculated automatically by the function (time consuming). The sequences are extracted from the alignemnt object.
#'
#' @return A vector of length equal to the length of aligned sequences
#'
#' @export
#'
#' @keywords conservation_metrics
#'
#' @examples
#' #No example due to external software requirements
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
                leftover_prob = .preprocess_hmm_output(leftover_hmm)$probabilities
                leftover_pos = .preprocess_hmm_output(leftover_hmm)$alignment_positions
                sub_prob = .preprocess_hmm_output(sub_hmm)$probabilities
                sub_pos = .preprocess_hmm_output(sub_hmm)$alignment_positions
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


#' Read a substitution matrix
#'
#' This function facilitates reading of substitution matrices for further use
#'
#' @param matrix_name A string with path to the substitution matrix in a text file to be read
#'
#' @return
#' \item{names}{A vector of characters with amino acid names included in the matrix}
#' \item{matrix}{A numeric matrix with values}
#'
#' @export
#'
#' @examples
#' path = system.file("extdata", "GONNET.txt", package = "BALCONY")
#' sub_mat = substitution_mtx(path)
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

#' Calculate substitution rate matrix between two amino acids
#'
#' This function is used to calculate Landgraf conservation metric. D_matrix contains substitution rates between two amino acids in the alignment, according to the following formula:\cr \cr
#' \deqn{D(a,b)= (d(a,a)-d(a,b))/d(a,a)}\cr
#' where:\cr
#' \eqn{d(a,a)} is a probability of AA substitution by itself\cr
#' \eqn{d(a,b)} is a probability of substitution of amino acid a with other amino acid.
#'
#' @param substitution_matrix A matrix with probablity of substitutions, e.g.  Gonnet substitution matrix
#'
#' @return A matrix of substitution probablities for all amino acids
#' @export
#'
#' @examples
#' data("gonnet")
#' distance=D_matrix(gonnet)
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

#' Calculate Landgraf conservation score
#'
#' This function calculates Landgraf conservarion score
#'
#' @param matrix_name A string with path to the file woith substitution matrix to be used to calculate the Landgraf conservation score. Optional parameter, if not provided the Gonnet substitution matrix is used (according to author's suggestion)
#' @param alignment An alignment object read with \code{\link[seqinr]{read.alignment}} function
#' @param weights A vector with weight for each sequence in the alignment  to be used to calculate the Landgraf conservation score e.g. each sequence similarity to the consensus sequence from the alignment - output from \code{\link{cons2seqs_ident}} fuction
#'
#' @return A vector of length equal to the length of aligned sequences
#'
#' @note Please note that the Shannon matric formula can be found in the paper listed in "See Also" section below.
#' Also, this function originally calculates the entropy values which can be used to estimate the conservativity score according to the following formula:
#' \deqn{conservation = 1 - entropy}
#'
#' @references http://onlinelibrary.wiley.com/doi/10.1002/prot.10146/abstract
#'
#' @export
#'
#' @examples
#' data("small_alignment")
#' alignment = small_alignment
#' threshold_consensus = 30
#' consensus_seq=consensus(alignment, threshold_consensus);
#' consensus_sequences_identity=cons2seqs_ident(alignment, consensus_seq)
#' score = landgraf_conservativity(alignment = alignment, weights = consensus_sequences_identity)
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
        pb <- progress::progress_bar$new(
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

#' Calculate pseudo counts for alignment
#'
#' This function calculates pseudo-counts (as shown in \href{https://doi.org/10.1093/bioinformatics/12.2.135}{Henikoff et al. (1996)}) for an alignment with the use of substitution matrices. It is recommended to estimate amino acid frequencies for alignments with small number of sequences (in order to calculate reliable entropy scores)
#'
#' @param alignment data loaded with \code{\link[seqinr]{read.alignment}}
#' @param substitution_mtx Matrix with amino acids substitution frequencies. Default: GONNET
#'
#' @return Matrix with pseudo counts of size 21x number of alignment columns
#'
#' @export
#'
#' @references Henikoff et al.(1996) Using substitution probabilities to improve position-specific scoring matrices, Bioinformatics, 12, 135–143\cr
#' Claverie (1994) Some useful statistical properties of position-weight matrices.
#' Comput. Chem., 18, 287-293
#'
#' @note Please note that when using other scoring matrix user needs to make sure that all alignment symbols are present there. Missing symbol will issue an error.
#'
#' @examples
#' data("alignment")
#' PC <- calculate_pseudo_counts(alignment)
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
        # note that other substitution matrices may have different symbols (like '*', or other letters)
        mtx_alignment = alignment2matrix(alignment)
        pseudoCounts = matrix(NA,
                              nrow = AA_COUNT,
                              ncol = dim(mtx_alignment)[2])
        B = apply(
            X = mtx_alignment,
            MARGIN = 2,
            FUN = function(X)
                (5 * length(unique(X)))
        )
        calc_ba <- function(column, B, substitution_mtx) {
            pseudocounts = matrix(NA, nrow = AA_COUNT, ncol = 1)
            names(pseudocounts) <- append(Biostrings::AA_STANDARD, "-")
            N = sum(table(column))
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

