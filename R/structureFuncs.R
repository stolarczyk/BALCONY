#' Superimpose structural data of interest on sequence after the alignmment
#'
#' Create sequence of a protein structure model based on numbers of amino acids given in a text file (list of IDs and numbers in protein)
#'
#' This function is useful to create sequence covered with structural data provided in a .txt file. This sequence can be compared with alignment to check the conservation for interesting amino acid(s). Additionally, if path to the PDB file is provided the function corrects the output accordingly to the information in REMARK465 on missing amino acids.
#'
#' @param structure_list A list of structure data used for further evolutionary analysis. It can be text file(s) read by the \code{\link{read_structure}} function (text file with 2 columns: numbers of amino acids and 3-letters codes of AA; First row needs to contain markers)
#' @param sequence_id The id/name of the target sequence in alignment which will be a base of structure sequence
#' @param alignment An alignment object read with \code{\link[seqinr]{read.alignment}} function, must contain the target sequence
#' @param pdb_path A string specifying the path to the PDB file with structural information. Optional parameter, required if the structure is incomplete e.g. fragments such as loops are missing
#' @param chain_identifier A character specifying the chain of interest e.g. "A" or "B"
#' @param shift A numeric value. In case there is a need to adjust the amino acids numeration due to missing amino acids at the beginning of the structure (that are not considered in the PDB file REMARK465 section)
#'
#' @return \item{structure_matrix}{A matrix of characters "S" and "N" marking on sequence the structural element; "S" - amino acid forms the analyzed structure, "N" - amino acid which does not form the structure. Number of rows of the matrix corresponds to the number of structures analyzed}
#' \item{structure_numbers}{A vector containing the numbers of the amino acids in the sequence of interest (no gaps)}
#' \item{structure_probabilities}{A matrix of numeric values: probabilities of corresponding to the structural information from first element of the output, which helps to reduce the effect of non-consistent structural amino acids on the conservativity analysis of the structure of interest}
#' @export
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
            length_alignment = get_align_params(alignment)$col_no
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

#' Get "REMARK 465" data from PDB file
#'
#' This function extracts the data concerning missing amino acids in PDB protein structure from the PDB file
#'
#' @param pdb_path A string specifying the path tp the PDB file
#' @param chain_identifier A character specifying the chain to be considered
#'
#' @return \item{aa_numbers}{A numeric vector of indices of missing amino acids}
#' \item{chain}{A character specifying the chain which was considered in remark 465 data extraction}
#' @export
#'
#' @examples
#' require(Rpdb)
#' chain_identifier = "A"
#' pdb_path = system.file("extdata", "4jnc.pdb", package = "BALCONY")
#' print(pdb_path)
#' remark465_data = get_remarks465_pdb(pdb_path,chain_identifier)
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

#' Get IDs of structure(s) elements from aligned sequences (MSA)
#'
#' This function allows to obtain positions in aligned sequences for analyzed structure (e.g. functionally related amino acids dispersed in sequence) based on sequence corresponding to the crystal structure.
#'
#' It facilitates the management and oparation on the entropy values calculated for given MSA.
#'
#' @param structure The output of create_structure_seq function
#'
#' @return
#' \item{proteinIndices}{A sorted vector of amino acids of analyzed sequence in MSA}
#' \item{strucureIndices}{A list of sorted vectors of amino acids indices in aligned sequence for each structure}
#'
#' @export
#'
#' @examples
#' data("structure")
#' #creating library uniprot - PDB
#' lib=list(c("Q84HB8","4I19","4QA9"),
#'          c("P34913","4JNC"),
#'          c("P34914","1EK2","1CR6","1EK1","1CQZ"))
#' pdb_name = "1CQZ" #A string with path to PDB file
#' uniprot=find_seqid(pdb_name,lib)
#' tunnel=create_structure_seq(structure,uniprot,alignment)
#' structure_index=get_structures_idx(tunnel)
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

#' Get MSA-based calculated entropy for chosen protein
#'
#' This function allows to obtain vector of entropies for one complete protein sequence from MSA (gaps introduced in alignment are omitted)
#'
#' This function can be used on list of entropies or list with one element for one entropy score.
#'
#' @param protein_index Indices of given protein aminoacids in aligned sequence
#' @param score_list A list of entropy scores calculated for MSA
#'
#' @return A list where each element is a vector of entropy values provided in entropy_scores_list
#' @export
#'
#' @examples
#' data("structure")
#' data("alignment")
#' pdb_name = "1CQZ" #A string with path to PDB file
#' uniprot="P43914"
#' chain_identifier = "B"
#' structure_index=get_structures_idx(structure)
#' \donttest{entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
#'                                    Escore_entropy = Escore_conservativity(alignment))
#' prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)
#'
#' # In case of one entropy score
#' entropy_scores_list = list()
#' entropy_scores_list[[1]] = Schneider_entropy = schneider_conservativity(alignment)
#' prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)}
get_prot_entropy <- function(protein_index, score_list) {
    #documentation get_prot_entropy.Rd
    #allows to get idx of whole protein in alignment
    #returns list of entropy for protein
    if (is.list(score_list)) {
        prot_cons = list()
        for (i in seq(1, length(score_list))) {
            prot_cons[[i]] = score_list[[i]][protein_index]
        }
        names(prot_cons) <- names(score_list)
    } else{
        print("score_list is not a list!")
        stop()
    }
    return(prot_cons)
}

#' Plot entropies for protein
#'
#' This function plots entropies of protein. Plots might be superimposed or not.
#'
#' This function produces plots for given values, on X axis are amino acids, on Y axis are values of entropy/conservation. Legend contains score names for description values.
#'
#' @param protein_conservation A list or a vectors of protein conservation/entropies. The output of \code{\link{get_prot_entropy}} function
#' @param colors (optional) A vector of colors for each plot, default: rainbow
#' @param impose (optional) A boolean, if True/T plots are superimposed, if False/F plots are printed separately, default: TRUE
#' @param prot_name (optional) A string with structure name (to be used in the tile of the plot), default: NULL
#' @param legend_pos (optional) A string witch legend position - one of following: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default: "bottomleft"
#'
#' @return This function produces plots
#'
#' @keywords plot
#'
#' @export
#'
#' @importFrom grDevices dev.new rainbow
#'
#' @examples
#' data("alignment")
#' data("structure")
#' uniprot="P34914"
#' structure_index=get_structures_idx(structure)
#' \donttest{entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
#'                                    Escore_entropy = Escore_conservativity(alignment))
#' prot_entropy=get_prot_entropy(structure_index$proteinIndices, entropy_scores_list)
#'
#' plot_entropy(prot_entropy, colors = c("red","green","blue"),
#'              impose = TRUE, prot_name = "Murine Epoxide Hydrolase",
#'           legend_pos = "bottomright")}
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
        } else {
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

#' Get entropy of amino acids (for region of interest) in given protein
#'
#' This function allows to get values of entropy/conservation for amino acids dispersed in sequence of given protein. It works well with a list of dispersed amino acids in one protein.
#'
#' This function allows to obtain entropy (calculated on MSA) for dispersed amino acids in protein e.g. surface, binding site, tunnels etc. The input is a list of few structure indices in given protein sequence. Function calculates position of those in aligned sequence and returns a vector/matrix or a list of matrices with entropy values.
#'
#' @param structure_index A is a list of indices in alignment of protein and structures. Output output of \code{\link{get_structures_idx}} function
#' @param score_list A list of entropies for whole alignment
#'
#' @return A list of matrices. Rows are entropy scores, columns are
#'
#' @keywords structure
#'
#' @export
#'
#' @examples
#' data("structure")
#' data("alignment")
#'
#' #creating library uniprot - PDB
#' uniprot="P34914"
#' tunnel=create_structure_seq(structure,uniprot,alignment)
#' indices=get_structures_idx(structure)
#' protein_index = indices$proteinIndices
#' structure_index = indices$structureIndices
#' \donttest{entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
#'                                    Escore_entropy = Escore_conservativity(alignment))
#' structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)}
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

#' Combine the entropy data for structure building amino acids with their indices
#'
#' This function combines the entropy data for structure building amino acids with its indices. It prepares the data for convenient visualization or processing.
#'
#' @param structure A structure data, the output of \code{\link{create_structure_seq}} function
#' @param structure_entropy The entropy values for the structure building residues, the output of \code{\link{get_structures_entropy}} function
#'
#' @return Each element is a list of entropy values (matrix of entropy scores) and indices of residues building structure in protein of interest.
#'
#' @keywords structure
#'
#' @export
#'
#' @examples
#' data("alignment")
#' data("structure")
#' uniprot="P34914"
#' indices=get_structures_idx(structure)
#' protein_index = indices$proteinIndices
#' structure_index = indices$structureIndices
#' \donttest{entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
#'                                    Escore_entropy = Escore_conservativity(alignment))
#' structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)
#' structure_profile = prepare_structure_profile(structure, structure_entropy)}
prepare_structure_profile <-
    function(structure, structure_entropy) {
        #structure-> output of create_structure_seq
        #structure_entropy -> list of entropy scores for alignment
        megalist = list()
        score_count = dim(structure_entropy[[1]])[1] # ????
        names <- rownames(structure[[1]])
        for (i in seq(1, length(structure[[1]][, 1]))) {
            StruEnt = list(entropy = c(), idx = c())
            StruEnt$idx = as.numeric(structure[[2]][which(structure[[1]][i, ] == "S")])
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

#' Plot structure entropy on protein background
#'
#' This function enables to visually asses the stucture(s) entropy in comparison with protein's entropy
#'
#' For each entropy score from structure_profiles (these must correspond to prot_entropy) this function plots separate plots. Each plot presents entropy score for whole protein each structure is marked as one of 21 symbols available in generic \code{\link[graphics]{plot}} function.
#'
#' @param protein_entropy A list of entropy values for protein of interest. Output of \code{\link{get_prot_entropy}} function
#' @param structure_profiles Output of \code{\link{prepare_structure_profile}} function
#' @param pdb_name (optional) A string with protein's name e.g. PDB ID.
#' @param colors (optional) A vector of strings with colors to be used to plot the stucture markers of length equal to number of structures in structure profiles, default: rainbow()
#' @param structure_names (optional) A vector of strings to be displayed as names in the legend, default: "stru <no>"
#' @param legend_pos (optional) A string witch legend position - one of following: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default: "bottomleft"
#'
#' @return This function produces plot
#'
#' @keywords plot
#'
#' @export
#'
#' @examples
#' data("alignment")
#' data("structure")
#' indices=get_structures_idx(structure)
#' protein_index = indices$proteinIndices
#' structure_index = indices$structureIndices
#' \donttest{entropy_scores_list=list(Schneider_entropy = schneider_conservativity(alignment),
#'                                    Escore_entropy = Escore_conservativity(alignment))
#' structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)
#' structure_profile = prepare_structure_profile(structure, structure_entropy)
#' prot_entropy=get_prot_entropy(protein_index, entropy_scores_list)
#'
#' plot_structure_on_protein(prot_entropy, structure_profile)}
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

#' Compare conservation metrics
#'
#' This function is designed to compare the conservation metrics used in the analysis. This way the user can notice the significant correlation or differences between these to evaluate their performance in a specific case.
#'
#' This function allows to show the scatterplots of an entropy scores. The protein is marked as gray points, the structures are marked with symbols. It is useful to visualise differences between entropy scores, and choose the best one for further analysis.
#'
#' @param protein_entropy List of entropy scores values for a whole protein (output of \code{\link{get_prot_entropy}})
#' @param structure_profile Each element is a list of entropy values (matrix of entropy scores) and indices of residues building structure in protein of interest  (output of \code{\link{prepare_structure_profile}})
#' @param pdb_name The name of the analyzed protein
#'
#' @return This function produces a set of scatter plots facilitating the visual inspection of entropy metrics dependancies
#'
#' @export
#'
#' @keywords conservativity_metrics
#'
#' @examples
#' data("alignment")
#' alignment = delete_isoforms(alignment)
#' data("structure")
#' uniprot="P34913"
#' indices=get_structures_idx(structure)
#' protein_index = indices$proteinIndices
#' structure_index = indices$structureIndices
#' \donttest{entropy_scores_list=list(
#'   Schneider_entropy = schneider_conservativity(alignment),
#'   Escore_entropy = Escore_conservativity(alignment)
#' )
#' structure_entropy=get_structures_entropy(structure_index, entropy_scores_list)
#'
#' structure_profile = prepare_structure_profile(structure, structure_entropy)
#' protein_entropy=get_prot_entropy(protein_index, entropy_scores_list)
#' compare_cons_metrics(protein_entropy, structure_profile, "1CQZ")
#' }
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

#' Perform Kolmogorov-Smirnov test for structural data
#'
#' This function facilitates the comparison of conservativity of structure of interest with the rest of the protein. For example comparison of tunnel conservativity with overall protein conservativity.
#'
#' @param protein_entropy A list of calculated entropy scores (vectors of numeric values). The output of \code{\link{get_prot_entropy}} function
#' @param structure_entropy A list where each element is a list of structure indices in the protein and matrix with corresponding entropy values (each row is a separate score metric)
#' @param alternative A numeric value indicating the character of alternative hypothesis of the test to be performed: 1 - two sided test, 2 - less, 3 - greater, following the generic \code{\link[stats]{ks.test}} function.
#' @param pdb_name (optional) A string with name of the reference protein, default: "Reference"
#' @param range (optional) A numeric vector with region of reference protein to be excluded from the data set. Useful when protein consists of additional chains with outstandingly low/high entropy values which may distort result of the test, default: NULL
#' @param make_plot (optional) A logical indicating if cumulative distribution functions should be displayed, default: TRUE
#'
#' @return A matrix of p-values for each entropy metric (rows) and structure (columns)
#' @export
#' @importFrom stats ecdf
#'
#' @examples
#' data("alignment")
#' data("structure")
#' \donttest{entropy_data=list(Schneider.entropy=schneider_conservativity(alignment),
#'                             Escore.entropy = Escore_conservativity(alignment),
#'                             Kabat.entropy =  kabat_conservativity(alignment))
#' indices=get_structures_idx(structure)
#' protein_index = indices$proteinIndices
#' structure_index = indices$structureIndices
#' prot_cons=get_prot_entropy(protein_index,entropy_data)
#' stru_entropy=get_structures_entropy(structure_index,entropy_data)
#' profiles_for_structure=prepare_structure_profile(structure, stru_entropy)
#' EQUAL=kolmogorov_smirnov_test(protein_entropy = prot_cons,
#'                               structure_entropy = profiles_for_structure,
#'                               alternative = 1,range = c(1:233),make_plot = TRUE)}
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

#' Calculate real-value Evolutionary Trace (ET)
#'
#' This function allows to calculate real-valued ET for MSA.
#'
#' Here, the real-valued ET is calculated using an evolutionary tree calculated for given alignment. The tree is calculated using UPGMA method. Real-valued ET score can be used as complimentary analysis of evolutionary entropy measures.
#'
#' @param alignment Data loaded with \code{\link[seqinr]{read.alignment}} function
#'
#' @return A vector of real valued ET score corresponding to each MSA column
#'
#' @references Mihalek, Res, Lichtarge, 2004
#'
#' @keywords conservation_metrics
#' @export
#'
#' @examples
#' data("small_alignment")
#' alignment = small_alignment
#' weights = get_seq_weights(alignment)
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
    pb <- progress::progress_bar$new(
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
