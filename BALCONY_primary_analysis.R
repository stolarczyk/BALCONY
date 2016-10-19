# Load libraries
library(seqinr)
library(Biostrings)
# Source of BALCONY functions
source("~/Uczelnia/PROJEKT/BALCONY/BALCONYfunctions.R")
# Set working directory
setwd("~/Uczelnia/PROJEKT/BALCONY")
getwd()
####################### Read data and analysis parameters
# Alignment data (fasta format)
#file=read.alignment(file="aln2_312_pro.fasta", format="fasta", forceToLower=F) #Read alignment to variable
file=read.alignment(file="aln2_312_pro_unix.fasta", format="fasta", forceToLower=F) #Read alignment to variable UNIX
# Read structure data files 
myFiles <- list.files(pattern = "*.txt");
structure_names = c();
for (i in seq(1,length(myFiles),by = 1)){
  structure_names= append(structure_names,strsplit(strsplit(myFiles[i],"_")[[1]][2],"[.]")[[1]][1]);
}
structure_names_list = as.list(structure_names)
structure_list = list();
i=1;
for (structure in myFiles){
  temp = read.table(structure)
  structure_list[[i]]= temp;
  i = i + 1;
}
# Set PDB name
pdb_name = "4JNC"; 
# Specify alignment position to examine
aa_position = 420;
# Set the threshold for consensus calculation
threshold_consensus= 30;
# Set the grouping method for calculation of consensus to the peptide sequence similarity: general, hydrophobicity, size, aromaticity
grouping_method = 'general';
# Set the threshold for detecting key amino acids (the percentage of all at the given position)
threshold_variations = 1;
# Magical shift ...
shift=170;
# Library mapping Uniprot names to PDB
lib = list(
  c("Q6U6J0","4QLA"), #creating library uniprot - PDB
  c("Q9UR30","3G0I","3G02","1QO7"),
  c("Q84HB8","4I19","4QA9"),
  c("P34913","4JNC"),
  c("P34914","1EK2","1CR6","1EK1","1CQZ"),
  c("G9BEX6","4O08","4IO0","2E3J","4INZ"),
  c("Q41413","3CXU","2CJP")
)
# Calculating consensus sequence 
consensus_seq=consensus(file, threshold_consensus);
# Extract just the sequences
aligned_sequences=file[[3]];
# Extract just the number of sequences
ilosc_seq=file[[1]];
# Get alignment parameters (#rows,#columns)
parameters=alignment_parameters(aligned_sequences);
# Convert aligned sequences to matrix
aligned_sequences_matrix=alignment2matrix(parameters,aligned_sequences);
# Calculate the identidy of consensus sequence to each sequence in the dataset
consensus_sequences_identity=cons2seqs_ident(aligned_sequences,ilosc_seq, consensus_seq)
# Calculating group consensus sequence to AA identity (instead of amino acids their group representatives are taken into consideration. Groups are established according to various AA properties - defined by the user)
group_consensus=cons2seqs_sim(parameters,aligned_sequences_matrix,consensus_seq,grouping_method);
# Following lines find the most similar and the least similar sequences to the consensus (detecting outliers, which can be excluded from the analysis)
best_consensus_true = cons_best_for(consensus_sequences_identity, file); 
worst_consensus_true = worst_cons_for(consensus_sequences_identity, file); 
best_consensus_group = cons_best_for(group_consensus, file); 
worst_consensus_group = worst_cons_for(group_consensus, file); 
list_most_common = most_common(consensus_sequences_identity, file);
#### ERROR, look into! calculate_AA_variation & calculate_GROUP_variation
# Calculating amino acids variations on each alignment (protein) position
var_aa = calculate_AA_variation(parameters,aligned_sequences,threshold_variations);
# Calculating amino acids groups variations on each alignment (protein) position
var_group = calculate_GROUP_variation(parameters,aligned_sequences,threshold_variations);
variations_matrix = display_AA_variation(var_aa);
#find reference sequence
uniprot=find_seqid(pdb_name,lib);
my_seq=find_seq(uniprot, file,1);
# add tunnels
structure=create_structure_seq(structure_list,uniprot,file,3);
structure_matrix=display_structure(structure,structure_list);
# set residue indexes
structure_numbers=show_numbers(structure);  
############# CALCULATE CONSERVATION

final_output=rbind(variations_matrix,structure_matrix,structure_numbers,uniprot);
stats=TG_conservativity1(final_output,var_aa);

how_long = parameters[[2]]
iteracja = seq(1:how_long)
scoreK = rep(1, each = how_long);
scoreS = rep(1, each = how_long);
for (i in iteracja) {
  if (i %% 100 == 0) {
    print(i)
  }
  column = aligned_sequences_matrix[,i]
  scoreK[i] = conservativity(column)$Kabat;
  scoreS[i] = conservativity(column)$Schneider
}

####FINAL CSV
SCORE_LIST=list(scoreS,rep(1, each = how_long),stats$relative_conservativity,scoreK)
final_CSV=create_final_CSV("NaszPlik.csv",variations_matrix, structure_matrix,structure_numbers,uniprot,file,SCORE_LIST)



