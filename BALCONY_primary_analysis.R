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
# file=read.alignment(file="aln2_312_pro.fasta", format="fasta", forceToLower=F) #Read alignment to variable
file=read.alignment(file="aln2_312_pro_unix.fasta", format="fasta", forceToLower=F) #Read alignment to variable UNIX
# Structure data files
### ZROBIC AUTOMATYCZNE WCZYTYWANIE PLIKOW ZE STRUKTURA PO NAZWIE!
myFiles <- list.files(pattern = "*.txt");
structure_names = c();
for (i in seq(1,length(myFiles),by = 1)){
  structure_names= append(structure_names,strsplit(strsplit(myFiles[i],"_")[[1]][2],"[.]")[[1]][1]);
}
tunnel_file = as.list(structure_names)
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
best_consensus_true = cons_best_for(consensus_sequences_identity, file); write.table(
  best_consensus_true,"best_consensus_true.txt",col.names = F, row.names = F
)
worst_consensus_true = worst_cons_for(consensus_sequences_identity, file); write.table(
  worst_consensus_true,"worst_consensus_true.txt",col.names = F, row.names = F
)
best_consensus_group = cons_best_for(group_consensus, file); write.table(
  best_consensus_group,"best_consensus_group.txt",col.names = F, row.names = F
)
worst_consensus_group = worst_cons_for(group_consensus, file); write.table(
  worst_consensus_group,"worst_consensus_group.txt",col.names = F, row.names = F
)
list_most_common = most_common(consensus_sequences_identity, file);
#### ERROR, look into! calculate_AA_variation & calculate_GROUP_variation
# Calculating amino acids variations on each alignment (protein) position
var_aa = calculate_AA_variation(parameters,aligned_sequences,threshold_variations);
# Calculating amino acids groups variations on each alignment (protein) position
var_group = calculate_GROUP_variation(parameters,aligned_sequences,threshold_variations);
variations_matrix = display_AA_variation(var_aa);
co

