# Load libraries
library(seqinr)
library(Biostrings)
# Source of BALCONY functions
source("~/Uczelnia/PROJEKT/BALCONY/BALCONYfunctions.R")
# Set working directory
#setwd("~/Uczelnia/PROJEKT/BALCONY/analiza_manual/")
#getwd()
####################### Read data and analysis parameters
# Alignment data (fasta format)
file=read.alignment(file="aln2_312_pro.fasta", format="fasta", forceToLower=F) #Read alignment to variable
#file=read.alignment(file="aln2_312_pro_unix.fasta", format="fasta", forceToLower=F) #Read alignment to variable UNIX
file = delete_isoforms(file);

# Read structure data files 
myFiles <- list.files(pattern = "*.txt");
structure_names = c();
for (i in seq(1,length(myFiles),by = 1)){
  structure_names= append(structure_names,strsplit(myFiles,"[.]")[[i]][1]);
}
structure_names_list = as.list(structure_names)
structure_list = list();
i=1;
for (structure in myFiles){
  temp = read.table(structure)
  structure_list[[i]]= temp;
  i = i + 1;
}
names(structure_list)=structure_names
# Set PDB name
pdb_name = "1CQZ"; 
# Specify alignment position to examine
alignmnent_position = 925;
# Set the threshold for consensus calculation
threshold_consensus= 30;
# Set the grouping method for calculation of consensus to the peptide sequence similarity: general, hydrophobicity, size, aromaticity
grouping_method = 'general';
# Set the threshold for detecting key amino acids (the percentage of all at the given position)
threshold_variations = 0.01;
# Magical shift ...
shift=3;
# Substitution matrix name for Landgraf conservation
matrix_name="GONNET.txt";
# Library mapping Uniprot names to PDB
lib = list(
  c("Q6U6J0","4QLA"), 
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
parameters=alignment_parameters(alignment);
# Convert aligned sequences to matrix
aligned_sequences_matrix=alignment2matrix(alignment);
# Calculate the identidy of consensus sequence to each sequence in the dataset
consensus_sequences_identity=cons2seqs_ident(alignment,consensus_seq)
# Calculating group consensus sequence to AA identity (instead of amino acids their group representatives are taken into consideration. Groups are established according to various AA properties - defined by the user)
group_consensus=cons2seqs_sim(alignment,consensus_seq,grouping_method);
# Following line find the most similar and the least similar sequences to the consensus (detecting outliers, which can be excluded from the analysis)
outliers=noteworthy_sequences(consensus_sequences_identity, alignment);
# Calculating amino acids variations on each alignment (protein) position
var_aa = calculate_AA_variation(alignment, threshold_variations);
# Calculating amino acids groups variations on each alignment (protein) position
var_group = calculate_AA_variation(alignment,threshold_variations, grouped = T, grouping_method = "general");
variations_matrix = var_aa$matrix;
variations_matrix_group = var_group$matrix;
#find reference sequence
uniprot=find_seqid(pdb_name,lib);
my_seq=find_seq(uniprot, file);
# add structure and name the rows
structure=create_structure_seq(structure_list,uniprot,file,shift);
structure_matrix=structure$structure_matrix; 
rownames(structure_matrix) = structure_names
# set residue indexes
structure_numbers=structure$structure_numbers;
# bind the results into one table 
final_output=rbind(variations_matrix,structure_matrix,structure_numbers);
# Calculate TG entropy score for all alignment positions
TG_entropy=TG_conservativity(var_aa);
#lub
TG_entropy_group=TG_conservativity(var_group)
# Calculate Schneider, Kabat & Landgraf entropy scores for chosen alignmnet position
conservativity = conservativity(aligned_sequences_matrix)
#lub
aligned_sequences_matrix_groups = aligned_sequences_matrix2groups(aligned_sequences_matrix,'general')
conservativity_groups = conservativity(aligned_sequences_matrix_groups)
Landgraf = landgraf_conservation(matrix_name,alignment,weights = consensus_sequences_identity)
# Write final output - amino acid variations, structure data, sequence numbers and conservation scores combined
# Need to calculate scores for all the positions to combine them with the output table!
entropy_data=list(Schneider.entropy=conservativity$Schneider,Landgraf.entropy = Landgraf,TG.entropy = TG_entropy,Kabat.entropy = conservativity$Kabat)
entropy_data=list(Schneider.entropy=conservativity$Schneider,TG.entropy = TG_entropy,Kabat.entropy = conservativity$Kabat)
final_CSV=create_final_CSV("BALCONY_OUTPUT",var_aa, structure,uniprot,file,entropy_data)

