# Load libraries
require(seqinr)
require(Biostrings)
# Source of BALCONY functions
source("~/Uczelnia/PROJEKT/BALCONY/BALCONYfunctions.R")
# Set working directory
setwd("~/Uczelnia/PROJEKT/BALCONY/analiza_manual/")
getwd()
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
#pdb_name = "4JNC"; 
pdb_name = "1CQZ"; 
# Specify alignment position to examine
alignmnent_position = 925;
# Set the threshold for consensus calculation
threshold_consensus= 30;
# Set the grouping method for calculation of consensus to the peptide sequence similarity: general, hydrophobicity, size, aromaticity
grouping_method = 'general';
# Set the threshold for detecting key amino acids (the percentage of all at the given position)
threshold_variations = 0.01;
# Substitution matrix name for Landgraf conservation
matrix_name="GONNET.txt";
# Library mapping Uniprot names to PDB
dictionary = list(
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
# Convert aligned sequences to matrix
aligned_sequences_matrix=alignment2matrix(file);
# Calculate the identidy of consensus sequence to each sequence in the dataset
consensus_sequences_identity=cons2seqs_ident(file, consensus_seq)
# Calculating group consensus sequence to AA identity (instead of amino acids their group representatives are taken into consideration. Groups are established according to various AA properties - defined by the user)
group_consensus=cons2seqs_sim(file,consensus_seq,grouping_method);
# Following line find the most similar and the least similar sequences to the consensus (detecting outliers, which can be excluded from the analysis)
noteworthy_seqs= noteworthy_sequences(consensus_sequences_identity, file);
# Calculating amino acids variations on each alignment (protein) position
var_aa = calculate_AA_variation(file,threshold_variations);
# Calculating amino acids groups variations on each alignment (protein) position
var_aa_group = calculate_AA_variation(file,threshold_variations, grouped = T,grouping_method = "general");
variations_matrix = var_aa$matrix;
variations_matrix_grouped=var_aa_group$matrix
#find reference sequence
uniprot=find_seqid(pdb_name,dictionary);
my_seq=find_seq(uniprot, file);
# add structure and name the rows
structure=create_structure_seq(structure_list,uniprot,file,pdb_path = "1cqz.pdb",chain_identifier = "B");
# Calculate TG entropy score for all alignment positions
TG_entropy=TG_conservativity(file);
# Calculate Schneider, Kabat & Landgraf entropy scores for chosen alignmnet position
conservativity = kabat_conservativity(file)
Landgraf = landgraf_conservativity(alignment_file = file,weights = consensus_sequences_identity)
# Write final output - amino acid variations, structure data, sequence numbers and conservation scores combined
# Need to calculate scores for all the positions to combine them with the output table!
entropy_data=list(Schneider.entropy=schneider_conservativity(file),TG.entropy = TG_conservativity(file),Kabat.entropy =  kabat_conservativity(file))
final_CSV=create_final_CSV("BALCONY_OUTPUT",var_aa, structure,uniprot,file,entropy_data)

barplotshow(position = 34, AA_variation = var_aa)
