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
file=read.alignment(file="~/Uczelnia/PROJEKT/BALCONY/analiza_manual/aln2_312_pro.fasta", format="fasta", forceToLower=F) #Read alignment to variable
#file=read.alignment(file="aln2_312_pro_unix.fasta", format="fasta", forceToLower=F) #Read alignment to variable UNIX
file = delete_isoforms(file);

myFiles <- list.files(pattern = "*.txt");
structure_list = read_structure(myFiles)
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
  c("Q41413","3CXU","2CJP"),
  c("Q9ZAG3","1NWW")
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
structure=create_structure_seq(structure_list = structure_list,sequence_id = uniprot,alignment = file,pdb_path ="../package/BALCONY/inst/extdata/4jnc.pdb",chain_identifier = "B",shift = 236);
# Calculate TG entropy score for all alignment positions
# Write final output - amino acid variations, structure data, sequence numbers and conservation scores combined
# Need to calculate scores for all the positions to combine them with the output table!
entropy_data=list(Schneider.entropy=schneider_conservativity(file),TG.entropy = TG_conservativity(file),Kabat.entropy =  kabat_conservativity(file))
final_CSV=create_final_CSV("BALCONY_OUTPUT",var_aa, structure,uniprot,file,entropy_data)

# structure_analysis ------------------------------------------------------

# get indices of tunnels in alignment
structure=create_structure_seq(structure_list = structure_list,sequence_id = uniprot,alignment = file,pdb_path ="../package/BALCONY/inst/extdata/4jnc.pdb",chain_identifier = "B",shift = 236);
structure = exclude_low_probability_structures(structure = structure,threshold = 0.4)
tunnel_index=get_structures_idx(structure)
stru_entropy=get_structures_entropy(tunnel_index,entropy_data)
#find indices of protein in MSA
prot_cons=get_prot_entropy(tunnel_index[[1]],entropy_data)
#plot entropy
plot_entropy(prot_cons, impose=T,legend_pos = "topleft")
#get entropy for your structure

profiles_for_structure=prepare_structure_profile(structure, stru_entropy)

### PLOT STRUCTURE
plot_structure_on_protein(prot_cons, structure_profiles = profiles_for_structure,pdb_name, legend_pos = "topleft")
#######TUNNELS ENTROPY
plot(structure$structure_probabilities[4,][which(!is.nan(structure$structure_probabilities[4,]))],y = profiles_for_structure$T4_4JNC$entropy["TG.entropy",],ylim = c(0,1),xlim = c(0,1))
compare_cons_metrics(prot_cons,profiles_for_structure, pdb_name)
EQUAL=smirnof_kolmogorov_test(protein_cons = prot_cons, structure_cons = profiles_for_structure,alternative = 1, pdb_name = pdb_name,range = c(1:233),make_plot = T)
write(EQUAL, file = "SK_twosided.csv", sep = ",")
LESS=smirnof_kolmogorov_test(protein_cons = prot_cons, structure_cons = profiles_for_structure,alternative = 2, pdb_name = pdb_name,range = c(1:233),make_plot = T)
write(LESS, file = "less.csv", sep = ",")

GREATER=smirnof_kolmogorov_test(prot_cons, profiles_for_structure,3, pdb_name,F)
write(GREATER, file = "GREATER.txt", sep = " ")




