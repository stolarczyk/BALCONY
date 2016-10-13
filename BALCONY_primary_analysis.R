# Load libraries
library(seqinr)
library(Biostrings)
# Source of BALCONY functions
source("~/Uczelnia/PROJEKT/BALCONY/BALCONY/BALCONYfunctions.R")
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
threshold = 30;
# Set the grouping method for calculation of consensus to the peptide sequence similarity: general, hydrophobicity, size, aromaticity
grouping_method = 'general';
# Set the threshold for detecting key amino acids (the percentage of all at the given position)
proog = 0.1;
# Magical shift ...
shift=170;
# Calculating consensus sequence 
consensus_seq=consensus(file, threshold);
# Extract just the sequences
aligned_sequences=file[[3]];
# Extract just the number of sequences
ilosc_seq=file[[1]];
# Get alignment parameters (#rows,#columns)
parametry=alignment_parameters(aligned_sequences);
# Convert aligned sequences to matrix
aligned_sequences_matrix=alignment2matrix(parametry,aligned_sequences);

