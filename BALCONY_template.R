library(seqinr)
library(Biostrings)
options(warn = -1)
##############################################################################################################
###############################################################################################################
###############################################################################################################
######ADD YOUR WORKING DIRECTORY AND SOURCE FILE LOCATION AND COMMENT OTHERS. COMMENT ALL AFTER YOU ARE DONE!######
###############################################################################################################
###############################################################################################################
###############################################################################################################

#setwd("~/Dokumenty/GIT/BALCONY/")
setwd("~/Uczelnia/PROJEKT/BALCONY/BALCONY/")
getwd()
source("~/Uczelnia/PROJEKT/BALCONY/BALCONY/BALCONYfunctions.R")

#source("~/Dokumenty/GIT/BALCONY/BALCONYfunctions.R")

#source("~/Dokumenty/GIT/BALCONY/BALCONYfunctions.R")

############################################################################### Wczytaj pliki z tunelami do wyświetlenia
#reading file(s) containing information on tunnels and alignment file

file=read.alignment(file="aln2_312_pro.fasta", format="fasta", forceToLower=F) #Read alignment to variable
file=read.alignment(file="aln2_312_pro_unix.fasta", format="fasta", forceToLower=F) #Read alignment to variable UNIX
#tunnel_file1=read.table("t1_4jnc.txt");#tunnel_file2=read.table("t2_4jnc.txt");#tunnel_file3=read.table("t3_4jnc.txt");
tunnel_file1 = read.table("T1_mysz_1cqz.txt");#Set tunnel/structure files
tunnel_file2 = read.table("T2_mysz_1cqz.txt");
tunnel_file3 = read.table("T3_mysz_1cqz.txt");
tunnel_file = list(tunnel_file1,tunnel_file2,tunnel_file3)
pdb_name = "4JNC"; #set pdb name here
aa_position = 420; #SET Which peptide position AA variations you want to explore!!
threshold = 30; #SET THE THRESHOLD HERE!
grouping_method = 'general';
proog = 0.1; #threshold for detecting key amino acids (the percentage of all at the given position)

####################################################### calculate consensus
consensus_seq=consensus(file, threshold);#function invocation
aligned_sequences=file[[3]];#aligned seqs as strings
ilosc_seq=file[[1]];
parametry=alignment_parameters(aligned_sequences);
aligned_sequences_matrix=alignment2matrix(parametry,aligned_sequences);
####################################################### policzenie dopasowania consensusu
true_consensus=cons2seqs_ident(aligned_sequences,ilosc_seq, consensus_seq)
mean_tp=mean(true_consensus);
hist(true_consensus,breaks=15, col="firebrick3",cex.main=1.8,cex.lab=1.4,main="",xlim=c(0.3,1),ylim=c(0,250),xlab="Sequence to consensus identity score"); 
######################################################## liczy grup. konsensus
konsensus_grupowy=cons2seqs_sim(parametry,aligned_sequences_matrix,consensus_seq,grouping_method);
mean_tpG=mean(konsensus_grupowy);
hist(konsensus_grupowy,breaks=15, col="firebrick3",cex.main=1.8,cex.lab=1.4,main="",ylim=c(0,250),xlim=c(0.3,1),xlab="Sequence to consensus identity score");

################################################################################# analiza outputu
best_consensus_true = cons_best_for(true_consensus, file); write.table(
  best_consensus_true,"best_consensus_true.txt",col.names = F, row.names = F
)
worst_consensus_true = worst_cons_for(true_consensus, file); write.table(
  worst_consensus_true,"worst_consensus_true.txt",col.names = F, row.names = F
)
best_consensus_group = cons_best_for(konsensus_grupowy, file); write.table(
  best_consensus_group,"best_consensus_group.txt",col.names = F, row.names = F
)
worst_consensus_group = worst_cons_for(konsensus_grupowy, file); write.table(
  worst_consensus_group,"worst_consensus_group.txt",col.names = F, row.names = F
)
list_most_common = most_common(true_consensus, file);
write.table(list_most_common,"most_common.txt",col.names = F, row.names = F);
################################################################################# AA variations
var_aa = calculate_AA_variation(parametry,aligned_sequences,proog);
var_group = calculate_GROUP_variation(parametry,aligned_sequences,proog);
variations_matrix = display_AA_variation(var_aa);
################################################################################# plik z tunelami
lib = list(
  c("Q6U6J0","4QLA"), #creating library uniprot - PDB
  c("Q9UR30","3G0I","3G02","1QO7"),
  c("Q84HB8","4I19","4QA9"),
  c("P34913","4JNC"),
  c("P34914","1EK2","1CR6","1EK1","1CQZ"),
  c("G9BEX6","4O08","4IO0","2E3J","4INZ"),
  c("Q41413","3CXU","2CJP")
)

uniprot = find_seqid(pdb_name,lib);
my_seq = find_seq(uniprot, file,1);
shift=170#for 4JNC structure
structure = create_structure_seq(tunnel_file,uniprot,file,shift);
structure_matrix = display_structure(structure,tunnel_file);
structure_numbers = show_numbers(structure);

final_csv = create_final_csv(
  variations_matrix,structure_matrix,structure_numbers,landgraf_normalized,score3,TG_s,uniprot,file
);
barplotshow(aa_position,var_aa);
final_output=rbind(variations_matrix,structure_matrix,structure_numbers);

stats=TG_conservativity1(final_output,var_aa);
conservativity_new=conservativity(aligned_sequences_matrix);



#do testów
how_long = parametry[[2]]
iteracja = seq(1:how_long)
score = rep(1, each = how_long);
for (i in iteracja) {
  if (i %% 100 == 0) {
    print(i)
  }
  column = aligned_sequences_matrix[,i]
  score[i] = conservativity(column)$Kabat;
}
plot(score)
