library(seqinr)
library(Biostrings)

#read alignment
file=read.alignment(file="aln2_312_pro.fasta", format="fasta", forceToLower=F) #Read alignment to variable

#read in tunnels
tunnel_file1=read.table("T1_mysz_1cqz.txt");#Set tunnel/structure files
tunnel_file2=read.table("T2_mysz_1cqz.txt");
tunnel_file3=read.table("T3_mysz_1cqz.txt");
tunnel_file=list(tunnel_file1,tunnel_file2,tunnel_file3)

# the structure
pdb_name="1CQZ";


threshold=30; # general sequence homology level
grouping_method = 'general';
proog=0.1;

# calculate consensus
consensus_seq=consensus(file, threshold);#function invocation
aligned_sequences=file[[3]];#aligned seqs as strings
ilosc_seq=file[[1]];
parametry=alignment_parameters(aligned_sequences);
aligned_sequences_matrix=alignment2matrix(parametry,aligned_sequences);

# sequence to consensus identity
true_consensus=cons2seqs_ident(aligned_sequences,ilosc_seq, consensus_seq)
mean_tp=mean(true_consensus);
hist(true_consensus,breaks=15, col="firebrick3",cex.main=1.8,cex.lab=1.4,main="",xlim=c(0.3,1),ylim=c(0,250),xlab="Sequence to consensus identity score"); 

#amino acid variation for a whole alignment
var_aa=calculate_AA_variation(parametry,aligned_sequences,proog);  #specific AA
var_group=calculate_GROUP_variation(parametry,aligned_sequences,proog); #group
variations_matrix=display_AA_variation(var_aa);

#Visualisation: PREPARE TUNNELS
lib=list(c("Q6U6J0","4QLA"), #creating library uniprot - PDB
         c("Q9UR30","3G0I","3G02","1QO7"),
         c("Q84HB8","4I19","4QA9"),
         c("P34913","4JNC"),
         c("P34914","1EK2","1CR6","1EK1","1CQZ"),
         c("G9BEX6","4O08","4IO0","2E3J","4INZ"),
         c("Q41413","3CXU","2CJP"))

uniprot=find_seqid(pdb_name,lib); #find reference sequence
my_seq=find_seq(uniprot, file,1);
structure=create_structure_seq(tunnel_file,uniprot,file,3); # add tunnels

structure_matrix=display_structure(structure,tunnel_file);
structure_numbers=show_numbers(structure);  # set residue indexes

final_csv = create_final_csv(
  variations_matrix,structure_matrix,structure_numbers,landgraf_normalized,score3,TG_s,uniprot,file
);
#final_csv=create_final_csv(variations_matrix,structure_matrix,structure_numbers);
final_output=rbind(variations_matrix,structure_matrix,structure_numbers,landgraf_normalized,score3,uniprot);

stats=TG_conservativity1(final_output,var_aa);
stats
# get indexes of tunnels 
tunnel_index=list()
for (i in seq(1:length(structure))){
  tunnel_index[[i]]=which(structure[[i]]=="T");
}

#caclucate onservation scores for alignment
#corr.R

shanon=score2
schneider=score3
landgraf=landgraf_normalized
max_freq=maxx

#whole protein conservation
whole_prot=which(structure[[1]]!="-");
b1=shanon[whole_prot]
b2=schneider[whole_prot]
b3=landgraf[whole_prot]
b4=TG_s[whole_prot]
b5=max_freq[whole_prot]
b=rbind(b1,b2,b3,b4,b5)
plot(b1, col ="violet",type="l", main=paste("Entropy scores for 1CQZ"),pch = 20, 
     ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
par(new=T)
plot(b4, col="magenta",type="l", pch=20, 
     ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
par(new=T)
plot(b5, col="blue",type="l",pch=20, 
     ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
legend(10,20,c("Schneider","TG","max_freq"),col=c("red","magenta","orange"))
hist(b2,main="Schneider entropy for 4JNC",xlim=c(0.0,1.0), xlab="entropy")
hist(b4,main="TG entropy for 4JNC",xlim=c(0.0,1.0), xlab="entropy")
hist(b5,main="Maximum  AA occurence for 4JNC",xlim=c(0.0,1.0), xlab="entropy")



result=list()
lengths=list()
for (i in seq(1:length(structure))){
  lengths[[i]]=length(tunnel_index[[i]])
  a1=shanon[tunnel_index[[i]]]
  a2=schneider[tunnel_index[[i]]]
  a3=landgraf[tunnel_index[[i]]]
  a4=TG_s[tunnel_index[[i]]]
  a5=max_freq[tunnel_index[[i]]]
  a=rbind(a1,a2,a3,a4,a5)
  #plot
  #plot(tunnel_index[[3]],a2, col ="red",type="l", main=paste("Entropy scores for tunnel",i),pch = 20,xlim=c(498,1375), 
  #    ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
  plot(a2, col ="red",type="l", main=paste("Entropy scores for tunnel",i),pch = 20,xlim=c(0,lengths[[i]]), 
       ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
  par(new=T)
  plot(a4, col="magenta",type="l", pch=20,xlim=c(0,lengths[[i]]), 
       ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
  par(new=T)
  plot(a5, col="orange",type="l",pch=20,xlim=c(0,lengths[[i]]), 
       ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
  legend(10,20,c("Schneider","TG","max_freq"),col=c("red","magenta","orange"))
  hist(a2,main=paste("Schnieder entropy for tunnel",i),xlim=c(0.0,1.0), xlab="entropy")
  hist(a4,main=paste("TG entropy for tunnel",i),xlim=c(0.0,1.0), xlab="entropy")
  hist(a5,main=paste("max_freq for tunnel",i),xlim=c(0.0,1.0), xlab="entropy")
  rownames(a) <- c("Shanon","Schneider","Landgraf","TG","max_freq") 
  
  result[[i]]=a
}



#profilet1=entropy_profile(tunnel_file, sequence_id, alignment_file, shift,prot_entropy, 1) 
shift=3
#schneider_entropy
tunnel_idx=profilet1[[2]]
plot(b2[tunnel_idx], type="l",xlim=c(0,80),ylim=c(0,1))
par(new=T)
profilet1=entropy_profile(tunnel_file, uniprot, file, shift,b2, 1) 
profilet2=entropy_profile(tunnel_file, uniprot, file, shift,b2, 2) 
profilet3=entropy_profile(tunnel_file, uniprot, file, shift,b2, 3)
######## ważne ;D 
#SCHNEIDER
#te same rozkłady 
k_smirnof_test_equal=list()
for( i in seq(1:length(structure))){
  k_smirnof_test_equal[[i]]=ks.test(result[[i]][2,],b2)
}
k_smirnof_test_equal
# różne rozklady
k_smirnof_test_less=list()
for( i in seq(1:length(structure))){
  k_smirnof_test_less[[i]]=ks.test(result[[i]][2,],b2, alternative = "less")
}
k_smirnof_test_less

k_smirnof_test_grater=list()
for( i in seq(1:length(structure))){
  k_smirnof_test_grater[[i]]=ks.test(result[[i]][2,],b2, alternative = "greater")
}
k_smirnof_test_grater

############TG

ks_TG_score_equal=list()
for( i in seq(1:length(structure))){
  ks_TG_score_equal[[i]]=ks.test(result[[i]][3,],b3)
}
ks_TG_score_equal

ks_TG_score_less=list()
for( i in seq(1:length(structure))){
  ks_TG_score_less[[i]]=ks.test(result[[i]][3,],b3, alternative = "less")
}
ks_TG_score_less


ks_TG_score_greater=list()
for( i in seq(1:length(structure))){
  ks_TG_score_greater[[i]]=ks.test(result[[i]][3,],b3, alternative = "greater")
}
ks_TG_score_greater


# plot values of ODF :D
dystrybuanta_TG=ecdf(result[[1]][4,])
plot(dystrybuanta_TG, xlim=c(0,1), col="blue", xlab="entropy",main="CDF of entropy for 4JNC's tunnels")
par(new=T)

dystrybuanta_PROT=ecdf(b4)
plot(dystrybuanta_PROT,xlim=c(0,1),xlab="", col="grey", main="")
par(new=T)
dystrybuanta_TG_t2=ecdf(result[[3]][4,])
plot(dystrybuanta_TG_t2,xlim=c(0,1),xlab="",col="red", main="")
par(new=T)
dystrybuanta_TG_tn=ecdf(result[[2]][4,])
plot(dystrybuanta_TG_tn,,xlim=c(0,1),xlab="",col="green", main="")
par(new=T)
legend('topleft',c("protein","T1","T2","T3"),lty=c(1,1,1,1),lwd=c(2.5,2.5),col=c("grey","blue","green","red")) 
######porownanie tunel=tunel
## korelacja tuneli 
plot(result[[1]][2,], result[[1]][3,], col="blue")
par(new=T)
plot(b2,b3, col="red")

plot(result[[2]][2,], result[[2]][3,], col="blue")
par(new=T)
plot(b2,b3, col="red")

plot(result[[3]][2,], result[[3]][3,], col="blue")
par(new=T)
plot(b2,b3, col="red")

plot(result[[1]][2,], col="blue")
par(new=T)
plot(result[[2]][2,], col="magenta")
par(new=T)
plot(result[[3]][2,], col="orange")
## konserwatywność wzdłuż tunelu 
plot(b2,b4,xlim=c(0,1),, col="grey",ylim=c(0,1),xlab="Schneider entropy",ylab="our score",main=paste("Conservation of amino acids in",pdb_name),pch=19)
par(new=T)
plot(profilet1[[1]],profilet1_TG[[1]], col="blue",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
par(new=T)
plot(profilet2[[1]],profilet2_TG[[1]], col="green",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
par(new=T)
plot(profilet3[[1]],profilet3_TG[[1]], col="red",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
legend('bottomright',c("protein","T1","T2","T3"),lty=c(1,1,1,1),lwd=c(2.5,2.5),col=c("grey","blue","green","red")) 

#konserwatywność całego białka
plot(b2, col ="black",type="l", main=paste("Entropy scores for murine soluble epoxide hydrolase"),pch = 20, xlim=c(0,560),
     ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
par(new=T)
plot(profilet1[[2]],b2[profilet1[[2]]],col="blue",xlim=c(0,560),ylim=c(0,1),xlab="",ylab="",main="",pch=21)
par(new=T)
plot(profilet2[[2]],b2[profilet2[[2]]], col="green",xlim=c(0,560),ylim=c(0,1),xlab="",ylab="",main="",pch=20)
par(new=T)
plot(profilet3[[2]],b2[profilet3[[2]]],col="red",xlim=c(0,560),ylim=c(0,1),xlab="",ylab="",main="",pch=22)
legend('topleft',c("protein","T1","T2","T3"),lty=c(1,0,0,0),pch=c(-1,21,20,22),lwd=c(2.5,2.5),col=c("black","blue","green","red")) 

