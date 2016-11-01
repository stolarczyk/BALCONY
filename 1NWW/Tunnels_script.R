library(seqinr)
library(Biostrings)

##############      UPLOAD DATA
#read alignment
file=read.alignment(file="limoEH.fasta", format="fasta", forceToLower=F) #Read alignment to variable

tunnel_file1=read.table("T1_lim.txt");#Set tunnel/structure files
tunnel_file2=read.table("T2_lim.txt");
tunnel_file3=read.table("T3_lim.txt");
tunnel_file=list(tunnel_file1,tunnel_file2,tunnel_file3)
# the structure
pdb_name="1NWW";

prot_descript="R.erytropolis limonene-1,2-EH"

# general sequence homology level
threshold=30;
grouping_method = 'general';
proog=0.1;
shift=3

#creating library uniprot - PDB
lib=list(c("Q9ZAG3","1NWW"))
##############    CALCULATE ENTROPIES
# calculate consensus
consensus_seq=consensus(file, threshold);#function invocation
aligned_sequences=file[[3]];#aligned seqs as strings
no_seq=file[[1]];
parametry=alignment_parameters(aligned_sequences);
aligned_sequences_matrix=alignment2matrix(parametry,aligned_sequences);

# optional: you can perform quick quality check
true_consensus=cons2seqs_ident(aligned_sequences,no_seq, consensus_seq)

#amino acid variation for a whole alignment
var_aa=calculate_AA_variation(parametry,aligned_sequences,proog);  #specific AA
variations_matrix=display_AA_variation(var_aa);


uniprot=find_seqid(pdb_name,lib); #find reference sequence
my_seq=find_seq(uniprot, file,1);
structure=create_structure_seq(tunnel_file,uniprot,file,shift); # add tunnels

structure_matrix=display_structure(structure,tunnel_file);
structure_numbers=show_numbers(structure);  # set residue indexes

####CALCULATE CONSERVATION SCORE
final_output=rbind(variations_matrix,structure_matrix,structure_numbers,uniprot);
stats=TG_conservativity(final_output,var_aa);
conservativity_new=conservativity(aligned_sequences_matrix);
source('calculate_landgraf.txt')

SCORE_LIST=list()

SCORE_LIST[[1]]<- conservativity_new[[1]]
SCORE_LIST[[2]]<- conservativity_new[[2]]
SCORE_LIST[[3]]<-conservativity_new[[3]]
SCORE_LIST[[4]]<- landgraf_normalized
SCORE_LIST[[5]]<-stats
names(SCORE_LIST)<-c("Shannon","Schneider", "Kabat", "Landgraf", "TG")
final_CSV=create_final_CSV("1nww",variations_matrix, structure_matrix,structure_numbers,uniprot,file,SCORE_LIST)

# get indexes of tunnels 
tunnel_index=list()
for (i in seq(1:length(structure))){
  tunnel_index[[i]]=which(structure[[i]]=="T");
}

#FIND INDEXES OF PROTEIN AA IN ALIGNMENT
whole_prot=which(structure[[1]]!="-");
prot_cons=list()
for(i in seq(1, length(SCORE_LIST))){
  prot_cons[[i]]=SCORE_LIST[[i]][whole_prot]
}
names(prot_cons)<-c("Shannon","Schneider", "Kabat", "Landgraf", "TG")

#####PLOT PROTEIN ENTROPY
colors<- rainbow(length(prot_cons))
plot(prot_cons[[1]],ylim=c(0,1), col=colors[1],xlab="amino acid",ylab="entropy score", main=paste("entropy score for ",pdb_name), type="l")
for(i in seq(2, length(prot_cons))){
  par(new=T)
  plot(prot_cons[[i]],ylim=c(0,1), col=colors[i],xlab="",ylab="", main="", type="l")
}
legend("topleft", c("Shannon","Schneider", "Kabat", "Landgraf", "TG"), col=colors, lty =c(4))


#######TUNNELS ENTROPY 
EntropyTunel=list()
lengths=list()
for (i in seq(1:length(structure))){
  lengths[[i]]=length(tunnel_index[[i]])
  Shannon_tunel=SCORE_LIST[[1]][tunnel_index[[i]]]
  Schneider_tunel=SCORE_LIST[[2]][tunnel_index[[i]]]
  Kabat_tunel=SCORE_LIST[[3]][tunnel_index[[i]]]
  Landgraf_tunel=SCORE_LIST[[4]][tunnel_index[[i]]]
  TG_tunel=SCORE_LIST[[5]][tunnel_index[[i]]]
  
  a=rbind(Shannon_tunel,Schneider_tunel,Kabat_tunel,Landgraf_tunel,TG_tunel)
  EntropyTunel[[i]]=a
}
Names=c("Shannon","Schneider", "Kabat", "Landgraf", "TG")
for(i in seq(1:length(prot_cons))){
  profilet1=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 1) 
  profilet2=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 2) 
  profilet3=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 3)
  
  plot(prot_cons[[i]], col ="black",type="l", main=paste(Names[i],"entropy score for ",prot_descript, pdb_name),pch = 20, xlim=c(0,my_seq$len),
       ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
  par(new=T)
  plot(profilet1[[2]],prot_cons[[i]][profilet1[[2]]],col="blue",xlim=c(0,my_seq$len),ylim=c(0,1),xlab="",ylab="",main="",pch=21)
  par(new=T)
  plot(profilet2[[2]],prot_cons[[i]][profilet2[[2]]], col="green",xlim=c(0,my_seq$len),ylim=c(0,1),xlab="",ylab="",main="",pch=20)
  par(new=T)
  plot(profilet3[[2]],prot_cons[[i]][profilet3[[2]]],col="red",xlim=c(0,my_seq$len),ylim=c(0,1),xlab="",ylab="",main="",pch=22)
  legend('topleft',c("protein","T1","T2","T3"),lty=c(1,0,0,0),pch=c(-1,21,20,22),lwd=c(2.5,2.5),col=c("black","blue","green","red")) 
  
}
## SCATTERPLOTS FOR PROTEINS AND TUNNELS
for(i in seq(1,length(prot_cons))){
  p1x=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 1)
  p2x=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 2)
  p3x=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 3)
  
  for( j in seq(1, length(prot_cons))){
    if (j>i){
      p1y=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[j]], 1)
      p2y=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[j]], 2)
      p3y=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[j]], 3)
      
      plot(prot_cons[[i]],prot_cons[[j]],xlim=c(0,1),ylim=c(0,1),col="grey",xlab=paste(Names[i]),ylab=paste(Names[j]),main=paste("Entropy of amino acids in",prot_descript),pch=19)
      par(new=T)
      plot(p1x[[1]],p1y[[1]], col="blue",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
      par(new=T)
      plot(p2x[[1]],p2y[[1]], col="green",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
      par(new=T)
      plot(p3x[[1]],p3y[[1]], col="red",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
      legend('bottomright',c("protein","T1","T2","T3"),pch=c(19,19,19,19),col=c("grey","blue","green","red")) 
      
    }
  }
}


#KOLMOGOROV-SMIRNOF TEST FOR QUALITY OF DISTRIBURIONS
EqualKS=list()
for(k in seq(1, length(prot_cons))){
  k_smirnof_test_equal=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_equal[[i]]=ks.test(EntropyTunel[[i]][k,],prot_cons[[k]])
  }
  EqualKS[[k]]=k_smirnof_test_equal
}
names(EqualKS)<-Names


GreaterKS=list()
for(k in seq(1, length(prot_cons))){
  k_smirnof_test_greater=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_greater[[i]]=ks.test(EntropyTunel[[i]][k,],prot_cons[[k]],alternative = "greater")
  }
  GreaterKS[[k]]=k_smirnof_test_greater
}
names(GreaterKS)<-Names

SmallerKS=list()
for(k in seq(1, length(prot_cons))){
  k_smirnof_test_smaller=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_smaller[[i]]=ks.test(EntropyTunel[[i]][k,],prot_cons[[k]],alternative = "less")
  }
  SmallerKS[[k]]=k_smirnof_test_smaller
}
names(SmallerKS)<-Names
#########DISTRIBUTION PLOT
for( i in seq(1:length(prot_cons))){
  dystrybuanta_PROT=ecdf(prot_cons[[i]])
  plot(dystrybuanta_PROT,xlim=c(0,1),xlab="", col="grey", main=paste("CDF of",Names[i], "entropy for",pdb_name," tunnels"))
  par(new=T)
  dystrybuanta_t1=ecdf(EntropyTunel[[1]][i,])
  plot(dystrybuanta_t1, xlim=c(0,1), col="blue", xlab="entropy",main="")
  par(new=T)
  dystrybuanta_t2=ecdf(EntropyTunel[[2]][i,])
  plot(dystrybuanta_t2,xlim=c(0,1),xlab="",col="red", main="")
  par(new=T)
  dystrybuanta_t3=ecdf(EntropyTunel[[3]][i,])
  plot(dystrybuanta_t3,xlim=c(0,1),xlab="",col="green", main="")
  legend('topleft',c("protein","T1","T2","T3"),lty=c(1,1,1,1),lwd=c(2.5,2.5),col=c("grey","blue","green","red")) 
  
}
