entropy_data

# get indices of tunnels in alignment
tunnel_index=get_structure_idx(structure)

#find indices of protein in MSA 
prot_cons=get_prot_entropy(tunnel_index[[1]],entropy_data)
#plot entropy
plot_entropy(prot_cons)
#get entropy for your structure
stru_entropy=get_structures_entropy(tunnel_index,entropy_data)
profils_for_structure=entropy_for_all(structure_list,uniprot,file,3,prot_cons)

### PLOT STRUCTURE
plot_structure_on_protein(prot_cons, profils_for_structure,pdb_name)
#######TUNNELS ENTROPY 

## SCATTERPLOTS FOR PROTEINS AND TUNNELS
compare_cons_metrics(prot_cons,profils_for_structure, pdb_name)


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
