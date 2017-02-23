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
#######TUNNELS ENTROPY 
EntropyTunel=list()
lengths=list()
for (i in seq(1:length(structure))){
  lengths[[i]]=length(tunnel_index[[i]])
  Shannon_tunel=entropy_data[[1]][tunnel_index[[i]]]
  Schneider_tunel=entropy_data[[2]][tunnel_index[[i]]]
  Kabat_tunel=entropy_data[[3]][tunnel_index[[i]]]
  Landgraf_tunel=entropy_data[[4]][tunnel_index[[i]]]
  TG_tunel=entropy_data[[5]][tunnel_index[[i]]]
  
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
