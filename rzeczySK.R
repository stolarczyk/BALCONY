### sprawdzenie SK bez tuneli
EqualKS=list()
for(k in seq(1, length(bezTuneli))){
  k_smirnof_test_equal=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_equal[[i]]=ks.test(EntropyTunel[[i]][k,],bezTuneli[[k]])
  }
  EqualKS[[k]]=k_smirnof_test_equal
}
names(EqualKS)<-Names


GreaterKS=list()
for(k in seq(1, length(bezTuneli))){
  k_smirnof_test_greater=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_greater[[i]]=ks.test(EntropyTunel[[i]][k,],bezTuneli[[k]],alternative = "greater")
  }
  GreaterKS[[k]]=k_smirnof_test_greater
}
names(GreaterKS)<-Names
ListaPvalue=list()
vector1=list()
for (i in seq(1,length(bezTuneli))){
  for(j in seq(1, length(structure))){
    vector1[[j]]=EqualKS[[i]][[j]]$p.value
  }
  ListaPvalue[[i]]=vector1
}
SmallerKS=list()
for(k in seq(1, length(bezTuneli))){
  k_smirnof_test_smaller=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_smaller[[i]]=ks.test(EntropyTunel[[i]][k,],bezTuneli[[k]],alternative = "less")
  }
  SmallerKS[[k]]=k_smirnof_test_smaller
}
names(SmallerKS)<-Names

for( i in seq(1:length(bezTuneli))){
  dystrybuanta_PROT=ecdf(bezTuneli[[i]])
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

i=5
j=3
p1x=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 1)
p2x=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 2)
p3x=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 3)

p1y=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[j]], 1)
p2y=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[j]], 2)
p3y=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[j]], 3)

plot(prot_cons[[5]],abs(log(prot_cons[[3]]))/max(abs(log(prot_cons[[3]]))),xlim=c(0,1),ylim=c(0,1),col="grey",xlab=paste(Names[i]),ylab=paste(Names[j]),main=paste("Entropy of amino acids in",prot_descript),pch=19)
par(new=T)
plot(p1x[[1]],abs(log(p1y[[1]]))/max(abs(log(p1y[[1]]))), col="blue",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
par(new=T)
plot(p2x[[1]],abs(log(p2y[[1]]))/max(abs(log(p2y[[1]]))), col="green",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
par(new=T)
plot(p3x[[1]],abs(log(p3y[[1]]))/max(abs(log(p3y[[1]]))), col="red",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",pch=19)
legend('topright',c("protein","T1","T2","T3"),pch=c(19,19,19,19),col=c("grey","blue","green","red")) 

