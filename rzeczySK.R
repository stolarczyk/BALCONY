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