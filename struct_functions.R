#Tunels functions
#checked only for windows
get_structure_idx<- function(structure){
  #get idx of structure in alignemnt
  #return sorted list of indices in MSA, one list element is one structure
  #the first element is indices of whole protein in alignment
  stru_index=list()
  for (i in seq(1:length(structure))){
    stru_index[[i]]=which(structure[[i]]=="T");
  }
  whole_prot=which(structure[[1]]!="-");
  out=list(proteinIndices=whole_prot,structureIndices=stru_index)
  return(out)
}
get_prot_entropy<- function(whole_prot, SCORE_LIST,score_descript){
  #allows to get idx of whole protein in alignment
  #returns list of entropy for protein
 
  prot_cons=list()
  for(i in seq(1, length(SCORE_LIST))){
    prot_cons[[i]]=SCORE_LIST[[i]][whole_prot]
  }
  names(prot_cons)<-score_descript
  return(prot_cons)
}

plot_entropy<- function(prot_cons, description){
  #plots scores, sutomagically uses name of pdb FIXME
  colors<- rainbow(length(prot_cons))
  plot(prot_cons[[1]],ylim=c(0,1), col=colors[1],xlab="amino acid",ylab="entropy score", main=paste("entropy score for ",pdb_name), type="l")
  for(i in seq(2, length(prot_cons))){
    par(new=T)
    plot(prot_cons[[i]],ylim=c(0,1), col=colors[i],xlab="",ylab="", main="", type="l")
  }
  legend("topleft",description, col=colors, lty =c(4))
}