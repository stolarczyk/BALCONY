#Tunels functions
#checked only for windows
get_structure_idx<- function(structure){
  #documentation get_structure_idx.Rd
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
  #documentation get_prot_entropy.Rd
  #allows to get idx of whole protein in alignment
  #returns list of entropy for protein
 
  prot_cons=list()
  for(i in seq(1, length(SCORE_LIST))){
    prot_cons[[i]]=SCORE_LIST[[i]][whole_prot]
  }
  names(prot_cons)<-score_descript
  return(prot_cons)
}

plot_entropy<- function(prot_cons, description, colors,impose){
  #plots scores on one plot, if colors are not specified plot as rainbow
  #automagically uses name of pdb FIXME
  # recive list of entropy scores and list of Names in this list
  if(missing(colors)){
    colors<- rainbow(length(prot_cons))
  }
  if(missing(impose)){
    stacked<-T
  }
 
  plot(prot_cons[[1]],ylim=c(0,1), col=colors[1],xlab="amino acid",ylab="entropy score", main=paste("entropy score for ",pdb_name), type="l")
  for(i in seq(2, length(prot_cons))){
    par(new=impose)
    plot(prot_cons[[i]],ylim=c(0,1), col=colors[i],xlab="",ylab="", main="", type="l")
  }
  legend("topleft",description, col=colors, lty =c(4))
}

### entropy for stuff
get_structures_entropy<- function(structure_index, SCORE_LIST, NAMES){
  #structure_index is a list of indexes in alignment of protein and structures 
  #SCORE_LIST list of entropies for whole alignment
  #NAMES names of entropy scores
  #output is a list of matrixes where each row contains values of entropy for AA in structure
  t_index=structure_index$structureIndices
  Entropy=list()
  lengths=list()
  for (i in seq(1:length(t_index))){
   lengths[[i]]=length(t_index[[i]])
   output=matrix(NA,nrow=length(SCORE_LIST),ncol=lengths[[i]])
    for (j in seq(1:length(SCORE_LIST))){
      output[j,]=SCORE_LIST[[j]][t_index[[i]]]
    }
    rownames(output)<-NAMES
    Entropy[[i]]=output
   
  }
  return(Entropy)
}

profil_for_all<- function(tunnel_file, uniprot, file, shift,prot_cons){
  profilet=list()
  stru_numb=length(tunnel_file)
  names<- paste(names(prot_cons))
  #k- possible entropy score iterator
  megalist=list()
  for (k in seq(1,length(prot_cons))){
    
    for(i in seq(1,stru_numb)){
      profilet[[i]]=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[k]], i)
    #profilet[[i]]<- setNames(mapply(paste("stru",i))
    }
    megalist[[names[k]]]<-profilet
    
  }
  return(megalist)
}

profils_for_structure=profil_for_all(tunnel_file,uniprot,file,3,prot_cons)
length(profils_for_structure$s)
##### FINISH ME! :D 
plot_structure_entropy(profils_for_structure, prot_cons,prot_descript, pdb_name)
plot_structure_entropy<-function(global_profil, prot_cons, prot_descript, pdb_name, strucure_names, colors ){
  #dodać wejście funkcji - listę entropii
  #global_profil-lista indeksów AAi ich entropii dla wszystkich struktur
  #wyświetli wykresy max dla 25 struktur
  #global_profil=profils_for_structure
  #global_profil[[i=1]]
  num_stru=length(global_profil[[1]])
  if(missing(strucure_names)) {
    strucure_names=c()
    for(f in seq(1,num_stru)){
      strucure_names[f]=paste("stru",f)
    }
  if(missing(colors)) {
    colors=rainbow(num_stru)
  } 
  NAMES<- paste(names(prot_cons))
  for( i in seq(1, length(global_profil))){
    len=c()
    TMPlist=global_profil[[i]]
    for (k in seq(1, length(TMPlist))){
      len[k]=length(TMPlist[[k]][[1]])
    }
    StruLen=max(len)
    plot(prot_cons[[i]], col ="black",type="l", main=paste(NAMES[i],"entropy score for ",prot_descript, pdb_name),pch = 20, xlim=c(0,length(prot_cons[[1]])),
         ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
    for (j in seq(1,length(len))){

      par(new=T)
      plot(TMPlist[[i]][[2]],TMPlist[[i]][[1]],col=colors[j],xlim=c(0,length(prot_cons[[1]])),ylim=c(0,1),xlab="",ylab="",main="",pch=j)
      
    }
    legend('topleft',c("protein",structure_names),lty=c(1,rep(0,j)),pch=c(-1,seq(1,j)),lwd=c(2.5,2.5),col=colors)
  }
  }
}
Names=c("Shannon","Schneider", "Kabat", "Landgraf", "TG")
bezTuneli=list()
profilet=list()
prot_cons=2
i=2
for(i in seq(1:length(prot_cons))){
  profilet[[i]]=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], [[i]]) 
  profilet2=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 2) 
  profilet3=entropy_profile(tunnel_file, uniprot, file, shift,prot_cons[[i]], 3)
  bezTuneli[[i]]=prot_cons[[i]][-intersect(intersect(profilet1[[2]],profilet3[[2]]),profilet3[[2]])]
  profilet3
  StruLen=max(length(profilet3[[1]]),length(profilet2[[1]]),length(profilet1[[1]]))
  ###################### STRUCTURES ON PROTEIN
  plot(prot_cons[[i]], col ="black",type="l", main=paste(Names[i],"entropy score for ",prot_descript, pdb_name),pch = 20, xlim=c(0,my_seq$len),
       ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
  par(new=T)
  plot(profilet1[[2]],prot_cons[[i]][profilet1[[2]]],col="blue",xlim=c(0,my_seq$len),ylim=c(0,1),xlab="",ylab="",main="",pch=21)
  par(new=T)
  plot(profilet2[[2]],prot_cons[[i]][profilet2[[2]]], col="green",xlim=c(0,my_seq$len),ylim=c(0,1),xlab="",ylab="",main="",pch=20)
  par(new=T)
  plot(profilet3[[2]],prot_cons[[i]][profilet3[[2]]],col="red",xlim=c(0,my_seq$len),ylim=c(0,1),xlab="",ylab="",main="",pch=22)
  legend('topleft',c("protein","T1","T2","T3"),lty=c(1,0,0,0),pch=c(-1,21,20,22),lwd=c(2.5,2.5),col=c("black","blue","green","red")) 
  ##############################3 PROFILES OF TUNNEL
  plot(profilet1[[1]], type="l", main=paste(Names[i]," entropy profile for tunnels in", pdb_name),xlim=c(0,StruLen),ylim=c(0,1), col="blue", xlab="distance from acitve site [AA]",ylab=paste(Names[i]," entropy score"))
  par(new=T)
  plot(profilet2[[1]], type="l", main="",xlim=c(0,StruLen),ylim=c(0,1), col="red", xlab="",ylab="")
  par(new=T)
  plot(profilet3[[1]], type="l", main="",xlim=c(0,StruLen),ylim=c(0,1), col="green", xlab="",ylab="")
  par(new=T)
  legend('bottomright',c("T1","T2","T3"),lty=c(1,1,1),lwd=c(2.5,2.5),col=c("blue","green","red")) 
  
}

nowa_struktura=read_structure("1QCZ_wt_petla.txt",uniprot,file, 3)
read_structure<- function(structure_file, sequence_id, alignment_file, shift){
  #structure_file- string with a name of structure file eg. "stru1.txt"
  #structutre_file="1QCZ_wt_petla.txt"
  #shift=3
  #sequence_id=uniprot
  #alignment_file= read.alignment(file="aln2_312_pro.fasta", format="fasta", forceToLower=F)
  file=as.data.frame(read.table(structutre_file))
  DAT <- data.frame(lapply(file, as.character), stringsAsFactors=FALSE)
  numb=as.numeric(DAT[1,])
  AA=DAT[2,]
  if(shift!=0){
    struct=numb+rep(shift,length(numb))
  }
  base_seq = find_seq(sequence_id, alignment_file,1)
  seq = rep("N",each = base_seq$len);
  
  seq[struct] = "T"

  aa_positions = which(s2c(base_seq$sequence) != "-")
  just_align = alignment_file[[3]]
  paramet = alignment_parameters(just_align)
  length_alignment = dim(alignment2matrix(paramet, just_align))[2]
  structures = rep("-",each = length_alignment);
  j = 1;
  for (a in aa_positions) {
    #aligning the tunnels information with the alignment sequence
    structures[a] = seq[j];
    j = j + 1;
  }
  structures=c(paste(sequence_id,"STRU"),structure_loop)
  return(structures)
  }