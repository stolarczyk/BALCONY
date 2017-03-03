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
  names(out[[2]])= names(structure)
  return(out)
}
get_prot_entropy<- function(whole_prot,score_list){
  #documentation get_prot_entropy.Rd
  #allows to get idx of whole protein in alignment
  #returns list of entropy for protein
 
  prot_cons=list()
  for(i in seq(1, length(score_list))){
    prot_cons[[i]]=score_list[[i]][whole_prot]
  }
  names(prot_cons)<-names(score_list)
  return(prot_cons)
}

plot_entropy<- function(prot_cons, colors,impose=NULL){
  #plots scores on one plot, if colors are not specified plot as rainbow
  #automagically uses name of pdb FIXME
  # recive list of entropy scores and list of Names in this list
  if(missing(colors)){
    colors<- rainbow(length(prot_cons))
  }
  if(is.null(impose)){
    impose<-T
  }
 
  plot(prot_cons[[1]],ylim=c(0,1), col=colors[1],xlab="amino acid",ylab="entropy score", main=paste("entropy score for ",pdb_name), type="l")
  for(i in seq(2, length(prot_cons))){
    par(new=impose)
    plot(prot_cons[[i]],ylim=c(0,1), col=colors[i],xlab="",ylab="", main="", type="l")
  }
  legend("topleft",names(prot_cons), col=colors, lty =c(1))
}

### entropy for stuff
get_structures_entropy<- function(structure_index, score_list){
  #structure_index is a list of indexes in alignment of protein and structures 
  #SCORE_LIST list of entropies for whole alignment
  #output is a list of matrixes where each row contains values of entropy for AA in structure
  t_index=structure_index$structureIndices
  Entropy=list()
  lengths=list()
  for (i in seq(1:length(t_index))){
   lengths[[i]]=length(t_index[[i]])
   output=matrix(NA,nrow=length(score_list),ncol=lengths[[i]])
    for (j in seq(1:length(score_list))){
      output[j,]=score_list[[j]][t_index[[i]]]
    }
    rownames(output)<-names(score_list)
    Entropy[[i]]=output
   
  }
  return(Entropy)
}

entropy_for_all<- function(tunnel_file, uniprot, file, shift,prot_cons){
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

profils_for_structure=profile_for_all(tunnel_file,uniprot,file,3,prot_cons)
length(profils_for_structure$s)
##### FINISH ME! :D 
plot_structure_on_protein<- function(protein_entropy, structure_profils, pdb_name, colors, structure_names=NULL){
  #protein entropy- list of an entropy values (with different scores)
  #structure_profils - list of entropy scores for structures (as list) where [[1]] entropy value
    # and [[2]] index in protein
  protein_entropy=prot_cons
  structure_profils=profils_for_structure
  
  if(length(protein_entropy) == length(structure_profils)){
    score_names=names(protein_entropy)
    prot_lenght=length(protein_entropy[[1]])
    StruLen=length(structure_profils[[1]])
    if(missing(colors)) {
      colors=rainbow(StruLen)
    } 
    if(is.null(structure_names)){
      structure_names=c()
      for(i in seq(1,StruLen)){
        structure_names[i]=paste("stru", i)
      }
    }
    for(i in seq(1, length(protein_entropy))){
      
      plot(protein_entropy[[i]], col ="black",type="l", main=paste(score_names[i]," score for ", pdb_name),pch = 20, xlim=c(0,prot_lenght),
           ylim=c(0.0,1.0), xlab='Amino Acid', ylab='Entropy')
      for(j in seq(1,StruLen)){
        par(new=T)
        plot(structure_profils[[i]][[j]][[2]],structure_profils[[i]][[j]][[1]], col=colors[[j]], main="",pch = j, xlim=c(0,prot_lenght),
             ylim=c(0.0,1.0), xlab='', ylab='')
      }
      legend('topleft',c(pdb_name,structure_names),lty=c(1,rep(0,j)),pch=c(-1,seq(1,j)),lwd=c(2.5,2.5),col=c("black",colors))
      
      
    }
  }
  else print("The lists contain different number of conservation/entropy scores!")
}

# STATISTICS --------------------------------------------------------------
protein=prot_cons
structure_cons=profils_for_structure
compare_cons_metrics(prot_cons,profils_for_structure, pdb_name)
compare_cons_metrics<- function(protein, structures_cons, pdb_name){
  metrics_count=length(protein)
  structures_count=length(structure_cons[[1]])
  colors=rainbow(structures_count)
  structure_names=c()
  for(i in seq(1,structures_count)){
    structure_names[i]=paste("stru", i)}
  for(i in seq(1, metrics_count)){
    for(j in seq(1, metrics_count)){
      if(i!=j){
        plot(protein[[i]], protein[[j]],main=paste('Scatterplot of', names(protein)[i],"vs. ", names(protein)[j]), xlim=c(0,1), ylim=c(0,1), xlab=names(protein)[i],ylab=names(protein)[j],pch=20, col="black")
        for(k in seq(1,structures_count)){
          par(new=T)
          plot(structure_cons[[i]][[k]][[1]], structure_cons[[j]][[k]][[1]], col=colors[k], pch = k, main="", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="")
        }
        legend('bottomright',c(pdb_name,structure_names),pch=c(20,seq(1,k)),col=c("black",colors))
      }
    }
  }
  
}
protein_cons=protein
smirnof_kolmogorov_test<-function(protein_cons, structure_cons,alternative,plot=T){
  metrics_count=length(protein_cons)
  structures_count=length(structure_cons[[1]])
  alt_hip=c("two.sided","less", "greater")[alternative]
  for(i in seq(1,metrics_count)){
    for(j in seq(1, structures_count)){
      i=1; j=1
      reference=protein_cons[[i]][-structure_cons[[i]][[j]][[2]]]
      ks.test(reference, structure_cons[[1]][[1]][[1]])$p.value
      }
  }
  if(plot=T){
    colors=rainbow(structures_count)
  }
}

# READ NEW STRUCTURE FILE -------------------------------------------------


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