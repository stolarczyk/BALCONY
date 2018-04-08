pseudo_counts_sub<- function(alignment, substitution_mtx =NULL){
  if (is.null(substitution_mtx)) {
    load("./documentation/data/Gonnet_mtx.rd") # JEDNA KROPKA! MOŻE TRZEBA TO POPRAWIĆ W LANDGRAFIE 
    #substitution_mtx <- as.numeric(gonnet_mtx[[2]])
    substitution_mtx<- apply(X = gonnet_mtx[[2]],MARGIN = 2,FUN = function(X) as.numeric(as.character(X)))
    substitution_mtx=exp(substitution_mtx)
    colnames(substitution_mtx)<- gonnet_mtx[[1]]
    rownames(substitution_mtx)<- gonnet_mtx[[1]]
  }
  #warning - other substitution matrices may have different symbols (like '*', or other letters)
  #substitution_mtx<- exp(PAM250)
    #probab = create_probab()
  mtx_alignment = alignment2matrix(alignment)
  pseudoCounts= matrix(NA,nrow = length(append(AA_STANDARD,"-")),ncol =dim(mtx_alignment)[2])
  B = apply(X=mtx_alignment,MARGIN =2, FUN = function(X) (5*length(unique(X))))
  calc_ba<-function(column,B,substitution_mtx){
    pseudocounts = matrix(NA,nrow = length(append(AA_STANDARD,"-")),ncol =1)
    names(pseudocounts) <- append(AA_STANDARD,"-")
    N = sum((table(column)))
    AA=table(column)
    for(i in names(pseudocounts)){
    to_sum <- rep(NA, length(AA))
     for(j in seq_along(names(AA))){
       to_sum[j]<-(AA[names(AA)[j]]/N)*substitution_mtx[j,i]
     }
    pseudocounts[i]<- B*sum(to_sum)
    }
    return(t(pseudocounts))
  }
    
  for(column in seq_len(ncol(mtx_alignment))){
    pseudoCounts[,column]<- calc_ba(mtx_alignment[,column],B[column], substitution_mtx)
  }
  return(t(pseudoCounts))
}
  