RealValET_conservativity <- function(alignment){
  ET_groups = seq(1,length(alignment$nam))
  dist <- seqinr::dist.alignment(alignment)
  clust <- hclust(dist, method = "average") #upgma
  groups <- list() # groups for all nodes
  for (i in ET_groups) {
    groups[[i]] <- cutree(clust,i)
  }
  
  aligned_matrix <- alignment2matrix(alignment = alignment)
  
  
  shannon_group <- function(group){
    counts = table(group)
    len = length(group)
    freq = counts/len
    Shannon = -sum(freq * log(freq))
    
  }
  #test= aligned_matrix[,754]
  realValET <- function(column, group_tree){
    nodes = length(group_tree) - 1
    node_group = rep(0,nodes)
    for (i in 1:nodes) {
      temp = split(column, group_tree[[i]])
      group = rep(0,length(temp))
      for (j in 1:length(temp)) {
        group[j] = shannon_group(temp[[j]])
      }
      node_group[i] = (1/i)*sum(group)
    }
    RO = 1 + sum(node_group)
    return(RO)
  }
  
  x = rep(0,dim(aligned_matrix)[2])
  for (i in 1:dim(aligned_matrix)[2]) {
    x[i] = realValET(aligned_matrix[,i], groups)
  }
  
  return(x)
  #ET = lapply(X = aligned_matrix[,700-705],FUN = realValET,group_tree=groups)
}