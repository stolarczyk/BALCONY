library(Biostrings)
library(BALCONY)
library(readr)

preprocess_hmm_output <- function(hmm_out){
  size = dim(hmm_out)
  hmm_out = hmm_out[seq(from = 1,to = size[1]-1,by = 3),-1]
  sites = as.numeric(unlist(hmm_out[,21]))
  probs = hmm_out[,1:20]
  probs = exp(-probs)
  return(list(probabilities=probs,alignment_positions=sites))
}

##get the system hmmbuild binary absolute path
hmmbuild_path = "/home/mstolarczyk/Downloads/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild"

pairwiseAlignemnt_scores = readRDS("/home/mstolarczyk/Uczelnia/PROJEKT/BALCONY_git/BALCONY/package/BALCONY/dev/pairwiseAlignment_mtx.rds")

#perform hierarchical clustering on the distance matrix obtained from pairwise alignemnt matrix
dendro = hclust(dist(pairwiseAlignemnt_scores),method = "average")
#cut the tree and get clusters
clusters = cutree(dendro,h = 8000)
#get the number of clusters
no_clusters = max(clusters)
RE=matrix(data = NA,nrow = no_clusters,ncol = length(s2c(alignment$seq[1])))
CRE = c()
Z=c()
#iterate over each clusters
for(i in seq_len(no_clusters)){
  #create copies of the alignemnt object to work on 
  sub_alignment = alignment
  leftover_alignment = alignment
  #get the indices of the sequences in the cluster of interest
  to_delete = which(clusters == clusters[i])
  #delete approproate sequences and their names in the objects and adjust the sequeces counts
  sub_alignment$seq = sub_alignment$seq[to_delete]
  leftover_alignment$seq = leftover_alignment$seq[-to_delete]
  sub_converted_seq = lapply(sub_alignment$seq,s2c)
  leftover_converted_seq = lapply(leftover_alignment$seq,s2c)
  sub_alignment$nam = sub_alignment$nam[to_delete]
  leftover_alignment$nam = leftover_alignment$nam[-to_delete]
  leftover_alignment$nb = leftover_alignment$nb - length(to_delete)
  sub_alignment$nb = length(to_delete)
  #write updated objects to files
  leftover_fasta_name = "leftover_MSA.fasta"
  sub_fasta_name = "sub_MSA.fasta"
  write.fasta(sequences = leftover_converted_seq,names = leftover_alignment$nam,file.out = leftover_fasta_name)
  write.fasta(sequences = sub_converted_seq,names = sub_alignment$nam,file.out = sub_fasta_name)
  #prepare commands to run HMMER
  sub_hmm_command = paste(hmmbuild_path,"sub_hmm.out",sub_fasta_name)
  leftover_hmm_command = paste(hmmbuild_path,"leftover_hmm.out",leftover_fasta_name)
  system(command = sub_hmm_command,wait = T)
  system(command = leftover_hmm_command,wait = T)
  leftover_hmm <- read_table2("leftover_hmm.out", col_names = FALSE, skip = 21)
  sub_hmm <- read_table2("sub_hmm.out", col_names = FALSE, skip = 21)
  leftover_prob = preprocess_hmm_output(leftover_hmm)$probabilities
  leftover_pos = preprocess_hmm_output(leftover_hmm)$alignment_positions
  sub_prob = preprocess_hmm_output(sub_hmm)$probabilities
  sub_pos = preprocess_hmm_output(sub_hmm)$alignment_positions
  intersection_pos = intersect(sub_pos,leftover_pos)

  for(pos in seq_len(length(s2c(alignment$seq[1])))){
    which_pos = which(intersection_pos == pos)
    if(length(which_pos)>0){
      pre_RE = c()
      for(aa in seq_len(dim(sub_prob)[2])){
        sub_prob_pos = which(sub_pos == pos)
        leftover_prob_pos = which(leftover_pos == pos)
        pre_RE[aa] = sub_prob[sub_prob_pos,aa]*log(sub_prob[sub_prob_pos,aa]/leftover_prob[leftover_prob_pos,aa])
      }
      RE[i,pos] = sum(pre_RE,na.rm = T)
    } else{
      RE[i,pos]=0
    }
    CRE[pos] = sum(RE[,pos])
  }
}

for(pos in seq_len(length(s2c(alignment$seq[1])))){
  Z[pos] = (CRE[pos]-mean(CRE))/sd(CRE)
}


