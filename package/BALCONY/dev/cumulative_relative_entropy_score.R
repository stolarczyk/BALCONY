library(Biostrings)
library(BALCONY)

hmmbuild_path = "/home/mstolarczyk/Downloads/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild"

pairwiseAlignemnt_scores = readRDS("/home/mstolarczyk/Uczelnia/PROJEKT/BALCONY_git/BALCONY/package/BALCONY/pairwiseAlignment_mtx.rds")
dendro = hclust(pairwiseAlignemnt_scores,method = "average") 
clusters = cutree(dendro,h = 8000)
no_clusters = max(clusters)

for(i in seq_len(no_clusters)){
  sub_alignment = alignment
  leftover_alignment = alignment
  to_delete = which(clusters == clusters[i])
  sub_alignment$seq = sub_alignment$seq[to_delete]
  leftover_alignment$seq = leftover_alignment$seq[-to_delete]
  sub_converted_seq = lapply(sub_alignment$seq,s2c)
  leftover_converted_seq = lapply(leftover_alignment$seq,s2c)
  sub_alignment$nam = sub_alignment$nam[to_delete]
  leftover_alignment$nam = leftover_alignment$nam[-to_delete]
  leftover_alignment$nb = leftover_alignment$nb - length(to_delete)
  sub_alignment$nb = length(to_delete)
  leftover_fasta_name = paste("leftover_MSA_cluster_",i,".fasta",sep = "")
  sub_fasta_name = paste("sub_MSA_cluster_",i,".fasta",sep = "")
  write.fasta(sequences = converted_seq,names = sub_alignment$nam,file.out = fasta_name)
  
  hmm_command = paste(hmmbuild_path,"hmm.out",fasta_name)
  system(command = hmm_command)
}

