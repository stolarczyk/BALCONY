library(Biostrings)
library(BALCONY)
data("alignment")
convert_then_align <- function(seqj,seqi){
  seqi = s2c(seqi)
  seqi = c2s(seqi[seqi!= "-"])
  seqj = s2c(seqj)
  seqj = c2s(seqj[seqj!= "-"])
  return(pairwiseAlignment(seqi, seqj,scoreOnly=T))
}

score_mtx = matrix(NA,nrow = alignment$nb,ncol = alignment$nb)
for(i in seq_len(alignment$nb)){
  seqi = alignment$seq[i]
  print(i)
  score_mtx[i,] = sapply(alignment$seq, convert_then_align, seqi)
}
