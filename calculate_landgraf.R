consensus=consensus_seq
consensusG=consensus;
#Dividing AAs into groups of similar AAs in consensus
for (i in seq(1,length(consensus))){ 
  if (consensus[i] =="I" || consensus[i] =="L" || consensus[i] == "V")
    consensusG[i] = "G1"
  if (consensus[i] =="F" || consensus[i] =="Y" || consensus[i] == "W")
    consensusG[i] = "G2"
  if (consensus[i] =="A" || consensus[i] =="G")
    consensusG[i] = "G3"
  if (consensus[i] =="P")
    consensusG[i] = "G4"
  if (consensus[i] =="N" || consensus[i] =="Q")
    consensusG[i] ="G5"
  if (consensus[i] =="R" || consensus[i] =="H" || consensus[i] == "K")
    consensusG[i] = "G6"
  if (consensus[i] =="D" || consensus[i] =="E")
    consensusG[i] = "G7"
  if (consensus[i] =="M" || consensus[i] =="C" )
    consensusG[i] = "G8"
}
parametry=alignment_parameters(aligned_sequences);
weights=cons2seqs_sim(parametry,aligned_sequences_matrix,consensus_seq,grouping_method);
w2=abs(weights-1)
landgraf_normalized=Landgraf_conservation("GONNET.txt",aligned_sequences_matrix,weights = w2)

