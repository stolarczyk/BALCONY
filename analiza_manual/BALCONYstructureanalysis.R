
# Analysis of protein's structure part ------------------------------------


# get indices of tunnels in alignment
tunnel_index=get_structures_idx(structure)
vex=list()
vex[[1]]=entropy_data[[1]]
#find indices of protein in MSA 
prot_cons=get_prot_entropy(tunnel_index[[1]],entropy_data)
#plot entropy
plot_entropy(prot_cons, legend_pos = "topleft")
#get entropy for your structure
stru_entropy=get_structures_entropy(tunnel_index,entropy_data)
profiles_for_structure=prepare_structure_profile(structure, stru_entropy)

### PLOT STRUCTURE
plot_structure_on_protein(prot_cons, profils_for_structure,pdb_name,structure_names = c("T1", "T2", "T3"), legend_pos = "topleft")
#######TUNNELS ENTROPY 

## SCATTERPLOTS FOR PROTEINS AND TUNNELS
library(scales)
compare_cons_metrics(prot_cons,profils_for_structure, pdb_name)

#KOLMOGOROV-SMIRNOF TEST FOR QUALITY OF DISTRIBURIONS
EQUAL=smirnof_kolmogorov_test(prot_cons, profils_for_structure,1, pdb_name,T)
write(EQUAL, file = "SK_twosided.txt", sep = " ")
LESS=smirnof_kolmogorov_test(prot_cons, profils_for_structure,2, pdb_name,F)
write(LESS, file = "LESS.txt", sep = " ")

GREATER=smirnof_kolmogorov_test(prot_cons, profils_for_structure,3, pdb_name,F)
write(GREATER, file = "GREATER.txt", sep = " ")



EqualKS=list()
for(k in seq(1, length(prot_cons))){
  k_smirnof_test_equal=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_equal[[i]]=ks.test(EntropyTunel[[i]][k,],prot_cons[[k]])
  }
  EqualKS[[k]]=k_smirnof_test_equal
}
names(EqualKS)<-Names


GreaterKS=list()
for(k in seq(1, length(prot_cons))){
  k_smirnof_test_greater=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_greater[[i]]=ks.test(EntropyTunel[[i]][k,],prot_cons[[k]],alternative = "greater")
  }
  GreaterKS[[k]]=k_smirnof_test_greater
}
names(GreaterKS)<-Names

SmallerKS=list()
for(k in seq(1, length(prot_cons))){
  k_smirnof_test_smaller=list()
  for( i in seq(1:length(structure))){
    k_smirnof_test_smaller[[i]]=ks.test(EntropyTunel[[i]][k,],prot_cons[[k]],alternative = "less")
  }
  SmallerKS[[k]]=k_smirnof_test_smaller
}
names(SmallerKS)<-Names

