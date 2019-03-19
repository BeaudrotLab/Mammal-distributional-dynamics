# Create phylogenetic and trail matrices to use as input for FPDist

library(picante)
library(ecodist)
library(vegan)
library(FD)

#input data
mat<-read.delim("mammal_com.csv", sep=",", row.names=1)
com<-read.delim("species_list.csv", sep=",", row.names=1)
traits<-read.delim("mammal_traits.csv", sep=",", row.names=1)

#load unresolved phylo
phylo<-read.nexus("mammalST_MSW05_all.tre")
phylo<-phylo[[1]]

#load resolved phylo
phylo_res<-read.nexus("FritzTree.rs200k.100trees.tre")
phylo_names<-read.csv("global_tip_labels_mod.csv")
phylo_names1<-as.character(phylo_names$sp)

### match phylo data with incidence
phylomat<-list()
for (i in 1:100)
{
  my_phylo_res<-phylo_res[[i]]
  #modify species names in phylo to match TEAM data
  #write.table(phylo_res$tip.label,file="globa_tip_labels.txt")
  my_phylo_res$tip.label<-phylo_names1
  
### prune phylo to match camera trap data
myphylosp<-match.phylo.data(my_phylo_res, com)
phy<-myphylosp$phy
phydist<-as.data.frame(cophenetic(phy))
phylomat[[i]]<-as.matrix(phydist)
}

mat.new=matrix(nrow=126,ncol=126)
for (k in 1:126){
  for (l in 1:126){
    
    cm=numeric()
    for(h in 1:100){
      cm=c(cm,phylomat[[h]][k,l])
    }
    mat.new[k,l]=mean(cm)
  }
  
}

mean_phylo_mat<-data.frame(mat.new)
row.names(mean_phylo_mat)<-row.names(phydist)
names(mean_phylo_mat)<-names(phydist)

#write.csv(mean_phylo_mat, file ="mean_phylomat.csv")

#Sort mean_phylo_mat alphabetically by species name
mean_phylo_mat2 <- (mean_phylo_mat[order(rownames(mean_phylo_mat)),order(colnames(mean_phylo_mat))])
write.csv(mean_phylo_mat2, file="mean_phylomat2.csv")


#create trait matrix
threetraits<-traits[1:3]
trait_mat_equal<-as.matrix(gowdis(threetraits))
write.matrix(trait_mat_equal, file="mammal_trait_matrix.csv", sep=",")

