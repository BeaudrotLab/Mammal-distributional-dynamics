

# Create dendrograms of ecological similarity for BIF, RNF, UDZ, VB


mat<-read.delim("mammal_com.csv", sep=",", row.names=1)

mat_BIF <- mat[,1:60]
mat_BIF <- mat_BIF[rowSums(mat_BIF)>0,]
com_BIF <- data.frame(as.character(rownames(mat_BIF)))
rownames(com_BIF) <- com_BIF[,1]
com_BIF <- com_BIF[,-1]

mat_RNF <- mat[,209:268]
mat_RNF <- mat_RNF[rowSums(mat_RNF)>0,]
com_RNF <- data.frame(as.character(rownames(mat_RNF)))
rownames(com_RNF) <- com_RNF[,1]
com_RNF <- com_RNF[,-1]

mat_UDZ <- mat[,269:329]
mat_UDZ <- mat_UDZ[rowSums(mat_UDZ)>0,]
com_UDZ <- data.frame(as.character(rownames(mat_UDZ)))
rownames(com_UDZ) <- com_UDZ[,1]
com_UDZ <- com_UDZ[,-1]

mat_VB <- mat[,330:389]
mat_VB <- mat_VB[rowSums(mat_VB)>0,]
com_VB <- data.frame(as.character(rownames(mat_VB)))
rownames(com_VB) <- com_VB[,1]
com_VB <- com_VB[,-1]

mat_YAN <- mat[,390:450]
mat_YAN <- mat_YAN[rowSums(mat_YAN)>0,]
com_YAN <- data.frame(as.character(rownames(mat_YAN)))
rownames(com_YAN) <- com_YAN[,1]
com_YAN <- com_YAN[,-1]


mat <- mat_YAN
com <- com_YAN

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


plot(phy, main="")

# Really need a dendrogram of the EDI by site with species names
library(cluster)
library(blockmodeling)
source("funphylocom.R")
source("EDI_func_2.R")
YAN
################################################
#Combines Functional and Phylogenetic matrices
################################################
phylomat_mod <- read.csv(file="mammal_trait_matrix.csv")
traitmat_mod <- read.csv(file="mammal_trait_matrix.csv")

spp=phylomat_mod[,1]

phylomat=as.matrix(phylomat_mod[,-1])
traitmat=as.matrix(traitmat_mod[,-1])

rownames(phylomat)=spp
colnames(phylomat)=spp

rownames(traitmat)=spp
colnames(traitmat)=spp

FPDist.global=FPDist(PDist=phylomat,FDist=traitmat,a=0.7,p=2)

EDI_com <- FPDist.global[which(rownames(FPDist.global) %in% rownames(com)), which(rownames(FPDist.global) %in% rownames(com))]
clustered <- agnes(EDI_com, diss=TRUE, stand=TRUE, method="average")
plot(clustered, which.plots=c(2), cex=1, main="")



#EDI_BIF <- FPDist.global[which(rownames(FPDist.global) %in% rownames(com)), which(rownames(FPDist.global) %in% rownames(com))]
#clustered_BIF <- agnes(EDI_BIF, diss=TRUE, stand=TRUE, method="average")
#plot(clustered_BIF, which.plots=c(2), cex=1, main="")
