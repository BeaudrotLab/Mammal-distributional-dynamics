# Create dendrograms of ecological similarity for each TEAM site
library(cluster)
library(fields)

# Load files and functions
phylomat_mod <- read.csv(file="mean_phylomat2.csv")
traitmat_mod <- read.csv(file="mammal_trait_matrix.csv")
load("SitesBinary.RData")
source("funphylocom.R")
source("EDI_func_2.R")

################################################
#Combines Functional and Phylogenetic matrices
################################################

spp=phylomat_mod[,1]

phylomat=as.matrix(phylomat_mod[,-1])
traitmat=as.matrix(traitmat_mod[,-1])

rownames(phylomat)=spp
colnames(phylomat)=spp

rownames(traitmat)=spp
colnames(traitmat)=spp

FPDist.global=FPDist(PDist=phylomat,FDist=traitmat,a=0.7,p=2)

############################################################
########## Extract matrix for each TEAM site ###############
############################################################

VBcom <- gsub(" ", "_", rownames(SitesBinary[[1]]))
UDZcom <- gsub(" ", "_", rownames(SitesBinary[[2]]))
BIFcom <- gsub(" ", "_", rownames(SitesBinary[[3]]))
PSHcom <- gsub(" ", "_", rownames(SitesBinary[[4]]))
YANcom <- gsub(" ", "_", rownames(SitesBinary[[5]]))
NAKcom <- gsub(" ", "_", rownames(SitesBinary[[6]]))
RNFcom <- gsub(" ", "_", rownames(SitesBinary[[7]]))

VB.FPDist <- FPDist.global[charmatch(VBcom, rownames(FPDist.global)), charmatch(VBcom, rownames(FPDist.global))]
VB.FPDist <- VB.FPDist[-23, -23] # Remove Tapirus terrestris (no observations)
UDZ.FPDist <- FPDist.global[charmatch(UDZcom, rownames(FPDist.global)), charmatch(UDZcom, rownames(FPDist.global))]
BIF.FPDist <- FPDist.global[charmatch(BIFcom, rownames(FPDist.global)), charmatch(BIFcom, rownames(FPDist.global))]
BIF.FPDist <- BIF.FPDist[-11,-11] # Remove Cricetomys emini
PSH.FPDist <- FPDist.global[charmatch(PSHcom, rownames(FPDist.global)), charmatch(PSHcom, rownames(FPDist.global))]
PSH.FPDist <- PSH.FPDist[-17, -17] #Remove NA
YAN.FPDist <- FPDist.global[charmatch(YANcom, rownames(FPDist.global)), charmatch(YANcom, rownames(FPDist.global))]
NAK.FPDist <- FPDist.global[charmatch(NAKcom, rownames(FPDist.global)), charmatch(NAKcom, rownames(FPDist.global))]
RNF.FPDist <- FPDist.global[charmatch(RNFcom, rownames(FPDist.global)), charmatch(RNFcom, rownames(FPDist.global))]

################################################################################
########## Plot TEAM site mammal community similarity as a dendrograms ######### 
################################################################################
par(mar=c(1,2,0,1), oma=c(1, 0, 1, 0))
set.panel(2, 2)

pdf(file="Dendrograms1-4.pdf", onefile=TRUE, width=8, height=11)
par(mar=c(1,2,0,1), oma=c(1, 0, 1, 0))
set.panel(2, 2)
plot(agnes(VB.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=3, y=0.85, labels="VB")
plot(agnes(UDZ.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=3, y=0.7, labels="UDZ")
plot(agnes(BIF.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=3, y=0.7, labels="BIF")
plot(agnes(PSH.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=4, y=0.6, labels="PSH")
dev.off()

pdf(file="Dendrograms5-7.pdf", onefile=TRUE, width=8, height=11)
par(mar=c(1,2,0,1), oma=c(1, 0, 1, 0))
set.panel(2, 2)
plot(agnes(YAN.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=3, y=0.85, labels="YAN")
plot(agnes(NAK.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=3, y=0.58, labels="NAK")
plot(agnes(RNF.FPDist, diss=TRUE, stand=TRUE, method="average"), main="", xlab="", ylab="", sub="", which.plots=2, cex=0.8, hang=-1)
text(x=2, y=0.58, labels="RNF")
dev.off()

