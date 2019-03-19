library(blockmodeling)

phylomat_mod=read.csv(file="mean_phylomat2.csv")

traitmat_mod=read.csv(file="mammal_trait_matrix.csv")

load('SitesBinary.RData')

source('funphylocom.R') #includes function FPDist()

source('EDI_func_2.R')

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

################################################
#Plot phylogenetic, trait and combined matrices
#################################################

#Uncomment the following lines if you want to plot the matrices.

#plot.mat(phylomat) #This is a cool function of blockmodeling to plot matrices
#dev.copy2pdf(file="phylomat.pdf")

#plot.mat(traitmat)
#dev.copy2pdf(file="traitmat.pdf")

#plot.mat(FPDist.global)
#dev.copy2pdf(file="FPDist.global.pdf")

#######################################################
#All observed
#Loops through all sites calculating observed EDIs
########################################################
Obs.edi=vector("list",7) 
# (1) VB; (2) UDZ; (3) BIF; (4) PSH; (5) YAN; (6) NAK; (7) RNF

for (s in 1:length(SitesBinary)){

	CT=as.data.frame(SitesBinary[[s]])
	Obs.edi[[s]]=EDI(CT)
	
	}

##################
# Nulls
###################
#z=1
Z=vector("list",7)

for (z in 1:length(SitesBinary)){
Site=SitesBinary[[z]] #Original site species table
Site.edi=Obs.edi[[z]] #EDI values for each spp in the camera trap


ls=lapply(Site.edi, function(x) length(x)) #Counts # of species in each camera trap
spp=as.numeric(Site[,1]) #saves the species number unique ID
weights=as.numeric(rowSums(Site[,-1])/length(Site.edi)) #Calculates the proportion of cameras in which the species is present in this site to later weight the resampling.
vls=numeric(10000) 

W<-vector("list",length(Site.edi)) #empty vector for standardized re-sampled ecological distance matrix
for (j in 1:length(Site.edi)){ #Camera trap
	
	to=as.numeric(ls[j])
	for(h in 1:to){ #Species in each camera trap
		
		n=as.numeric(names(Site.edi[[j]][h])) #chooses a focal species in a camera trap
		
		for(u in 1:10000){
			
			if (ls[j]>0){#number of species co-occuring with focal species
			hm=as.numeric(ls[j])-1 
			spp2=spp[-which(spp==n)] #species co-ocurring with species n 
			weights2=weights[-which(spp==n)] #weight of coocurring species
			
			#The next line of code samples from the pool of species present in the site (region) random species to co-occur.  This is weighted by the proportion of camera traps in which each species is present.
			rnd=sample(spp2,hm,replace=FALSE,prob=weights2)
			vls[u]=sum(as.numeric(FPDist.global[n,rnd]))/length(rnd) #This calculates the average relationship between the focal species n and others sampled from the regional pool.		
					}else{vls[h]=0}
			
		}
		
		W[[j]][h]=-(as.numeric(Site.edi[[j]][h])-mean(vls))/sd(vls)
	}
}
cat("The value of z is",z)
Z[[z]]=W
rm(W)
}

# (1) VB; (2) UDZ; (3) BIF; (4) PSH; (5) YAN; (6) NAK; (7) RNF

names(Z)=c("VB","UDZ","BIF","PSH","YAN","NAK","RNF")

############################################################
#Identify each element in the list by the species id number
############################################################
#There might be a more elegant way of doing this...but this one works

for (q in 1:7){#Number of sites
	for (w in 1:length(Obs.edi[[q]])){
		names(Z[[q]][[w]])=names(Obs.edi[[q]][[w]])
	}
}	
	
save(Z,file="Scaled_FPDist.RData")
save(spp,file="spp.RData")





