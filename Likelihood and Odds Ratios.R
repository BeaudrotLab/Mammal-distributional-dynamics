# Purpose: Run likelihood ratio tests to compare top model to null model for species with 
#          top models containing covariates and delta AIC > 2 from null model. For GEB revisions.
# Author: Lydia Beaudrot
# Date: July 27, 2018; modified for GEB proofs on March 18, 2019

# NOTE BENE: For reporting odds ratios less than one, use 1/odds_ratios calculated with get.se (e.g., provides the number of times more likely to decline)


rm(list=ls())
library(unmarked)
library(plyr)
library(AICcmodavg)


######################################
# Define helper functions
######################################

isEmpty <- function(x) {
  return(length(x)==0)
}


CondNum <- function(model){
  max(eigen(hessian(model))$values)/min(eigen(hessian(model))$values)
}

get.se=function(model){
  
  c=exp(coef(model))
  var.diag=diag(vcov(model))
  or.se=sqrt(c^2*var.diag)
  df=data.frame(odds_ratios=c,se=or.se)
  return(df)
  
}


######################################
# Load formatted data
######################################


load("Occupancy Modelling/All_covs.RData")
load("Occupancy Modelling/Species7sites_Include.RData")
load("Occupancy Modelling/BIOTIC_Include.RData")

Species_data <- Species7sites_Include
BIOTIC <- BIOTIC_Include

nms=names(Species_data)

# Empty objects for loop
results.all=list()
mods.all=list()
results.table.ma=list()
results.table.aic=list()
colext.transformed=list()


# Function for formatting data for individual species (fm.data())

fm.data <- function(k){
  
sp.name <- names(Species_data)[k]
sp.name
species <- Species_data[[k]]
site <- substr(names(Species_data)[[k]],1,3)
site

# Define covariates

#siteCovs
site_covs <- paste(site, "covs", sep="_")
covs <- All_covs[names(All_covs)==site_covs]
Elevation <- unlist(as.matrix(sapply(covs, "[", 1)))
ForestLossCT <- unlist(as.matrix(sapply(covs, "[", 2)))
ForestGainCT <- unlist(as.matrix(sapply(covs, "[", 3)))
Biotic <- BIOTIC[[k]]
site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT, Biotic)

#yearlySiteCovs
Tmin <- as.data.frame(sapply(covs, "[", 4))
Tmax <- as.data.frame(sapply(covs, "[", 5))
Tvar <- as.data.frame(sapply(covs, "[", 6))
Tsd <- as.data.frame(sapply(covs, "[", 7))
Tmean <- as.data.frame(sapply(covs, "[", 8))

# Define number of primary periods
to=dim(Tmin)[2]
years=as.character(1:to)
years=matrix(years,nrow(species),to,byrow=TRUE)

# Create object with data formatted for unmarked

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin, Tmax=Tmax, Tvar=Tvar, Tsd=Tsd, Tmean=Tmean), siteCovs=site.covs, numPrimary=dim(Tmin)[2])
umf
}

############## BEGIN LIKELIHOOD RATIO TESTS #################

k <- 50 # NAK.Atherurus macrourus; ~1 ~ 1 ~ Tmin ~ Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm2.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Tmin,
             pformula=~Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm2.2)
get.se(fm2.2)

k <- 30 # PSH.Atherurus macrourus; ~1 ~Biotic ~1 ~Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm6.1=colext(psiformula=~1,
             gammaformula=~Biotic,
             epsilonformula=~1,
             pformula=~Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm6.1)
get.se(fm6.1)

k <- 8 # UDZ.Bdeogale crassicauda; ~1 ~ 1 ~ Tmax ~ Tmax
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm3.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Tmax,
             pformula=~Tmax,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm3.2)
get.se(fm3.2)

k <- 21 # BIF.Caracal aurata; ~1 ~ Tvar ~ 1 ~ Tvar
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm5.1=colext(psiformula=~1,
             gammaformula=~Tvar,
             epsilonformula=~1,
             pformula=~Tvar,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm5.1)
get.se(fm5.1)

k <- 22 # BIF.Cephalophus nigrifrons; ~1 ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm7.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Biotic + Tmin,
             pformula=~Biotic + Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm7.2)
get.se(fm7.2)

k <- 23 # BIF.Cephalophus silvicultor; ~1 ~ 1 ~ Biotic * Tmin ~ Biotic * Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm8.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Biotic * Tmin,
             pformula=~Biotic * Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm8.2)
get.se(fm8.2)

k <- 10 # UDZ.Cephalophus spadix; ~1 ~ 1 ~ Biotic * Tmin ~ Biotic * Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm8.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Biotic * Tmin,
             pformula=~Biotic * Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm8.2)
get.se(fm8.2)

k <- 11 # UDZ.Cercocebus sanjei; ~1 ~ Tvar * Biotic ~ 1 ~ Tvar * Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm12.1=colext(psiformula=~1,
              gammaformula=~Tvar * Biotic,
              epsilonformula=~1,
              pformula=~Tvar * Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm12.1)
get.se(fm12.1)


k <- 24 # BIF.Cercopithecus lhoesti; ~1 ~ Tmax + Biotic ~ 1 ~ Tmax + Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm9.1=colext(psiformula=~1,
             gammaformula=~Tmax + Biotic,
             epsilonformula=~1,
             pformula=~Tmax + Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm9.1)
get.se(fm9.1)

k <- 13 # UDZ.Cricetomys gambianus; ~1 ~ 1 ~ Tmax ~ Tmax
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm3.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Tmax,
             pformula=~Tmax,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm3.2)
get.se(fm3.2)

k <- 1; # VB_.Cuniculus paca; ~1 ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm7.2=colext(psiformula=~1,
             gammaformula=~1,
             epsilonformula=~Biotic + Tmin,
             pformula=~Biotic + Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm7.2)
get.se(fm7.2)

k <- 40 # YAN.Cuniculus paca; ~1 ~ Biotic + Tmin ~ 1 ~ Biotic + Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm7.1=colext(psiformula=~1,
             gammaformula=~Biotic + Tmin,
             epsilonformula=~1,
             pformula=~Biotic + Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm7.1)
get.se(fm7.1)

k <- 41 # YAN.Dasyprocta fuliginosa; ~1 ~ Tmin ~ 1 ~ Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm2.1=colext(psiformula=~1,
             gammaformula=~Tmin,
             epsilonformula=~1,
             pformula=~Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm2.1)
get.se(fm2.1)

k <- 2 # VB_.Dasyprocta punctata; ~1 ~ Tmin ~ 1 ~ Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm2.1=colext(psiformula=~1,
             gammaformula=~Tmin,
             epsilonformula=~1,
             pformula=~Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm2.1)
get.se(fm2.1)

k <- 3 # VB_.Dasypus novemcinctus; ~1 ~ Tvar * Biotic ~ Tvar * Biotic ~ Tvar * Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm12=colext(psiformula=~1,
            gammaformula=~Tvar * Biotic,
            epsilonformula=~Tvar * Biotic,
            pformula=~Tvar * Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm12)
get.se(fm12)

k <- 43 # YAN.Eira barbara; ~1 ~ Tvar ~ 1 ~ Tvar
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm5.1=colext(psiformula=~1,
             gammaformula=~Tvar,
             epsilonformula=~1,
             pformula=~Tvar,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm5.1)
get.se(fm5.1)

k <- 56 # RNF.Fossa fossana; ~Elevation ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm17.2=colext(psiformula=~Elevation,
             gammaformula=~1,
             epsilonformula=~Biotic + Tmin,
             pformula=~Biotic + Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm17.2)
get.se(fm17.2)

k <- 14 # UDZ.Genetta servalina; ~Elevation ~ Tmax * Biotic ~ 1 ~ Tmax * Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm20.1=colext(psiformula=~1,
              gammaformula=~Tmax * Biotic,
              epsilonformula=~1,
              pformula=~Tmax * Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm20.1)
get.se(fm20.1)

k <- 45 # YAN.Mazama americana; ~Elevation ~ Tmax ~ 1 ~ Tmax
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm13.1=colext(psiformula=~Elevation,
             gammaformula=~Tmax,
             epsilonformula=~1,
             pformula=~Tmax,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm13.1)
get.se(fm13.1)

k <- 5 # VB_.Mazama temama; ~1 ~ Tmax * Biotic ~ 1 ~ Tmax * Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm10.1=colext(psiformula=~1,
              gammaformula=~Tmax * Biotic,
              epsilonformula=~1,
              pformula=~Tmax * Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm10.1)
get.se(fm10.1)

k <- 46 # YAN.Nasua nasua; ~Elevation ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm17.2=colext(psiformula=~Elevation,
              gammaformula=~1,
              epsilonformula=~Biotic + Tmin,
              pformula=~Biotic + Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm17.2)
get.se(fm17.2)

k <- 27 # BIF.Pan troglodytes; ~Elevation ~ Tmax + Biotic ~ Tmax + Biotic ~ Tmax + Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm90.1=colext(psiformula=~Elevation,
              gammaformula=~Tmax + Biotic,
              epsilonformula=~Tmax + Biotic,
              pformula=~Tmax + Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm90.1)
get.se(fm90.1)



k <- 18 # UDZ.Paraxerus vexillarius; ~Elevation ~ Tmin ~ 1 ~ Tmin
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm12.1=colext(psiformula=~Elevation,
             gammaformula=~Tmin,
             epsilonformula=~1,
             pformula=~Tmin,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm12.1)
get.se(fm12.1)
# NB use 1/col(Tmin) for inverse value

k <- 6 # VB_.Pecari tajacu; ~1 ~ Tvar ~ 1 ~ Tvar
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm5.1=colext(psiformula=~1,
             gammaformula=~Tvar,
             epsilonformula=~1,
             pformula=~Tvar,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm5.1)
get.se(fm5.1)

k <- 29 # BIF.Potamochoerus larvatus; ~Elevation ~ 1 ~ Biotic ~ Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm95.1=colext(psiformula=~Elevation,
             gammaformula=~1,
             epsilonformula=~Biotic,
             pformula=~Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm95.1)
get.se(fm95.1)

k <- 19 # UDZ.Potamochoerus larvatus; ~1 ~ Tmax ~ 1 ~ Tmax
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm96.1=colext(psiformula=~1,
              gammaformula=~Tmax,
              epsilonformula=~1,
              pformula=~Tmax,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm96.1)
get.se(fm96.1)

k <- 53 # NAK.Tragulus kanchil; ~Elevation ~ Tmax + Biotic ~ 1 ~ Tmax + Biotic
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm90.1=colext(psiformula=~Elevation,
              gammaformula=~Tmax + Biotic,
              epsilonformula=~1,
              pformula=~Tmax + Biotic,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm90.1)
get.se(fm90.1)

k <- 38 # PSH.Tragulus kanchil; ~1 ~ 1 ~ Tvar ~ Tvar
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm97.2=colext(psiformula=~1,
              gammaformula=~1,
              epsilonformula=~Tvar,
              pformula=~Tvar,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm97.2)
get.se(fm97.2)

k <- 39 # PSH.Tupaia glis; ~1 ~ Tvar ~ 1 ~ Tvar
fm0=colext(psiformula=~1,
           gammaformula=~1,
           epsilonformula=~1,
           pformula=~1,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
fm5.1=colext(psiformula=~1,
             gammaformula=~Tvar,
             epsilonformula=~1,
             pformula=~Tvar,data=fm.data(k),method="L-BFGS-B",control=list(maxit=20000))
LRT(fm0, fm5.1)
get.se(fm5.1)
