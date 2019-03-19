# Purpose: R Code for Dynamic Occupancy Modeling (Abiotic & Biotic Drivers)
# Publication Title: Local temperature and ecological similarity drive distributional dynamics of tropical mammals worldwide
# Publication Authors: Beaudrot, Acevedo, Lessard et al.
# Journal: Global Ecology & Biogeography


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


############################################################################
################ Begin loop to run all models and extract results ##########
############################################################################

for(k in 1:length(nms)){
  print(k)
  
  # Define species for analysis and site using index value for list of all species
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
  
  # Create list to hold model set for model selection
  
  mods=list()
  
  # Null ##################################################################################
  try((fm0=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm0")) {
    if(CondNum(fm0)<5000){
      if(CondNum(fm0)>0){mods=c(mods,fm0)}
    } 
  }
  
##### COVARIATE MODELS ########################################################################
  
  # Tmin only ##################################################################################
  try((fm2=colext(psiformula=~1,
                  gammaformula=~Tmin,
                  epsilonformula=~Tmin,
                  pformula=~Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2")) {
    if(CondNum(fm2)<5000){
      if(CondNum(fm2)>0){mods=c(mods,fm2)}
    } 
  }
  
  try((fm2.1=colext(psiformula=~1,
                    gammaformula=~Tmin,
                    epsilonformula=~1,
                    pformula=~Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.1")) {
    if(CondNum(fm2.1)<5000){
      if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
    } 
  }
  
  try((fm2.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmin,
                    pformula=~Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.2")) {
    if(CondNum(fm2.2)<5000){
      if(CondNum(fm2.2)>0){mods=c(mods,fm2.2)}
    } 
  }
  
  try((fm2.3=colext(psiformula=~Elevation,
                    gammaformula=~Tmin,
                    epsilonformula=~Tmin,
                    pformula=~Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.3")) {
    if(CondNum(fm2.3)<5000){
      if(CondNum(fm2.3)>0){mods=c(mods,fm2.3)}
    } 
  }
  
  try((fm2.4=colext(psiformula=~Elevation,
                    gammaformula=~Tmin,
                    epsilonformula=~1,
                    pformula=~Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.4")) {
    if(CondNum(fm2.4)<5000){
      if(CondNum(fm2.4)>0){mods=c(mods,fm2.4)}
    } 
  }
  
  try((fm2.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmin,
                    pformula=~Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.5")) {
    if(CondNum(fm2.5)<5000){
      if(CondNum(fm2.5)>0){mods=c(mods,fm2.5)}
    } 
  }
  
  
  # Tmax ##################################################################################
  try((fm3=colext(psiformula=~1,
                  gammaformula=~Tmax,
                  epsilonformula=~Tmax,
                  pformula=~Tmax,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3")) {
    if(CondNum(fm3)<5000){
      if(CondNum(fm3)>0){mods=c(mods,fm3)}
    } 
  }
  
  try((fm3.1=colext(psiformula=~1,
                    gammaformula=~Tmax,
                    epsilonformula=~1,
                    pformula=~Tmax,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.1")) {
    if(CondNum(fm3.1)<5000){
      if(CondNum(fm3.1)>0){mods=c(mods,fm3.1)}
    }
  }
  
  
  try((fm3.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmax,
                    pformula=~Tmax,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.2")) {
    if(CondNum(fm3.2)<5000){
      if(CondNum(fm3.2)>0){mods=c(mods,fm3.2)}
    } 
  }
  
  try((fm3.3=colext(psiformula=~Elevation,
                    gammaformula=~Tmax,
                    epsilonformula=~Tmax,
                    pformula=~Tmax,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.3")) {
    if(CondNum(fm3.3)<5000){
      if(CondNum(fm3.3)>0){mods=c(mods,fm3.3)}
    } 
  }
  
  try((fm3.4=colext(psiformula=~Elevation,
                    gammaformula=~Tmax,
                    epsilonformula=~1,
                    pformula=~Tmax,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.4")) {
    if(CondNum(fm3.4)<5000){
      if(CondNum(fm3.4)>0){mods=c(mods,fm3.4)}
    }
  }
  
  
  try((fm3.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmax,
                    pformula=~Tmax,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.5")) {
    if(CondNum(fm3.5)<5000){
      if(CondNum(fm3.5)>0){mods=c(mods,fm3.5)}
    } 
  }
  
  # Tvar #############################################################################
  
  try((fm5=colext(psiformula=~1,
                  gammaformula=~Tvar,
                  epsilonformula=~Tvar,
                  pformula=~Tvar,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5")) {
    if(CondNum(fm5)<5000){
      if(CondNum(fm5)>0){mods=c(mods,fm5)}
    } 
  }
  
  try((fm5.1=colext(psiformula=~1,
                    gammaformula=~Tvar,
                    epsilonformula=~1,
                    pformula=~Tvar,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.1")) {
    if(CondNum(fm5.1)<5000){
      if(CondNum(fm5.1)>0){mods=c(mods,fm5.1)}
    } 
  }
  
  try((fm5.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tvar,
                    pformula=~Tvar,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.2")) {
    if(CondNum(fm5.2)<5000){
      if(CondNum(fm5.2)>0){mods=c(mods,fm5.2)}
    } 
  }
  
  try((fm5.3=colext(psiformula=~Elevation,
                    gammaformula=~Tvar,
                    epsilonformula=~Tvar,
                    pformula=~Tvar,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.3")) {
    if(CondNum(fm5.3)<5000){
      if(CondNum(fm5.3)>0){mods=c(mods,fm5.3)}
    } 
  }
  
  try((fm5.4=colext(psiformula=~Elevation,
                    gammaformula=~Tvar,
                    epsilonformula=~1,
                    pformula=~Tvar,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.4")) {
    if(CondNum(fm5.4)<5000){
      if(CondNum(fm5.4)>0){mods=c(mods,fm5.4)}
    } 
  }
  
  try((fm5.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tvar,
                    pformula=~Tvar,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.5")) {
    if(CondNum(fm5.5)<5000){
      if(CondNum(fm5.5)>0){mods=c(mods,fm5.5)}
    } 
  } 
  
  
  
  # Biotic only ##################################################################################
  try((fm6=colext(psiformula=~1,
                  gammaformula=~Biotic,
                  epsilonformula=~Biotic,
                  pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6")) {
    if(CondNum(fm6)<5000){
      if(CondNum(fm6)>0){mods=c(mods,fm6)}
    }
  }
  
  try((fm6.1=colext(psiformula=~1,
                    gammaformula=~Biotic,
                    epsilonformula=~1,
                    pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.1")) {
    if(CondNum(fm6.1)<5000){
      if(CondNum(fm6.1)>0){mods=c(mods,fm6.1)}
    }
  }
  
  try((fm6.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic,
                    pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.2")) {
    if(CondNum(fm6.2)<5000){
      if(CondNum(fm6.2)>0){mods=c(mods,fm6.2)}
    }
  }
  
  try((fm6.3=colext(psiformula=~Elevation,
                    gammaformula=~Biotic,
                    epsilonformula=~Biotic,
                    pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.3")) {
    if(CondNum(fm6.3)<5000){
      if(CondNum(fm6.3)>0){mods=c(mods,fm6.3)}
    }
  }
  
  try((fm6.4=colext(psiformula=~Elevation,
                    gammaformula=~Biotic,
                    epsilonformula=~1,
                    pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.4")) {
    if(CondNum(fm6.4)<5000){
      if(CondNum(fm6.4)>0){mods=c(mods,fm6.4)}
    }
  }
  
  try((fm6.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Biotic,
                    pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.5")) {
    if(CondNum(fm6.5)<5000){
      if(CondNum(fm6.5)>0){mods=c(mods,fm6.5)}
    }
  } 
  
  
  
  # Biotic + Tmin ##########################################################################
  try((fm7=colext(psiformula=~1,
                  gammaformula=~Biotic + Tmin,
                  epsilonformula=~Biotic + Tmin,
                  pformula=~Biotic + Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7")) {
    if(CondNum(fm7)<5000){
      if(CondNum(fm7)>0){mods=c(mods,fm7)}
    } 
  }
  
  try((fm7.1=colext(psiformula=~1,
                    gammaformula=~Biotic + Tmin,
                    epsilonformula=~1,
                    pformula=~Biotic + Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.1")) {
    if(CondNum(fm7.1)<5000){
      if(CondNum(fm7.1)>0){mods=c(mods,fm7.1)}
    } 
  }
  
  try((fm7.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~Biotic + Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.2")) {
    if(CondNum(fm7.2)<5000){
      if(CondNum(fm7.2)>0){mods=c(mods,fm7.2)}
    } 
  }
  
  try((fm7.3=colext(psiformula=~Elevation,
                    gammaformula=~Biotic + Tmin,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~Biotic + Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.3")) {
    if(CondNum(fm7.3)<5000){
      if(CondNum(fm7.3)>0){mods=c(mods,fm7.3)}
    } 
  }
  
  try((fm7.4=colext(psiformula=~Elevation,
                    gammaformula=~Biotic + Tmin,
                    epsilonformula=~1,
                    pformula=~Biotic + Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.4")) {
    if(CondNum(fm7.4)<5000){
      if(CondNum(fm7.4)>0){mods=c(mods,fm7.4)}
    } 
  }
  
  try((fm7.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~Biotic + Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.5")) {
    if(CondNum(fm7.5)<5000){
      if(CondNum(fm7.5)>0){mods=c(mods,fm7.5)}
    } 
  }
  
  
  #Biotic * Tmin ##########################################################################
  try((fm8=colext(psiformula=~1,
                  gammaformula=~Biotic * Tmin,
                  epsilonformula=~Biotic * Tmin,
                  pformula=~Biotic * Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8")) {
    if(CondNum(fm8)<5000){
      if(CondNum(fm8)>0){mods=c(mods,fm8)}
    } 
  }
  
  try((fm8.1=colext(psiformula=~1,
                    gammaformula=~Biotic * Tmin,
                    epsilonformula=~1,
                    pformula=~Biotic * Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.1")) {
    if(CondNum(fm8.1)<5000){
      if(CondNum(fm8.1)>0){mods=c(mods,fm8.1)}
    } 
  }
  
  try((fm8.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic * Tmin,
                    pformula=~Biotic * Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.2")) {
    if(CondNum(fm8.2)<5000){
      if(CondNum(fm8.2)>0){mods=c(mods,fm8.2)}
    } 
  }
  
  try((fm8.3=colext(psiformula=~Elevation,
                    gammaformula=~Biotic * Tmin,
                    epsilonformula=~Biotic * Tmin,
                    pformula=~Biotic * Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.3")) {
    if(CondNum(fm8.3)<5000){
      if(CondNum(fm8.3)>0){mods=c(mods,fm8.3)}
    } 
  }
  
  try((fm8.4=colext(psiformula=~Elevation,
                    gammaformula=~Biotic * Tmin,
                    epsilonformula=~1,
                    pformula=~Biotic * Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.4")) {
    if(CondNum(fm8.4)<5000){
      if(CondNum(fm8.4)>0){mods=c(mods,fm8.4)}
    } 
  }
  
  try((fm8.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Biotic * Tmin,
                    pformula=~Biotic * Tmin,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.5")) {
    if(CondNum(fm8.5)<5000){
      if(CondNum(fm8.5)>0){mods=c(mods,fm8.5)}
    } 
  }
  
  
  # Tmax + Biotic ##########################################################################
  try((fm9=colext(psiformula=~1,
                  gammaformula=~Tmax + Biotic,
                  epsilonformula=~Tmax + Biotic,
                  pformula=~Tmax + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9")) {
    if(CondNum(fm9)<5000){
      if(CondNum(fm9)>0){mods=c(mods,fm9)}
    } 
  }
  
  try((fm9.1=colext(psiformula=~1,
                    gammaformula=~Tmax + Biotic,
                    epsilonformula=~1,
                    pformula=~Tmax + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.1")){
    if(CondNum(fm9.1)<5000){
      if(CondNum(fm9.1)>0){mods=c(mods,fm9.1)}
    } 
  }
  
  try((fm9.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmax + Biotic,
                    pformula=~Tmax + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.2")) {
    if(CondNum(fm9.2)<5000){
      if(CondNum(fm9.2)>0){mods=c(mods,fm9.2)}
    } 
  }
  
  try((fm9.3=colext(psiformula=~Elevation,
                    gammaformula=~Tmax + Biotic,
                    epsilonformula=~Tmax + Biotic,
                    pformula=~Tmax + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.3")) {
    if(CondNum(fm9.3)<5000){
      if(CondNum(fm9.3)>0){mods=c(mods,fm9.3)}
    } 
  }
  
  try((fm9.4=colext(psiformula=~Elevation,
                    gammaformula=~Tmax + Biotic,
                    epsilonformula=~1,
                    pformula=~Tmax + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.4")){
    if(CondNum(fm9.4)<5000){
      if(CondNum(fm9.4)>0){mods=c(mods,fm9.4)}
    } 
  }
  
  try((fm9.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmax + Biotic,
                    pformula=~Tmax + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.5")) {
    if(CondNum(fm9.5)<5000){
      if(CondNum(fm9.5)>0){mods=c(mods,fm9.5)}
    } 
  }
  
  
  
  #Tmax * Biotic ##########################################################################
  try((fm10=colext(psiformula=~1,
                   gammaformula=~Tmax * Biotic,
                   epsilonformula=~Tmax * Biotic,
                   pformula=~Tmax * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10")) {
    if(CondNum(fm10)<5000){
      if(CondNum(fm10)>0){mods=c(mods,fm10)}
    } 
  }
  
  try((fm10.1=colext(psiformula=~1,
                     gammaformula=~Tmax * Biotic,
                     epsilonformula=~1,
                     pformula=~Tmax * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.1")) {
    if(CondNum(fm10.1)<5000){
      if(CondNum(fm10.1)>0){mods=c(mods,fm10.1)}
    } 
  }
  
  try((fm10.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmax * Biotic,
                     pformula=~Tmax * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.2")) {
    if(CondNum(fm10.2)<5000){
      if(CondNum(fm10.2)>0){mods=c(mods,fm10.2)}
    } 
  }
  
  try((fm10.3=colext(psiformula=~Elevation,
                     gammaformula=~Tmax * Biotic,
                     epsilonformula=~Tmax * Biotic,
                     pformula=~Tmax * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.3")) {
    if(CondNum(fm10.3)<5000){
      if(CondNum(fm10.3)>0){mods=c(mods,fm10.3)}
    } 
  }
  
  try((fm10.4=colext(psiformula=~Elevation,
                     gammaformula=~Tmax * Biotic,
                     epsilonformula=~1,
                     pformula=~Tmax * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.4")) {
    if(CondNum(fm10.4)<5000){
      if(CondNum(fm10.4)>0){mods=c(mods,fm10.4)}
    } 
  }
  
  try((fm10.5=colext(psiformula=~Elevation,
                     gammaformula=~1,
                     epsilonformula=~Tmax * Biotic,
                     pformula=~Tmax * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.5")) {
    if(CondNum(fm10.5)<5000){
      if(CondNum(fm10.5)>0){mods=c(mods,fm10.5)}
    } 
  }
  
  
  # Tvar * Biotic ##########################################################################
  try((fm12=colext(psiformula=~1,
                   gammaformula=~Tvar * Biotic,
                   epsilonformula=~Tvar * Biotic,
                   pformula=~Tvar * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12")) {
    if(CondNum(fm12)<5000){
      if(CondNum(fm12)>0){mods=c(mods,fm12)}
    } 
  }
  
  try((fm12.1=colext(psiformula=~1,
                     gammaformula=~Tvar * Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.1")) {
    if(CondNum(fm12.1)<5000){
      if(CondNum(fm12.1)>0){mods=c(mods,fm12.1)}
    } 
  }
  
  try((fm12.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar * Biotic,
                     pformula=~Tvar * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.2")) {
    if(CondNum(fm12.2)<5000){
      if(CondNum(fm12.2)>0){mods=c(mods,fm12.2)}
    } 
  }
  
  
  try((fm12.3=colext(psiformula=~Elevation,
                     gammaformula=~Tvar * Biotic,
                     epsilonformula=~Tvar * Biotic,
                     pformula=~Tvar * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.3")) {
    if(CondNum(fm12.3)<5000){
      if(CondNum(fm12.3)>0){mods=c(mods,fm12.3)}
    } 
  }
  
  try((fm12.4=colext(psiformula=~Elevation,
                     gammaformula=~Tvar * Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.4")) {
    if(CondNum(fm12.4)<5000){
      if(CondNum(fm12.4)>0){mods=c(mods,fm12.4)}
    } 
  }
  
  try((fm12.5=colext(psiformula=~Elevation,
                     gammaformula=~1,
                     epsilonformula=~Tvar * Biotic,
                     pformula=~Tvar * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.5")) {
    if(CondNum(fm12.5)<5000){
      if(CondNum(fm12.5)>0){mods=c(mods,fm12.5)}
    } 
  }
  
  
  #Tvar + Biotic ##########################################################################
  try((fm13=colext(psiformula=~1,
                   gammaformula=~Tvar + Biotic,
                   epsilonformula=~Tvar + Biotic,
                   pformula=~Tvar + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13")) {
    if(CondNum(fm13)<5000){
      if(CondNum(fm13)>0){mods=c(mods,fm13)}
    } 
  }
  
  try((fm13.1=colext(psiformula=~1,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.1")) {
    if(CondNum(fm13.1)<5000){
      if(CondNum(fm13.1)>0){mods=c(mods,fm13.1)}
    } 
  }
  
  try((fm13.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~Tvar + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.2")) {
    if(CondNum(fm13.2)<5000){
      if(CondNum(fm13.2)>0){mods=c(mods,fm13.2)}
    } 
  }
  
  try((fm13.3=colext(psiformula=~Elevation,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~Tvar + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.3")) {
    if(CondNum(fm13.3)<5000){
      if(CondNum(fm13.3)>0){mods=c(mods,fm13.3)}
    } 
  }
  
  try((fm13.4=colext(psiformula=~Elevation,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.4")) {
    if(CondNum(fm13.4)<5000){
      if(CondNum(fm13.4)>0){mods=c(mods,fm13.4)}
    } 
  }
  
  try((fm13.5=colext(psiformula=~Elevation,
                     gammaformula=~1,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~Tvar + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.5")) {
    if(CondNum(fm13.5)<5000){
      if(CondNum(fm13.5)>0){mods=c(mods,fm13.5)}
    } 
  }
  
  
  
  ############# INCLUDE QUADRATICS #####################
  
  # Tmin^2 only ##################################################################################
  try((fm22=colext(psiformula=~1,
                   gammaformula=~Tmin^2,
                   epsilonformula=~Tmin^2,
                   pformula=~Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22")) {
    if(CondNum(fm22)<5000){
      if(CondNum(fm22)>0){mods=c(mods,fm22)}
    } 
  }
  
  try((fm22.1=colext(psiformula=~1,
                     gammaformula=~Tmin^2,
                     epsilonformula=~1,
                     pformula=~Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.1")) {
    if(CondNum(fm22.1)<5000){
      if(CondNum(fm22.1)>0){mods=c(mods,fm22.1)}
    } 
  }
  
  try((fm22.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmin^2,
                     pformula=~Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.2")) {
    if(CondNum(fm22.2)<5000){
      if(CondNum(fm22.2)>0){mods=c(mods,fm22.2)}
    } 
  }
  
  try((fm22.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Tmin^2,
                     epsilonformula=~Tmin^2,
                     pformula=~Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.3")) {
    if(CondNum(fm22.3)<5000){
      if(CondNum(fm22.3)>0){mods=c(mods,fm22.3)}
    } 
  }
  
  try((fm22.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Tmin^2,
                     epsilonformula=~1,
                     pformula=~Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.4")) {
    if(CondNum(fm22.4)<5000){
      if(CondNum(fm22.4)>0){mods=c(mods,fm22.4)}
    } 
  }
  
  try((fm22.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tmin^2,
                     pformula=~Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.5")) {
    if(CondNum(fm22.5)<5000){
      if(CondNum(fm22.5)>0){mods=c(mods,fm22.5)}
    } 
  }
  
  
  
  # Tmax^2 ##################################################################################
  try((fm23=colext(psiformula=~1,
                   gammaformula=~Tmax^2,
                   epsilonformula=~Tmax^2,
                   pformula=~Tmax^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm23")) {
    if(CondNum(fm23)<5000){
      if(CondNum(fm23)>0){mods=c(mods,fm23)}
    } 
  }
  
  try((fm23.1=colext(psiformula=~1,
                     gammaformula=~Tmax^2,
                     epsilonformula=~1,
                     pformula=~Tmax^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm23.1")) {
    if(CondNum(fm23.1)<5000){
      if(CondNum(fm23.1)>0){mods=c(mods,fm23.1)}
    }
  }
  
  
  try((fm23.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmax^2,
                     pformula=~Tmax^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm23.2")) {
    if(CondNum(fm23.2)<5000){
      if(CondNum(fm23.2)>0){mods=c(mods,fm23.2)}
    } 
  }
  
  
  try((fm23.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Tmax^2,
                     epsilonformula=~Tmax^2,
                     pformula=~Tmax^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm23.3")) {
    if(CondNum(fm23.3)<5000){
      if(CondNum(fm23.3)>0){mods=c(mods,fm23.3)}
    } 
  }
  
  try((fm23.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Tmax^2,
                     epsilonformula=~1,
                     pformula=~Tmax^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm23.4")) {
    if(CondNum(fm23.4)<5000){
      if(CondNum(fm23.4)>0){mods=c(mods,fm23.4)}
    }
  }
  
  
  try((fm23.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tmax^2,
                     pformula=~Tmax^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm23.5")) {
    if(CondNum(fm23.5)<5000){
      if(CondNum(fm23.5)>0){mods=c(mods,fm23.5)}
    } 
  }
  
  
  
  # Tvar^2 #############################################################################
  
  try((fm25=colext(psiformula=~1,
                   gammaformula=~Tvar^2,
                   epsilonformula=~Tvar^2,
                   pformula=~Tvar^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm25")) {
    if(CondNum(fm25)<5000){
      if(CondNum(fm25)>0){mods=c(mods,fm25)}
    } 
  }
  
  try((fm25.1=colext(psiformula=~1,
                     gammaformula=~Tvar^2,
                     epsilonformula=~1,
                     pformula=~Tvar^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm25.1")) {
    if(CondNum(fm25.1)<5000){
      if(CondNum(fm25.1)>0){mods=c(mods,fm25.1)}
    } 
  }
  
  try((fm25.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar^2,
                     pformula=~Tvar^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm25.2")) {
    if(CondNum(fm25.2)<5000){
      if(CondNum(fm25.2)>0){mods=c(mods,fm25.2)}
    } 
  }
  
  try((fm25.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Tvar^2,
                     epsilonformula=~Tvar^2,
                     pformula=~Tvar^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm25.3")) {
    if(CondNum(fm25.3)<5000){
      if(CondNum(fm25.3)>0){mods=c(mods,fm25.3)}
    } 
  }
  
  try((fm25.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Tvar^2,
                     epsilonformula=~1,
                     pformula=~Tvar^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm25.4")) {
    if(CondNum(fm25.4)<5000){
      if(CondNum(fm25.4)>0){mods=c(mods,fm25.4)}
    } 
  }
  
  try((fm25.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tvar^2,
                     pformula=~Tvar^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm25.5")) {
    if(CondNum(fm25.5)<5000){
      if(CondNum(fm25.5)>0){mods=c(mods,fm25.5)}
    } 
  }
  
  
  # Biotic #############################################################################
  
  try((fm26.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Biotic,
                     epsilonformula=~Biotic,
                     pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm26.3")) {
    if(CondNum(fm26.3)<5000){
      if(CondNum(fm26.3)>0){mods=c(mods,fm26.3)}
    } 
  }
  
  try((fm26.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Biotic,
                     epsilonformula=~1,
                     pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm26.4")) {
    if(CondNum(fm26.4)<5000){
      if(CondNum(fm26.4)>0){mods=c(mods,fm26.4)}
    } 
  }
  
  try((fm26.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Biotic,
                     pformula=~Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm26.5")) {
    if(CondNum(fm26.5)<5000){
      if(CondNum(fm26.5)>0){mods=c(mods,fm26.5)}
    } 
  }
  
  
  
  # Biotic + Tmin^2 ##########################################################################
  try((fm27=colext(psiformula=~1,
                   gammaformula=~Biotic + Tmin^2,
                   epsilonformula=~Biotic + Tmin^2,
                   pformula=~Biotic + Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm27")) {
    if(CondNum(fm27)<5000){
      if(CondNum(fm27)>0){mods=c(mods,fm27)}
    } 
  }
  
  try((fm27.1=colext(psiformula=~1,
                     gammaformula=~Biotic + Tmin^2,
                     epsilonformula=~1,
                     pformula=~Biotic + Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm27.1")) {
    if(CondNum(fm27.1)<5000){
      if(CondNum(fm27.1)>0){mods=c(mods,fm27.1)}
    } 
  }
  
  try((fm27.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Biotic + Tmin^2,
                     pformula=~Biotic + Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm27.2")) {
    if(CondNum(fm27.2)<5000){
      if(CondNum(fm27.2)>0){mods=c(mods,fm27.2)}
    } 
  }
  
  try((fm27.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Biotic + Tmin^2,
                     epsilonformula=~Biotic + Tmin^2,
                     pformula=~Biotic + Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm27.3")) {
    if(CondNum(fm27.3)<5000){
      if(CondNum(fm27.3)>0){mods=c(mods,fm27.3)}
    } 
  }
  
  try((fm27.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Biotic + Tmin^2,
                     epsilonformula=~1,
                     pformula=~Biotic + Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm27.4")) {
    if(CondNum(fm27.4)<5000){
      if(CondNum(fm27.4)>0){mods=c(mods,fm27.4)}
    } 
  }
  
  try((fm27.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Biotic + Tmin^2,
                     pformula=~Biotic + Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm27.5")) {
    if(CondNum(fm27.5)<5000){
      if(CondNum(fm27.5)>0){mods=c(mods,fm27.5)}
    } 
  }
  
  
#  #Biotic * Tmin^2 ##########################################################################
  try((fm28=colext(psiformula=~1,
                   gammaformula=~Biotic * Tmin^2,
                   epsilonformula=~Biotic * Tmin^2,
                   pformula=~Biotic * Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm28")) {
    if(CondNum(fm28)<5000){
      if(CondNum(fm28)>0){mods=c(mods,fm28)}
    } 
  }
  
  try((fm28.1=colext(psiformula=~1,
                     gammaformula=~Biotic * Tmin^2,
                     epsilonformula=~1,
                     pformula=~Biotic * Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm28.1")) {
    if(CondNum(fm28.1)<5000){
      if(CondNum(fm28.1)>0){mods=c(mods,fm28.1)}
    } 
  }
  
  try((fm28.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Biotic * Tmin^2,
                     pformula=~Biotic * Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm28.2")) {
    if(CondNum(fm28.2)<5000){
      if(CondNum(fm28.2)>0){mods=c(mods,fm28.2)}
    } 
  }
  
  try((fm28.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Biotic * Tmin^2,
                     epsilonformula=~Biotic * Tmin^2,
                     pformula=~Biotic * Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm28.3")) {
    if(CondNum(fm28.3)<5000){
      if(CondNum(fm28.3)>0){mods=c(mods,fm28.3)}
    } 
  }
  
  try((fm28.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Biotic * Tmin^2,
                     epsilonformula=~1,
                     pformula=~Biotic * Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm28.4")) {
    if(CondNum(fm28.4)<5000){
      if(CondNum(fm28.4)>0){mods=c(mods,fm28.4)}
    } 
  }
  
  try((fm28.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Biotic * Tmin^2,
                     pformula=~Biotic * Tmin^2,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm28.5")) {
    if(CondNum(fm28.5)<5000){
      if(CondNum(fm28.5)>0){mods=c(mods,fm28.5)}
    } 
  }
  
  
  
  # Tmax^2 + Biotic ##########################################################################
  try((fm29=colext(psiformula=~1,
                   gammaformula=~Tmax^2 + Biotic,
                   epsilonformula=~Tmax^2 + Biotic,
                   pformula=~Tmax^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm29")) {
    if(CondNum(fm29)<5000){
      if(CondNum(fm29)>0){mods=c(mods,fm29)}
    } 
  }
  
  try((fm29.1=colext(psiformula=~1,
                     gammaformula=~Tmax^2 + Biotic,
                     epsilonformula=~1,
                     pformula=~Tmax^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm29.1")){
    if(CondNum(fm29.1)<5000){
      if(CondNum(fm29.1)>0){mods=c(mods,fm29.1)}
    } 
  }
  
  try((fm29.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmax^2 + Biotic,
                     pformula=~Tmax^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm29.2")) {
    if(CondNum(fm29.2)<5000){
      if(CondNum(fm29.2)>0){mods=c(mods,fm29.2)}
    } 
  }
  
  
  try((fm29.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Tmax^2 + Biotic,
                     epsilonformula=~Tmax^2 + Biotic,
                     pformula=~Tmax^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm29.3")) {
    if(CondNum(fm29.3)<5000){
      if(CondNum(fm29.3)>0){mods=c(mods,fm29)}
    } 
  }
  
  try((fm29.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Tmax^2 + Biotic,
                     epsilonformula=~1,
                     pformula=~Tmax^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm29.4")){
    if(CondNum(fm29.4)<5000){
      if(CondNum(fm29.4)>0){mods=c(mods,fm29.4)}
    } 
  }
  
  try((fm29.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tmax^2 + Biotic,
                     pformula=~Tmax^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm29.5")) {
    if(CondNum(fm29.5)<5000){
      if(CondNum(fm29.5)>0){mods=c(mods,fm29.5)}
    } 
  }
  
  
  
  #Tmax^2 * Biotic ##########################################################################
  try((fm30=colext(psiformula=~1,
                   gammaformula=~Tmax^2 * Biotic,
                   epsilonformula=~Tmax^2 * Biotic,
                   pformula=~Tmax^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("30")) {
    if(CondNum(30)<5000){
      if(CondNum(30)>0){mods=c(mods,30)}
    } 
  }
  
  try((30.1=colext(psiformula=~1,
                   gammaformula=~Tmax^2 * Biotic,
                   epsilonformula=~1,
                   pformula=~Tmax^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("30.1")) {
    if(CondNum(30.1)<5000){
      if(CondNum(fm30.1)>0){mods=c(mods,fm30.1)}
    } 
  }
  
  try((fm30.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmax^2 * Biotic,
                     pformula=~Tmax^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm30.2")) {
    if(CondNum(fm30.2)<5000){
      if(CondNum(fm30.2)>0){mods=c(mods,fm30.2)}
    } 
  }
  
  
  try((fm30=colext(psiformula=~Elevation^2,
                   gammaformula=~Tmax^2 * Biotic,
                   epsilonformula=~Tmax^2 * Biotic,
                   pformula=~Tmax^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("30.3")) {
    if(CondNum(30.3)<5000){
      if(CondNum(30.3)>0){mods=c(mods,30.3)}
    } 
  }
  
  try((30.4=colext(psiformula=~Elevation^2,
                   gammaformula=~Tmax^2 * Biotic,
                   epsilonformula=~1,
                   pformula=~Tmax^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("30.4")) {
    if(CondNum(30.4)<5000){
      if(CondNum(fm30.4)>0){mods=c(mods,fm30.4)}
    } 
  }
  
  try((fm30.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tmax^2 * Biotic,
                     pformula=~Tmax^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm30.5")) {
    if(CondNum(fm30.5)<5000){
      if(CondNum(fm30.5)>0){mods=c(mods,fm30.5)}
    } 
  }
  
  
  # Tvar^2 * Biotic ##########################################################################
  try((fm32=colext(psiformula=~1,
                   gammaformula=~Tvar^2 * Biotic,
                   epsilonformula=~Tvar^2 * Biotic,
                   pformula=~Tvar^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32")) {
    if(CondNum(fm32)<5000){
      if(CondNum(fm32)>0){mods=c(mods,fm32)}
    } 
  }
  
  try((fm32.1=colext(psiformula=~1,
                     gammaformula=~Tvar^2 * Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.1")) {
    if(CondNum(fm32.1)<5000){
      if(CondNum(fm32.1)>0){mods=c(mods,fm32.1)}
    } 
  }
  
  try((fm32.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar^2 * Biotic,
                     pformula=~Tvar^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.2")) {
    if(CondNum(fm32.2)<5000){
      if(CondNum(fm32.2)>0){mods=c(mods,fm32.2)}
    } 
  }
  
  
  try((fm32.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Tvar^2 * Biotic,
                     epsilonformula=~Tvar^2 * Biotic,
                     pformula=~Tvar^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.3")) {
    if(CondNum(fm32.3)<5000){
      if(CondNum(fm32.3)>0){mods=c(mods,fm32.3)}
    } 
  }
  
  try((fm32.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Tvar^2 * Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.4")) {
    if(CondNum(fm32.4)<5000){
      if(CondNum(fm32.4)>0){mods=c(mods,fm32.4)}
    } 
  }
  
  try((fm32.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tvar^2 * Biotic,
                     pformula=~Tvar^2 * Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.5")) {
    if(CondNum(fm32.5)<5000){
      if(CondNum(fm32.5)>0){mods=c(mods,fm32.5)}
    } 
  }
  
  
  #Tvar^2 + Biotic ##########################################################################
  try((fm33=colext(psiformula=~1,
                   gammaformula=~Tvar^2 + Biotic,
                   epsilonformula=~Tvar^2 + Biotic,
                   pformula=~Tvar^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm33")) {
    if(CondNum(fm33)<5000){
      if(CondNum(fm33)>0){mods=c(mods,fm33)}
    } 
  }
  
  try((fm33.1=colext(psiformula=~1,
                     gammaformula=~Tvar^2 + Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm33.1")) {
    if(CondNum(fm33.1)<5000){
      if(CondNum(fm33.1)>0){mods=c(mods,fm33.1)}
    } 
  }
  
  try((fm33.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar^2 + Biotic,
                     pformula=~Tvar^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm33.2")) {
    if(CondNum(fm33.2)<5000){
      if(CondNum(fm33.2)>0){mods=c(mods,fm33.2)}
    } 
  }
  
  
  try((fm33.3=colext(psiformula=~Elevation^2,
                     gammaformula=~Tvar^2 + Biotic,
                     epsilonformula=~Tvar^2 + Biotic,
                     pformula=~Tvar^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm33.3")) {
    if(CondNum(fm33.3)<5000){
      if(CondNum(fm33.3)>0){mods=c(mods,fm33.3)}
    } 
  }
  
  try((fm33.4=colext(psiformula=~Elevation^2,
                     gammaformula=~Tvar^2 + Biotic,
                     epsilonformula=~1,
                     pformula=~Tvar^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm33.4")) {
    if(CondNum(fm33.4)<5000){
      if(CondNum(fm33.4)>0){mods=c(mods,fm33.4)}
    } 
  }
  
  try((fm33.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~Tvar^2 + Biotic,
                     pformula=~Tvar^2 + Biotic,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm33.5")) {
    if(CondNum(fm33.5)<5000){
      if(CondNum(fm33.5)>0){mods=c(mods,fm33.5)}
    } 
  }
  


 #ForestLossCT ##########################################################################
  try((fm34=colext(psiformula=~1,
                   gammaformula=~ForestLossCT,
                   epsilonformula=~ForestLossCT,
                   pformula=~ForestLossCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm34")) {
    if(CondNum(fm34)<5000){
      if(CondNum(fm34)>0){mods=c(mods,fm34)}
    } 
  }
  
  try((fm34.1=colext(psiformula=~1,
                     gammaformula=~ForestLossCT,
                     epsilonformula=~1,
                     pformula=~ForestLossCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm34.1")) {
    if(CondNum(fm34.1)<5000){
      if(CondNum(fm34.1)>0){mods=c(mods,fm34.1)}
    } 
  }
  
  try((fm34.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~ForestLossCT,
                     pformula=~ForestLossCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm34.2")) {
    if(CondNum(fm34.2)<5000){
      if(CondNum(fm34.2)>0){mods=c(mods,fm34.2)}
    } 
  }
  
  
  try((fm34.3=colext(psiformula=~Elevation^2,
                     gammaformula=~ForestLossCT,
                     epsilonformula=~ForestLossCT,
                     pformula=~ForestLossCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm34.3")) {
    if(CondNum(fm34.3)<5000){
      if(CondNum(fm34.3)>0){mods=c(mods,fm34.3)}
    } 
  }
  
  try((fm34.4=colext(psiformula=~Elevation^2,
                     gammaformula=~ForestLossCT,
                     epsilonformula=~1,
                     pformula=~ForestLossCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm34.4")) {
    if(CondNum(fm34.4)<5000){
      if(CondNum(fm34.4)>0){mods=c(mods,fm34.4)}
    } 
  }
  
  try((fm34.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~ForestLossCT,
                     pformula=~ForestLossCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm34.5")) {
    if(CondNum(fm34.5)<5000){
      if(CondNum(fm34.5)>0){mods=c(mods,fm34.5)}
    } 
  }
 




 #ForestGainCT ##########################################################################
  try((fm35=colext(psiformula=~1,
                   gammaformula=~ForestGainCT,
                   epsilonformula=~ForestGainCT,
                   pformula=~ForestGainCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm35")) {
    if(CondNum(fm35)<5000){
      if(CondNum(fm35)>0){mods=c(mods,fm35)}
    } 
  }
  
  try((fm35.1=colext(psiformula=~1,
                     gammaformula=~ForestGainCT,
                     epsilonformula=~1,
                     pformula=~ForestGainCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm35.1")) {
    if(CondNum(fm35.1)<5000){
      if(CondNum(fm35.1)>0){mods=c(mods,fm35.1)}
    } 
  }
  
  try((fm35.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~ForestGainCT,
                     pformula=~ForestGainCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm35.2")) {
    if(CondNum(fm35.2)<5000){
      if(CondNum(fm35.2)>0){mods=c(mods,fm35.2)}
    } 
  }
  
  
  try((fm35.3=colext(psiformula=~Elevation^2,
                     gammaformula=~ForestGainCT,
                     epsilonformula=~ForestGainCT,
                     pformula=~ForestGainCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm35.3")) {
    if(CondNum(fm35.3)<5000){
      if(CondNum(fm35.3)>0){mods=c(mods,fm35.3)}
    } 
  }
  
  try((fm35.4=colext(psiformula=~Elevation^2,
                     gammaformula=~ForestGainCT,
                     epsilonformula=~1,
                     pformula=~ForestGainCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm35.4")) {
    if(CondNum(fm35.4)<5000){
      if(CondNum(fm35.4)>0){mods=c(mods,fm35.4)}
    } 
  }
  
  try((fm35.5=colext(psiformula=~Elevation^2,
                     gammaformula=~1,
                     epsilonformula=~ForestGainCT,
                     pformula=~ForestGainCT,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm35.5")) {
    if(CondNum(fm35.5)<5000){
      if(CondNum(fm35.5)>0){mods=c(mods,fm35.5)}
    } 
  }


  
  ######################################
  # Run Model Selection
  ######################################
  
  # If zero models for a species converged:
  
  if(isEmpty(mods)==TRUE){
    results.all[[k]]=NA  
    results.table.ma[[k]]=NA
    names(results.table.ma)[k] <- nms[k]
    results.table.aic[[k]]=NA
    names(results.table.aic)[k] <- nms[k]
    colext.transformed[[k]]=NA
    names(colext.transformed)[k] <- nms[k]
  }else{
    models=fitList(fits=mods)
    (ms <- modSel(models))
    results.all[[k]]=ms  	
    
    
    # For species with models that converged, run model selection:
    
    models=fitList(fits=mods)
    
    (ms <- modSel(models))
    
    results.all[[k]]=ms
    mods.all[[k]]=mods
    
    
    ######################################
    # Extract Results
    ######################################
    
    toExport <- as(ms,"data.frame")
    
    null0.aic <- toExport$delta[toExport$formula=="~1 ~ 1 ~ 1 ~ 1"]
    
    #if null didn't converge
    if(isEmpty(null0.aic)==TRUE){
      nulls <- NA
    }else{
      nulls <- rbind(toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",])
                        }
    
    if((null0.aic==0) || (isEmpty(null0.aic)==TRUE)){
      results.table.ma[[k]] <- rbind(nulls)
      results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
      names(results.table.ma)[k] <- nms[k]
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
      names(results.table.aic)[k] <- nms[k]
      
    }else{
      results.table.ma[[k]] <- rbind(toExport[1,], nulls)
      results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
      names(results.table.ma)[k] <- nms[k]
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
      results.table.aic[[k]] <- cbind(nms[k], results.table.aic[[k]])
      names(results.table.aic)[k] <- nms[k]
    }
  }    
  
  # Remove all models and results
  
  rm(fm0, 
     fm2,fm2.1,fm2.2,fm2.3,fm2.4,fm2.5,
     fm3,fm3.1,fm3.2,fm3.3,fm3.4,fm3.5,
     fm5,fm5.1,fm5.2,fm5.3,fm5.4,fm5.5,
     fm6,fm6.1,fm6.2,fm6.3,fm6.4,fm6.5,
     fm7,fm7.1,fm7.2,fm7.3,fm7.4,fm7.5,
     fm8,fm8.1,fm8.2,fm8.3,fm8.4,fm8.5,
     fm9,fm9.1,fm9.2,fm9.3,fm9.4,fm9.5,
     fm10,fm10.1,fm10.2,fm10.3,fm10.4,fm10.5,
     fm12,fm12.1,fm12.2,fm12.3,fm12.4,fm12.5,
     fm13,fm13.1,fm13.2,fm13.3,fm13.4,fm13.5,
     fm22,fm22.1,fm22.2,fm22.3,fm22.4,fm22.5,
     fm23,fm23.1,fm23.2,fm23.3,fm23.4,fm23.5,
     fm25,fm25.1,fm25.2,fm25.3,fm25.4,fm25.5,
     fm26.3,fm26.4,fm26.5, 
     fm27,fm27.1,fm27.2,fm27.3,fm27.4,fm27.5,
     fm28,fm28.1,fm28.2,fm28.3,fm28.4,fm28.5,
     fm29,fm29.1,fm29.2,fm29.3,fm29.4,fm29.5,
     fm30,fm30.1,fm30.2,fm30.3,fm30.4,fm30.5,
     fm31,fm31.1,fm31.2,fm31.3,fm31.4,fm31.5,
     fm32,fm32.1,fm32.2,fm32.3,fm32.4,fm32.5,
     fm33,fm33.1,fm33.2,fm33.3,fm33.4,fm33.5,
     fm34,fm34.1,fm34.2,fm34.3,fm34.4,fm34.5,
     fm35,fm35.1,fm35.2,fm35.3,fm35.4,fm35.5,
     mods, ms, temp, toExport,
     null0.aic, nulls)
}

############################################################################
# End loop
############################################################################




######################################
# Write Results
######################################

# Coerce the lists of results into dataframes and write to files

results.table.ma.df <- ldply(results.table.ma, data.frame)
results.table.aic.df <- ldply(results.table.aic, data.frame)
colext.transformed.df <- ldply(colext.transformed, data.frame)

write.csv(results.table.ma.df, file="results.table.ma_Covariate.csv")
write.csv(results.table.aic.df, file="results.table.aic_Covariate.csv")
write.csv(colext.transformed.df, file="colext.transformed.csv")

for(i in 1:length(nms)) {
  outputname <- paste(nms[i], "colextAIC.covariates", "csv", sep=".")
  if(is.na(results.all[[i]])==TRUE){
    output <- NA
  }else{
    output <- results.all[[i]]@Full
  } 
  write.csv(output, file=outputname)
}