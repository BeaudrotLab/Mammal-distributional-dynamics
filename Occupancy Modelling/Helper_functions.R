# PURPOSE:FUNCTIONS FOR PLOTTING GRAPHS IN FIGURES for covariates of colonization and extinction IN GLOBAL ECOLOGY & BIOGEOGRAPHY SUBMISSION
# AUTHOR: Lydia Beaudrot
# Date: February 27, 2018

# Functions for extracting data for plots

create.umf <- function(k){
  load("Occupancy Modelling/All_covs.RData")
  load("Occupancy Modelling/Species7sites_Include.RData")
  load("Occupancy Modelling/BIOTIC_Include.RData")
  
  Species_data <- Species7sites_Include
  BIOTIC <- BIOTIC_Include
  
  nms=names(Species_data)
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
return(umf)
}

extract.Tmin <- function(k){
  load("Occupancy Modelling/All_covs.RData")
  load("Occupancy Modelling/Species7sites_Include.RData")
  load("Occupancy Modelling/BIOTIC_Include.RData")
  
  Species_data <- Species7sites_Include
  BIOTIC <- BIOTIC_Include
  
  nms=names(Species_data)
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
  
  return(Tmin)
  
}

extract.Tmax <- function(k){
  load("Occupancy Modelling/All_covs.RData")
  load("Occupancy Modelling/Species7sites_Include.RData")
  load("Occupancy Modelling/BIOTIC_Include.RData")
  
  Species_data <- Species7sites_Include
  BIOTIC <- BIOTIC_Include
  
  nms=names(Species_data)
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
  Tmax <- as.data.frame(sapply(covs, "[", 5))
  return(Tmax)
}

extract.Tvar <- function(k){
  load("Occupancy Modelling/All_covs.RData")
  load("Occupancy Modelling/Species7sites_Include.RData")
  load("Occupancy Modelling/BIOTIC_Include.RData")
  
  Species_data <- Species7sites_Include
  BIOTIC <- BIOTIC_Include
  
  nms=names(Species_data)
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
  Tvar <- as.data.frame(sapply(covs, "[", 6))
  return(Tvar)
}

extract.Biotic <- function(k){
  load("Occupancy Modelling/All_covs.RData")
  load("Occupancy Modelling/Species7sites_Include.RData")
  load("Occupancy Modelling/BIOTIC_Include.RData")
  
  Species_data <- Species7sites_Include
  BIOTIC <- BIOTIC_Include
  
  nms=names(Species_data)
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
  return(Biotic)
  
}



# Functions for plotting for temperature covariates of colonization and extinction

# Create translucent colors
redtrans <- rgb(238, 0, 0, 150, maxColorValue=255)
goldtrans <- rgb(255, 242, 0, 150, maxColorValue=255)
orangetrans <- rgb(255, 165, 0, 200, maxColorValue=255)
bluetrans <- rgb(30, 144, 255, 150, maxColorValue=255)
  
plot.colTmax <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="", xlab="Tmax", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=redtrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmax, type="line", lwd=2)
}

plot.colTmin <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, xlab="Tmin", ylab="", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=orangetrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

plot.colTvar <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="", xlab="Tvar", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=goldtrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}


plot.extTmax <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="", xlab="Tmax", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=redtrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmax, type="line", lwd=2)
}

plot.extTmin <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, ylab="", xlab="Tmin", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=orangetrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

plot.extTvar <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="", xlab="Tvar", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=goldtrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}

# Functions to plot biotic estimates + standard errors for colonization and extinction 

plot.colBiotic <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Biotic, pch=19, las=1, ylab="", xlab="Biotic", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Biotic, rev(outputdata$Biotic)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=bluetrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Biotic, type="line", lwd=2)
}

plot.extBiotic <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Biotic, pch=19, las=1, ylab="", xlab="Biotic", bty="n", xlim=c(-3,3), ylim=c(0,1), cex=0.8, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Biotic, rev(outputdata$Biotic)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=bluetrans, border=NA)
  lines(outputdata$Predicted ~ outputdata$Biotic, type="line", lwd=2)
}


# 16 colext top models used for 33 populations (excludes models with interactions)

f.c.Tmax       <-  function(umf){colext(psiformula=~1, gammaformula=~Tmax, epsilonformula=~1, pformula=~Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.cE.Tmax      <-  function(umf){colext(psiformula=~Elevation, gammaformula=~Tmax, epsilonformula=~1, pformula=~Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.c.Tmin       <-  function(umf){colext(psiformula=~1, gammaformula=~Tmin, epsilonformula=~1, pformula=~Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.cE.Tmin      <-  function(umf){colext(psiformula=~Elevation, gammaformula=~Tmin, epsilonformula=~1, pformula=~Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.c.Tvar       <-  function(umf){colext(psiformula=~1, gammaformula=~Tvar, epsilonformula=~1, pformula=~Tvar, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.c.Biotic     <-  function(umf){colext(psiformula=~1, gammaformula=~Biotic, epsilonformula=~1, pformula=~Biotic, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.c.BioticTmin <-  function(umf){colext(psiformula=~1, gammaformula=~Biotic + Tmin, epsilonformula=~1, pformula=~Biotic + Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.c.BioticTmax <-  function(umf){colext(psiformula=~1, gammaformula=~Biotic + Tmax, epsilonformula=~1, pformula=~Biotic + Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))}

f.cE.BioticTmax <- function(umf){colext(psiformula=~Elevation, gammaformula=~Biotic + Tmax, epsilonformula=~1, pformula=~Biotic + Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))}

f.ceE.BioticTmax <-function(umf){colext(psiformula=~Elevation, gammaformula=~Biotic + Tmax, epsilonformula=~Biotic + Tmax, pformula=~Biotic + Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))}    

f.e.Tmax       <-  function(umf){colext(psiformula=~1, gammaformula=~1, epsilonformula=~Tmax, pformula=~Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.e.Tmin       <-  function(umf){colext(psiformula=~1, gammaformula=~1, epsilonformula=~Tmin, pformula=~Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.e.Tvar       <-  function(umf){colext(psiformula=~1, gammaformula=~1, epsilonformula=~Tvar, pformula=~Tvar, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.eE.Biotic    <-  function(umf){colext(psiformula=~Elevation, gammaformula=~1, epsilonformula=~Biotic, pformula=~Biotic, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.e.BioticTmin <-  function(umf){colext(psiformula=~1, gammaformula=~1, epsilonformula=~Biotic + Tmin, pformula=~Biotic + Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
f.eE.BioticTmin <- function(umf){colext(psiformula=~Elevation, gammaformula=~1, epsilonformula=~Biotic + Tmin, pformula=~Biotic + Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))}
