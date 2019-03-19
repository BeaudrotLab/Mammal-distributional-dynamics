# Plots for species with interaction terms in top TB model without consideration of elevation.

# Population & top model run on 12-21-2017    Occ       Col               Ext                   Detection
# BIF.species.BIF.Cephalophus silvicultor	    ~1          ~ 1               ~ Biotic * Tmin       ~ Biotic * Tmin
# UDZ.species.UDZ.Cephalophus spadix	        ~1          ~ 1               ~ Biotic * Tmin       ~ Biotic * Tmin
# UDZ.species.UDZ.Cercocebus sanjei	          ~1          ~ Tvar * Biotic   ~ 1                   ~ Tvar * Biotic
# UDZ.species.UDZ.Genetta servalina	          ~Elevation  ~ Tmax * Biotic   ~ 1                   ~ Tmax * Biotic
# VB_.species.VB_.Mazama temama	              ~1          ~ Tmax * Biotic   ~ 1                   ~ Tmax * Biotic
# VB_.species.VB_.Dasypus novemcinctus	      ~1          ~ Tvar * Biotic   ~ Tvar * Biotic       ~ Tvar * Biotic

rm(list=ls())
library(unmarked)
library(fields)

# Functions for plotting interactions for covariates of colonization and detection

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

plot.colTmax <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="", xlab="Tmax", bty="n", xlim=c(-3,3), ylim=c(0,1), cex.axis=1.5, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmax, lty=1, lwd=1)
}

extra.colTmax <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmax, lty=2, lwd=1)
}

plot.colTvar <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="", xlab="Tvar", bty="n", xlim=c(-3,3), ylim=c(0,1), cex.axis=1.5, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tvar, lty=1, lwd=1)
}

extra.colTvar <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tvar, lty=2, lwd=1)
}


plot.extTvar <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="", xlab="Tvar", bty="n", xlim=c(-3,3), ylim=c(0,1), cex.axis=1.5, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tvar, lty=1, lwd=1)
}

extra.extTvar <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tvar, lty=2, lwd=1)
}


plot.colTmin <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, ylab="", xlab="Tmin", bty="n", xlim=c(-3,3), ylim=c(0,1), cex.axis=1.5, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmin, lty=1, lwd=1)
}

extra.colTmin <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
  polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmin, lty=2, lwd=1)
}

plot.extTmin <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, ylab="", xlab="Tmin", bty="n", xlim=c(-3,3), ylim=c(0,1), cex.axis=1.5, col="white", axes=FALSE)
  Axis(side = 1, labels=FALSE, lwd.ticks=0)
  Axis(side = 2, labels=FALSE, lwd.ticks=0)
  polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmin, lty=1, lwd=1)
}

extra.extTmin <- function(inputmodel, inputdata){
  outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
  polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
          c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
            ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255), border=NA)
  lines(outputdata$Predicted ~ outputdata$Tmin, lty=2, lwd=1)
}

# light green: col=rgb(0,200,100,75,maxColorValue=255)


load("Occupancy Modelling/All_covs.RData")
load("Occupancy Modelling/Species7sites_Include.RData")
load("Occupancy Modelling/BIOTIC_Include.RData")

Species_data <- Species7sites_Include
BIOTIC <- BIOTIC_Include

nms=names(Species_data)



############### BEGIN PLOTTING #######################
par(mar=c(2,1,1,1), mfrow=c(3,2), mgp=c(2,1,0))
#par(mar=c(3,1,0,1), mfrow=c(3,2))


# UDZ.species.UDZ.Cercocebus sanjei k=11
k=11
umf <- create.umf(k)
Tvar <- extract.Tvar(k)
fm.c.BTvar <- colext(psiformula=~1, gammaformula=~Biotic * Tvar, epsilonformula=~1, pformula=~Biotic * Tvar, data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm.c.BTvar
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTvar(inputmodel, inputdata)
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTvar(inputmodel, inputdata)
mtext("Cercocebus sanjei (UDZ)", side=3, line=0, cex=0.5)



# BIF.Cephalophus silvicultor k=23
k=23
umf <- create.umf(k)
Tmin <- extract.Tmin(k)
fm.e.BTmin <- colext(psiformula=~1, gammaformula=~1, epsilonformula=~Biotic * Tmin, pformula=~Biotic * Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm.e.BTmin
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.extTmin(inputmodel, inputdata)
mtext("Cephalophus silvicultor (BIF)", side=3, line=0, cex=0.5, font=3)

# UDZ.species.UDZ.Genetta servalina k=14
k=14
umf <- create.umf(k)
Tmax <- extract.Tmax(k)
fm.c.BTmax.elev <- colext(psiformula=~Elevation, gammaformula=~Biotic * Tmax, epsilonformula=~1, pformula=~Biotic * Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm.c.BTmax.elev
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTmax(inputmodel, inputdata)
mtext("Genetta servalina (UDZ)", side=3, line=0, cex=0.5, font=3)

# UDZ.species.UDZ.Cephalophus spadix k=10
k=10
umf <- create.umf(k)
Tmin <- extract.Tmin(k)
fm.e.BTmin <- colext(psiformula=~1, gammaformula=~1, epsilonformula=~Biotic * Tmin, pformula=~Biotic * Tmin, data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm.e.BTmin
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.extTmin(inputmodel, inputdata)
mtext("Cephalophus spadix (UDZ)", side=3, line=0, cex=0.5, font=3)


# VB_.species.VB_.Mazama temama	k=5
k=5
umf <- create.umf(k)
Tmax <- extract.Tmax(k)
fm.c.BTmax <- colext(psiformula=~1, gammaformula=~Biotic * Tmax, epsilonformula=~1, pformula=~Biotic * Tmax, data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm.c.BTmax
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTmax(inputmodel, inputdata)
mtext("Mazama temama (VB)", side=3, line=0, cex=0.5, font=3)


# VB_.species.VB_.Dasypus novemcinctus k=3
k=3
umf <- create.umf(k)
Tvar <- extract.Tvar(k)
fm.ce.BTvar <- colext(psiformula=~1, gammaformula=~Biotic * Tvar, epsilonformula=~Biotic * Tvar, pformula=~Biotic * Tvar, data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm.ce.BTvar
#inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
#plot.colTvar(inputmodel, inputdata)
#inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
#extra.colTvar(inputmodel, inputdata)
#mtext("Dasypus novemcinctus (VB)", side=3, line=0, cex=0.5, font=3)


inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.extTvar(inputmodel, inputdata)
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.extTvar(inputmodel, inputdata)
mtext("Dasypus novemcinctus (VB)", side=3, line=0, cex=0.5, font=3)
