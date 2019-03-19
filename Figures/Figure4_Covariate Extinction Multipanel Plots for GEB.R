# PURPOSE: Code for figure with covariates of extinction for GEB
# AUTHOR: Lydia Beaudrot
# Date: March 2, 2018; edited March 18, 2019
# Extinction covariate plots (excludes species with interaction term in top model)
# Populations wtih top models run on 12-21-2017; reanalyzed (with same results for proofs)    

rm(list=ls())
library(unmarked)
source(file="Occupancy Modelling/Helper_functions.R")
load("Occupancy Modelling/Species7sites_Include.RData")
load("Occupancy Modelling/BIOTIC_Include.RData")

Species_data <- Species7sites_Include
BIOTIC <- BIOTIC_Include
nms=names(Species_data)

# Set population using corresponding row number (k) in nms

############### BEGIN PLOTTING #######################
######### EXTINCTION (5, 3) ##########
par(mar=c(1,0.5,2,0.5), mfcol=c(5,3), mgp=c(2,1,0))


########### TMIN (5) ###########


#56	Fossa fossana (RNF)	f.eE.BioticTmin	~Elevation ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 56
Tmin <- extract.Tmin(k)
inputmodel <- f.eE.BioticTmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Fossa fossana (RNF)", side=3, cex=0.55, font=3)


#46	Nasua nasua (YAN)	f.eE.BioticTmin	~Elevation ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 46
Tmin <- extract.Tmin(k)
inputmodel <- f.eE.BioticTmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Nasua nasua (YAN)", side=3, cex=0.55, font=3)


#50	Atherurus macrourus (NAK)	f.e.Tmin	~1 ~ 1 ~ Tmin ~ Tmin
k <- 50
Tmin <- extract.Tmin(k)
inputmodel <- f.e.Tmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Atherurus macrourus (NAK)", side=3, cex=0.55, font=3)


#1	Cuniculus paca (VB)	f.e.BioticTmin	~1 ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 1
Tmin <- extract.Tmin(k)
inputmodel <- f.e.BioticTmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Cuniculus paca (VB)", side=3, cex=0.55, font=3)

#22	Cephalophus nigrifrons (BIF)	f.e.BioticTmin	~1 ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 22
Tmin <- extract.Tmin(k)
inputmodel <- f.e.BioticTmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Cephalophus nigrifrons (BIF)", side=3, cex=0.55, font=3)

########### TVAR (1) ###########

#38	Tragulus kanchil (PSH)	f.e.Tvar	~1 ~ 1 ~ Tvar ~ Tvar
k <- 38
Tvar <- extract.Tvar(k)
inputmodel <- f.e.Tvar(create.umf(k))
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60))
plot.extTvar(inputmodel, inputdata)
mtext("Tragulus kanchil (PSH)", side=3, cex=0.55, font=3)



########### TMAX (3) ###########
#13	Cricetomys gambianus (UDZ)	f.e.Tmax	~1 ~ 1 ~ Tmax ~ Tmax
k <- 13
Tmax <- extract.Tmax(k)
inputmodel <- f.e.Tmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Cricetomys gambianus (UDZ)", side=3, cex=0.55, font=3)


#8	Bdeogale crassicauda (UDZ)	f.e.Tmax	~1 ~ 1 ~ Tmax ~ Tmax
k <- 8
Tmax <- extract.Tmax(k)
inputmodel <- f.e.Tmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Bdeogale crassicauda (UDZ)", side=3, cex=0.55, font=3)


#27	Pan troglodytes (BIF)	f.ceE.BioticTmax	~Elevation ~ Biotic + Tmax ~ Biotic + Tmax ~ Biotic + Tmax
k <- 27
Tmax <- extract.Tmax(k)
inputmodel <- f.ceE.BioticTmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Pan troglodytes (BIF)", side=3, cex=0.55, font=3)


# BIOTIC (6)	

#46	Nasua nasua (YAN)	f.eE.BioticTmin	~Elevation ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 46
Biotic <- extract.Biotic(k)
inputmodel <- f.eE.BioticTmin(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmin=rep(0, length=60), Int=rep(0, length=60))
plot.extBiotic(inputmodel, inputdata)
mtext("Nasua nasua (YAN)", side=3, cex=0.55, font=3)

#27	Pan troglodytes (BIF)	f.ceE.BioticTmax	~Elevation ~ Biotic + Tmax ~ Biotic + Tmax ~ Biotic + Tmax
k <- 27
Biotic <- extract.Biotic(k)
inputmodel <- f.ceE.BioticTmax(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmax=rep(0, length=60), Int=rep(0, length=60))
plot.extBiotic(inputmodel, inputdata)
mtext("Pan troglodytes (BIF)", side=3, cex=0.55, font=3)

#29	Potamochoerus larvatus (BIF)	f.eE.Biotic	~Elevation ~ 1 ~ Biotic ~ Biotic
k <- 29
Biotic <- extract.Biotic(k)
inputmodel <- f.eE.Biotic(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60))
plot.extBiotic(inputmodel, inputdata)
mtext("Potamochoerus larvatus (BIF)", side=3, cex=0.55, font=3)


#1	Cuniculus paca (VB)	f.e.BioticTmin	~1 ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 1
Biotic <- extract.Biotic(k)
inputmodel <- f.e.BioticTmin(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmin=rep(0, length=60), Int=rep(0, length=60))
plot.extBiotic(inputmodel, inputdata)
mtext("Cuniculus paca (VB)", side=3, cex=0.55, font=3)


#22	Cephalophus nigrifrons (BIF)	f.e.BioticTmin	~1 ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 22
Biotic <- extract.Biotic(k)
inputmodel <- f.e.BioticTmin(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmin=rep(0, length=60), Int=rep(0, length=60))
plot.extBiotic(inputmodel, inputdata)
mtext("Cephalophus nigrifrons (BIF)", side=3, cex=0.55, font=3)


#56	Fossa fossana (RNF)	f.eE.BioticTmin	~Elevation ~ 1 ~ Biotic + Tmin ~ Biotic + Tmin
k <- 56
Biotic <- extract.Biotic(k)
inputmodel <- f.eE.BioticTmin(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmin=rep(0, length=60), Int=rep(0, length=60))
plot.extBiotic(inputmodel, inputdata)
mtext("Fossa fossana (RNF)", side=3, cex=0.55, font=3)




# ADD AXES AND OTHER ANNOTATIONS IN KEYNOTE FILE "Extinction Risk Figures Global Ecology & Biogeography"

