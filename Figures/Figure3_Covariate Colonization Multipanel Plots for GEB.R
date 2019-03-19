# PURPOSE: Code for multipanel figure with covariates of colonization for GEB
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
######### COLONIZATION (5, 4) ##########
par(mar=c(1,0.5,2,0.5), mfcol=c(5,4), mgp=c(2,1,0))




############# TMIN #############



#41	Dasyprocta fuliginosa (YAN)	f.c.Tmin	~1 ~ Tmin ~ 1 ~ Tmin
k <- 41
Tmin <- extract.Tmin(k)
inputmodel <- f.c.Tmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60))
plot.colTmin(inputmodel, inputdata)
mtext("ERASE ME!!!!", side=3, cex=0.55, font=3)


## REPEAT CODE FOR FIRST MODEL AS PLACE HOLDER, THEN MANUALLY BLOCK OUT EXTRA GRAPH ##
plot.colTmin(inputmodel, inputdata)
mtext("Dasyprocta fuliginosa (YAN)", side=3, cex=0.55, font=3)


#2	Dasyprocta punctata (VB)	f.c.Tmin	~1 ~ Tmin ~ 1 ~ Tmin
k <- 2
Tmin <- extract.Tmin(k)
inputmodel <- f.c.Tmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Dasyprocta punctata (VB)", side=3, cex=0.55, font=3)


#18	Paraxerus vexillarius (UDZ)	f.cE.Tmin	~Elevation ~ Tmin ~ 1 ~ Tmin
k <- 18
Tmin <- extract.Tmin(k)
inputmodel <- f.cE.Tmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Paraxerus vexillarius (UDZ)", side=3, cex=0.55, font=3)

#40	Cuniculus paca (YAN)	f.c.BioticTmin	~1 ~ Biotic + Tmin ~ 1 ~ Biotic + Tmin
k <- 40
Tmin <- extract.Tmin(k)
inputmodel <- f.c.BioticTmin(create.umf(k))
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Cuniculus paca (YAN)", side=3, cex=0.55, font=3)


## REPEAT CODE FOR LAST MODEL AS PLACE HOLDER, THEN MANUALLY BLOCK OUT EXTRA GRAPH ##
plot.colTmin(inputmodel, inputdata)
mtext("ERASE ME!!!!", side=3, cex=0.55, font=3)



####### TVAR ###########

#21	Caracal aurata (BIF)	f.c.Tvar	~1 ~ Tvar ~ 1 ~ Tvar
k <- 21
Tvar <- extract.Tvar(k)
inputmodel <- f.c.Tvar(create.umf(k))
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Caracal aurata (BIF)", side=3, cex=0.55, font=3)


#39	Tupaia glis (PSH)	f.c.Tvar	~1 ~ Tvar ~ 1 ~ Tvar
k <- 39
Tvar <- extract.Tvar(k)
inputmodel <- f.c.Tvar(create.umf(k))
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Tupaia glis (PSH)", side=3, cex=0.55, font=3)

#43	Eira barbara (YAN)	f.c.Tvar	~1 ~ Tvar ~ 1 ~ Tvar
k <- 43
Tvar <- extract.Tvar(k)
inputmodel <- f.c.Tvar(create.umf(k))
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Eira barbara (YAN)", side=3, cex=0.55, font=3)

#6	Pecari tajacu (VB)	f.c.Tvar	~1 ~ Tvar ~ 1 ~ Tvar
k <- 6
Tvar <- extract.Tvar(k)
inputmodel <- f.c.Tvar(create.umf(k))
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Pecari tajacu (VB)", side=3, cex=0.55, font=3)



######## Tmax ##########
#54	Tragulus kanchil (NAK)	f.cE.BioticTmax	~Elevation ~ Biotic + Tmax ~ 1 ~ Biotic + Tmax  
k <- 54
Tmax <- extract.Tmax(k)
inputmodel <- f.cE.BioticTmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
mtext("Tragulus kanchil (NAK)", side=3, cex=0.55, font=3)


#45	Mazama americana (YAN)	f.cE.Tmax	~Elevation ~ Tmax ~ 1 ~ Tmax
k <- 45
Tmax <- extract.Tmax(k)
inputmodel <- f.cE.Tmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60))
plot.colTmax(inputmodel, inputdata)
mtext("Mazama americana (YAN)", side=3, cex=0.55, font=3)


#19	Potamochoerus larvatus (UDZ)	f.c.Tmax	~1 ~ Tmax ~ 1 ~ Tmax
k <- 19
Tmax <- extract.Tmax(k)
inputmodel <- f.c.Tmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60))
plot.colTmax(inputmodel, inputdata)
mtext("Potamochoerus larvatus (UDZ)", side=3, cex=0.55, font=3)


#27	Pan troglodytes (BIF)	f.ceE.BioticTmax	~Elevation ~ Biotic + Tmax ~ Biotic + Tmax ~ Biotic + Tmax
k <- 27
Tmax <- extract.Tmax(k)
inputmodel <- f.ceE.BioticTmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
mtext("Pan troglodytes (BIF)", side=3, cex=0.55, font=3)

#24	Cercopithecus lhoesti (BIF)	f.c.BioticTmax	~1 ~ Biotic + Tmax ~ 1 ~ Biotic + Tmax
k <- 24
Tmax <- extract.Tmax(k)
inputmodel <- f.c.BioticTmax(create.umf(k))
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
mtext("Cercopithecus lhoesti (BIF)", side=3, cex=0.55, font=3)




######### BIOTIC ##########
#nms	label	FUNCTION TO USE	formula
#40	Cuniculus paca (YAN)	f.c.BioticTmin	~1 ~ Biotic + Tmin ~ 1 ~ Biotic + Tmin
k <- 40
Biotic <- extract.Biotic(k)
inputmodel <- f.c.BioticTmin(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmin=rep(0, length=60), Int=rep(0, length=60))
plot.colBiotic(inputmodel, inputdata)
mtext("Cuniculus paca (YAN)", side=3, cex=0.55, font=3)


#27	Pan troglodytes (BIF)	f.ceE.BioticTmax	~Elevation ~ Biotic + Tmax ~ Biotic + Tmax ~ Biotic + Tmax
k <- 27
Biotic <- extract.Biotic(k)
inputmodel <- f.ceE.BioticTmax(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmax=rep(0, length=60), Int=rep(0, length=60))
plot.colBiotic(inputmodel, inputdata)
mtext("Pan troglodytes (BIF)", side=3, cex=0.55, font=3)


#30	Atherurus macrourus (PSH)	f.c.Biotic	~1 ~ Biotic ~ 1 ~ Biotic
k <- 30
Biotic <- extract.Biotic(k)
inputmodel <- f.c.Biotic(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60))
plot.colBiotic(inputmodel, inputdata)
mtext("Atherurus macrourus (PSH)", side=3, cex=0.55, font=3)

#24	Cercopithecus lhoesti (BIF)	f.c.BioticTmax	~1 ~ Biotic + Tmax ~ 1 ~ Biotic + Tmax
k <- 24
Biotic <- extract.Biotic(k)
inputmodel <- f.c.BioticTmax(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmax=rep(0, length=60), Int=rep(0, length=60))
plot.colBiotic(inputmodel, inputdata)
mtext("Cercopithecus lhoesti (BIF)", side=3, cex=0.55, font=3)

#54	Tragulus kanchil (NAK)	f.cE.BioticTmax	~Elevation ~ Biotic + Tmax ~ 1 ~ Biotic + Tmax  
k <- 54
Biotic <- extract.Biotic(k)
inputmodel <- f.cE.BioticTmax(create.umf(k))
inputdata <- data.frame(Biotic=seq(min(Biotic),max(Biotic),length=60), Tmax=rep(0, length=60), Int=rep(0, length=60))
plot.colBiotic(inputmodel, inputdata)
mtext("Tragulus kanchil (NAK)", side=3, cex=0.55, font=3)

#par(mar=c(3,0.5,1,0))


# ADD AXES AND OTHER ANNOTATIONS IN KEYNOTE FILE "Extinction Risk Figures Global Ecology & Biogeography"

