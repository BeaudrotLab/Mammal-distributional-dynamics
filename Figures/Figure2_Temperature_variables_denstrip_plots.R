# Create figure showing distribution of temperature trends from individual camera trap data

y <- (-1.3):2
par(mar=c(3,3,1,1))
plot(y, type="n", 
     xlim=c(0.5,7), 
     bty="n", las=1, ylab="Temperature (C)", 
     cex.axis=1.5, bg="transparent", mgp=c(3,1,0), xaxt="n")
axis(1, at=c(1,2,3,4,5,6,7), 
     labels=c("BIF", "NAK", "PSH", "RNF", "UDZ", "VB", "YAN"), 
     cex.axis=1.5)

sites <- list(BIF.Tmax, NAK.Tmax, PSH.Tmax, RNF.Tmax, UDZ.Tmax, VB.Tmax, YAN.Tmax)
sites <- list(BIF.Tmin, NAK.Tmin, PSH.Tmin, RNF.Tmin, UDZ.Tmin, VB.Tmin, YAN.Tmin)
sites <- list(BIF.Tvar, NAK.Tvar, PSH.Tvar, RNF.Tvar, UDZ.Tvar, VB.Tvar, YAN.Tvar)


hold <- list()

for(i in 1:length(sites)){
  ys <- as.numeric(colnames(sites[[i]]))
  mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
  hold[[i]] <- apply(as.matrix(sites[[i]]), MARGIN=1, FUN=mod.fn)
  denstrip(x=hold[[i]], na.rm=TRUE, at=i, horiz=FALSE, width=0.5,
         ticks=c(min(hold[[i]], na.rm=TRUE), max(hold[[i]], na.rm=TRUE)), 
         mticks=median(hold[[i]], na.rm=TRUE))
}
