#' Functional traits, phylogenies, communities, simulations
#'
#' Simulates the influence of functional traits and phylogenies on
#' ecological communities, via the \code{\link{fpcomSims}} function. 
#' 
#' @docType package
#' @name funphylocom
#' @import ape picante
#' @aliases funphylocom package-funphylocom
NULL

#' should i use this ??
#' http://stackoverflow.com/a/12429344/2047693
if(getRversion() >= "2.15.1")  utils::globalVariables("pic3")

#' Functional-phylogenetic community simulations
#'
#' Simulations for assessing the relative importance of phylogenetic 
#' versus functional information for understanding variation in
#' community composition.
#'
#' Simulations are based on the following procedure:
#' \describe{
#'
#'	\item{Phylogeny}{A tree is simulated using \code{rcoal} with
#'	default settings.}
#'
#'	\item{Trait evolution}{Two traits are simulated under Brownian motion. 
#'	One trait is called 'observed' and the other 'unknown'. Both traits
#'  influence species' probabilies of occurrence, but only the 'observed'
#'	trait is included in the output. The idea behind this distinction is that
#'	information about the 'unknown' trait will possibly be present in
#'	the phylogeny, which is assumed known.}
#'
#'	\item{Probabilities of occurrence}{For each species at each site, these
#'	probabilities depend logistically on two gradients; one is observed
#'	and the other is unknown, each corresponding to the observed and unknown
#'	traits (see the \code{sim.envObserved} and \code{sim.envUnknown} arguments).
#'	The corresponding traits are the logit-scale slopes of these logistic
#'	curves -- in other words each species depends on the gradient differently
#'	depending on their traits. When \code{p} equals one (zero), only the observed
#'	(unknown) gradient and trait determine probability of occurrence. When 
#'	\code{p} is between zero and one, both the observed and unknown
#'	trait-gradient combinations have some influence.}
#'
#'	\item{Occurrence}{Each species is present at each site with probability
#'	given by the corresponding probability of occurrence.}
#'
#' }
#'
#' @param n Number of sites to simulate.
#' @param m Number of species to simulate.
#' @param p Number between 0 and 1 giving the relative importance of
#'	observed versus unknown traits in determining community structure.
#' @param diverg.obs A vector with m elements giving the divergence
#'	effects on the observed trait for each species.
#' @param diverg.unk A vector with m elements giving the divergence
#'	effects on the unknown trait for each species.
#' @param sim.envObserved A numeric vector of simulated values for the 
#'	observed environmental variables.
#' @param sim.envUnknown A numeric vector of simulated values for the 
#'	unknown environmental variables.
#' @param site.names Character vector of site names.
#' @param spp.names Character vector of species names.
#' @return An object of class \code{fpcomSims} with components:
#'	\item{comm}{An n-by-m matrix of presences and absences}
#'	\item{probs}{An n-by-m matrix of probabilities of occurrence}
#'	\item{traits}{A length-m vector of values for the observed trait}
#'	\item{env}{A length-n vector of values for the observed gradient}
#'	\item{tree}{A phylo object with the phylogenetic tree relating
#'	the m species}
#' @export
fpcomSims <- function(n, m, p = 0.5,
	diverg.obs = rep(0, m),
	diverg.unk = rep(0, m),
	sim.envObserved = as.vector(scale(rnorm(n))),
	sim.envUnknown = as.vector(scale(rnorm(n))),
	site.names = numnames(n,"site"), 
	spp.names = numnames(m,"sp")
){
	
	names(sim.envObserved) <- names(sim.envUnknown) <- site.names
	
	# simulate evolutionary history
	sim.tree <- rcoal(m, tip.label = spp.names)
	sim.tree$tip.label <- spp.names
	
	# store the resulting traits for the species
	sim.traits.obs <- as.vector(t(chol(vcv(sim.tree))) %*% rnorm(m, sd = 1))
	sim.traits.unk <- as.vector(t(chol(vcv(sim.tree))) %*% rnorm(m, sd = 1))
	sim.traits.obs <- sim.traits.obs + diverg.obs
	sim.traits.unk <- sim.traits.unk + diverg.unk
	names(sim.traits.obs) <- spp.names

	# simulate the contemporary communities
	sim.eta <- (p*outer(sim.envObserved, sim.traits.obs)) + 
			   ((1-p)*outer(sim.envUnknown, sim.traits.unk))
	sim.eta <- scale(as.vector(sim.eta))
	dim(sim.eta) <- c(n, m)
	sim.p <- exp(sim.eta)/(1+exp(sim.eta))
	sim.comm <- matrix(rbinom(n*m, 1, sim.p), n, m)
	dimnames(sim.comm) <- dimnames(sim.p) <- list(site.names, spp.names)

	ecosysfunc <- (p*sim.envObserved) + ((1-p)*sim.envUnknown)

	out <- list(comm = sim.comm, 
		probs = sim.p,
		traits = structure(sim.traits.obs, unknown = sim.traits.unk),
		env = sim.envObserved,
		ecosysfunc = ecosysfunc,
		tree = sim.tree,
		tuning = p)
	class(out) <- "fpcomSims"
	return(out)
}

update.fpcomSims <- function(object, nnew = 1, 
	sim.envObserved = rnorm(nnew),
	sim.envUnknown = rnorm(nnew),
	site.names = numnames(length(object$env) + nnew, "site"),
	spp.names = numnames(length(object$traits),"sp"), ...){

	m <- ncol(object$probs)
	
	sim.traits.obs <- structure(object$traits, unknown = NULL)
	sim.traits.unk <- attr(object$traits, 'unknown')
	p <- object$tuning

	# simulate the new contemporary communities
	sim.eta <- (p*outer(sim.envObserved, sim.traits.obs)) + 
			   ((1-p)*outer(sim.envUnknown, sim.traits.unk))
	sim.eta <- scale(as.vector(sim.eta))
	dim(sim.eta) <- c(nnew, m)
	sim.p <- exp(sim.eta)/(1+exp(sim.eta))
	sim.comm <- matrix(rbinom(nnew*m, 1, sim.p), nnew, m)
	
	sim.p <- rbind(object$probs, sim.p)
	sim.comm <- rbind(object$comm, sim.comm)
	sim.envObserved <- c(object$env, sim.envObserved)
		
	dimnames(sim.comm) <- dimnames(sim.p) <- list(site.names, spp.names)
	names(sim.envObserved) <- site.names
	
	ecosysfunc <- sim.envObserved + sim.envUnknown
	
	out <- list(comm = sim.comm, 
		probs = sim.p,
		traits = object$traits,
		env = sim.envObserved,
		ecosysfunc = ecosysfunc,
		tree = object$tree,
		tuning = p)
	class(out) <- "fpcomSims"
	return(out)
}

#' Plot fpcomSims objects
#'
#' This plot method can produce four different graphical summaries of the
#' output of the \code{\link{fpcomSims}} function.
#'
#' \describe{
#'
#'	\item{\code{plottype} equals \code{"distance"}}{A scatterplot
#'	of the phylogenetic versus functional distances between the species 
#'	pairs are produced. Species numbers are used to label the points. The
#'	\code{\link{cophenetic}} function is used to compute the phylogenetic
#'	distances and the \code{\link{dist}} function is used to compute the
#'	functional Euclidean distances. Both distances are standardised such
#'	that the maximum distance is one.}
#'	
#'	\item{\code{plottype} equals \code{"traitgram"}}{A \code{\link{traitgram}}
#'	is produced, which combines both the phylogeny and the observed trait.}
#'
#'	\item{\code{plottype} equals \code{"gradient"}}{A scatterplot of the 
#'	probabilities of occurrence versus the observed gradient is produced.
#'	The species are identified by their numbers. Note that these probabilities
#'	of occurrence depend only on the two gradients and the two traits, and
#'	therefore are not subject to variation among sites with identical gradient
#'	values. Such variation is expressed in the \code{comm} element of
#'	\code{fpcomSims} objects, which gives not probability of occurrence but
#'	occurrence itself. Note that if the \code{p} argument to \code{fpcomSims}
#'	is one, then only the observed trait and gradient determine probability of
#'	occurrence and the plot consists of perfect sigmoid curves.  But if \code{p}
#'	is zero then only the unknown gradient and trait determine probability of
#'	occurrence and the plot is just noise. In this latter case we would expect
#'	phylogenetic distance to provide more important information about communities.}
#'
#'	\item{\code{plottype} equals \code{"ordination"}}{A special type of 
#'	ordination of the species is produced. This ordination is based on the
#'	probabilities of occurrence and not on occurrence itself, and so it is able
#'	to fully explain all variation on two axes (if \code{p} is not either zero
#'	of one). Therefore, such an ordination is not possible in practice with real
#'	data but it is useful in this case for fully representing distances between
#'	species in species distribution space. Technically it is a singular value
#'	decomposition of the logit transformed probabilities of occurrence. Note
#'	well the percentages of variation explained by the two axes, and that this
#'	ordination is gauranteed to have 100 percent of the variation explained by
#'	the first two axes.}
#' }
#'
#' @param x An \code{fpcomSims} object.
#' @param y Not used at the moment.
#' @param cex.add The cex graphics parameter for additional graphical. 
#'	elements not produced by the \code{\link{plot}} command.
#' @param plottype The type of plot to produce -- see details.
#' @param ... Additional parameters to be passed to \code{\link{plot}}.
#' @method plot fpcomSims
#' @return No return value, called for its side-effect of producing a plot.
#' @export
plot.fpcomSims <- function(x,y,cex.add=0.8,
	plottype=c("distance","traitgram","gradient","ordination"),...){
	
	if(plottype[1]=="distance"){
		PDist <- cophenetic(x$tree)/max(cophenetic(x$tree))
		PDist <- PDist[order(row.names(PDist)),order(row.names(PDist))]
		FDist <- as.matrix(dist(x$traits))/max(as.matrix(dist(x$traits)))
		plot(PDist,FDist,type="n",
			xlab="phylogenetic distance",
			ylab="functional distance",
			...)
		text(as.dist(PDist),as.dist(FDist),
			paste(as.dist(col(PDist)),as.dist(row(FDist)),sep=","),
			cex=cex.add)
	}
	
	if(plottype[1]=="traitgram"){
		traitgram(x$traits,x$tree,...)
	}
	
	if(plottype[1]=="gradient"){
		plot(x$env,x$probs[,1],type="n",ylim=c(0,1),
			xlab="gradient",
			ylab="probability of occurrence",
			...)
		for(i in 1:ncol(x$comm)){
			text(x$env,x$probs[,i],i,cex=cex.add)
		}
	}
	
	if(plottype[1]=="ordination"){
		x.ord <- svd(log(x$probs/(1-x$probs)))
		plot(x.ord$v[,1:2]%*%diag(x.ord$d[1:2]),type="n",asp=1,
			xlab=paste("Ordination axis I (",
				100*round(x.ord$d[1]/sum(x.ord$d),2),"%)",sep=""),
			ylab=paste("Ordination axis II (",
				100*round(x.ord$d[2]/sum(x.ord$d),2),"%)",sep=""),
			...)
		text(x.ord$v[,1:2]%*%diag(x.ord$d[1:2]),labels=names(x$traits),cex=cex.add)
	}
}



#' Generates numbered names
#'
#' Generates numbered names based on a character prefix that will
#' be alphanumerically sorted as expected.
#'
#' @param n Number of names.
#' @param prefix Character prefix to go in front of the numbers.
#' @return A character vector of numbered names.
#' @export
numnames <- function(n, prefix = "name"){
	n <- as.integer(n)
	if(n < 1) stop("number of names must be one or more")
	zeropad.code <- paste("%0", nchar(n), ".0f", sep = "")
	numpart <- sprintf(zeropad.code, 1:n)
	paste(prefix, numpart, sep = "")
}

#' Tail of a list
#'
#' utility function to make extractions of the 
#' trait values at the tips easier
#'
#' @param x A list.
#' @param n Number of lines.
#' @param ... Not used.
#' @return The first \code{n} elements of each list element.
#' @method tail list
tail.list <- function(x,n=6L,...){
	sapply(x,function(xi)rev(rev(xi)[seq_len(n)]))
}

#' Functional phylogenetic distance
#'
#' Computates a functional phylogenetic distance matrix from 
#' phylogenetic and functional distance matrices and two
#' weighting parameters.
#'
#' @param PDist A phylogenetic distance matrix.
#' @param FDist A functional distance matrix.
#' @param a A number between 0 and 1 giving the amount of weight.
#'	to put on \code{PDist} relative to \code{FDist}.
#' @param p A number giving the \code{p}-norm.
#' @param ord Order rows and columns of the matrices?  (defaults
#'  to \code{TRUE}).
#' @return A distance matrix.
#' @note This function is not very user friendly yet.  There are
#'	no doubt many use cases that I've ignored.
#' @export
FPDist <- function(PDist, FDist, a, p, ord = TRUE){
	PDist <- as.matrix(PDist, rownames.force = TRUE)
	FDist <- as.matrix(FDist, rownames.force = TRUE)
	if(
		is.null(rownames(FDist)) ||
		is.null(colnames(FDist)) ||
		is.null(rownames(PDist)) ||
		is.null(colnames(PDist))
	) stop('distance matrices must have row and column names')
	if(
		any(rownames(FDist) != colnames(FDist)) ||
		any(rownames(PDist) != colnames(PDist))
	) stop('row and column names must match for distance matrices')
	if(ord){
		FDist <- FDist[order(rownames(FDist)), order(rownames(FDist))]
		PDist <- PDist[order(rownames(PDist)), order(rownames(PDist))]
	}
	if(
		any(rownames(FDist) != rownames(PDist)) ||
		any(colnames(FDist) != colnames(PDist))
	) stop('FDist and PDist must have same row and column names')
	
	FDist <- FDist/max(FDist)
	PDist <- PDist/max(PDist)
	
	((a*(PDist^p)) + ((1-a)*(FDist^p)))^(1/p)
}

#' Rao's quadratic entropy
#'
#' Computes Rao's quadratic entropy from a distance matrix and
#' community matrix.
#'
#' @param D A species by species distance matrix.
#' @param X A sites by species community matrix.
#' @param ord Order rows and columns of the matrices? (defaults
#'  to \code{TRUE}).  Note that sites are never ordered.
#' @return A vector of diversity indices (one for each site).
#' @export
rao <- function(X, D, ord = TRUE){
	if(
		is.null(rownames(D)) ||
		is.null(colnames(D))
	) stop('distance matrix must have row and column names')
	if(is.null(colnames(X)))
		stop('community matrix must have column (i.e. species) names')
	if(any(rownames(D) != colnames(D)))
		stop('row and column names must match for distance matrix')
	if(ord){
		X <- X[, order(colnames(X))]
		D <- D[order(rownames(D)), order(rownames(D))]
	}
	if(any(colnames(X) != colnames(D)))
		stop('species names must match')
	
	# row sums
	rs <- apply(X, 1, sum)
	
	# get relative abundances by dividing by row sums
	X.rel <- sweep(X, 1, rs, FUN = '/')
	
	# make relative abundances be zero for all species
	# at sites for which all species have zero abundance.
	# this is needed b/c of division by zero.
	X.rel[is.nan(X.rel)] <- 0
	
	# combine rao index, (x %*% D %*% x)/2, for each site
	out <- apply(X.rel, 1, function(x) x %*% D %*% x)/2

	names(out) <- rownames(X)
	return(out)
}


#' Mean pairwise distance
#' 
#' Calculates mean pairwise distance separating taxa in a community
#' (slightly modified from the \code{picante} function: 
#' \code{\link{mpd}})
#'
#' Two changes from original \code{\link{mpd}} function:  
#' (1) avoid NA values (TODO: better description here).
#' (2) give useful error messages for poorly named \code{samp} and \code{dis}.
#' Also please be aware that if \code{abundance.weighted = FALSE}, then
#' \code{samp} should be explicitly a presence-absence matrix (i.e. with
#' zeros and ones only).
#'
#' @param samp Community data matrix.
#' @param dis Interspecific distance matrix.
#' @param abundance.weighted Should mean pairwise diatnces be weighted
#'  by species abundance? (default = FALSE).
#' @param ord Order rows and columns of the matrices? (defaults
#'  to \code{TRUE}).  Note that sites are never ordered.
#' @return Vector of MPD values for each community.
#' @author Code modified from original version by Steven Kembel.
#' @export
mpd. <- function (samp, dis, abundance.weighted = FALSE, ord = TRUE) 
{

	if(
		is.null(rownames(dis)) ||
		is.null(colnames(dis))
	) stop('distance matrix must have row and column names')
	if(is.null(colnames(samp)))
		stop('community matrix must have column (i.e. species) names')
	if(any(rownames(dis) != colnames(dis)))
		stop('row and column names must match for distance matrix')
	if(ord){
		samp <- samp[, order(colnames(samp))]
		dis <- dis[order(rownames(dis)), order(rownames(dis))]
	}
	if(any(colnames(samp) != colnames(dis)))
		stop('species names must match')

	
    N <- dim(samp)[1]
    mpd <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            if (abundance.weighted) {
                sample.weights <- t(as.matrix(samp[i, sppInSample, 
                  drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                  drop = FALSE])
                mpd[i] <- weighted.mean(sample.dis, sample.weights)
            }
            else {
                mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
            }
        }
        else {
            mpd[i] <- 0  # used to be NA in picante.  only change here.
        }
    }
    mpd
}


#' Null models for functional-phylogenetic diversity
#'
#' Simulate expectations (under a null model) of mean pairwise distance for 
#' a set of communities with different species richness.
#'
#' If \code{plot == TRUE}, then a surface is drawn giving the
#' null distribution.  Lighter shades of gray give larger intervals 
#' with categories: 0.005-0.995 = 99\%, 0.025-0.975 = 95\%, 0.05-0.95 = 90\%,
#' 0.25-0.75 = 50\%.
#'
#' @param tab The community by species presence table.
#' @param Dist The distance (functional, phylogenetic, FPDist) on which the 
#'  mean pairwise distance is based.
#' @param Sim The number of permutations of the presence vector used to 
#'  make the estimations.
#' @param Plot TRUE or FALSE to make the plot of the expected average 
#'  mean pairwise distance, and the 5-95\% confidence interval.
#' @param ord Order rows and columns of the matrices? (defaults
#'  to \code{TRUE}).  Note that sites are never ordered.
#' @param disp99 Display the 99\% interval?
#' @return TODO
#' @export
ConDivSim<-function(tab, Dist, Sim, Plot = TRUE, ord = TRUE, disp99 = FALSE){
	if(
		is.null(rownames(Dist)) ||
		is.null(colnames(Dist))
	) stop('distance matrix must have row and column names')
	if(is.null(colnames(tab)))
		stop('community matrix must have column (i.e. species) names')
	if(any(rownames(Dist) != colnames(Dist)))
		stop('row and column names must match for distance matrix')
	if(ord){
		tab <- tab[, order(colnames(tab))]
		Dist <- Dist[order(rownames(Dist)), order(rownames(Dist))]
	}
	if(any(colnames(tab) != colnames(Dist)))
		stop('species names must match')

	## Calculate the observed species richness 
	SpeRich.obs <- rowSums(tab)
	Occur.obs <- colSums(tab)
	
	## Check that numbers are high enough for randomisation
	## (e.g. no communities with one species only)
	if(any(SpeRich.obs < 2)) 
		stop('all communities must have more than one species present')
	if(any(Occur.obs < 1)) 
		stop('all species must be present at more than zero communities')

	## Calculate the range of species richness in the communities
	if(min(SpeRich.obs)>2){
		SpRich <- (min(SpeRich.obs)-1):(max(SpeRich.obs)+1)
	}
	else {
		SpRich <- min(SpeRich.obs):(max(SpeRich.obs)+1)
	}

	## Calculate the regional species richness(gamma)
	NbSp <- ncol(tab)
	
	## Create the matrix in which to store the results
	Randomiz <- matrix(0, length(SpRich), 11)
	colnames(Randomiz)<-c("SpRich","ExpMeanMPD",
                              "ExpQ0.005MPD","ExpQ0.025MPD",
                              "ExpQ0.05MPD","ExpQ0.25MPD",
                              "ExpQ0.75MPD","ExpQ0.95MPD",
                              "ExpQ0.975MPD","ExpQ0.995MPD",
                              "ExpSdMPD")
	row.names(Randomiz) <- 1:length(SpRich)
	
	## calculate the observed mean pairwise distance for the community
	## (with the unweighted method from mpd in picante, it means only 
	## the lower triangle!)
	MPD.obs <- mpd.(tab, Dist, abundance.weighted = FALSE)
  
	## Loop on the species richness range to calculate the random 
	## expectations of mean pairwise distance
	for (k in 1:length(SpRich)){
    
    	# Vector of local richness within the possible range SpRich
    	locrich <- SpRich[k]
    
    	# Create the basic vector of presence
	    Vec.Pres <- c(rep(1, locrich), rep(0, NbSp-locrich))
    
		# Estimate the random expectations of mean pairwise distances 
		# for the random community with locrich as the species richness
		Sim.Div <- NULL

		for (i in 1:Sim) {
        	Vec.Pres.Sim <- sample(Vec.Pres)
	        sample.dis <- Dist[
	        	as.logical(Vec.Pres.Sim), 
	        	as.logical(Vec.Pres.Sim)
	        ]
	        
	        # line in mpd from picante with abundance.weighted = FALSE
        	Sim.Div[i] <- mean(sample.dis[lower.tri(sample.dis)])
      	}

		# Store the results in Randomiz
	    Randomiz[k,1]<-locrich
	    Randomiz[k,2]<-mean(Sim.Div)
	    Randomiz[k, 3:10] <- quantile(Sim.Div, 
	    	c(0.005, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975, 0.995))
    	Randomiz[k,11]<-sd(Sim.Div)
	}
	
	## Make the graph if required
	if(Plot == TRUE){
    
   		plot(Randomiz[,1:2],
    		ylim = c(min(Randomiz[, 3], MPD.obs), max(Randomiz[, 10], MPD.obs)),
			type="n", xlab="Species richness",ylab="Mean pairwise distance")
		
		## Surface for each of the quantile pair:
		if(disp99){
			polygon( ## 0.005-0.995 = 99%
				c(Randomiz[,1], Randomiz[length(SpRich):1, 1]), 
				c(Randomiz[,3], Randomiz[length(SpRich):1, 10]),
				col=grey(0.9), border=NA)
    	}
    	
		polygon( ## 0.025-0.975 = 95%
			c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
			c(Randomiz[,4], Randomiz[length(SpRich):1,9]),
			col=grey(0.8),border=NA)
    	
		polygon( ## 0.05-0.95 = 90%
			c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
			c(Randomiz[,5], Randomiz[length(SpRich):1,8]),
			col=grey(0.7),border=NA)
    
		polygon( ## 0.25-0.75 = 50%
			c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
			c(Randomiz[,6], Randomiz[length(SpRich):1,7]),
			col=grey(0.6),border=NA)
    	
		## Average of the expected mean pairwise distance
		lines(Randomiz[,1:2], col = "black", lwd = 2)
    
		## Labelling the communities
		if(is.null(row.names(tab))){
			points(SpeRich.obs, MPD.obs)
		}
		else{
			text(SpeRich.obs, MPD.obs, labels = row.names(tab))
		}
    
    }
	
	return(Randomiz)
}




#' Null models for (TODO:  better title)
#'
#' Check for a given community how the mean pairwise distance changes with 
#' the parameter a in FPDist
#'
#' If \code{plot == TRUE}, then a surface is drawn giving the
#' null distribution.  Lighter shades of gray give larger intervals 
#' with categories: 0.005-0.995 = 99\%, 0.025-0.975 = 95\%, 0.05-0.95 = 90\%,
#' 0.25-0.75 = 50\%.
#'
#' @param comm A vector over the species in the regional pool, giving the
#'  presence or absence of the species occurring in the community
#'  (1 for present species, 0 for absent species).
#' @param PDist Phylogenetic distance matrix over the species in the regional 
#'  pool.
#' @param FDist Functional distance matrix over the species in the regional 
#'  pool.
#' @param Sim The number of permutations of the presence vector used to 
#'  make the estimations.
#' @param p Parameter passed to \code{\link{FPDist}} function.
#' @param Plot TRUE or FALSE to make the plot of the expected average 
#'  mean pairwise distance over a range of a-values (\code{a} is a parameter
#'  in the \code{\link{FPDist}} function), and the 5-95\% confidence interval.
#' @param ord Order rows and columns of the matrices? (defaults
#'  to \code{TRUE}).  Note that sites are never ordered.
#' @param disp99 Display the 99\% interval?
#' @param num.a Number of a-values to compute.
#' @return TODO
#' @export
ConDivSimComm<-function(comm, PDist, FDist, Sim, p = 2,
	Plot = TRUE, ord = TRUE, disp99 = FALSE, num.a = 20){
	
	if(!is.data.frame(comm)) stop('comm must be a one-row data frame')
	if(nrow(comm) != 1) stop('comm must be a one-row data frame')
	if(is.null(names(comm)))
		stop('community vector must have names')
	if(
		(ncol(comm) != nrow(PDist)) ||
		(ncol(comm) != ncol(PDist)) ||
		(ncol(comm) != nrow(FDist)) ||
		(ncol(comm) != ncol(FDist))
	) stop('number of species must match between comm, PDist, and FDist')
	if(
		is.null(rownames(PDist)) ||
		is.null(colnames(PDist)) ||
		is.null(rownames(FDist)) ||
		is.null(colnames(FDist))
	) stop('distance matrices must have row and column names')
	if(
		any(rownames(PDist) != colnames(PDist)) ||
		any(rownames(FDist) != colnames(FDist))
	) stop('row and column names must match for distance matrices')
	if(ord){
		comm <- comm[, order(colnames(comm))]
		PDist <- PDist[order(rownames(PDist)), order(rownames(PDist))]
		FDist <- FDist[order(rownames(FDist)), order(rownames(FDist))]
	}
	if(
		any(colnames(comm) != colnames(PDist)) ||
		any(colnames(comm) != colnames(FDist))
	) stop('species names must match between comm, PDist, and FDist')
	
	
	## Create the vector of a values
	a <- seq(0, 1, length.out = num.a)
  
  	## Calculate the regional species richness (gamma)
	NbSp<-length(comm)
  
	## Calculate the local species richness (alpha)
	locrich<-sum(comm)
	
	## Create the basic vector of presence
	Vec.Pres<-c(rep(1,locrich),rep(0,NbSp-locrich))
	
	## Create the matrix in which to store the results
	Randomiz<-matrix(0,length(a),12)
	colnames(Randomiz)<-c("SpRich","ExpMeanMPD",
                              "ExpQ0.01MPD","ExpQ0.05MPD",
                              "ExpQ0.1MPD","ExpQ0.25MPD",
                              "ExpQ0.75MPD","ExpQ0.9MPD",
                              "ExpQ0.95MPD","ExpQ0.99MPD",
                              "ExpSdMPD","Obs,MPD")
	row.names(Randomiz)<-1:length(a)
  
  
	## Loop on the a values: for each one calculate the random 
	## expectations and observed values of mean pairwise distance
	for(k in 1:length(a)){

		# Calculate the FPDistance based on the current a value
		FPDist.temp <- FPDist(PDist, FDist, a[k], p = p)
      
		# Calculate the observed mean pairwise distance for the 
		# community (with the unweighted method from mpd in picante)
		obs.dis <- FPDist.temp[as.logical(comm), as.logical(comm)]
		
		# line in mpd from picante with abundance.weighted = FALSE
		Randomiz[k,12] <- mean(obs.dis[lower.tri(obs.dis)])

		# declare Sim.Div to be filled in the loop below
		Sim.Div<-NULL

		# Estimate the random expectations of mean pairwise distances 
		# for the community		
		for (i in 1:Sim) {
			Vec.Pres.Sim<-sample(Vec.Pres)
			sample.dis <- FPDist.temp[
				as.logical(Vec.Pres.Sim), 
				as.logical(Vec.Pres.Sim)
			]
			# line in mpd from picante with abundance.weighted = FALSE
			Sim.Div[i]<-mean(sample.dis[lower.tri(sample.dis)])
		}
      
		# Store the results in Randomiz
	    Randomiz[k,1] <- a[k]
	    Randomiz[k,2]<-mean(Sim.Div)
	    Randomiz[k, 3:10] <- quantile(Sim.Div, 
	    	c(0.005, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975, 0.995))
    	Randomiz[k,11]<-sd(Sim.Div)    	
	}


	if(Plot == TRUE){
    	plot(
    		Randomiz[,1:2], 
    		ylim = c(min(Randomiz[,c(3,12)]), max(Randomiz[,c(10,12)])), 
    		type="l", 
    		main = paste(
    			"Community= ",row.names(comm)," - Sp Richness =",locrich
    		), 
    		xlab = "a", ylab = "Mean pairwise distance ")
    		
    	## Surface for each of the quantile pair:
		if(disp99) {
			polygon( ## 0.005-0.995 = 99%
				c(Randomiz[,1], Randomiz[length(a):1, 1]), 
				c(Randomiz[,3], Randomiz[length(a):1, 10]),
				col=grey(0.9), border=NA)
		}
    	
		polygon( ## 0.025-0.975 = 95%
			c(Randomiz[,1], Randomiz[length(a):1,1]),
			c(Randomiz[,4], Randomiz[length(a):1,9]),
			col=grey(0.8),border=NA)
    	
		polygon( ## 0.05-0.95 = 90%
			c(Randomiz[,1], Randomiz[length(a):1,1]),
			c(Randomiz[,5], Randomiz[length(a):1,8]),
			col=grey(0.7),border=NA)
    
		polygon( ## 0.25-0.75 = 50%
			c(Randomiz[,1], Randomiz[length(a):1,1]),
			c(Randomiz[,6], Randomiz[length(a):1,7]),
			col=grey(0.6),border=NA)

 
		lines(Randomiz[,1:2])
		
		lines(Randomiz[,c(1,12)], lty = 2)
	}
	
	out <- list(CommName=row.names(comm),CommSpRich=locrich,Simul=Randomiz)
	return(out)
}


#' Community traitgram
#' 
#' Displays a subset of species from a species pool within an Ackerly
#' traitgram (\code{\link{traitgram}}).
#'
#' This function is a modification of the \code{\link{traitgram}} function
#' in the \code{picante} package.
#' 
#' @param x See \code{\link{traitgram}}.
#' @param phy See \code{\link{traitgram}}.
#' @param xaxt See \code{\link{traitgram}}.
#' @param underscore See \code{\link{traitgram}}.
#' @param show.names See \code{\link{traitgram}}.
#' @param show.xaxis.values See \code{\link{traitgram}}.
#' @param method See \code{\link{traitgram}}.
#' @param edge.color A single color or a vector of two colors.
#' @param edge.lwd.in A single width value for species in the subset.
#' @param edge.lwd.out A single width value for species out of the subset 
#'  (ignored if \code{is.null(Group)}).
#' @param tip.color A single color or a vector of two colors.
#' @param tip.cex A single size value.
#' @param Group Either a subset of the tip labels in \code{phy} or
#'  the indices associated with the species in the community.
#' @param ... Additional arguments to \code{\link{traitgram}}.
#' @return No return value, but a traitgram is plotted.
#' @author David Ackerly, modified by Cecile Albert.
#' @export
comm.traitgram <- function (x, phy, xaxt = "s", underscore = FALSE, 
	show.names = TRUE, show.xaxis.values = TRUE, 
	method = c("ace", "pic"), edge.color = c('black', 'light grey'), 
	edge.lwd.in = 2, edge.lwd.out = 1, 
	tip.color = c('black', 'light grey'), 
	tip.cex = par("cex"), Group = NULL,
        name.jitter = 0, ...)
{
	
	# comments indicate modifications to original traitgram function
	
	## edge.color and tip.color must be either of length 1 or 2
	if(is.null(Group) && !((length(edge.color) %in% 1:2)))
		stop('edge.color must be either of length 1 or 2')
	if(is.null(Group) && !((length(tip.color) %in% 1:2)))
		stop('tip.color must be either of length 1 or 2')
	if(!is.null(Group) && (length(edge.color) != 2))
		stop('with non-null Group, edge.color must be of length 2')
	if(!is.null(Group) && (length(tip.color) != 2))
		stop('with non-null Group, tip.color must be of length 2')	
	
	
	## I added that to be sure that both trait and tree tip labels 
	## are in the same order
	x <- x[phy$tip.label]
    
    ## Identify species present in the community
    if(is.character(Group)) Group <- which(phy$tip.label %in% Group)

    Ntaxa = length(phy$tip.label)
    Ntot = Ntaxa + phy$Nnode
    phy = node.age(phy)
    ages = phy$ages[match(1:Ntot, phy$edge[, 2])]
    ages[Ntaxa + 1] = 0
    if (class(x) %in% c("matrix", "array")) {
        xx = as.numeric(x)
        names(xx) = row.names(x)
    }
    else xx = x
    if (!is.null(names(xx))) {
        umar = 0.1
        if (!all(names(xx) %in% phy$tip.label)) {
            print("trait and phy names do not match")
            return()
        }
        xx = xx[match(phy$tip.label, names(xx))]
    }
    else umar = 0.1
    lmar = 0.2
    if (xaxt == "s")
        if (show.xaxis.values)
            lmar = 1
        else lmar = 0.5
    if (method[1] == "ace")
        xanc = ace(xx, phy)$ace
    else xanc = pic3(xx, phy)[, 3]
    xall = c(xx, xanc)
    a0 = ages[phy$edge[, 1]]
    a1 = ages[phy$edge[, 2]]
    x0 = xall[phy$edge[, 1]]
    x1 = xall[phy$edge[, 2]]
    
	## Create color and width vectors for the edges
    #edge.lwd <- rep(edge.lwd, length.out = nrow(phy$edge))
    #tip.cex <- rep(tip.cex, length.out = Ntaxa)
    
    #if(is.null(Group)){
	#	edge.color <- rep(edge.color[1], length.out = nrow(phy$edge))
	#	tip.color <- rep(tip.color[1], length.out = Ntaxa)
   	#}
   	#else{
   		
   	#}
    #edge.color<-rep("light grey",nrow(phy$edge))
    #edge.color[edge.Group]<-"black"
    #tip.color<-rep("light grey",Ntaxa)
    #tip.color[Group]<-"black"
   #}
####
    tg = par(bty = "n", mai = c(lmar, 0.1, umar, 0.1))
    if (show.names) {
        maxNameLength = max(nchar(names(xx)))
        ylim = c(min(ages), max(ages) * (1 + maxNameLength/50))
        if (!underscore)
            names(xx) = gsub("_", " ", names(xx))
    }
    else ylim = range(ages)
    plot(range(c(x0, x1)), range(c(a0, a1)), type = "n", xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", bty = "n", ylim = ylim,
        cex.axis = 0.8)
    if (xaxt == "s")
        if (show.xaxis.values)
            axis(1, labels = TRUE)
        else axis(1, labels = FALSE)
    
    ## if no grouping, plot normal traitgram but with edge color or width
    ## as specified in the arguments
    if(is.null(Group)){
		segments(x0, a0, x1, a1,col = edge.color[1], lwd = edge.lwd.in)
    }
    ## else if grouping, plot traitgram with different colors for edges
    ## corresponding to tips belonging or not to the community
    else{
		edge.Group <- which.edge(phy,Group)
		segments(
			x0[-edge.Group], 
			a0[-edge.Group], 
			x1[-edge.Group], 
			a1[-edge.Group],
			col=edge.color[2], lwd = edge.lwd.out)
		segments(
			x0[edge.Group],
			a0[edge.Group],
			x1[edge.Group],
			a1[edge.Group],
			col=edge.color[1], lwd = edge.lwd.in)
    }
    if (show.names) {
    
    	## if no grouping, display normal tip labels
    	if(is.null(Group)){
        	text(
        		sort(xx+name.jitter), max(ages), labels = names(xx)[order(xx)],
            	adj = -0, srt = 90, cex = tip.cex, col = tip.color[1])
		}
		## else if grouping, display tips with different colors 
		## corresponding to tips belonging or not to the community
		else{
			tip.Group <- numeric(Ntaxa)
			tip.Group[Group] <- 1
			tip.Group <- as.logical(tip.Group)
			text(
				sort(xx+name.jitter)[!tip.Group[order(xx)]],
				max(ages), 
				labels = names(xx)[order(xx)][!tip.Group[order(xx)]], 
				adj = -0, srt = 90,cex = tip.cex, col = tip.color[2])
			text(
				sort(xx+name.jitter)[tip.Group[order(xx)]],
				max(ages),
				labels = names(xx)[order(xx)][tip.Group[order(xx)]],
				adj = -0, srt = 90, cex = tip.cex, col = tip.color[1])                            
    	}
    }
    on.exit(par(tg))
}


#' Functional phylogenetic diversity generalised linear model
#'
#' Calculate a generalised linear model using functional-phylogenetic 
#' diversity indices as predictors, given a value for the phylogenetic 
#' weighting parameter, \code{a}.
#'
#' @param tab The community by species presence table.
#' @param y An ecosystem function (or other site characteristic being
#'	used as a response variable).
#' @param PDist A phylogenetic distance matrix.
#' @param FDist A functional distance matrix.
#' @param a A number between 0 and 1 giving the amount of weight
#'	to put on \code{PDist} relative to \code{FDist}.
#' @param p A number giving the \code{p}-norm.
#' @param index Character string indicating the diversity index used 
#'	to combine distances and community data. Currently, either \code{'mpd'}
#'	or \code{'rao'}.
#' @param abundance.weighted Should mean pairwise distances be weighted
#'  by species abundance? (default = FALSE).  Only relevant if 
#'  \code{index = 'mpd'} and \code{tab} is an abundance matrix.
#' @param ... Additional arguments to pass to \code{\link{glm}}.
#' @return A \code{\link{glm}} object.
#' @export
FPDglm <- function(tab, y, PDist, FDist, a, p = 2, index = 'mpd', 
	abundance.weighted = FALSE, ...){
	# TODO: make the 'index' argument more generalizable by
	# allowing user supplied functions
	FPDist.temp <- FPDist(PDist, FDist, a, p)
	if(index == 'rao') fpd <- rao(tab, FPDist.temp)
	else if(index == 'mpd') fpd <- mpd.(tab, FPDist.temp, abundance.weighted = abundance.weighted)
	else stop('index not recognised')
	glm(y ~ fpd, ...)
}


FPDglm_ap <- function(ap, abundance.weighted = FALSE, ...){
	FPDglm(a = ap[1], p = ap[2], abundance.weighted = abundance.weighted, ...)
}


#' Calculate  Grid search functional phylogenetic diversity generalised linear model
#'
#' Calculate deviances over a grid of the phylogenetic weighting parameter, a,
#' for generalised linear models using functional-phylogenetic diversity 
#' indices as predictors.
#'
#' @param tab The community by species presence table.
#' @param y An ecosystem function (or other site characteristic being
#'	used as a response variable).
#' @param PDist A phylogenetic distance matrix.
#' @param FDist A functional distance matrix.
#' @param a A vector of numbers between 0 and 1 giving the amount of weight
#'	to put on \code{PDist} relative to \code{FDist}.  It is best if a is an
#'  evenly-spaced grid.
#' @param p A number giving the \code{p}-norm.
#' @param index Character string indicating the diversity index used 
#'	to combine distances and community data. Currently, either \code{'mpd'}
#'	or \code{'rao'}.
#' @param abundance.weighted Should mean pairwise distances be weighted
#'  by species abundance? (default = FALSE).  Only relevant if 
#'  \code{index = 'mpd'} and \code{tab} is an abundance matrix.
#' @param saveGLMs Should GLMs for every point along the grid be saved as
#'      an attribute, \code{glms}?
#' @param ... Additional arguments to pass to \code{\link{FPDglm}}.
#' @return TODO.
#' @export
FPDglm_grid <- function(tab, y, PDist, FDist, a = seq(0, 1, 0.01), 
                        p = 2, index = 'mpd', abundance.weighted = FALSE,
                        saveGLMs = FALSE, ...){
	
    aps <- merge(a, p)
    glms <- lapply(as.data.frame(t(aps)),
                   FPDglm_ap,
                   tab = tab, y = y, 
                   PDist = PDist, FDist = FDist,
                   index = index, 
                   abundance.weighted = abundance.weighted, ...)
    loglikes <- sapply(glms, logLik)
    slopes <- sapply(glms, coef)[2]
    likelihood <- exp(loglikes - max(loglikes))
    delta <- mean(diff(a))
    posterior <- likelihood/(sum(delta*likelihood))
    mean.a <- delta * sum(posterior * a)
    mode.a <- a[which.max(posterior)]
    var.a <- delta * sum(posterior * ((a - mean.a)^2))
    surf <- data.frame(
        a = aps$x, p = aps$y, # sorry for the x and y -- artifact of default merge behaviour
        loglikes, slopes, posterior)
    out <- list(surf = surf, delta = delta,
                mean.a = mean.a, mode.a = mode.a, var.a = var.a)
    if(saveGLMs) attr(out, "glms") <- glms
    class(out) <- 'FPDglm_grid'
    return(out)
}

# TODO:  document this method!
#'@export
plot.FPDglm_grid <- function(x, y, ...){
	## xx <- x$surf$a
	## yy <- x$surf$posterior

	## plot(xx, yy, type = 'l', las = 1,
	## 	ylab = 'Posterior density',
	## 	xlab = 'Phylogenetic weighting parameter, a',
	## 	...)
	## rug(xx)
    xlab <- expression(paste('Phylogenetic weighting parameter, ', italic(a)))
    ylab <- 'Approximate posterior density'
    ylim <- c(0, max(x$surf$posterior))
    hpd <- with(x$surf, a.hpd(a, posterior))
    par(mar = c(5, 5, 1, 1))
    plot(posterior ~ a, data = x$surf,
         las = 1, type = 'n',
         xlab = xlab, ylab = ylab, ylim = ylim)
    with(hpd, polygon(c(a, a[length(a):1]), c(posterior, rep(0, length(a))),
                      col = grey(0.9), border = NA))
    lines(posterior ~ a, data = x$surf)
}


#' Highest posterior density region for a
#'
#' Find points on a grid within a 100\code{p}\% highest posterior 
#' density region for the tuning parameter, a
#'
#' The 100\code{p}\% highest posterior density region for 'a' is 
#' the subset of the interval between 0 and 1, which contains 
#' 100\code{p}\% of the probability.
#'
#' @param a_grid A vector of (preferably evenly spaced) 'a' 
#'	values (between 0 and 1).
#' @param posterior The values of the posterior density at
#'	each point in \code{a_grid}.
#' @param level Size of the highest posterior density region.
#' @return A data frame with two columns:  the values of the
#'	grid within the hpd region and the value of the posterior
#'	at each point in this grid.
#' @export
a.hpd <- function(a_grid, posterior, level = 0.95){
	out <- data.frame(a = a_grid, posterior = posterior)
	if(level == 1L) return(out)
	if(level == 0L) return(out[-(1:length(posterior)),])
	dscrt.post <- posterior/sum(posterior)
	post.ord <- order(-posterior)
	hpd.levels <- rep(0, length(posterior))
	for(n in 1:length(posterior)){
		hpd.levels[n] <- sum(dscrt.post[post.ord[1:n]])
		if(hpd.levels[n] > level) break
	}
	n <- n - 1
	post.ord <- sort(post.ord[1:n])
	hpd.points <- a_grid[post.ord]
	out <- out[post.ord, ]
	print(hpd.levels[n])
	return(out)
}


a.postsim <- function(a_grid, posterior, n = 1){
	dscrt.post <- posterior/sum(posterior)
	replicate(n, a_grid[min(which(runif(1) < cumsum(dscrt.post)))])
}

#' Simulate phylogenetically correlated trait data
#' 
#' @param m Number of species.
#' @param q Number of relevant traits.
#' @param brown.sd Standard deviation of Brownian motion.
#' @param diverge.sd Standard deviation of noise due to divergence from 
#'  Brownian motion.
#' @export
#' @return A list with the following components:
#'  \item{tree}{phylogenetic tree}
#'  \item{traits}{species-by-traits matrix}
#'  \item{FDist}{functional distance matrix}
#'  \item{PDist}{phylogenetic distance matrix}
fd.pd.sim <- function(m, q, brown.sd, diverge.sd){

	# simulate tree
	tr <- rcoal(m)
	
	# simulate traits 
	# (1st step = brownian motion)
	sim.traits <- t(chol(vcv(tr))) %*% matrix(rnorm(q*m, sd = brown.sd), m, q)
	# (2nd step = noise due to divergence)
	sim.traits <- sim.traits + matrix(rnorm(m*q, sd = diverge.sd), m, q)
	
	# calculate distance matrices
	FDist <- as.matrix(dist(sim.traits))
	PDist <- as.matrix(as.dist(cophenetic(tr)))
	
	# standardize distance matrices
	FDist <- FDist / max(FDist)
	PDist <- PDist / max(PDist)
	
	out <- list(
		tree = tr,
		traits = sim.traits,
		FDist = FDist,
		PDist = PDist
	)
	return(out)
}


#' Performance of FPDist in simulations
#'
#' @param a Phylogenetic weighting parameter.
#' @param q.obs Number of observed traits.
#' @param object Output from \code{fd.pd.sim}.
#' @param plot Should plot results?
#' @param ... Additional objects to pass to \code{\link{cor}}.
#' @export
#' @return correlation between true and estimated functional distances.
fd.est.sim <- function(a, q.obs, object, plot = FALSE,...){
	
	# calculate FDist matrix for a reduced set of 'observed' traits
	FDist.obs <- as.matrix(dist(object$traits[, 1:q.obs]))
	
	# standardize
	FDist.obs <- FDist.obs/max(FDist.obs)

	# calculate FPDist matrix
	# (ord = FALSE is important!  otherwise row and col names don't
	# match between distance matrices, thereby throwing off the 
	# correlation computed below.)
	FDist.est <- FPDist(object$PDist, FDist.obs, a, 2, ord = FALSE)
	
	# get one triangle of the distance matrices as a numeric vector,
	# so that true distances (object$FDist) can be correlated with
	# the estimated distances (FDist.est)
	y <- as.numeric(as.dist(object$FDist))
	x <- as.numeric(as.dist(FDist.est))
	
	if(plot) plot(x, y, xlab = '', ylab = '', las = 1, tck = 0.02, bty = 'n', mgp = c(3, 0.5, 0))
	
	# return the correlation
	return(cor(x, y, ...))
}


#' Performance of FPDist over a range of simulations
#'
#' @param m Number of species.
#' @param q Number of traits.
#' @param q.step.size Coarseness of the number of observed traits grid.
#' @param a.grid.size Coarseness of the a-value grid.
#' @export
#' @return A 3d array with the following dimensions:  
#'  (grid of a-values)-by-(num obs traits)-by-(trait divergence)
fd.est.sim.cor <- function(m, q, q.step.size = 2, a.grid.size = 101){
	
	# create a grid of a-values
	a.grid <- seq(0, 1, 1/(a.grid.size-1))

	# proportion of variance attributable to divergence, relative to
	# brownian motion.  create a grid of such proportions.
	varprop <- seq(0, 1, 0.1)

	# create a grid of the number of observed traits.
	q.obs <- seq(1, q-1, 2)

	# create a data frame of all possible crosses of a.grid and q.obs.
	df <- merge(a.grid, q.obs)
	names(df) <- c('a', 'q.obs')

	# allocate a list to store correlations
	cors <- list()
	
	# loop over each varprop
	for(i in seq_along(varprop)){
		
		# assign variances to brownian motion and divergence
		brown.sd <- sqrt(1-varprop[i])
		diverge.sd <- sqrt(varprop[i])
		
		# simulate data
		sims <- fd.pd.sim(m, q, brown.sd, diverge.sd)
		
		# calculate the correlations by applying fd.est.sim
		cors[[i]] <- mapply(fd.est.sim, df$a, df$q.obs,
			MoreArgs = list(object = sims))
		
		# correct the dimensions of the cors object
		dim(cors[[i]]) <- c(a.grid.size, length(q.obs))
	}
	
	return(simplify2array(cors))
}

rposta <- function(n, surf)
  sample(surf$a, n, TRUE, surf$posterior)
