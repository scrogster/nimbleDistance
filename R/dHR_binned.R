#' Hazard rate distribution for use in \code{nimble} models
#'
#' \code{dHR_binned} and \code{dHR_binned_V} provide hazard-rate detection
#' distributions for binned data that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dHR_binned dHR_binned_V
#'
#' @name dHR_binned
#' @param x distance data. either a single value (dHR) or a vector of values (dHR_V). As data are binned these
#' can be any value within the bounds of each bin (e.g. a midpoint).
#' @param b shape of the hazard rate function
#' @param sigma scale of the hazard rate function
#' @param Xmax right truncation distance for integration of the likelihood function
#' @param breaks vector of break points for binning. Typically between zero and Xmax.
#' @param point logical value indicating that point surveys were used rather than line-transects. Defaults to zero.
#' @param Xmin left truncation distance for integration of the likelihood function
#' @param log if TRUE, return the log-likelihood
#'
#' @author Michael Scroggie
#' @examples
#'
#' #direct invocation of functions from R to evaluate likelihoods
#' dHR_binned(x=20, b=1, sigma=40, Xmax=100, breaks=seq(0, 100, by=20),  point=0)
#' dHR_binned_V(x=c(20, 21, 41), b=1, sigma=40, Xmax=100, breaks=seq(0, 100, by=20),  point=0)
#'
#'N=500
#'true_y<-runif(N, 0, 100)
#'#hazard rate truth
#'sigma_true<-40
#'b_true<- 5
#'#detection
#'p_detect<- 1-exp(-(true_y/sigma_true)^(-b_true))
#'plot(p_detect~true_y, type="p")
#'detect<-rbinom(N, 1, p_detect)
#'#observations
#'y<-true_y[detect==1]
#'rug(y, side=1)
#'nind<-length(y)
#'hist(y, breaks=20)
#'
#'distCode_binned<-nimbleCode({
#'	y[1:nind]~dHR_binned_V(b=b, sigma=sigma, Xmax=100, breaks=br[], point=0)
#'	sigma~dunif(10, 70)
#'	b~dunif(2, 20)
#'})
#'#inits and monitors
#'inits <- function() list(sigma=50, b=5)
#'params <- c("sigma", "b")
#'
#'#generate some MCMC samples
#'samples <- nimbleMCMC(
#' code = distCode_binned,
#'	constants = list(y=y, nind=nind, br=c(0, 20, 40, 60, 80, 100)), ## the distance data
#'	inits = inits,
#'	monitors = params,
#'	niter = 1000,
#'	nburnin = 500)
#'
#'hist(samples[,"sigma"])
#'hist(samples[,"b"])
#' @rdname dHR_binned
#' @export
dHR_binned<- nimbleFunction(
	#just use lines rather than transect as POC.
	run = function(x = double(0),  #single distance observation
								 b = double(0), #scalar b
								 sigma = double(0), #scalar sigma
								 Xmax = double(0, default=100),  #right truncation distance
								 point=double(0, default = 0),
								 breaks=double(1),  #vector of breaks
								 log = double(0, default = 0)) {
		returnType(double(0))
		#check consistency of breakpoints
		if(min(x)<0){nimStop("Error: Distance less than zero")}
		if(max(breaks)<Xmax){nimStop("Error: break greater than Xmax")}
		if(min(breaks)!=0){nimStop("Error: expect minimum break at zero")}
		#compute bin totals
		nbreaks=length(breaks)
		nbins=nbreaks-1
		xtots=numeric(nbins)
		xtots[1:nbins]<-0
		#compute bin totals
		for(i in 1:nbins){
			xtots[i]<-x>=breaks[i] & x<breaks[i+1]
		}
		#compute the bin probabilities conditional on sigma, b and breaks
		#integralHR <- function(b=5, sigma=30,Xmin,Xmax, point=0){
		p=numeric(nbins) #store probabilities
		tot_area=RintegralHR(b=b, sigma=sigma, 0, Xmax, point=point)
		for(i in 1:nbins){
			p[i]=RintegralHR(b=b, sigma=sigma, breaks[i], breaks[i+1], point=point)/tot_area
		}
		# #multinomial log-likelihood for the observed bin counts
		LL=dmulti(x=xtots, size=1, prob=p, log=TRUE)
		#hack to deal with really small likelihood values
		if(LL==-Inf){LL<- -1000}
		if(log){ return(LL) } else return(exp(LL))
	}
)

#' #' @rdname dHR_binned
#' #' @export
 rHR_binned<- nimbleFunction(
 	run = function(n = integer(),
 								 b = double(0),
 								 sigma = double(0),
 								 Xmax  = double(0, default=100),
 								 point = double(0, default = 0),
 								 breaks=double(1) #vector of breaks
 								 ) {
 		returnType(double(0))
 		k<-0
 		while(k==0){
 			xrand<-runif(1, 0, Xmax)
 			p <- 1-exp(-(xrand/sigma)^(-b))
 			k=rbinom(1, 1, p)
 		}
 		return(xrand)
 	}
)


#'#' @rdname dHR_binned
#' @export
 dHR_binned_V<- nimbleFunction(
 	#just use lines rather than transect as POC.
 	run = function(x = double(1),  #vector of distance observation
 								 b = double(0),
 								 sigma = double(0), #scalar sigma
 								 Xmax = double(0, default=100),  #right truncation distance
 								 point=double(0, default=0),
 								 breaks=double(1),  #vector of breaks
 								 log = double(0, default = 0)) {
 		returnType(double(0))
 		#check consistency of breakpoints
 		if(min(x)<0){nimStop("Error: Distance less than zero")}
 		if(max(breaks)<Xmax){nimStop("Error: break greater than Xmax")}
 		if(min(breaks)!=0){nimStop("Error: expect minimum break at zero")}
 		#compute bin totals
 		nbreaks=length(breaks)
 		nbins=nbreaks-1
 		xtots=numeric(nbins)
 		xtots[1:nbins]<-0
 		nind=length(x)
 		#compute bin totals
 		for(i in 1:nbins){
 			xtots[i]<-sum(x>=breaks[i] & x<breaks[i+1])
 		}
 		#compute the bin probabilities conditional on sigma, b and breaks
 		#integralHR <- function(b=5, sigma=30,Xmin,Xmax, point=0){
 		p=numeric(nbins) #store probabilities
 		tot_area=RintegralHR(b=b, sigma=sigma, 0, Xmax, point=point)
 		for(i in 1:nbins){
 			p[i]=RintegralHR(b=b, sigma=sigma, breaks[i], breaks[i+1], point=point)/tot_area
 		}
 		# #multinomial log-likelihood for the observed bin counts
 		LL=dmulti(x=xtots, size=nind, prob=p, log=TRUE)
 		#hack to deal with really small likelihood values
 		if(LL==-Inf){LL<- -1000}
 		if(log){ return(LL) } else return(exp(LL))
 	}
 )
#' @rdname dHR_binned
#' @export
rHR_binned_V<- nimbleFunction(
	run = function(n = integer(),
								 b = double(0),
								 sigma = double(0), #scalar sigma
								 Xmax = double(0, default=100),  #right truncation distance
								 point=double(0, default=0),
								 breaks=double(1) #vector of breaks
								 ) {
		returnType(double(0))
		k<-0
		while(k==0){
			xrand<-runif(1, 0, Xmax)
			p <- 1-exp(-(xrand/sigma)^(-b))
			k=rbinom(1, 1, p)
		}
		return(xrand)
	}
)


