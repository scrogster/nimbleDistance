#' Half-normal distribution for use in \code{nimble} models
#'
#' \code{dHN_binned} and \code{dHN_binned_V} provide half-normal detection
#' distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dHN_binned dHN_binned_V
#'
#' @name dHN_binned
#' @param x distance observations. either a single value (dHN_binned) or a vector of values (dHN_binned_V)
#' @param sigma scale of the half-normal distribution
#' @param Xmax right truncation distance for integration of the likelihood function
#' @param Xmin left truncation distance for integration of the likelihood function
#' @param log if TRUE, return the log-likelihood
#'
#' @author Michael Scroggie
#'
#' @keywords hazard rate
#' @examples
#'
#' #direct invocation of functions from R to evaluate the likelihood
#' dHN_binned(x=20,  sigma=40, breaks=seq(0, 100, by=25), Xmax=100)
#'
#'N=1000
#'true_y<-runif(N, 0, 100)
#'#hazard rate truth
#'sigma_true<-40
#'#half-normal detection
#'p_detect<- exp(-true_y^2/sigma_true^2)
#'plot(p_detect~true_y, type="p")
#'detect<-rbinom(N, 1, p_detect)
#'#observations
#'y<-true_y[detect==1]
#'nind<-length(y)
#'hist(y, breaks=20)
#'	distCodeV<-nimbleCode({
#'  y[1:nind]~dHN_binned_V(sigma, 100, breaks=seq(0, 100, by=20))
#'  sigma~dunif(10, 70)
#'  })
#'#inits and monitors
#'inits <- function() list(sigma=50)
#'params <- c("sigma")
#'
#'#generate some MCMC samples
#'samples <- nimbleMCMC(
#'	code = distCodeV,
#'	constants = list(y=y, nind=nind), ## the distance data
#'	inits = inits,
#'	monitors = params,
#'	niter = 2000,
#'	nburnin = 1000)
#'
#'layout(matrix(1:2, nrow=2))
#'hist(y, freq=FALSE)
#'lines(x=1:100, y=dHN(1:100, mean(samples[,"sigma"])), col="red")
#'hist(samples[,"sigma"], col="red")
#'
#' @rdname dHN_binned
#' @export
dHN_binned_V<- nimbleFunction(
	#just use lines rather than transect as POC.
	run = function(x = double(1),  #vector of distance observation
								 sigma = double(0), #scalar sigma
								 Xmax = double(0, default=100),  #right truncation distance
								 breaks=double(1),  #vector of breaks
								 log = logical(0, default = 0)) {
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
		 	xtots[i]<-sum(x>breaks[i] & x<breaks[i+1])
		 }
		 #precompute the bin probabilities conditional on sigma and breaks
		   p=numeric(nbins) #store probabilities
		   tot_area=pnorm(Xmax, 0, sigma)-pnorm(breaks[1], 0, sigma)
		 for(i in 1:nbins){
		 	p[i]=(pnorm(breaks[i+1], 0, sigma)-pnorm(breaks[i], 0, sigma))/tot_area
		 }
		# #multinomial log-likelihood for the observed bin counts
		 LL=dmulti(x=xtots, size=nind, prob=p, log=TRUE)
		 #hack to deal with really small likelihood values
     if(LL==-Inf){LL<- -745}
		 if(log){ return(LL) } else return(exp(LL))
	}
)
#' @rdname dHN_binned
#' @export
dHN_binned<- nimbleFunction(
	#just use lines rather than transect as POC.
	run = function(x = double(0),  #single distance observation
								 sigma = double(0), #scalar sigma
								 Xmax = double(0, default=100),  #right truncation distance
								 breaks=double(1),  #vector of breaks
								 log = logical(0, default = 0)) {
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
			xtots[i]<-sum(x>breaks[i] & x<breaks[i+1])
		}
		#precompute the bin probabilities conditional on sigma and breaks
		p=numeric(nbins) #store probabilities
		tot_area=pnorm(Xmax, 0, sigma)-pnorm(breaks[1], 0, sigma)
		for(i in 1:nbins){
			p[i]=(pnorm(breaks[i+1], 0, sigma)-pnorm(breaks[i], 0, sigma))/tot_area
		}
		# #multinomial log-likelihood for the observed bin counts
		LL=dmulti(x=xtots, size=nind, prob=p, log=TRUE)
		#hack to deal with really small likelihood values
		if(LL==-Inf){LL<- -745}
		if(log){ return(LL) } else return(exp(LL))
	}
)
