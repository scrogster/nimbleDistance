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
#' dHN_binned_V(x=c(20, 21, 41, 66),  sigma=40, breaks=seq(0, 100, by=25), Xmax=100)
#'N=1000
#'true_y<-runif(N, 0, 100)
#'#half normal truth
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
#'  y[1:nind]~dHN_binned_V(sigma, 100, breaks=br[1:6])
#'  sigma~dunif(10, 70)
#'  })
#'#inits and monitors
#'inits <- function() list(sigma=50)
#'params <- c("sigma")
#'
#'#generate some MCMC samples
#'samples <- nimbleMCMC(
#'	code = distCodeV,
#'	constants = list(y=y, nind=nind, br=seq(0, 100, by=20)), ## the distance data
#'	inits = inits,
#'	monitors = params,
#'	niter = 2000,
#'	nburnin = 1000)
#'
#'layout(matrix(1:2, nrow=2))
#'hist(y, freq=FALSE)
#'hist(samples[,"sigma"], col="red")
#'
#'
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

#' @rdname dHN_binned
#' @export
rHN_binned<- nimbleFunction(
	run = function(n = integer(),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 breaks=double(1)) #vector of breaks)
		{
		returnType(double(0))
		out<-numeric(n)
		point=0 #temporary hard coded until point added to main likelihood
		k<-0
		while(k==0){
			if(point==0) {xrand<-runif(1, 0, Xmax)} else
			{randpoint<-runif(2, -Xmax, Xmax)
			xrand<-sqrt(sum(randpoint^2))}
			p=exp(-(xrand^2)/(2*sigma^2))
			detect<-rbinom(1, 1, p)
			k= xrand*detect <Xmax & xrand*detect >0
		}
		return(xrand)
	}
)


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

#' @rdname dHN
#' @export
rHN_binned_V<- nimbleFunction(
	run = function(n = integer(),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 breaks=double(1)  #vector of breaks
								 ) {
		returnType(double(1))
		out<-numeric(n)
		point=0 #temporary hard coded until point added to main likelihood
		for(i in 1:n) {
			k<-0
			while(k==0){
				if(point==0) {xrand<-runif(1, 0, Xmax)} else
				{randpoint<-runif(2, -Xmax, Xmax)
				xrand<-sqrt(sum(randpoint^2))}
				p=exp(-(xrand^2)/(2*sigma^2))
				detect<-rbinom(1, 1, p)
				k= xrand*detect <Xmax & xrand*detect >0
			}
			out[i]<-xrand
		}
		return(out)
	}
)

