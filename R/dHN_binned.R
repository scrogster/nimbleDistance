#' Half-normal distribution for use in \code{nimble} models
#'
#' \code{dHN_binned_V} provide a half-normal detection
#' distribution that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dHN_binned dHN_binned_V
#'
#' @name dHN_binned
#' @param x distance observations. either a single value (dHN_binned) or a vector of values (dHN_binned_V)
#' @param sigma scale of the half-normal distribution
#' @param Xmax right truncation distance for integration of the likelihood function
#' @param breaks vector of breakpoints for binning the data
#' @param log if TRUE, return the log-likelihood
#' @param n   number of random values to generate
#'
#' @author Michael Scroggie
#'
#' @keywords half normal.
#' @examples
#' #direct invocation of functions from R to evaluate the likelihood
#' dHN_binned_V(x=c(20, 21, 41, 66),  sigma=40, breaks=seq(0, 100, by=25), Xmax=100)
#'N=1000
#'true_y<-runif(N, 0, 100)
#'#half normal truth
#'sigma_true<-40
#'Xmax=60
#'#half-normal detection
#'p_detect<- exp(-true_y^2/(2*sigma_true^2))
#'plot(p_detect~true_y, type="p")
#'detect<-rbinom(N, 1, p_detect)
#'#observations
#'y<-true_y[detect==1 & true_y<Xmax]
#'nind<-length(y)
#'hist(y, breaks=20)
#'	distCodeV<-nimbleCode({
#'  y[1:nind]~dHN_V(sigma, Xmax=100, point=0)
#'  sigma~dunif(10, 70)
#'  })
#'
#'#inits and monitors
#'inits <- function() list(sigma=50)
#'params <- c("sigma")
#'
#'#generate some MCMC samples
#'samples <- nimbleMCMC(
#'	code = distCodeV,
#'	data=list(y=y, br=seq(0, 100, by=20)),
#'	constants = list(nind=nind),
#'	inits = inits,
#'	monitors = params,
#'	niter = 5000,
#'	nburnin = 0)
#'
#'plot(samples[,"sigma"], col="red", type="l")
#'abline(h=sigma_true)
#'
#' @rdname dHN_binned
#' @export
dHN_binned<- nimbleFunction(
	#just use lines rather than transect as POC.
	run = function(x = double(0),  #scalar distance observation
								 sigma = double(0), #scalar sigma
								 Xmax = double(0, default=100),  #right truncation distance
								 breaks=double(1),  #vector of breaks
								 log = logical(0, default = 0)) {
		returnType(double(0))
		#check consistency of breakpoints
		if(x<0){nimStop("Error: Distance less than zero")}
		if(max(breaks)>Xmax){nimStop("Error: break greater than Xmax")}
		if(min(breaks)!=0){nimStop("Error: expect minimum break at zero")}
		#compute bin totals
		nbreaks=length(breaks)
		nbins=nbreaks-1
		xtots=numeric(nbins)
		xtots[1:nbins]<-0
		nind=1
		#compute bin totals
		for(i in 1:nbins){
			xtots[i]<-x>breaks[i] & x<breaks[i+1]
		}
		#precompute the bin probabilities conditional on sigma and breaks
		p=numeric(nbins) #store probabilities
		integral= 2*(pnorm(Xmax, 0, sigma)-0.5) * sqrt(pi/2) * sigma #total area between 0 and Xmax
		for(i in 1:nbins){
			A1<-2*(pnorm(breaks[i], 0, sigma)-0.5) * sqrt(pi/2) * sigma
			A2<-2*(pnorm(breaks[i+1], 0, sigma)-0.5) * sqrt(pi/2) * sigma
			p[i]<-(A2-A1)/integral
		}
		# #multinomial log-likelihood for the observed bin counts
		LL=dmulti(x=xtots, size=1, prob=p, log=TRUE)
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
		if(max(breaks)>Xmax){nimStop("Error: break greater than Xmax")}
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
		integral= 2*(pnorm(Xmax, 0, sigma)-0.5) * sqrt(pi/2) * sigma #total area between 0 and Xmax
		for(i in 1:nbins){
			A1<-2*(pnorm(breaks[i], 0, sigma)-0.5) * sqrt(pi/2) * sigma
			A2<-2*(pnorm(breaks[i+1], 0, sigma)-0.5) * sqrt(pi/2) * sigma
			p[i]<-(A2-A1)/integral
		}
		# #multinomial log-likelihood for the observed bin counts
		LL=dmulti(x=xtots, size=sum(xtots), prob=p, log=TRUE)
		if(log){ return(LL) } else return(exp(LL))
	}
)

#' @rdname dHN_binned
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

