#' Hazard rate distribution for use in \code{nimble} models
#'
#' \code{dHR} and \code{dHR_V} provide hazard-rate detection
#' distributions that can be used directly from R or in \code{nimble}
#' models. The corresponding randomisation routines \code{rHR} and \code{rHR_V} generate random distances
#' from the corresponding hazard-rate distribution.
#' \code{integralHR} is a helper function for computing the integral for the hazard-rate function with the specified parameters.
#'
#' @aliases dHR dHR_V rHR rHR_V
#'
#' @name dHR
#' @param x distance data. either a single value (dHR) or a vector of values (dHR_V)
#' @param b shape of the hazard rate function
#' @param sigma scale of the hazard rate function
#' @param Xmax right truncation distance for integration of the likelihood function
#' @param point logical value indicating that point surveys were used rather than line-transects. Defaults to zero.
#' @param Xmin left truncation distance for integration of the likelihood function
#' @param log if TRUE, return the log-likelihood
#' @param n   number of random values to generate
#'
#' @author Michael Scroggie
#'
#' @keywords hazard rate
#' @examples
#'
#' #direct invocation of functions from R to evaluate the likelihood
#' dHR(x=20, b=1, sigma=40, Xmax=100)
#' dHR(x=20, b=1, sigma=40, Xmax=100, point=1)
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
#'	distCodeV<-nimbleCode({
#'  y[1:nind]~dHR_V(b, sigma, 100, point=0)
#'  sigma~dunif(10, 70)
#'  b~dunif(2, 20)
#'  p_mean<-RintegralHR(b, sigma, 0, 100, point=0)/100
#'  })
#'#inits and monitors
#'inits <- function() list(sigma=50, b=5)
#'params <- c("sigma", "b", "p_mean")
#'
#'#generate some MCMC samples
#'samples <- nimbleMCMC(
#'	code = distCodeV,
#'	constants = list(y=y, nind=nind), ## the distance data
#'	inits = inits,
#'	monitors = params,
#'	niter = 1000,
#'	nburnin = 500)
#'
#'hist(samples[,"sigma"], col="red", xlab="sigma")
#'abline(v=sigma_true, col="blue", lwd=2)
#'hist(samples[,"b"], col="red", xlab="b")
#'abline(v=b_true, col="blue", lwd=2)
#'hist(y, freq=FALSE)

#helper functions for computing the integrals of the Hazard Rate distance function
#' @export
integralHR <- function(b=5, sigma=30,Xmin,Xmax, point=0){
	Fline <- function(x){1-exp(-(x/sigma)^(-b))} #for line transects
	Fpoint <- function(x){x*(1-exp(-(x/sigma)^(-b)))} #for point transects
	if(point==1) {F=Fpoint} else {F=Fline}
	return(integrate(F,Xmin,Xmax, rel.tol = .Machine$double.eps^0.5)[[1]])
}
#' @export
RintegralHR <- nimbleRcall(function(b=double(0), sigma=double(0), Xmin=double(0),
																		Xmax=double(0), point=double(0)){}, Rfun='integralHR',
													 returnType = double(0))

#' @rdname dHR
#' @export
dHR <- nimbleFunction(
	run = function(x = double(0),
								 b = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 point = double(0, default = 0),
								 log = double(0, default = 0)) {
		returnType(double(0))
    integral = RintegralHR(b, sigma, 0, Xmax, point)
    if(point==0) {p<-1-exp(-(x/sigma)^(-b)) } else  #line transects
		             {p <- x*(1-exp(-(x/sigma)^(-b))) } #points
		L<-p/integral
		LL<-log(L)
		if(log) return(LL)
		else return(L)
	},
	buildDerivs=TRUE
)

#' @rdname dHR
#' @export
rHR<- nimbleFunction(
	run = function(n = integer(),
								 b = double(0),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 point = double(0, default = 0)) {
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

#########################################################################################
#vectorized versions of dHR computes likelihood for a set of distances with common parameters
# much faster
#' @rdname dHR
#' @export
dHR_V <- nimbleFunction(
	run = function(x = double(1),
								 b = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 point = double(0, default = 0),
								 log = double(0, default = 0)) {
		returnType(double(0))
    integral = RintegralHR(b, sigma, 0, Xmax, point)
		# evaluate hazard rate function at x
    if(point==0) {p<-1-exp(-(x/sigma)^(-b)) } else  #line transects
                 {p <- x*(1-exp(-(x/sigma)^(-b))) } #points
		L<-p/integral
		LL<-sum(log(L))
		if(log) return(LL)
		else return(exp(LL))
	},
	buildDerivs=TRUE
)
#' @rdname dHR
#' @export
rHR_V<- nimbleFunction(
	run = function(n = integer(),
								 b = double(0),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 point = logical(0, default = 0)) {
		returnType(double(1))
		out<-numeric(n)
		for(i in 1:n) {
		k<-0
		while(k==0){
			xrand<-runif(1, 0, Xmax)
			p <- 1-exp(-(xrand/sigma)^(-b))
			k=rbinom(1, 1, p)
		}
		out[i]<-xrand
		}
		return(out)
	}
)

