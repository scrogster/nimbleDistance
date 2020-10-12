#' Hazard rate distribution for use in \code{nimble} models
#'
#' \code{dHR} and \code{dHR_V} provide hazard-rate detection
#' distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dHR dHR_V rHR rHR_V
#'
#' @name dHR
#' @param x distance data. either a single value (dHR) or a vector of values (dHR_V)
#' @param b shape of the hazard rate function
#' @param sigma scale of the hazard rate function
#' @param Xmax truncation distance for integration of the likelihood function
#' @param log if TRUE, return the log-likelihood
#'
#' @author Michael Scroggie
#'
#' @keywords hazard rate
#' @examples
#'
#' #direct invocation of functions from R to evaluate the likelihood
#' dHR(x=20, b=1, sigma=40, Xmax=100, exactint = 0)
#' dHR(x=20, b=1, sigma=40, Xmax=100, exactint = 1)
#'
#'
#' #vectorized versions
#' dHR_V(x=c(20, 10), b=1, sigma=40, Xmax=100, exactint = 0)
#' dHR_V(x=c(20, 10), b=1, sigma=40, Xmax=100, exactint = 1)
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
#'nind<-length(y)
#'hist(y, breaks=20)
#'distCodeV<-nimbleCode({
#'	y[1:nind]~dHR_V(b, sigma, 100)
#'	sigma~dunif(10, 70)
#'	b~dunif(2, 20)
#'	p_mean<-RintegralHR(b, sigma, 0, 100)/100
#'	})
#'	distCodeV<-nimbleCode({
#'  y[]~dHR_V(b, sigma, 100)
#'  sigma~dunif(10, 70)
#'  b~dunif(2, 20)
#'  p_mean<-RintegralHR(b, sigma, 0, 100)/100
#'  })
#'#inits and monitors
#'inits <- function() list(sigma=50, b=5)
#'params <- c("sigma", "b", "p_mean)
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
#'layout(matrix(1:3, nrow=3))
#'hist(samples[,"sigma"], col="red", xlab="sigma")
#'abline(v=sigma_true, col="blue", lwd=2)
#'hist(samples[,"b"], col="red", xlab="b")
#'abline(v=b_true, col="blue", lwd=2)
#'hist(samples[,"p_mean"], col="red", xlab="b", xlim=c(0, 1))
#'abline(v=RintegralHR(b_true, sigma_true, 0, 100)/100, col="blue", lwd=2)

#helper functions for computing the integrals of the Hazard Rate distance function
#' @export
integralHR <- function(b=5, sigma=30,Xmin,Xmax){
	F <- function(x){1-exp(-(x/sigma)^(-b))} #Expression of params
	return(integrate(F,Xmin,Xmax, rel.tol = .Machine$double.eps^0.5)[[1]])
}
#' @export
RintegralHR <- nimbleRcall(function(b=double(0), sigma=double(0), Xmin=double(0),
																		Xmax=double(0)){}, Rfun='integralHR',
													 returnType = double(0))

#' @rdname dHR
#' @export
dHR <- nimbleFunction(
	run = function(x = double(0),
								 b = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 log = logical(0, default = 0)) {
		returnType(double(0))
    integral = RintegralHR(b, sigma, 0, Xmax)
		p <- 1-exp(-(x/sigma)^(-b))
		L<-p/integral
		LL<-log(L)
		if(log) return(LL)
		else return(L)
	}
)

#' @rdname dHR
#' @export
rHR<- nimbleFunction(
	run = function(n = integer(),
								 b = double(0),
								 sigma = double(0),
								 Xmax  = double(0, default=100)) {
		returnType(double(0))
		return(0)
	}
)

#########################################################################################
#vectorized versions of dHR computes likelihood for a set of distances with common parameters
# much faster
#' @export
dHR_V <- nimbleFunction(
	run = function(x = double(1),
								 b = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 log = integer(0, default = 0)) {
		returnType(double(0))
    integral = RintegralHR(b, sigma, 0, Xmax)
		# evaluate hazard rate function at x
		p <- 1-exp(-(x/sigma)^(-b))
		L<-p/integral
		LL<-sum(log(L))
		if(log) return(LL)
		else return(exp(LL))
	}
)
#' @export
rHR_V<- nimbleFunction(
	run = function(n = integer(),
								 b = double(0),
								 sigma = double(0),
								 Xmax  = double(0, default=100)) {
		returnType(double())
		return(0)
	}
)

