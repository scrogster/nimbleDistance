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
#' @param exactint Use the analytic integral for the likelihood function, or a numerical approximation (default).
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
#'#inits and monitors
#'inits <- function() list(sigma=50, b=5)
#'params <- c("sigma", "b")
#'
#'#generate some MCMC samples
#'samples <- nimbleMCMC(
#'	code = distCodeV,
#'	constants = list(y=y, nind=nind), ## the distance data
#'	inits = inits,
#'	monitors = params,
#'	niter = 10000,
#'	nburnin = 5000)
#'
#'layout(matrix(1:2, nrow=2))
#'hist(samples[,"sigma"], col="red", breaks=50, xlab="sigma")
#'abline(v=sigma_true, col="blue", lwd=2); mtext(side=1, at=sigma_true, text="Truth")
#'hist(samples[,"b"], col="red", xlab="b", breaks=50)
#'abline(v=b_true, col="blue", lwd=2); mtext(side=1, at=b_true, text="Truth")

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

#incomplete gamma function required for exact integral
#' @export
inc_gamma <- function(a=1, x=5){
	return(gamma_inc(a, x))
}
#' @export
Rinc_gamma <- nimbleRcall(function(a=double(0), x=double(0)){}, Rfun='inc_gamma',
													returnType = double(0))

#' @rdname dHR
#' @export
dHR <- nimbleFunction(
	run = function(x = double(0),
								 b = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 exactint=integer(0, default=0),
								 log = integer(0, default = 0)) {
		returnType(double(0))
		#integral to normalize the likelihood (two approaches)
		#1. analytic integral, using the incomplete gamma function
		#2. numeric integral using R's numerical integration function
		if(exactint==1) { integral = Xmax - (sigma*Rinc_gamma(  -b^(-1), (Xmax/sigma)^(-b)) /b)} else
		                {integral = RintegralHR(b, sigma, 0, Xmax) }
		# evaluate hazard rate function at x
		p <- 1-exp(-(x/sigma)^(-b))
		L<-p/integral
		LL<-log(L)
		if(log) return(LL)
		else return(L)
	}
)

rHR<- nimbleFunction(
	run = function(n = integer(),
								 b = double(0),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 exactint=integer(0, default=0),
								 log = integer(0, default = 0)) {
		returnType(double(0))
		return(0)
	}
)

#########################################################################################
#vectorized versions of the functions
#' @export
dHR_V <- nimbleFunction(
	run = function(x = double(1),
								 b = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 exactint=integer(0, default=0),
								 log = integer(0, default = 0)) {
		returnType(double(0))
		#integral to normalize the likelihood (two approaches)
		#1. analytic integral, using the incomplete gamma function
		#2. numeric integral using R's numerical integration function
		if(exactint==1) { integral = Xmax - (sigma*Rinc_gamma(  -b^(-1), (Xmax/sigma)^(-b)) /b)} else
	                	{integral = RintegralHR(b, sigma, 0, Xmax) }
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
								 Xmax  = double(0, default=100),
								 exactint=integer(0, default=0),
								 log = integer(0, default = 0)) {
		returnType(double())
		return(0)
	}
)

