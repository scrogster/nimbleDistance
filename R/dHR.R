#' Hazard rate distribution for use in \code{nimble} models
#'
#' \code{dHR} and \code{dHR_V} provide hazard-rate detection
#' distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dHR dHR_V
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
#' dHR(x=20, b=1, sigma=40, Xmax=100, exactint = 0)
#' dHR(x=20, b=1, sigma=40, Xmax=100, exactint = 1)
#'
#'
#' #vectorized versions
#' dHR_V(x=c(20, 10), b=1, sigma=40, Xmax=100, exactint = 0)
#' dHR_V(x=c(20, 10), b=1, sigma=40, Xmax=100, exactint = 1)
#
#
#helper functions for computing the integrals of the Hazard Rate distance function
integralHR <- function(b=5, sigma=30,Xmin,Xmax){
	F <- function(x){1-exp(-(x/sigma)^(-b))} #Expression of params
	return(integrate(F,Xmin,Xmax, rel.tol = .Machine$double.eps^0.5)[[1]])
}
RintegralHR <- nimbleRcall(function(b=double(0), sigma=double(0), Xmin=double(0),
																		Xmax=double(0)){}, Rfun='integralHR',
													 returnType = double(0))

#incomplete gamma function required for exact integral
inc_gamma <- function(a=1, x=5){
	return(gamma_inc(a, x))
}
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
		if(exactint==1) { integral = Xmax - (sigma*Rmy_gamma_inc(  -b^(-1), (Xmax/sigma)^(-b)) /b)} else
		                {integral = RintegralHR(b, sigma, 0, Xmax) }
		# evaluate hazard rate function at x
		p <- 1-exp(-(x/sigma)^(-b))
		L<-p/integral
		LL<-log(L)
		if(log) return(LL)
		else return(L)
	}
)

#########################################################################################
#vectorized versions of the functions
#' @export
dHR_V <- nimbleFunction(
	run = function(x = double(1), b = double(0),
								 sigma = double(0), Xmax = double(0, default=100),
								 exactint=integer(0, default=0), log = integer(0, default = 0)) {
		returnType(double(0))
		#integral to normalize the likelihood (two approaches)
		#1. analytic integral, using the incomplete gamma function
		#2. numeric integral using R's numerical integration function
		if(exactint==1) { integral = Xmax - (sigma*Rmy_gamma_inc(  -b^(-1), (Xmax/sigma)^(-b)) /b)} else
	                	{integral = RintegralHR(b, sigma, 0, Xmax) }
		# evaluate hazard rate function at x
		p <- 1-exp(-(x/sigma)^(-b))
		L<-p/integral
		LL<-sum(log(L))
		if(log) return(LL)
		else return(exp(LL))
	}
)


#rHR_V <- nimbleFunction(
#	run = function(n = integer(), b = double(0),
#								 sigma = double(0), Xmax = double(0)) {
#		returnType(double(0))
#		return(0)
#	}
#)
