#' Half-normal distribution for use in \code{nimble} models
#'
#' \code{dHN} and \code{dHN_V} provide half-normal detection
#' distributions that can be used directly from R or in \code{nimble}
#' models. The corresponding randomisation routines \code{rHN} and \code{rHN_V} generate random distances
#' from the corresponding distance-detection distribution.
#'
#' @aliases dHN dHN_V rHN rHN_V
#'
#' @name dHN
#' @param x distance observations. either a single value (dHN) or a vector of values (dHN_V)
#' @param sigma scale of the half-normal distribution
#' @param Xmax right truncation distance for integration of the likelihood function
#' @param Xmin left truncation distance for integration of the likelihood function
#' @param point logical, if 1 compute likelihood for point transects, if 0 (default) compute likelihood for line transects
#' @param log if TRUE, return the log-likelihood
#'
#' @author Michael Scroggie
#'
#' @keywords hazard rate
#' @examples
#'
#' #direct invocation of functions from R to evaluate the likelihood
#' dHN(x=20,  sigma=40, Xmax=100)
#' dHN(x=20,  sigma=40, Xmax=100, point=1)
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
#'  y[1:nind]~dHN_V(sigma, 100, point=0)
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
#' @rdname dHN
#' @export
dHN <- nimbleFunction(
	run = function(x = double(0),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 point = double(0, default=0),
								 log = double(0, default = 0)) {
		returnType(double(0))
		#line transect likelihood
		if(point==0) {integral=pnorm(Xmax, 0, sigma)-0.5
		              g = dnorm(x, 0, sigma)
		              L=g/integral} else
		#point transect likelihood
		if(point==1) {integral <- sigma^2 * (1 - exp(-(Xmax^2)/(2*sigma^2)))
		             g <- exp(-(x^2)/(2*sigma^2))
		             L <- (x * g)/integral }
		LL<-log(L)
		if(log) return(LL)
		else return(L)
	}
)

#' @rdname dHN
#' @export
rHN<- nimbleFunction(
	run = function(n = integer(),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 point = double(0, default=0)) {
		returnType(double(0))
		out<-numeric(n)
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

#' @rdname dHN
#' @export
dHN_V <- nimbleFunction(
	run = function(x = double(1),
								 sigma = double(0),
								 Xmax = double(0, default=100),
								 point = double(0, default=0),
								 log = double(0, default = 0)) {
		returnType(double(0))
		#line transect likelihood
		if(point==0) {integral <- 2*(pnorm(Xmax, 0, sigma)-0.5) * sqrt(pi/2) * sigma
		g = exp(-x^2/(2*sigma^2))
		L=g/integral} else
			#point transect likelihood
			if(point==1) {integral <- sigma^2 * (1 - exp(-(Xmax^2)/(2*sigma^2)))
			g <- exp(-(x^2)/(2*sigma^2))
			L <- (x * g)/integral }
		LL<-sum(log(L))
		if(log) return(LL)
		else return(exp(LL))
	}
)

#' @rdname dHN
#' @export
rHN_V<- nimbleFunction(
	run = function(n = integer(),
								 sigma = double(0),
								 Xmax  = double(0, default=100),
								 point = double(0, default=0)) {
		returnType(double(1))
		out<-numeric(n)
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
