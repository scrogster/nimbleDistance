#helper functions for computing the integral of the Hazard Rate distance function
integralHR <- function(b=5, sigma=30,Xmin,Xmax){
	F <- function(x){1-exp(-(x/sigma)^(-b))} #Expression of params
	return(integrate(F,Xmin,Xmax, rel.tol = .Machine$double.eps^0.5)[[1]])
}
RintegralHR <- nimbleRcall(function(b=double(0), sigma=double(0), Xmin=double(0),
																Xmax=double(0)){}, Rfun='integralHR',
											 returnType = double(0))

#incomplete gamma function required for exact integral
inc_gamma <- function(a=1, x=5){
	return(gsl::gamma_inc(a, x))
}
Rinc_gamma <- nimbleRcall(function(a=double(0), x=double(0)){}, Rfun='inc_gamma',
														 returnType = double(0))

dHR <- nimbleFunction(
	run = function(x = double(0), b = double(0),
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
		LL<-log(L)
		if(log) return(LL)
		else return(L)
	}
)

rHR <- nimbleFunction(
	run = function(n = integer(), b = double(0),
								 sigma = double(0), Xmax = double(0)) {
		returnType(double(0))
		return(0)
	}
)

#########################################################################################
#vectorized versions
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

rHR_V <- nimbleFunction(
	run = function(n = integer(), b = double(0),
								 sigma = double(0), Xmax = double(0)) {
		returnType(double(0))
		return(0)
	}
)



#test out both integration approach
dHR(20, b=1, sigma=40, Xmax=100, exactint = 0)
dHR(20, b=1, sigma=40, Xmax=100, exactint = 1)

#test vectorized version
dHR_V(c(20, 10), b=1, sigma=40, Xmax=100, exactint = 0)
dHR_V(c(20, 10), b=1, sigma=40, Xmax=100, exactint = 1)
