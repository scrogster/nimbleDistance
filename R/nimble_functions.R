# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

	packageStartupMessage("Registering the following user-defined functions:\n ",
												"dHR", " dDHR_V", " dHN", " dHN_V")

	# Register the custom distributions
	suppressMessages({

	registerDistributions(list(
		dHR = list(BUGSdist="dHR(b, sigma, Xmax, point)",
		#					 Rdist = "rHR(b, sigma, Xmax)",
								 pqAvail = FALSE,
								 range = c(0, Inf))), verbose=FALSE)

	registerDistributions(list(
		dHR_V = list(BUGSdist="dHR_V(b, sigma, Xmax, point)",
		#						 Rdist="rHR_V(b, sigma, Xmax)",
									pqAvail = FALSE,
									types = c('value = double(1)', 'b = double(0)', 'sigma = double(0)',
														'Xmax = double(0)', 'point = double(0)'),
									range = c(0, Inf))), verbose=FALSE)

	registerDistributions(list(
		dHN = list(BUGSdist="dHN(sigma, Xmax, point)",
		#					 Rdist = "rHN(sigma, Xmax)",
							 pqAvail = FALSE,
							 range = c(0, Inf))), verbose=FALSE)

	registerDistributions(list(
		dHN_V = list(BUGSdist="dHN_V(sigma, Xmax, point)",
		#						 Rhist="rHN_V(sigma, Xmax)",
								 pqAvail = FALSE,
								 types = c('value = double(1)', 'sigma = double(0)',
								 					'Xmax = double(0)', 'point = double(0)'),
								 range = c(0, Inf))), verbose=FALSE)
	})
}
