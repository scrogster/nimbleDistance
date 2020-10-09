# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

	packageStartupMessage("Loading nimbleDistance. \nRegistering the following user-defined functions:\n ",
												"dHR", " dDHR_V")

	# Register the custom distributions
	suppressMessages({

	registerDistributions(list(
		dHR = list(BUGSdist="dHR(b, sigma, Xmax)",
								 pqAvail = FALSE,
								 range = c(0, Inf))), verbose=FALSE)

	registerDistributions(list(
		dHR_V = list(BUGSdist="dHR_V(b, sigma, Xmax)",
									pqAvail = FALSE,
									types = c('value = double(1)', 'b = double(0)', 'sigma = double(0)', 'Xmax = double(0)'),
									range = c(0, Inf))), verbose=FALSE)

	})
}
