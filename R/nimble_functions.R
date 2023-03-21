# This file registers all distributions when the package is loaded.
.onAttach <- function(libname, pkgname) {

	packageStartupMessage("Registering user-defined distributions:\n
												dHR, dHN, dHR_binned, dHN_binned")
	#Hazard-rate-----------------------------------------------------------
	registerDistributions(list(
		dHR = list(BUGSdist="dHR(b, sigma, Xmax, point)",
							 pqAvail = FALSE,
							 range = c(0, Inf))))

	registerDistributions(list(
		dHR_V = list(BUGSdist="dHR_V(b, sigma, Xmax, point)",
									pqAvail = FALSE,
									types = c('value = double(1)', 'b = double(0)', 'sigma = double(0)',
														'Xmax = double(0)', 'point = double(0)'),
									range = c(0, Inf))))
	#Half-normal-----------------------------------------------------------
	registerDistributions(list(
			dHN = list(BUGSdist="dHN(sigma, Xmax, point)",
								 pqAvail = FALSE,
								 range = c(0, Inf))))

	registerDistributions(list(
		dHN_V = list(BUGSdist="dHN_V(sigma, Xmax, point)",
								 pqAvail = FALSE,
								 types = c('value = double(1)', 'sigma = double(0)',
								 					'Xmax = double(0)', 'point = double(0)'),
								 range = c(0, Inf))))
	#Hazard-rate binned-----------------------------------------------------
	registerDistributions(list(
		dHR_binned_V = list(
			BUGSdist = "dHR_binned_V(b, sigma, Xmax,  point, breaks)",
			types = c('value = double(1)', 'b = double(0)', 'sigma = double(0)',
								'Xmax = double(0)', 'point=double(0)', 'breaks=double(1)'),
			pqAvail = FALSE,
			range = c(0, Inf))
	))

	registerDistributions(list(
		dHR_binned = list(
			BUGSdist = "dHR_binned(b, sigma, Xmax,  point, breaks)",
			types = c('value = double(0)', 'b = double(0)', 'sigma = double(0)',
								'Xmax = double(0)', 'point=double(0)', 'breaks=double(1)'),
			pqAvail = FALSE,
			range = c(0, Inf))
	))

	#Half normal binned-----------------------------------------------------
	registerDistributions(list(
		dHN_binned_V = list(
			BUGSdist = "dHN_binned_V(sigma, Xmax, breaks)",
			types = c('value=double(1)',
								'sigma=double(0)',
								'b=double(0)',
								'Xmax = double(0)',
								'breaks = double(1)'),
			pqAvail = FALSE,
			range = c(0, Inf))
	))

}
