<img src="inst/nimbleDistance.png"  align="right"  width="150">

## nimbleDistance: functions for fitting distance sampling models in [nimble](https://https://r-nimble.org/)

This repository provides some user-defined distributions that can be used to implement distance sampling models in nimble:

* likelihoods for common distance-detection distributions (hazard-rate and half-normal).

* a vignette illustrating use of the package to undertake a simple analysis of distance-detection data collected on line transects.

The package is at a very early stage of development. Suggestions for improvement and contributions of code are very welcome.

### Installation:

```
devtools::install_github("scrogster/nimbleDistance", build_vignettes = TRUE, INSTALL_opts = "--no-multiarch")
```
