<img align="right" width="250" src="man/figures/logo.png">

# eDITH

An R package to model transport of environmental DNA in river networks

[![CRAN](http://www.r-pkg.org/badges/version/eDITH)](http://CRAN.R-project.org/package=eDITH)

## Overview

`eDITH` (**eD**NA **I**ntegrating **T**ransport and **H**ydrology) allows interpreting environmental DNA data collected from river networks. It implements the eDITH model, which couples a geomorphological and hydrological characterization of a catchment, eDNA transport and decay dynamics, and a species distribution model, and allows transforming pointwise eDNA
data collected at a catchment into predicted maps of taxon density. 

Features:

* It provides estimations of detection probability of a taxon across a whole catchment based on spatially replicated eDNA data 
* It can handle both eDNA concentration data (e.g., from qPCR) and metabarcoding (read counts) data
* Model fit can be performed via Bayesian techniques or optimization algorithms 
* Covariates can be specified by the user and/or selected by means of Asymmetric Eigenvector Maps (AEMs)
* An interface to the `DHARMa` package for residual diagnostics is provided

`eDITH` requires the use of river networks defined as `river` objects, which can be built via the `rivnet` package.

## Installing the package

`eDITH` can be installed from CRAN:

```
install.packages("eDITH")
```

The development version can be installed from Github:

```
devtools::install_github("lucarraro/eDITH")
```