\name{eval_posterior_eDITH}
\alias{eval_posterior_eDITH}

\title{Evaluate posterior predictions from an eDITH run}

\description{
Add description
}

\usage{
eval_posterior_eDITH(x, river, covariates, Z.normalize=TRUE, source.area="AG")
}

\arguments{
  \item{x}{List as produced by \code{\link{run_eDITH_BT}}}.
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{covariates}{Data frame containing covariate values for all \code{river} reaches.}
  \item{Z.normalize}{Logical. Should covariates be Z-normalized?}
  \item{source.area}{Defines the extent of the source area of a node. Possible values are \code{"AG"} (if the source
  area is the reach surface, i.e. length*width), \code{"SC"} (if the source area is the subcatchment area), or,
  alternatively, a vector with length \code{river$AG$nodes}.}
}
\details{
The arguments of the type \code{*.prior} consist in the lists  of arguments required by \code{\link{dtrunc}} 
(except the first argument \code{x}).
}
\value{
The output list copies all objects of the input \code{x} list. The following objects are added:
\item{p_median}{Median of the posterior distribution of production rates.}
\item{C_median}{Median of the posterior distribution of eDNA values (concentrations or read numbers).}
\item{probDet_median}{Median of the posterior distribution of detection probability.}
\item{p_mean}{Mean of the posterior distribution of production rates.}
\item{C_mean}{Mean of the posterior distribution of eDNA values (concentrations or read numbers).}
\item{probDet_mean}{Mean of the posterior distribution of detection probability.}
All of these objects are vectors of length \code{river$AG$nNodes}.
}



