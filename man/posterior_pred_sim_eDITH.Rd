\name{posterior_pred_sim_eDITH}
\alias{posterior_pred_sim_eDITH}

\title{Perform posterior prediction simulations from an eDITH run}

\description{
Add description
}

\usage{
posterior_pred_sim_eDITH(x, river, nParamSets = 10000, nDrawsPerParamSet = 10)
}

\arguments{
  \item{x}{List as produced by \code{\link{run_eDITH_BT}}}.
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{nParamSets}{Number of unique parameter sets sampled from the posterior distribution.}
  \item{nDrawsPerParamSet}{Number of simulations run per parameter set.}
}
\details{
\code{nParamSets} can be higher than the number of unique parameter sets in the posterior distribution, 
since the sampling of posterior parameter sets is operated with replacement.
}
\value{
A matrix with dimensions \code{nParamSets*nDrawsPerParamSet}, \code{length(x$data$ID)}. Each row is a 
predictive posterior simulation. Each column corresponds to a site where eDNA data were observed 
(corresponding to the entries of argument \code{data} in \code{\link{run_eDITH_BT}}.
}



