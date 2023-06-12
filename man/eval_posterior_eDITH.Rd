\name{eval_posterior_eDITH}
\alias{eval_posterior_eDITH}

\title{Evaluate posterior predictions from an eDITH run}

\description{
Add description
}

\usage{
eval_posterior_eDITH(x, river)
}

\arguments{
  \item{x}{List as produced by \code{\link{run_eDITH_BT}}}.
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
}
\details{
Add details.
}
\value{
The output list copies all objects of the input \code{x} list. The following objects are added:
\item{p_median}{Median of the posterior distribution of production rates.}
\item{C_median}{Median of the posterior distribution of eDNA values (concentrations or read numbers).}
\item{probDet_median}{Median of the posterior distribution of detection probability.}
\item{p_mean}{Mean of the posterior distribution of production rates.}
\item{C_mean}{Mean of the posterior distribution of eDNA values (concentrations or read numbers).}
\item{probDet_mean}{Mean of the posterior distribution of detection probability.}
All of these objects are vectors of length \code{river$AG$nNodes}. However, if a custom likelihood was used 
in \code{run_eDITH_BT}, then \code{probDet_median} and \code{probDet_mean} are not evaluated, and are 
replaced by a vector of zero length.
}



