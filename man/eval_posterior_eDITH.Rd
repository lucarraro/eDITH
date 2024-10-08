\name{eval_posterior_eDITH}
\alias{eval_posterior_eDITH}

\title{Evaluate posterior predictions from an eDITH run}

\description{
Function that evaluates relevant quantities from a posterior sample of the parameters of an eDITH model
}

\usage{
eval_posterior_eDITH(x, river, quant = 0.5)
}

\arguments{
  \item{x}{List as produced by \code{\link{run_eDITH_BT}}.}
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{quant}{Vector of quantiles.}
}
\value{
The output list copies all objects of the input \code{x} list. The following objects are added:
\item{param_quantile}{Selected quantiles (along rows) of the posterior distribution of model parameters.}
\item{p_quantile}{Selected quantiles (along rows) of the posterior distribution of production rates.}
\item{C_quantile}{Selected quantiles (along rows) of the posterior distribution of eDNA values 
(concentrations or read numbers).}
\item{probDet_quantile}{Selected quantiles (along rows) of the posterior distribution of detection probability.}
\item{param_mean}{Mean of the posterior distribution of model parameters.}
\item{p_mean}{Mean of the posterior distribution of production rates.}
\item{C_mean}{Mean of the posterior distribution of eDNA values (concentrations or read numbers).}
\item{probDet_mean}{Mean of the posterior distribution of detection probability.}
All of these objects are vectors of length \code{river$AG$nNodes}. However, if a custom likelihood was used 
in \code{run_eDITH_BT}, then \code{probDet_quantile} and \code{probDet_mean} are not evaluated, and are 
replaced by a vector of zero length.
}

\examples{
library(rivnet)
data(wigger)
data(outSample)
out <- eval_posterior_eDITH(outSample, wigger)
plot(wigger, out$p_mean)
}

