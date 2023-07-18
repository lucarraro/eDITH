\name{posterior_pred_sim_eDITH}
\alias{posterior_pred_sim_eDITH}

\title{Predictive posterior simulations from an eDITH run}

\description{
This function performs predictive posterior simulations from a run of the eDITH model (via \code{run_eDITH_BT}).
These can be used for diagnostics purposes, in particular to assess scaled (quantile) residuals via the 
\code{DHARMa} package.
}

\usage{
posterior_pred_sim_eDITH(x, river, nParamSets = 10000, nDrawsPerParamSet = 10, 
	verbose = FALSE)
}

\arguments{
  \item{x}{List as produced by \code{\link{run_eDITH_BT}}.}
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{nParamSets}{Number of unique parameter sets sampled from the posterior distribution.}
  \item{nDrawsPerParamSet}{Number of simulations run per parameter set.}
  \item{verbose}{Logical. Should updates be printed on the console?}
}
\details{
\code{nParamSets} can be higher than the number of unique parameter sets in the posterior distribution, 
since the sampling of posterior parameter sets is operated with replacement.
}
\value{
A matrix with dimensions \code{length(x$data$ID)}-by-\code{nParamSets*nDrawsPerParamSet}. Each column is a 
predictive posterior simulation. Each row corresponds to a site where eDNA data were observed 
(corresponding to the entries of argument \code{data} in \code{\link{run_eDITH_BT}}. Matrix entries are
eDNA values (either concentrations or read numbers) predicted by the model for a given predictive posterior
simulation at a given observational site.
}
\seealso{\code{\link{DHARMa}}.}

\examples{
library(DHARMa)
data(outSample)
data(wigger)
data(dataC)
pps <- posterior_pred_sim_eDITH(outSample, wigger, nParamSets = 1000) 
# reduced nParamSets for illustrative purposes

sim.out <- createDHARMa(pps, dataC$values)
plot(sim.out)
}
