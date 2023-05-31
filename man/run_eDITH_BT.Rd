\name{run_eDITH_BT}
\alias{run_eDITH_BT}

\title{Run eDITH}

\description{
Add description
}

\usage{
run_eDITH_BT(data, river, covariates, Z.normalize = TRUE, no.det = TRUE, ll.type = "norm", 
source.area="AG",
    mcmc.settings = NULL, likelihood = NULL, prior = NULL, sampler.type = "DREAMzs", 
	tau.prior = list(spec = "lnorm", a = 0, b = Inf, meanlog = log(5), 
		sd = sqrt(log(5) - log(4))), 
	log_p0.prior = list(spec = "unif", min = -20, max = 0), 
	beta.prior = list(spec = "norm", sd = 1), 
	sigma.prior = list(spec = "unif", min = 0, 
		max = max(data$values, na.rm = TRUE)), 
	omega.prior = list(spec = "unif", min = 1, 
		max = 10 * max(data$values, na.rm = TRUE)), 
	Cstar.prior = list(spec = "unif", min = 0, 
		max = max(data$values, na.rm = TRUE)))
}

\arguments{
  \item{data}{eDNA data. Data frame containing columns \code{ID} (index of the AG node/reach where 
  the eDNA sample was taken) and \code{values} (value of the eDNA measurement, expressed as concentration
  or number of reads). }
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{covariates}{Data frame containing covariate values for all \code{river} reaches.}
  \item{Z.normalize}{Logical. Should covariates be Z-normalized?}
  \item{no.det}{Logical. Should a probability of non-detection be included in the model?}
  \item{ll.type}{Character. String defining the error distribution used in the log-likelihood formulation. 
  Allowed values are \code{norm} (for normal distribution), \code{lnorm} (for lognormal distribution), 
  \code{nbinom} (for negative binomial distribution). The latter choice is suited when eDNA data are 
  expressed as read numbers, while \code{norm} and \code{lnorm} are better suited to eDNA concentrations. }
  \item{source.area}{Defines the extent of the source area of a node. Possible values are \code{"AG"} (if the source
  area is the reach surface, i.e. length*width), \code{"SC"} (if the source area is the subcatchment area), or,
  alternatively, a vector with length \code{river$AG$nodes}.}
  \item{mcmc.settings}{List. It is passed as argument \code{settings} in \code{\link{runMCMC}}. Default is 
  \code{settings <- list(iterations = 2.7e6, burnin=1.8e6, message = TRUE, thin = 10).}}
  \item{likelihood}{Likelihood function to be passed as \code{likelihood} argument to \code{\link{createBayesianSetup}}.
  If not specified, it is generated based on arguments \code{no.det} and \code{ll.type}.}
  \item{prior}{Prior function to be passed as \code{prior} argument to \code{\link{createBayesianSetup}}.
  If not specified, it is generated based on the \code{*.prior} arguments provided.}
  \item{sampler.type}{Character. It is passed as argument \code{sampler} in \code{\link{runMCMC}}.}
  \item{tau.prior}{List that defines the prior distribution for the decay time parameter \code{tau}. See details. }
  \item{log_p0.prior}{List that defines the prior distribution for the natural logarithm of the baseline production rate 
  \code{p0}.  See details.}
  \item{beta.prior}{List that defines the prior distribution for the covariate effects \code{beta}. See details. If a single 
  \code{spec} is provided, the same prior distribution is specified for all \code{beta} parameters. Alternatively, if 
  \code{spec} (and the other arguments, if provided) is a vector with length equal to the number of covariates included, 
  different prior distributions can be specified for the different \code{beta} parameters.}
  \item{sigma.prior}{List that defines the prior distribution for the standard deviation of the measurement error 
 when \code{ll.type}  is \code{"norm"} or \code{"lnorm"}. It is not used if \code{ll.type = "nbinom"}. See details.}
  \item{omega.prior}{List that defines the prior distribution for the overdispersion parameter \code{omega} of the 
  measurement error when \code{ll.type = "nbinom"}. It is not used if  \code{ll.type} is \code{"norm"} or \code{"lnorm"}.
  See details.}
  \item{Cstar.prior}{List that defines the prior distribution for the \code{Cstar} parameter controlling the probability
  of no detection. It is only used if \code{no.det = TRUE}. See details.}
}
\details{
The arguments of the type \code{*.prior} consist in the lists  of arguments required by \code{\link{dtrunc}} 
(except the first argument \code{x}).
}
\value{
A list with objects:
\item{param_map}{Vector of named parameters corresponding to the maximum a posteriori estimate. It is the
output of the call to \code{\link{MAP}}.}
\item{p_map}{Vector of best-fit eDNA production rates corresponding to the maximum a posteriori parameter 
estimate \code{param_map}. It has length equal to \code{river$AG$nNodes}.}
\item{C_map}{Vector of best-fit eDNA values (in the same unit as \code{data$values}, i.e. concentrations or read numbers) 
corresponding to the maximum a posteriori parameter  estimate \code{param_map}. It has length equal to \code{river$AG$nNodes}.}
\item{probDet_map}{Vector of best-fit detection probabilities corresponding to the maximum a posteriori
parameter estimate \code{param_map}. It has length equal to \code{river$AG$nNodes}}.
\item{cI}{Output of the call to \code{\link{getCredibleIntervals}}.}
\item{gD}{Output of the call to \code{\link{gelmanDiagnostics}}.}
\item{outMCMC}{Object of class \code{mcmcSampler} returned by the call to \code{\link{runMCMC}}.} 
Moreover, arguments \code{ll.type}, \code{no.det} and \code{data} are added to the list.
}



