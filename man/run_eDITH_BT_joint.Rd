\name{run_eDITH_BT_joint}
\alias{run_eDITH_BT_joint}

\title{Run eDITH with BayesianTools  based on joint eDNA and direct sampling data}

\description{
Function that runs a Bayesian sampler estimating parameters of an eDITH model fitted on both 
eDNA and direct sampling data.
}

\usage{
run_eDITH_BT_joint(data, river, covariates = NULL, Z.normalize = TRUE,
           use.AEM = FALSE, n.AEM = NULL, par.AEM = NULL,
           no.det = FALSE, ll.type = "norm", source.area = "AG",
           mcmc.settings = NULL, likelihood = NULL, prior = NULL, 
		   sampler.type = "DREAMzs",
           tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), 
		                    sd=sqrt(log(5)-log(4))),
           log_p0.prior = list(spec="unif",min=-20, max=0),
           beta.prior = list(spec="norm",sd=1),
           sigma.prior = list(spec="unif",min=0, 
			                  max=max(data$values[data$type=="e"], 
							  na.rm = TRUE)),
           omega.prior = list(spec="unif",min=1, 
			                  max=10*max(data$values[data$type=="e"], 
							  na.rm = TRUE)),
           Cstar.prior = list(spec="unif",min=0, 
			                  max=max(data$values[data$type=="e"], 
							  na.rm = TRUE)),
           omega_d.prior = list(spec="unif",min=1,
                                max=10*max(c(0.11, data$values[data$type=="d"]), 
								na.rm = TRUE)),
           alpha.prior = list(spec="unif", min=0, max=1e6),
           verbose = FALSE)
}

\arguments{
  \item{data}{eDNA and direct observation data. Data frame containing columns \code{ID} (index of the AG node/reach where 
  the sample was taken), \code{values} (value of the eDNA or direct measurement) and \code{type} 
  (equal to \code{"e"} for eDNA data and to \code{"d"} for direct observation data). eDNA values 
  are expressed as concentration or number of reads; direct observations are expressed as numbers of individuals. }
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{covariates}{Data frame containing covariate values for all \code{river} reaches. If \code{NULL} (default
  option), production rates are estimated via AEMs.}
  \item{Z.normalize}{Logical. Should covariates be Z-normalized?}
  \item{use.AEM}{Logical. Should eigenvectors based on AEMs be used as covariates? If \code{covariates = NULL}, it is 
  set to \code{TRUE}. If \code{TRUE} and \code{covariates} are provided, AEM eigenvectors are appended to the 
  \code{covariates} data frame.}
   \item{n.AEM}{Number of AEM eigenvectors (sorted by the decreasing respective eigenvalue) to be used as covariates. If 
   \code{par.AEM$moranI = TRUE}, this parameter is not used. Instead, the eigenvectors with significantly positive spatial
   autocorrelation are used as AEM covariates.}
  \item{par.AEM}{List of additional parameters that are passed to \code{\link{river_to_AEM}} for calculation of AEMs. 
  In particular, \code{par.AEM$moranI = TRUE} imposes the use of AEM covariates with significantly positive spatial 
  autocorrelation based on Moran's I statistic. }
  \item{no.det}{Logical. Should a probability of non-detection be included in the model?}
  \item{ll.type}{Character. String defining the error distribution used in the log-likelihood formulation. 
  Allowed values are \code{norm} (for normal distribution), \code{lnorm} (for lognormal distribution), 
  \code{nbinom} (for negative binomial distribution) and \code{geom} (for geometric distribution). The two latter choices 
  are suited when eDNA data are expressed as read numbers, while \code{norm} and \code{lnorm} are better suited 
  to eDNA concentrations. }
  \item{source.area}{Defines the extent of the source area of a node. Possible values are \code{"AG"} (if the source
  area is the reach surface, i.e. length*width), \code{"SC"} (if the source area is the subcatchment area), or,
  alternatively, a vector with length \code{river$AG$nodes}.}
  \item{mcmc.settings}{List. It is passed as argument \code{settings} in \code{\link{runMCMC}}. Default is 
  \code{list(iterations = 2.7e6, burnin=1.8e6, message = TRUE, thin = 10).}}
  \item{likelihood}{Likelihood function to be passed as \code{likelihood} argument to \code{\link{createBayesianSetup}}.
  If not specified, it is generated based on arguments \code{no.det} and \code{ll.type}. If a custom likelihood is specified, 
  a custom \code{prior} must also be specified.}
  \item{prior}{Prior function to be passed as \code{prior} argument to \code{\link{createBayesianSetup}}.
  If not specified, it is generated based on the \code{*.prior} arguments provided. If a user-defined prior is provided,
  parameter names must be included in \code{prior$lower}, \code{prior$upper} (see example).}
  \item{sampler.type}{Character. It is passed as argument \code{sampler} in \code{\link{runMCMC}}.}
  \item{tau.prior}{List that defines the prior distribution for the decay time parameter \code{tau}. See details. }
  \item{log_p0.prior}{List that defines the prior distribution for the logarithm (in base 10) of the baseline production rate 
  \code{p0}.  See details. If \code{covariates = NULL}, this defines the prior distribution for the logarithm (in base 10) of 
  production rates for all \code{river} reaches.}
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
  \item{omega_d.prior}{Prior distribution for the overdispersion parameter for direct sampling density observations.}
  \item{alpha.prior}{Prior distribution for the inverse DNA shedding rate (i.e., the organismal density 
  that sheds a unit eDNA value per unit time).}
  \item{verbose}{Logical. Should console output be displayed?}
}
\details{
The arguments of the type \code{*.prior} consist in the lists  of arguments required by \code{\link{dtrunc}} 
(except the first argument \code{x}).

By default, AEMs are computed without attributing weights to the edges  of the river network. 
Use e.g. \code{par.AEM = list(weight = "gravity")} to attribute weights. 
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
parameter estimate \code{param_map}. It has length equal to \code{river$AG$nNodes}. If a custom \code{likelihood} is provided,
this is a vector of null length (in which case the user should calculate the probability of detection independently, based on
the chosen likelihood).}
\item{cI}{Output of the call to \code{\link{getCredibleIntervals}}.}
\item{gD}{Output of the call to \code{\link{gelmanDiagnostics}}.}
\item{covariates}{Data frame containing input covariate values (possibly Z-normalized).}
\item{source.area}{Vector of source area values.}
\item{outMCMC}{Object of class \code{mcmcSampler} returned by the call to \code{\link{runMCMC}}.} 
Moreover, arguments \code{ll.type} (possibly changed to \code{"custom"} if a custom likelihood is specified), \code{no.det} 
and \code{data} are added to the list.
}

\examples{
data(wigger)
data(dataCD)

# reduce number of iterations for illustrative purposes 
# (use default mcmc.settings to ensure convergence)
settings.short <- list(iterations = 1e3, thin = 10)
set.seed(1)
out <- run_eDITH_BT_joint(dataCD, wigger, mcmc.settings = settings.short)

\donttest{
library(rivnet)
# best-fit (maximum a posteriori) map of eDNA production rates
plot(wigger, out$p_map)

# best-fit map (maximum a posteriori) of detection probability
plot(wigger, out$probDet_map)


# compare best-fit vs observed values
data.e <- which(dataCD$type=="e")
data.d <- which(dataCD$type=="d")

plot(out$C_map[dataCD$ID[data.e]], dataCD$values[data.e],
	xlab="Modelled (MAP) eDNA concentrations", ylab="Observed eDNA concentrations")
abline(a=0, b=1) 

plot(out$p_map[dataCD$ID[data.d]], dataCD$values[data.d],
	xlab="Modelled (MAP) eDNA production rate", ylab="Observed density data")


}
}

