\name{run_eDITH_single}
\alias{run_eDITH_single}

\title{Run eDITH for a single parameter set}

\description{
Function that runs the eDITH model for a given parameter set
}

\usage{
run_eDITH_single(param, river, covariates,  Z.normalize = TRUE, 
				no.det = FALSE, ll.type = NULL, 
				data = NULL, source.area = "AG",
                tau.prior = list(spec="lnorm",a=0,b=Inf, 
					meanlog=log(5), sd=sqrt(log(5)-log(4))),
                log_p0.prior = list(spec="unif",min=-20, max=0),
                    beta.prior = list(spec="norm",sd=1),
                sigma.prior = list(spec="unif",min=0, 
					max=max(data$values, na.rm = TRUE)),
                omega.prior = list(spec="unif",min=1, 
					max=10*max(data$values, na.rm = TRUE)),
                 Cstar.prior = list(spec="unif",min=0, 
					max=max(data$values, na.rm = TRUE)))
}

\arguments{
  \item{param}{Parameter set. It has to be a named vector, with names:
  \describe{
  \item{\code{tau}}{Decay time (expressed in h).}
  \item{\code{log_p0}}{Natural logarithm of the baseline production rate.}
  \item{\code{beta_X}}{Effect size of covariate \code{X}. There must be as many \code{beta_X} as columns in \code{covariates}. 
  \code{X} must be the name of the corresponding column in \code{covariates}.}
  \item{\code{omega}, \code{sigma}, \code{Cstar}}{Parameters for estimation of the log-likelihood and detection probability. 
  Only required if \code{ll.type} is provided.}
  }}
  \item{river}{A \code{river} object generated via \code{\link{aggregate_river}}.}
  \item{covariates}{Data frame containing covariate values for all \code{river} reaches.}
  \item{Z.normalize}{Logical. Should covariates be Z-normalized?}
  \item{no.det}{Logical. Should a probability of non-detection be included in the model?}
  \item{ll.type}{Character. String defining the error distribution used in the log-likelihood formulation. 
  Allowed values are \code{norm} (for normal distribution), \code{lnorm} (for lognormal distribution), 
  \code{nbinom} (for negative binomial distribution) and \code{geom} (for geometric distribution). The two latter choices 
  are suited when eDNA data are expressed as read numbers, while \code{norm} and \code{lnorm} are better suited 
  to eDNA concentrations. }
  \item{data}{eDNA data. Data frame containing columns \code{ID} (index of the AG node/reach where 
  the eDNA sample was taken) and \code{values} (value of the eDNA measurement, expressed as concentration
  or number of reads). }
  \item{source.area}{Defines the extent of the source area of a node. Possible values are \code{"AG"} (if the source
  area is the reach surface, i.e. length*width), \code{"SC"} (if the source area is the subcatchment area), or,
  alternatively, a vector with length \code{river$AG$nodes}.}
  \item{tau.prior, log_p0.prior,beta.prior,sigma.prior,omega.prior,Cstar.prior}{Prior distribution for the relevant
  parameters of the eDITH model. Only used if both \code{ll.type} and \code{data} are provided.}
}
\value{
A list with objects:
\item{p}{Vector of eDNA production rates corresponding to the parameter 
set \code{param}. It has length equal to \code{river$AG$nNodes}.}
\item{C}{Vector of eDNA values (concentrations or read numbers) corresponding to the parameter set \code{param}. 
It has length equal to \code{river$AG$nNodes}.}
\item{probDet}{Vector of detection probabilities corresponding to the parameter set \code{param}. It is only computed if 
\code{ll.type} is provided.  It has length equal to \code{river$AG$nNodes}. }
\item{logprior}{Value of the log-prior distribution (computed only if \code{ll.type} and \code{data} are provided).}
\item{loglik}{Value of the log-likelihood distribution (computed only if \code{ll.type} and \code{data} are provided).}
\item{logpost}{Value of the log-posterior distribution (computed only if \code{ll.type} and \code{data} are provided).}
}
\seealso{
See \code{\link{run_eDITH_BT}}, \code{\link{run_eDITH_optim}} for details on parameters names and log-likelihood specification.
}
\examples{
library(rivnet)
data(wigger)

# calculate AEMs and use the first 10 as covariates
ae <- river_to_AEM(wigger)
covariates <- data.frame(ae$vectors[,1:10])
names(covariates) <- paste0("AEM",1:10) 
# covariates names must correspond to param names
set.seed(1); param <- c(3,-15, runif(10,-1,1))
names(param) <- c("tau", "log_p0", paste0("beta_AEM",1:10))
# param names must correspond to covariates names

out <- run_eDITH_single(param, wigger, covariates)

# add parameter sigma and compute detection probability
param <- c(param, 5e-12) 
names(param)[length(param)] <- "sigma"
# note that the value of sigma has to be within the range indicated by sigma.prior
out2 <- run_eDITH_single(param, wigger, covariates, ll.type="norm")

# include data and compute logprior, loglikelihood, logposterior
data(dataC)
out3 <- run_eDITH_single(param, wigger, covariates, 
		ll.type="norm", data=dataC)

}

