\name{run_eDITH_optim}
\alias{run_eDITH_optim}

\title{Optimize eDITH}

\description{
Function that performs search of optimal parameters of an eDITH model 
}

\usage{
run_eDITH_optim(data, river, covariates = NULL, Z.normalize = TRUE,
           use.AEM = FALSE, n.AEM = NULL, par.AEM = NULL,
           no.det = FALSE, ll.type = "norm", source.area = "AG",
           likelihood = NULL, sampler = NULL, n.attempts = 100, 
		   par.optim = NULL, verbose = FALSE)
}

\arguments{
  \item{data}{eDNA data. Data frame containing columns \code{ID} (index of the AG node/reach where 
  the eDNA sample was taken) and \code{values} (value of the eDNA measurement, expressed as concentration
  or number of reads). }
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
  \item{likelihood}{Likelihood function to be passed as \code{likelihood} argument to \code{\link{createBayesianSetup}}.
  If not specified, it is generated based on arguments \code{no.det} and \code{ll.type}. }
  \item{sampler}{Function generating sets of initial parameter values for the optimization algorithm. If \code{NULL}, 
  initial parameter values are drawn from the default prior distributions of \code{\link{run_eDITH_BT}}. See details.}
  \item{n.attempts}{Number of times the optimizing function \code{optim} is executed. Each time, \code{sampler} is used to 
  draw an initial parameter set.}
  \item{par.optim}{List of parameters to be passed to \code{\link{optim}}. By default, the likelihood is maximized (i.e., 
  \code{control$fnscale = -1}), and the maximum number of iterations is set to 1e6. The default optimization method is 
  "Nelder-Mead" (same default as in \code{optim}).}
  \item{verbose}{Logical. Should console output be displayed?}
}
\details{
If specified by the user, \code{sampler} must be a function with no input argument that produces as output a "named num" 
vector of parameters. Parameter names must be same as in the \code{likelihood}. See example (missing!!).
}
\value{
A list with objects:
\item{p}{Vector of best-fit eDNA production rates corresponding to the optimum parameter 
estimates \code{param}. It has length equal to \code{river$AG$nNodes}.}
\item{C}{Vector of best-fit eDNA values (in the same unit as \code{data$values}, i.e. concentrations or read numbers) 
corresponding to the optimum parameter  estimates \code{param}. It has length equal to \code{river$AG$nNodes}.}
\item{probDet}{Vector of best-fit detection probabilities corresponding to the optimum
parameter estimate \code{param_map}. It has length equal to \code{river$AG$nNodes}. If a custom \code{likelihood} is provided,
this is a vector of null length (in which case the user should calculate the probability of detection independently, based on
the chosen likelihood).}
\item{param}{Vector of named parameters corresponding to the best-fit estimate.}
\item{covariates}{Data frame containing input covariate values (possibly Z-normalized).}
\item{source.area}{Vector of source area values.}
\item{out_optim}{List as provided by \code{optim}. Only the result of the call to \code{optim} (out of \code{n.attempts})
yielding the highest likelihood is exported.} 
Moreover, arguments \code{ll.type} (possibly changed to \code{"custom"} if a custom likelihood is specified), \code{no.det} 
and \code{data} are added to the list.
}

\examples{
data(wigger)
data(dataC)
data(dataRead)

## fit eDNA concentration data - use AEMs as covariates
set.seed(4)
out <- run_eDITH_optim(dataC, wigger, n.AEM = 10, n.attempts = 1) 
# reduced n.AEM and n.attempts for illustrative purposes
# it is recommended to attempt optimization several times to ensure convergence

\donttest{
library(rivnet)
# best-fit map of eDNA production rates
plot(wigger, out$p)

# best-fit map of detection probability
plot(wigger, out$probDet)

# compare best-fit vs observed eDNA concentrations
plot(out$C[dataC$ID], dataC$values,
	xlab="Modelled concentrations", ylab="Observed concentrations")
abline(a=0, b=1) 

## fit eDNA read number data - use AEMs as covariates
set.seed(4)
out <- run_eDITH_optim(dataRead, wigger, ll.type = "nbinom", n.attempts = 1) 
# reduced n.attempts for illustrative purposes

## use user-defined covariates
covariates <- data.frame(urban = wigger$SC$locCov$landcover_1,
                         agriculture = wigger$SC$locCov$landcover_2,
                         forest = wigger$SC$locCov$landcover_3,
                         elev = wigger$AG$Z,
                         log_drainageArea = log(wigger$AG$A))
						 
out.cov <- 	run_eDITH_optim(dataC, wigger, covariates, n.attempts = 1) 
# reduced n.attempts for illustrative purposes	

# use user-defined covariates and AEMs
out.covAEM <- 	run_eDITH_optim(dataC, wigger, covariates, use.AEM = TRUE, n.attempts = 1) 
# reduced n.attempts for illustrative purposes				 
}
\dontrun{ # this takes hours!!
# use AEMs with significantly positive spatial autocorrelation
out.AEM.moran <- run_eDITH_optim(dataC, wigger, par.AEM = list(moranI = TRUE)) 
# reduced n.attempts for illustrative purposes	

}


}


