\docType{data}
\name{outSample}
\alias{outSample}
\title{Posterior sample from fitted eDITH model}
\format{A list.}
\usage{
data(outSample)
}
\description{
It is produced via:

\code{covariates <- data.frame(urban=wigger$SC$locCov$landcover_1,
                         agriculture=wigger$SC$locCov$landcover_2,
                         forest=wigger$SC$locCov$landcover_3,
                         elev=wigger$AG$Z,
                         log_drainageArea=log(wigger$AG$A))}

\code{set.seed(1)}

\code{outSample <- run_eDITH_BT(dataC, wigger, covariates, 
                           mcmc.settings=list(iterations=9e5, burnin = 6e5, message = TRUE, thin = 30))}
						   
}
\keyword{datasets}
