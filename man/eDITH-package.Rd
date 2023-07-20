\name{eDITH-package}
\alias{eDITH-package}
\alias{eDITH}
\docType{package}
\title{
  Model Transport of Environmental DNA In River Networks
}
\description{
  Runs the eDITH (eDNA Integrating Transport and Hydrology) model,
  which implements a mass balance of eDNA transport at a river network scale 
  coupled with a species distribution model to obtain maps of species distribution. 
  eDITH can work with both eDNA concentration (e.g., obtained via qPCR) or 
  metabarcoding (read count) data. Parameter estimation can be performed 
  via Bayesian techniques (via the \code{BayesianTools} package) or optimization
  algorithms. An interface to the \code{DHARMa} package for posterior predictive
  checks is provided.
}
\author{
Luca Carraro (\email{luca.carraro@hotmail.it})
}
\references{
Carraro, L., Hartikainen, H., Jokela, J., Bertuzzo, E., and Rinaldo, A. (2018). Estimating species distribution and abundance in river networks using environmental DNA. Proceedings of the National Academy of Sciences of the United States of America, 115(46), 11724-11729. doi:10.1073/pnas.1813843115

Carraro, L., Maechler, E., Wuethrich, R., and Altermatt, F. (2020). Environmental DNA allows upscaling spatial patterns of biodiversity in freshwater ecosystems. Nature Communications, 11(1) doi:10.1038/s41467-020-17337-8

Carraro, L., Stauffer, J. B., and Altermatt, F. (2021). How to design optimal eDNA sampling strategies for biomonitoring in river networks. Environmental DNA, 3(1), 157-172. doi:10.1002/edn3.137

Carraro, L., Blackman, R. C., and Altermatt, F. (2023). Modelling environmental DNA transport in rivers reveals highly resolved spatio-temporal biodiversity patterns. Scientific Reports, 13(1) doi:10.1038/s41598-023-35614-6
}
