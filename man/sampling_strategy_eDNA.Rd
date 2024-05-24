\name{sampling_strategy_eDNA}
\alias{sampling_strategy_eDNA}

\title{Determine optimal spatial arrangement for eDNA sampling}

\description{
Function that determines the optimal spatial arrangement for eDNA sampling sites
}

\usage{
sampling_strategy_eDNA(river, nSites)
}

\arguments{
  \item{river}{A \code{river} object generated via \code{\link{paths_river}}.}
  \item{nSites}{Number of sites to be deployed. Cannot be higher than \code{river$AG$nNodes}.}
}

\value{
A vector containing the ID of the \code{nSites} selected sites according to this strategy. 
Sites are sorted according to their rank (i.e., the first site in the vector is the first one 
that has been selected).
}

\examples{
library(rivnet)
data(wigger)
wigger <- paths_river(wigger) 
sites <- sampling_strategy_eDNA(wigger, 20)
plot(wigger)
points_colorscale(wigger$AG$X[sites], wigger$AG$Y[sites], 1:20)
title("Rank of selected sites")
  
}

