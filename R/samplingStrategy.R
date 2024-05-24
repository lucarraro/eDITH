sampling_strategy_eDNA <- function(river, nSites){

  # control that river has paths
  if (length(river$AG$downstreamPathLength)==0){
    stop("Missing fields in river.
         You should run paths_river prior to this function.")
  }

  if (nSites==0){
    sites <- integer(0)

  } else {
    nNodes <- river$AG$nNodes
    # find confluents
    confluent <- 1+numeric(nNodes)
    hasConfluence <- logical(nNodes)
    for (i in 1:nNodes){
      j <- river$AG$downNode[i]
      k <- NULL
      if (j != 0){
        k <- which(river$AG$downNode==j)
        k <- setdiff(k, i)
        if (length(k)>1){
          k <- k[which.max(river$AG$A[k])]
        }
      }
      if (length(k)==1){
        hasConfluence[i] <- TRUE
        confluent[i] <- k}
    }

    # score for most "equal" confluences
    scoreEqual <- 4*hasConfluence*river$AG$nUpstream[confluent]*river$AG$nUpstream/(river$AG$nUpstream+river$AG$nUpstream[confluent])^2
    scoreEqual[scoreEqual==0] <- min(scoreEqual[scoreEqual>0]) # change 0's to min non-0
    scoreEqual[river$AG$outlet] <- max(scoreEqual) # rank outlet as best node

    # score for river size
    scoreArea <- (river$AG$nUpstream/max(river$AG$nUpstream))^0.5

    # score for distance from high-score upstream&downstream sites
    scoreTmp <- scoreEqual*scoreArea
    scaleLength <- 0.01*sum(river$AG$leng) #max(river$AG$A)/river$cellsize/1e3
    avgScoreDist <- numeric(nNodes)
    pathlengthMatrix <- as.matrix(river$AG$downstreamPathLength)
    for (i in 1:nNodes){
      dd <- as.numeric(pathlengthMatrix[,i]) + as.numeric(pathlengthMatrix[i,])
      dd[dd==0] <- Inf
      avgScoreDist[i] <- sum(scoreTmp*exp(-dd/scaleLength))/sum(exp(-dd/scaleLength))
    }
    avgScoreDist[is.na(avgScoreDist)] <- 0

    scoreDist <- (nNodes+1-rank(avgScoreDist))/(nNodes+1)

    scoreFinal <- scoreTmp*scoreDist

    tmp <- sort(scoreFinal,decreasing=T,index.return=T)

    sites <- tmp$ix[1:nSites]
  }
  invisible(sites)
}


sampling_strategy_direct <- function(river, nSites){

  # control that river is aggregated
  if (length(river$AG$A)==0){
    stop("Missing fields in river.
         You should run aggregate_river prior to this function.")
  }

  if (nSites==0){
    sites <- integer(0)
  } else {
    sites <- numeric(nSites)
    dm <- as.matrix(dist(cbind(river$AG$X, river$AG$Y)))
    sites[1] <- which.min(river$AG$X)
    for (i in 2:nSites){
      vecDist <- as.matrix(dm[,sites[1:(i-1)]])
      sites[i] <- which.max(apply(vecDist,1,min))[1]
    }
  }
  invisible(sites)
}
