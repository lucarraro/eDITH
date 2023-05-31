eval_posterior_eDITH <- function(x, river, covariates, Z.normalize=TRUE, source.area="AG"){

  ss <- sort(river$AG$A, index.return=T); ss <- ss$ix

  q <- numeric(river$AG$nNodes)
  for (i in 1:river$AG$nNodes) q[i]<-river$AG$discharge[i]-sum(river$AG$discharge[which(river$AG$downNode==i)])

  # Z-normalize covariates
  if (Z.normalize){
    for (i in 1:length(covariates)){
      covariates[,i] <- (covariates[,i]-mean(covariates[,i]))/sd(covariates[,i])
    }
  }

  if (source.area=="AG"){
    source.area <- river$AG$width*river$AG$leng
  } else if (source.area=="SC"){source.area <- river$SC$A}

  mcmc.sample <- rbind(x$outMCMC$chain[[1]][-1,1:x$outMCMC$setup$numPars],
                       x$outMCMC$chain[[2]][-1,1:x$outMCMC$setup$numPars],
                       x$outMCMC$chain[[3]][-1,1:x$outMCMC$setup$numPars])
  mcmc.sample <- unique(mcmc.sample)
  probDetection_mat <- matrix(0,dim(mcmc.sample)[1],river$AG$nNodes)
  p_mat <- matrix(0,dim(mcmc.sample)[1],river$AG$nNodes)
  C_mat <- matrix(0,dim(mcmc.sample)[1],river$AG$nNodes)

  for (r in 1:dim(mcmc.sample)[1]){
    tau <- mcmc.sample[r,"tau"]*3600
    p0 <- 10^mcmc.sample[r,"log_p0"]
    beta <- mcmc.sample[r,grep("beta_",names(x$param_map))]
    p <- p0*exp(as.numeric(as.matrix(covariates) %*% as.matrix(beta)))
    C <- evalConc2_cpp(river, ss, source.area, tau, p, "AG")

    p_mat[r, ] <- p
    C_mat[r, ] <- C

    local_expected_C <- p*source.area*exp(-river$AG$leng/river$AG$velocity/tau/3600)/q
    if (x$ll.type=="norm") {
      probDetection_mat[r, ] <- 1 - pnorm(0, mean = local_expected_C, sd = mcmc.sample[r,"sigma"])
    } else if (x$ll.type=="lnorm"){
      probDetection_mat[r, ] <- 1 - plnorm(0, meanlog = log(local_expected_C^2/sqrt(mcmc.sample[r,"sigma"]^2 + local_expected_C^2)),
                                           sdlog = sqrt(log(mcmc.sample[r,"sigma"]^2/local_expected_C^2 + 1)))
    } else if (x$ll.type=="nbinom"){
      probDetection_mat[r, ] <- 1 - pnbinom(0, size = local_expected_C/(mcmc.sample[r,"omega"]-1),
                                            prob = 1/mcmc.sample[r,"omega"])
    }
  }

  x[["p_median"]] <- apply(p_mat,2,median)
  x[["p_mean"]]   <- apply(p_mat,2,mean)
  x[["C_median"]] <- apply(C_mat,2,median)
  x[["C_mean"]]   <- apply(C_mat,2,mean)
  x[["probDetection_median"]] <- apply(probDetection_mat,2,median)
  x[["probDetection_mean"]]   <- apply(probDetection_mat,2,mean)

  invisible(x)
}
