eval_posterior_eDITH <- function(x, river){

  ss <- sort(river$AG$A, index.return=T); ss <- ss$ix

  q <- numeric(river$AG$nNodes)
  for (i in 1:river$AG$nNodes) q[i]<-river$AG$discharge[i]-sum(river$AG$discharge[which(river$AG$downNode==i)])

  nChains <- length(x$outMCMC$chain)
  mcmc.sample <- NULL
  for (i in 1:nChains){
    mcmc.sample <- rbind(mcmc.sample, x$outMCMC$chain[[i]][-1,1:x$outMCMC$setup$numPars])}

  mcmc.sample <- unique(mcmc.sample)
  probDetection_mat <- matrix(0,dim(mcmc.sample)[1],river$AG$nNodes)
  p_mat <- matrix(0,dim(mcmc.sample)[1],river$AG$nNodes)
  C_mat <- matrix(0,dim(mcmc.sample)[1],river$AG$nNodes)

  for (r in 1:dim(mcmc.sample)[1]){
    tau <- mcmc.sample[r,"tau"]*3600
    p0 <- 10^mcmc.sample[r,"log_p0"]
    beta <- mcmc.sample[r,grep("beta_",names(x$param_map))]
    p <- p0*exp(as.numeric(as.matrix(x$covariates) %*% as.matrix(beta)))
    C <- evalConc2_cpp(river, ss, x$source.area, tau, p, "AG")

    p_mat[r, ] <- p
    C_mat[r, ] <- C

    local_expected_C <- p*x$source.area*exp(-river$AG$leng/river$AG$velocity/tau)/q
    if (x$ll.type=="norm") {
      pD <- 1 - pnorm(0, mean = local_expected_C, sd = mcmc.sample[r,"sigma"])
    } else if (x$ll.type=="lnorm"){
      pD <- 1 - plnorm(0, meanlog = log(local_expected_C^2/sqrt(mcmc.sample[r,"sigma"]^2 + local_expected_C^2)),
                       sdlog = sqrt(log(mcmc.sample[r,"sigma"]^2/local_expected_C^2 + 1)))
    } else if (x$ll.type=="nbinom"){
      pD <- 1 - pnbinom(0, size = local_expected_C/(mcmc.sample[r,"omega"]-1),
                        prob = 1/mcmc.sample[r,"omega"])
    } else {pD <- numeric(river$AG$nNodes)}

    if (x$no.det){
      Cstar <- mcmc.sample[r, "Cstar"]
      pD <- pD*(1-exp(-local_expected_C/Cstar))}

    probDetection_mat[r, ] <- pD
  }

  x[["p_median"]] <- apply(p_mat,2,median)
  x[["p_mean"]]   <- apply(p_mat,2,mean)
  x[["C_median"]] <- apply(C_mat,2,median)
  x[["C_mean"]]   <- apply(C_mat,2,mean)

  if (x$ll.type=="custom"){
    x[["probDetection_median"]] <- x[["probDetection_mean"]] <- numeric(0)
  } else {
    x[["probDetection_median"]] <- apply(probDetection_mat,2,median)
    x[["probDetection_mean"]]   <- apply(probDetection_mat,2,mean)
  }

  invisible(x)
}
