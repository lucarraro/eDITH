posterior_pred_sim_eDITH <- function(x, river, nParamSets = 10000, nDrawsPerParamSet = 10){

  ss <- sort(river$AG$A,index.return=T); ss <- ss$ix

  nChains <- length(x$outMCMC$chain)
  mcmc.sample <- NULL
  for (i in 1:nChains){
    mcmc.sample <- rbind(mcmc.sample, x$outMCMC$chain[[i]][-1,1:x$outMCMC$setup$numPars])}

  mcmc.sample <- unique(mcmc.sample)
  sample.length <- dim(mcmc.sample)[1]

  pps <- matrix(0,  nParamSets*nDrawsPerParamSet, length(x$data$ID))
  k <- 1
  for (i in 1:nParamSets){
    s <- sample(sample.length, 1)
    param <- mcmc.sample[s, ]
    tau <- param["tau"]*3600
    p0 <- 10^param["log_p0"]
    beta <- param[grep("beta_",names(x$param_map))]
    p <- p0*exp(as.numeric(as.matrix(x$covariates) %*% as.matrix(beta)))
    C <- evalConc2_cpp(river, ss, x$source.area, tau, p, "AG")
    for (j in 1:nDrawsPerParamSet){
      if (x$ll.type=="norm"){
        tmp <- rnorm(length(x$data$ID), mean = C[x$data$ID], sd = param["sigma"])

      } else if (x$ll.type=="lnorm"){
        tmp <- rlnorm(length(x$data$ID),
                           meanlog = log(C[x$data$ID]^2/sqrt(param["sigma"]^2 + C[x$data$ID]^2)),
                           sdlog = sqrt(log(param["sigma"]^2/C[x$data$ID]^2 + 1)))

      } else if (x$ll.type=="nbinom"){
        tmp <- rnbinom(length(x$data$ID),
                            size = C[x$data$ID]/(param["omega"]-1),
                            prob = 1/param["omega"])
      }
      if (x$no.det){
        Cstar <- param["Cstar"]
        tmp <- tmp*rbinom(length(x$data$ID), 1, 1-exp(-tmp/Cstar))
      }

      pps[k, ] <- tmp
      k <- k + 1
    }
  }

  invisible(pps)

}
