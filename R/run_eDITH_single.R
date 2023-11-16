run_eDITH_single <- function(param, river, covariates,  Z.normalize = TRUE, no.det=FALSE,
                             ll.type=NULL, data=NULL, source.area="AG", #likelihood = NULL, prior = NULL,
                             tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4))),
                             log_p0.prior = list(spec="unif",min=-20, max=0),
                             beta.prior = list(spec="norm",sd=1),
                             sigma.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE)),
                             omega.prior = list(spec="unif",min=1, max=10*max(data$values, na.rm = TRUE)),
                             Cstar.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE))){

  if (isTRUE(ll.type=="norm") & isTRUE(is.na(param["sigma"]))) stop('Missing sigma parameter for ll.type = "norm".')
  if (isTRUE(ll.type=="lnorm") & isTRUE(is.na(param["sigma"]))) stop('Missing sigma parameter for ll.type = "lnorm".')
  if (isTRUE(ll.type=="nbinom") & isTRUE(is.na(param["omega"]))) stop('Missing sigma parameter for ll.type = "nbinom".')
  if (no.det & isTRUE(is.na(param["Cstar"]))) stop('Missing Cstar parameter for no.det = TRUE.')

  if (length(river$AG$A)==0) {
    stop('river is not aggregated. You should run rivnet::aggregate_river on river prior to run_eDITH_BT.')}

  if (length(river$AG$discharge)==0) {
    stop('Missing hydrological data. You should run rivnet::hydro_river on river prior to run_eDITH_BT.')}

  # calculate additional hydraulic variables
  ss <- sort(river$AG$A,index.return=T); ss <- ss$ix
  q <- numeric(river$AG$nNodes)
  for (i in 1:river$AG$nNodes) q[i] <- river$AG$discharge[i] - sum(river$AG$discharge[which(river$AG$downNode==i)])

  if (source.area=="AG"){
    source.area <- river$AG$width*river$AG$leng
  } else if (source.area=="SC"){source.area <- river$SC$A}

  # Z-normalize covariates
  if (Z.normalize & !is.null(covariates)){
    for (i in 1:length(covariates)){
      covariates[,i] <- (covariates[,i]-mean(covariates[,i]))/sd(covariates[,i])
    }
  }

  # run eDITH
  out <- eval.pC.pD(param, river, ss, covariates, source.area,
                    q, ll.type, no.det)

  out$tau <- NULL
  if (!is.null(data) & !is.null(ll.type)){
    pr <- prepare.prior(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                        beta.prior, sigma.prior, omega.prior, Cstar.prior, river$AG$nNodes)
    allPriors <- pr$allPriors

    logprior <- density_generic(param, allPriors)
    loglik <- likelihood_generic(param, river, ss, source.area, covariates,
                                 data, no.det, ll.type)
    logpost <- logprior + loglik

    out[["logprior"]] <- logprior
    out[["loglik"]] <- loglik
    out[["logpost"]] <- logpost
  }
  invisible(out)
}
