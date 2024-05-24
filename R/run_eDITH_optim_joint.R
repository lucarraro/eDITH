run_eDITH_optim_joint <-
  function(data, river, covariates = NULL, Z.normalize = TRUE,
           use.AEM = FALSE, n.AEM = NULL, par.AEM = NULL,
           no.det = FALSE, ll.type = "norm", source.area = "AG",
           likelihood = NULL, sampler = NULL, n.attempts = 100,
           n.restarts = round(n.attempts/10), par.optim = NULL,
           tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4))),
           log_p0.prior = list(spec="unif",min=-20, max=0),
           beta.prior = list(spec="norm",sd=1),
           sigma.prior = list(spec="unif",min=0,
                              max=1*max(c(1e-6, data$values[data$type=="e"]), na.rm = TRUE)),
           omega.prior = list(spec="unif",min=1,
                              max=10*max(c(0.11, data$values[data$type=="e"]), na.rm = TRUE)),
           Cstar.prior = list(spec="unif",min=0,
                              max=1*max(c(1e-6, data$values[data$type=="e"]), na.rm = TRUE)),
           omega_d.prior = list(spec="unif",min=1,
                                max=10*max(c(0.11, data$values[data$type=="d"]), na.rm = TRUE)),
           alpha.prior = list(spec="unif", min=0, max=1e6),
           verbose = FALSE){

    if (is.null(par.optim$control)) par.optim$control <- list(fnscale = -1, maxit = 1e6)
    if (is.null(par.optim$control$fnscale)) par.optim$control$fnscale <- -1

    if (!no.det & ll.type =="lnorm" & any(data$values[data$type=="e"]==0)){
      stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
    }

     if (!all(names(data) %in% c("ID","values","type"))){
       stop("data must contain fields named 'ID', 'values' and 'type'.")
     }

    if (length(river$AG$A)==0) {
      stop('river is not aggregated. You should run rivnet::aggregate_river on river prior to run_eDITH_optim.')}

    if (length(river$AG$discharge)==0) {
      stop('Missing hydrological data. You should run rivnet::hydro_river on river prior to run_eDITH_optim.')}


    if (!is.null(likelihood)){ll.type="custom"} # doesn't this give errors?

    if (is.null(n.AEM)){n.AEM <- round(0.1*river$AG$nNodes)}

    if (is.null(covariates)){
      use.AEM <- TRUE
      if (verbose){
        if (isTRUE(par.AEM$moranI)){
          message("Covariates not specified. Production rates will be estimated
                      based on AEMs with significantly positive spatial autocorrelation",
                  appendLF=F)
        } else {
          message(sprintf("Covariates not specified. Production rates will be estimated
                      based on the first n.AEM = %d AEMs. \n",n.AEM),appendLF=F)}
      }
    }

    if (use.AEM){
      par.AEM$river <- river
      out <- do.call(river_to_AEM, par.AEM)
      if (!is.null(out$moranI)){ select.AEM <- which(out$moranI$pvalue < 0.05)
      } else {select.AEM <- 1:n.AEM}
      cov.AEM <- data.frame(out$vectors[,select.AEM])
      names(cov.AEM) <- paste0("AEM", select.AEM)
      covariates <- data.frame(c(covariates, cov.AEM))
    }

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

    out <- prepare.prior.joint(covariates, no.det, ll.type,
                         tau.prior, log_p0.prior,
                         beta.prior, sigma.prior, omega.prior,
                         Cstar.prior, omega_d.prior,
                         alpha.prior,
                         river$AG$nNodes)
    names.par <- out$names.par; allPriors <- out$allPriors
    lb <- ub  <- numeric(0)
    for (nam in names.par){
      lb <- c(lb, allPriors[[nam]]$a)
      ub <- c(ub, allPriors[[nam]]$b)
    }
    names(lb) <- names(ub)  <- names.par

    if(is.null(likelihood)){ # likelihood_generic is taken from run_eDITH_BT
      likelihood <- function(param){likelihood_generic.joint(param, river,
                                                       ss, source.area,
                                                       covariates,
                                                       data, no.det,
                                                       ll.type)}}

    if (is.null(sampler)){
      sampler <- function(n=1){sampler_generic(n,  allPriors = allPriors)}}

    density <- function(x){
      if (any(x<lb) | any(x>ub) ){ -Inf
      } else {
        density_generic(x, allPriors = allPriors)}}
    logpost <- function(param) {likelihood(param) + density(param) }
    par.optim$fn <- logpost

    ll_end_vec <- counts <- conv <- tau_vec <-  numeric(n.attempts)
    for (ind in 1:n.attempts){
      if (ind %in% seq(1,n.attempts,ceiling(n.attempts/n.restarts))){ # start n.restarts different times
        par.optim$par <- sampler(1)
      } else {par.optim$par <- out$par}
      out <- suppressWarnings(do.call(optim, par.optim))
      tau_vec[ind] <- out$par["tau"]
      counts[ind] <- out$counts["function"]
      conv[ind] <- out$convergence
      ll_end_vec[ind] <- out$value
      if (ind > 1){
        if (ll_end_vec[ind] > max(ll_end_vec[1:(ind-1)])) out_optim <- out
      } else {out_optim <- out}
      if (verbose) message(sprintf("%.2f%% done \r", ind/n.attempts*100), appendLF = FALSE)
    }
    if (verbose) message("100% done    \n", appendLF = FALSE)

    attempts.stats <- list(lp=ll_end_vec, counts=counts, conv=conv, tau=tau_vec)

    tmp <- eval.pC.pD(out_optim$par, river, ss, covariates, source.area,
                      q, ll.type, no.det)

    param <- out_optim$par; param["tau"] <- tmp$tau # replace with actual tau (in hours)

    out <- list(p = tmp$p, C = tmp$C, probDet = tmp$probDetection,
                param = param,
                ll.type=ll.type, no.det=no.det, data=data,
                covariates = covariates, source.area = source.area,
                out_optim = out_optim, attempts.stats = attempts.stats)

    invisible(out)

  }


likelihood_generic.joint <- function(param, river, ss, source.area, covariates, data,
                               no.det=FALSE, ll.type="norm"){

  edata <- subset(data, data$type=="e")
  kdata <- subset(data, data$type=="d")

  tmp <- eval.pC.pD(param, river, ss, covariates, source.area)
  ConcMod <- tmp$C
  DensityMod <- tmp$p*param["alpha"]  # p ~ density (an additional coeff couldn't be identified)

  if (no.det){
    phi <- exp(-ConcMod/param["Cstar"])
    sites_nondet <- edata$ID[edata$values==0]
    sites_data_det <- edata$values!=0
    sites_det <- edata$ID[sites_data_det]
    list_density <- prepare.list.density(0, param, ConcMod[sites_nondet], ll.type) # problem with a
    density0 <- do.call(dtrunc, list_density)
    prob_0detection <- (1-phi[sites_nondet])*density0
    foo <- is.nan(prob_0detection) | prob_0detection=="Inf"
    prob_0detection[foo] <- 1-phi[sites_nondet[foo]] # so that overall prob.nondet=1
    loglik_nondet <- log(phi[sites_nondet] + prob_0detection)
    loglik_nondet[loglik_nondet==-Inf] <- -1e4
    loglik_det <- log(1-phi[sites_det])
    loglik_det[loglik_det==-Inf] <- -1e4
  } else {
    loglik_nondet <- loglik_det <- 0
    sites_data_det <- 1:length(edata$ID)
    sites_det <- edata$ID
  }

  if (ll.type=="nbinom" | ll.type=="geom"){
    y <- round(edata$values[sites_data_det])
  } else {y <- edata$values[sites_data_det]}

  list_density <- prepare.list.density(y, param, ConcMod[sites_det], ll.type)
  loglik_values <- log(do.call(dtrunc, list_density))
  loglik_values[loglik_values==-Inf] <- -1e4
  loglik_values[is.nan(loglik_values)] <- -1e4

  loglik_kicknet <- log(dnbinom(kdata$values,
                                size= DensityMod[kdata$ID]/(param["omega_d"]-1),
                                prob=1/param["omega_d"]))
  loglik_kicknet[loglik_kicknet==-Inf] <- -1e4
  loglik_kicknet[is.nan(loglik_kicknet)] <- -1e4

  loglik <- sum(loglik_nondet) + sum(loglik_det) +
    sum(loglik_values) + sum(loglik_kicknet)

  return(loglik)
}


prepare.prior.joint <- function(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                          beta.prior, sigma.prior, omega.prior, Cstar.prior, omega_d.prior,
                          alpha.prior, nNodes){

  if (!is.null(covariates)){
    names.beta <- paste0("beta_",names(covariates))
    if (ll.type=="geom"){
      names.par <- c("tau","log_p0",names.beta)
    } else if (ll.type=="nbinom"){
      names.par <- c("tau","log_p0",names.beta,"omega","omega_d","alpha")
    } else {
      names.par <- c("tau","log_p0",names.beta,"sigma","omega_d","alpha")
    }
    if (no.det){names.par <- c(names.par,"Cstar")}

    # split beta prior among covariates
    if(length(beta.prior$spec==1)){
      for (nam in names.beta){
        assign(paste0(nam,".prior"),beta.prior)
      }
    } else {
      if (length(beta.prior$spec)!=length(names.beta)){
        stop("Number of elements in beta.prior does not match number of covariates")}
      fieldnames <- names(beta.prior)
      for (ind in 1:length(names.beta)){
        nam <- names.beta[ind]
        foo <- list()
        for (fn in fieldnames){
          foo[[fn]] <- beta.prior[[fn]][ind]
        }
        assign(paste0(nam,".prior"),foo)
      }
    }

    # assign boundaries and copy to allPriors
    allPriors <- list()
    for (nam in names.par){
      eval(parse(text=paste0(nam,".prior <- set_boundaries(",nam,".prior)")))
      eval(parse(text=paste0('allPriors[["',nam,'"]] <- ',nam,'.prior')))
    }
  } else {
    names.p <- character(nNodes)
    for (ind in 1:nNodes) names.p[ind] <- paste0("log.p",ind)
    names.par <- c("tau", names.p)

    # assign boundaries and copy to allPriors
    allPriors <- list()
    tau.prior <- set_boundaries(tau.prior)
    allPriors[["tau"]] <- tau.prior
    for (nam in names.par[-1]){
      eval(parse(text=paste0(nam,".prior <- set_boundaries(log_p0.prior)")))
      eval(parse(text=paste0('allPriors[["',nam,'"]] <- ',nam,'.prior')))
    }
    if (ll.type=="nbinom"){
      names.par <- c(names.par,"omega")
      omega.prior <- set_boundaries(omega.prior)
      allPriors[["omega"]] <- omega.prior
    } else {
      names.par <- c(names.par,"sigma")
      sigma.prior <- set_boundaries(sigma.prior)
      allPriors[["sigma"]] <- sigma.prior
    }
  }

  out <- list(names.par=names.par, allPriors=allPriors)
  invisible(out)
}





