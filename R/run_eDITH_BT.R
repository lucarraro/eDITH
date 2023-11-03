run_eDITH_BT <-
  function(data, river, covariates = NULL, Z.normalize = TRUE,
           use.AEM = FALSE, n.AEM = NULL, par.AEM = NULL,
           no.det = FALSE, ll.type = "norm", source.area = "AG",
           mcmc.settings = NULL, likelihood = NULL, prior = NULL, sampler.type = "DREAMzs",
           tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4))),
           log_p0.prior = list(spec="unif",min=-20, max=0),
           beta.prior = list(spec="norm",sd=1),
           sigma.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE)),
           omega.prior = list(spec="unif",min=1, max=10*max(data$values, na.rm = TRUE)),
           Cstar.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE)),
           verbose = FALSE){

    if (is.null(prior) & !is.null(likelihood)){
      stop('If a custom likelihood is specified, a custom prior must also be specified.')
    }

    if (!no.det & ll.type =="lnorm" & any(data$values==0)){
      stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
    }

    if (length(river$AG$A)==0) {
      stop('river is not aggregated. You should run rivnet::aggregate_river on river prior to run_eDITH_BT.')}

    if (length(river$AG$discharge)==0) {
    stop('Missing hydrological data. You should run rivnet::hydro_river on river prior to run_eDITH_BT.')}

    if (!is.null(likelihood)){ll.type="custom"}

    if (is.null(n.AEM)){n.AEM <- round(0.1*river$AG$nNodes)}

    if (is.null(covariates)){
      use.AEM <- TRUE
      if (verbose){
        if (isTRUE(par.AEM$moranI)){
          message("Covariates not specified. Production rates will be estimated
                      based on AEMs with significantly positive spatial autocorrelation",
                  appendLF=F)}
      } else {
        message(sprintf("Covariates not specified. Production rates will be estimated
                      based on the first n.AEM = %d AEMs. \n",n.AEM),appendLF=F)}
    }

    if (use.AEM){
      par.AEM$river <- river
      out <- do.call(river_to_AEM, par.AEM)
      if (!is.null(out$moranI)){ select.AEM <- which(out$moranI$pvalue < 0.05)
      } else {select.AEM <- 1:n.AEM}
      cov.AEM <- data.frame(out$vectors[,select.AEM])
      names(cov.AEM) <- paste0("AEM",1:select.AEM)
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

    if (is.null(prior)){
      out <- prepare.prior(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                           beta.prior, sigma.prior, omega.prior, Cstar.prior, river$AG$nNodes)
      names.par <- out$names.par; allPriors <- out$allPriors
      lb <- ub <- numeric(0)
      for (nam in names.par){
        lb <- c(lb, allPriors[[nam]]$a)
        ub <- c(ub, allPriors[[nam]]$b)
      }
      names(lb) <- names(ub) <- names.par

      density <- function(x){density_generic(x, allPriors = allPriors)}
      sampler <- function(n=1){sampler_generic(n, allPriors = allPriors)}

      prior <- createPrior(density = density, sampler = sampler, lower = lb, upper = ub)

    } else {
      names.par <- names(prior$lower)
      if (is.null(names.par)){
        stop('Missing parameter names in user-defined prior.
             Ensure that objects "lower" and "upper" of the prior contain parameter names.')
      }

      # re-structure sampler so that it yields parameter names
      sampler_no.name <- prior$sampler
      sampler2 <- function(n=1){ ss <- sampler_no.name(n); colnames(ss) <- names.par; return(ss[1,])}
      prior <- createPrior(density = prior$density, sampler = sampler2, lower = prior$lower, upper = prior$upper)
    }

    if(is.null(likelihood)){
      likelihood <- function(param){likelihood_generic(param, river, ss, source.area, covariates,
                                                       data, no.det, ll.type)}}

    #options(warn=-1) # ignore warnings

    setUp <- createBayesianSetup(likelihood = likelihood, prior = prior, names = names.par)

    if (is.null(mcmc.settings)){
      settings <- list(iterations = 2.7e6, burnin = 1.8e6, message = verbose, thin = 10)
    } else {  settings <- mcmc.settings}

    outMCMC <- runMCMC(bayesianSetup = setUp, sampler = sampler.type, settings = settings)
    map <- MAP(outMCMC)
    chains <- outMCMC$chain[[1]]
    n.chains <- length(outMCMC$chain)
    if (n.chains > 1){
      for (ind in 2:n.chains){ chains <- rbind(chains, outMCMC$chain[[ind]])}
    }
    cI <- getCredibleIntervals(chains)
    colnames(cI) <- c(names.par,"logpost","loglik","prior")
    gD <- gelmanDiagnostics(outMCMC)

    tmp <- eval.pC.pD(map$parametersMAP, river, ss, covariates, source.area,
                      q, ll.type, no.det)

    out <- list(p_map = tmp$p, C_map = tmp$C,
                probDet_map = tmp$probDetection,
                param_map = map$parametersMAP,
                ll.type=ll.type, no.det=no.det, cI=cI, gD=gD, data=data,
                covariates = covariates, source.area = source.area,
                outMCMC = outMCMC) # this is biggish (OK with thinning)

    invisible(out)
  }


set_boundaries <- function(ll){
  if (ll$spec=="unif" & !is.null(ll$min) & is.null(ll$a)) {ll$a <- ll$min}
  if (ll$spec=="unif" & !is.null(ll$max) & is.null(ll$b)) {ll$b <- ll$max}
  if(is.null(ll$a)){ll$a <- -Inf}
  if(is.null(ll$b)){ll$b <- Inf}
  invisible(ll)
}

density_generic = function(param, allPriors){
  nPars <- length(allPriors)
  d_comp <- numeric(nPars)
  names(d_comp) <- names(allPriors)
  for (ind in 1:nPars){
    nam <- names(allPriors)[ind]
    pri <- allPriors[[nam]]
    pri[["x"]] <- param[nam]
    d_comp[nam] <-  log(do.call(dtrunc, pri))
  }
  d_out <- sum(d_comp)
  return(d_out)
}

sampler_generic = function(n=1, allPriors){
  nPars <- length(allPriors)
  r_comp <- numeric(nPars)
  names(r_comp) <- names(allPriors)
  for (ind in 1:nPars){
    nam <- names(allPriors)[ind]
    pri <- allPriors[[nam]]
    pri[["n"]] <- n
    r_comp[nam] <- do.call(rtrunc, pri)
  }
  return(r_comp)
}

likelihood_generic <- function(param, river, ss, source.area, covariates, data,
                               no.det=FALSE, ll.type="norm"){

  tmp <- eval.pC.pD(param, river, ss, covariates, source.area)
  ConcMod <- tmp$C
  # if (!is.null(tau_min)){
  # tau <- (tau_min + (tau_max-tau_min)*(param["tau"])/(1+exp(param["tau"])))*3600
  # } else {
  # tau <- param["tau"]*3600}
  #
  # p <- eval.p(param, covariates)
  # ConcMod <- evalConc2_cpp(river, ss, source.area, tau, p, "AG")

  if (no.det){
    phi <- exp(-ConcMod/param["Cstar"])
    sites_nondet <- data$ID[data$values==0]
    sites_data_det <- data$values!=0
    sites_det <- data$ID[sites_data_det]
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
    sites_data_det <- 1:length(data$ID)
    sites_det <- data$ID
  }

  if (ll.type=="nbinom" | ll.type=="geom"){
    y <- round(data$values[sites_data_det])
  } else {y <- data$values[sites_data_det]}

  list_density <- prepare.list.density(y, param, ConcMod[sites_det], ll.type)
  loglik_values <- log(do.call(dtrunc, list_density))
  loglik_values[loglik_values==-Inf] <- -1e4
  loglik_values[is.nan(loglik_values)] <- -1e4
  loglik <- sum(loglik_nondet) + sum(loglik_det) + sum(loglik_values)

  return(loglik)
}


prepare.prior <- function(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                          beta.prior, sigma.prior, omega.prior, Cstar.prior, nNodes){

  if (!is.null(covariates)){
    names.beta <- paste0("beta_",names(covariates))
    if (ll.type=="geom"){
      names.par <- c("tau","log_p0",names.beta)
    } else if (ll.type=="nbinom"){
      names.par <- c("tau","log_p0",names.beta,"omega")
    } else {
      names.par <- c("tau","log_p0",names.beta,"sigma")
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


prepare.list.density <- function(x, param, ConcMod, ll.type){
  list_density <- list(x=x, spec=ll.type)
  if (ll.type=="norm"){
    list_density[["mean"]] <- ConcMod
    list_density[["sd"]] <- param["sigma"]
    list_density[["a"]] <- 0
  } else if (ll.type=="lnorm"){
    list_density[["meanlog"]] <- log(ConcMod^2/sqrt(param["sigma"]^2 + ConcMod^2))
    list_density[["sdlog"]] <- sqrt(log(param["sigma"]^2/ConcMod^2 + 1))
    list_density[["a"]] <- 0
  } else if (ll.type=="nbinom"){
    list_density[["size"]] <- ConcMod/(param["omega"]-1)
    list_density[["prob"]] <- 1/param["omega"]
    list_density[["a"]] <- -Inf # otherwise 0 is excluded from the distribution
  } else if (ll.type=="geom"){
    list_density[["prob"]] <- 1/(1+ConcMod)
    list_density[["a"]] <- -Inf # otherwise 0 is excluded from the distribution
  }
  invisible(list_density)
}

eval.p <- function(param, covariates){
  if (!is.null(covariates)){
    p <- 10^param["log_p0"]*exp(as.numeric(as.matrix(covariates) %*% as.matrix(param[grep("beta_",names(param))])))
  } else {
    p <- 10^param[grep("log.p",names(param))]
  }
  invisible(p)
}

eval.pC.pD <- function(param, river, ss, covariates, source.area,
                       q=NULL,ll.type=NULL, no.det=NULL){

  tau <- param["tau"]*3600
  p <- eval.p(param, covariates)
  C <- evalConc2_cpp(river, ss, source.area, tau, p, "AG")

  if (!is.null(ll.type)){
    local_expected_C <- p*source.area*exp(-river$AG$leng/river$AG$velocity/tau)/q
    if (ll.type=="geom"){
      probDetection <- 1 - pgeom(0, prob = 1/(1+local_expected_C))
    } else if (ll.type=="norm") {
      probDetection <- 1 - pnorm(0, mean = local_expected_C, sd = param["sigma"])
    } else if (ll.type=="lnorm"){
      probDetection <- 1 - plnorm(0, meanlog =  log(local_expected_C^2/sqrt(param["sigma"]^2 + local_expected_C^2)),
                                  sdlog = sqrt(log(param["sigma"]^2/local_expected_C^2 + 1)))
    } else if (ll.type=="nbinom"){
      probDetection <- 1 - pnbinom(0, size = local_expected_C/(param["omega"]-1),
                                   prob = 1/param["omega"])
    } else {probDetection = numeric(0)}

    if (no.det) probDetection <- probDetection*(1-exp(-local_expected_C/param["Cstar"]))
  }

  out <- list(p=p, C=C, tau=tau/3600)
  if (!is.null(ll.type)){out[["probDetection"]] <- probDetection}

  invisible(out)
}
