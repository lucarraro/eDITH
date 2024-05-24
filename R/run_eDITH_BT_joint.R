run_eDITH_BT_joint <-
  function(data, river, covariates = NULL, Z.normalize = TRUE,
           use.AEM = FALSE, n.AEM = NULL, par.AEM = NULL,
           no.det = FALSE, ll.type = "norm", source.area = "AG",
           mcmc.settings = NULL, likelihood = NULL, prior = NULL, sampler.type = "DREAMzs",
           tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4))),
           log_p0.prior = list(spec="unif",min=-20, max=0),
           beta.prior = list(spec="norm",sd=1),
           sigma.prior = list(spec="unif",min=0, max=max(data$values[data$type=="e"], na.rm = TRUE)),
           omega.prior = list(spec="unif",min=1, max=10*max(data$values[data$type=="e"], na.rm = TRUE)),
           Cstar.prior = list(spec="unif",min=0, max=max(data$values[data$type=="e"], na.rm = TRUE)),
           omega_d.prior = list(spec="unif",min=1,
                                max=10*max(c(0.11, data$values[data$type=="d"]), na.rm = TRUE)),
           alpha.prior = list(spec="unif", min=0, max=1e6),
           verbose = FALSE){

    if (is.null(prior) & !is.null(likelihood)){
      stop('If a custom likelihood is specified, a custom prior must also be specified.')
    }

    if (!no.det & ll.type =="lnorm" & any(data$values==0)){
      stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
    }

    if (!all(names(data) %in% c("ID","values","type"))){
      stop("data must contain fields named 'ID', 'values' and 'type'.")
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

    if (is.null(prior)){
      out <- prepare.prior.joint(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                           beta.prior, sigma.prior, omega.prior, Cstar.prior, omega_d.prior,
                           alpha.prior,
                           river$AG$nNodes)
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
      likelihood <- function(param){likelihood_generic.joint(param, river, ss, source.area, covariates,
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



