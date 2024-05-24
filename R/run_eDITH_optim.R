run_eDITH_optim <-
  function(data, river, covariates = NULL, Z.normalize = TRUE,
           use.AEM = FALSE, n.AEM = NULL, par.AEM = NULL,
           no.det = FALSE, ll.type = "norm", source.area = "AG",
           likelihood = NULL, sampler = NULL, n.attempts = 100,
           n.restarts = round(n.attempts/10), par.optim = NULL,
           tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4))),
           log_p0.prior = list(spec="unif",min=-20, max=0),
           beta.prior = list(spec="norm",sd=1),
           sigma.prior = list(spec="unif",min=0, max=1*max(data$values, na.rm = TRUE)),
           omega.prior = list(spec="unif",min=1, max=10*max(data$values, na.rm = TRUE)),
           Cstar.prior = list(spec="unif",min=0, max=1*max(data$values, na.rm = TRUE)),
           verbose = FALSE){

    if (is.null(par.optim$control)) par.optim$control <- list(fnscale = -1, maxit = 1e6)
    if (is.null(par.optim$control$fnscale)) par.optim$control$fnscale <- -1

    if (!no.det & ll.type =="lnorm" & any(data$values==0)){
      stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
    }

    if (!all(names(data) %in% c("ID","values"))){
      stop("data must contain fields named 'ID' and 'values'.")
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

    out <- prepare.prior(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                         beta.prior, sigma.prior, omega.prior, Cstar.prior, river$AG$nNodes)
    names.par <- out$names.par; allPriors <- out$allPriors
    lb <- ub  <- numeric(0)
    for (nam in names.par){
      lb <- c(lb, allPriors[[nam]]$a)
      ub <- c(ub, allPriors[[nam]]$b)
    }
    names(lb) <- names(ub)  <- names.par

    if(is.null(likelihood)){ # likelihood_generic is taken from run_eDITH_BT
      likelihood <- function(param){likelihood_generic(param, river,
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


