run_eDITH_optim <-
  function(data, river, covariates, Z.normalize = TRUE, no.det = FALSE,
           ll.type = "norm", source.area = "AG",
           likelihood = NULL,  n.attempts = 100, ...){ # no parallel option for the moment
                                                       # no covariate-free option!
                                                       # maybe use AEM when no covariates are given

    dots <- list(...) # list of parameters to be passed to optim
    if (is.null(dots$control)) dots$control <- list(fnscale = -1, maxit = 1e6)
    if (is.null(dots$control$fnscale)) dots$control$fnscale <- -1

    if (!no.det & ll.type =="lnorm" & any(data$values==0)){
      stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
    }

    if (!is.null(likelihood)){ll.type="custom"} # doesn't this give errors?

    # if (is.null(covariates)){
    #   message(sprintf("Covariates not specified. Production rates will be estimated
    #                   independently for the %d reaches. \n",river$AG$nNodes),appendLF=F)}

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

    # default values
    tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4)))
    log_p0.prior = list(spec="unif",min=-20, max=0)
    beta.prior = list(spec="norm",sd=1)
    sigma.prior = list(spec="unif",min=0, max=1*max(data$values, na.rm = TRUE))
    omega.prior = list(spec="unif",min=1, max=10*max(data$values, na.rm = TRUE))
    Cstar.prior = list(spec="unif",min=0, max=1*max(data$values, na.rm = TRUE))

    out <- eDITH:::prepare.prior(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                         beta.prior, sigma.prior, omega.prior, Cstar.prior, river$AG$nNodes)
    names.par <- out$names.par; allPriors <- out$allPriors
    # lb <- ub  <- numeric(0)
    # for (nam in names.par){
    #   lb <- c(lb, allPriors[[nam]]$a)
    #   ub <- c(ub, allPriors[[nam]]$b)
    #}
    #names(lb) <- names(ub)  <- names.par

    if(is.null(likelihood)){ # likelihood_generic is taken from run_eDITH_BT
      likelihood <- function(param){eDITH:::likelihood_generic(param, river,
                                                               ss, source.area,
                                                               covariates,
                                                               data, no.det,
                                                               ll.type)}}
    dots$fn <- likelihood

    sampler <- function(n=1){eDITH:::sampler_generic(n,
                                                     no.det=no.det,
                                                     allPriors = allPriors)}
    ll_end_vec <- counts <- conv <- tau_vec <-  numeric(n.attempts)
    for (ind in 1:n.attempts){
      dots$par <- sampler(1)
      out <- suppressWarnings(do.call(optim, dots))
      tau_vec[ind] <- out$par["tau"]
      counts[ind] <- out$counts["function"]
      conv[ind] <- out$convergence
      if (out$par["tau"]<0) out$value = -Inf # discard solutions with negative tau (unnecessary with logit transf)
      ll_end_vec[ind] <- out$value
      if (ind > 1){
        if (ll_end_vec[ind] > max(ll_end_vec[1:(ind-1)])) out_optim <- out
      } else {out_optim <- out}
        message(sprintf("%.2f%% done \r", ind/n.attempts*100), appendLF = FALSE)
    }
    message("100% done    \n", appendLF = FALSE)


  tmp <- eval.pC.pD(out_optim$par, river, ss, covariates, source.area,
                    q, ll.type, no.det)

  param <- out_optim$par; param["tau"] <- tmp$tau # replace with actual tau (in hours)

  out <- list(p = tmp$p, C = tmp$C, probDetection = tmp$probDetection,
              param = param,
              ll.type=ll.type, no.det=no.det, data=data,
              covariates = covariates, source.area = source.area,
              out_optim = out_optim)

  invisible(out)

  }

# use logit link for tau
#tau <- tau_min + (tau_max*tau_min)*(par["tau"])/(1+exp(par["tau"]))
