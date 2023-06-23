run_eDITH_optim <-
  function(data, river, covariates = NULL, Z.normalize = TRUE, no.det = FALSE,
           ll.type = "norm", source.area = "AG",
           likelihood = NULL, par.info = list(),
           parallel = FALSE, ...){

    dots <- list(...) # list of parameters to be passed to optim
    if (is.null(dots$method)) dots$method <- "L-BFGS-B"
    if (is.null(dots$control)) dots$control <- list(fnscale = -1)

    # par.info = list(names=, lb=, ub=, init=)

    if (!no.det & ll.type =="lnorm" & any(data$values==0)){
      stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
    }

    if (!is.null(likelihood)){ll.type="custom"}

    if (is.null(covariates)){
      message(sprintf("Covariates not specified. Production rates will be estimated
                      independently for the %d reaches. \n",river$AG$nNodes),appendLF=F)}

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
    sigma.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE))
    omega.prior = list(spec="unif",min=1, max=10*max(data$values, na.rm = TRUE))
    Cstar.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE))

    out <- eDITH:::prepare.prior(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                         beta.prior, sigma.prior, omega.prior, Cstar.prior, river$AG$nNodes)
    names.par <- out$names.par; allPriors <- out$allPriors
    lb <- ub <- init <- numeric(0)
    for (nam in names.par){
      lb <- c(lb, allPriors[[nam]]$a)
      ub <- c(ub, allPriors[[nam]]$b)
      ii <- do.call(qtrunc, c(p=0.5, allPriors[[nam]])) # use median as initial value
      init <- c(init, ii)
    }
    names(lb) <- names(ub) <- names(init) <- names.par

    if (is.null(par.info$init)) par.info$init <- init
    if (is.null(par.info$lb)) par.info$lb <- lb
    if (is.null(par.info$ub)) par.info$ub <- ub
    if (is.null(par.info$names)) par.info$names.par <- names.par

    if(is.null(likelihood)){ # likelihood_generic is taken from run_eDITH_BT
      likelihood <- function(param){eDITH:::likelihood_generic(param, river, ss, source.area, covariates,
                                                       data, no.det, ll.type)}}
    par.info$init[["sigma"]] <- 1e-5 # temporary patch to see what happens
    dots$par <- par.info$init
    dots$fn <- likelihood
    if (is.null(dots$lower)) dots$lower <- par.info$lb
    if (is.null(dots$upper)) dots$upper <- par.info$ub

    if (parallel==FALSE){
    out_optim <- do.call(optim, dots)
    } else {
      out_optim <- do.call(optimParallel, dots)
    }

  invisible(out_optim)

  }
