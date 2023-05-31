run_eDITH_BT <-
function(data, river, covariates, Z.normalize=TRUE, no.det=TRUE, ll.type="norm",
         source.area = "AG",
         mcmc.settings = NULL, likelihood = NULL, prior = NULL, sampler.type = "DREAMzs",
         #tau.prior = list(spec="norm",a=0,b=Inf, mean=5, sd=2),
         tau.prior = list(spec="lnorm",a=0,b=Inf, meanlog=log(5), sd=sqrt(log(5)-log(4))),
         log_p0.prior = list(spec="unif",min=-20, max=0),
         beta.prior = list(spec="norm",sd=1),
         sigma.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE)),
         omega.prior = list(spec="unif",min=1, max=10*max(data$values, na.rm = TRUE)),
         Cstar.prior = list(spec="unif",min=0, max=max(data$values, na.rm = TRUE))){

  out <- prepare.prior(covariates, no.det, ll.type, tau.prior, log_p0.prior,
                       beta.prior, sigma.prior, omega.prior, Cstar.prior)
  names.par <- out$names.par; allPriors <- out$allPriors

  ss <- sort(river$AG$A,index.return=T); ss <- ss$ix

  q <- numeric(river$AG$nNodes)
  for (i in 1:river$AG$nNodes) q[i] <- river$AG$discharge[i] - sum(river$AG$discharge[which(river$AG$downNode==i)])

  if (source.area=="AG"){
    source.area <- river$AG$width*river$AG$leng
  } else if (source.area=="SC"){source.area <- river$SC$A}

  Z.normalize=TRUE
  # Z-normalize covariates
  if (Z.normalize){
    for (i in 1:length(covariates)){
      covariates[,i] <- (covariates[,i]-mean(covariates[,i]))/sd(covariates[,i])
    }
  }

  if (!no.det & ll.type =="lnorm" & any(data$values==0)){
    stop('If eDNA data are log-normally distributed, non-detections must be accounted for.')
  }

  lb <- ub <- numeric(0)
  for (nam in names.par){
    lb <- c(lb, allPriors[[nam]]$a)
    ub <- c(ub, allPriors[[nam]]$b)
  }
  names(lb) <- names(ub) <- names.par

  density <- function(x){density_generic(x, no.det=no.det, allPriors = allPriors)}
  sampler <- function(n=1){sampler_generic(n, no.det=no.det, allPriors = allPriors)}

  if (is.null(prior)){prior <- createPrior(density = density, sampler = sampler, lower = lb, upper = ub)}

  if(is.null(likelihood)){
    likelihood <- function(param){likelihood_generic(param, river, ss, source.area, covariates,
                                                     data, no.det, ll.type)}}

  #options(warn=-1) # ignore warnings

  setUp <- createBayesianSetup(likelihood=likelihood,prior=prior, names=names.par)

  if (is.null(mcmc.settings)){
  settings <- list(iterations = 2.7e6, burnin=1.8e6, message = TRUE, thin = 10)
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


  # covariates must be Z-normalized!
  tau_map <- map$parametersMAP["tau"]*3600
  p0_map <- 10^map$parametersMAP["log_p0"]
  beta_map <- map$parametersMAP[grep("beta_",names(map$parametersMAP))]
  p_map <- p0_map*exp(as.numeric(as.matrix(covariates) %*% as.matrix(beta_map)))
  C_map <- evalConc2_cpp(river,ss,source.area,tau_map,p_map,"AG")

  local_expected_C <- p_map*source.area*exp(-river$AG$leng/river$AG$velocity/tau_map/3600)/q

  if (ll.type=="norm") {
    probDetection <- 1 - pnorm(0, mean = local_expected_C, sd = map$parametersMAP[["sigma"]])
  } else if (ll.type=="lnorm"){
    probDetection <- 1 - plnorm(0, meanlog =  log(local_expected_C^2/sqrt(map$parametersMAP["sigma"]^2 + local_expected_C^2)),
                                sdlog = sqrt(log(map$parametersMAP["sigma"]^2/local_expected_C^2 + 1)))
  } else if (ll.type=="nbinom"){
  probDetection <- 1 - pnbinom(0, size = local_expected_C/(map$parametersMAP[["omega"]]-1),
                               prob = 1/map$parametersMAP[["omega"]])
  }

  out <- list(p_map=p_map, C_map=C_map,
              param_map = map$parametersMAP,
              probDet_map = probDetection,
              ll.type=ll.type, no.det=no.det, cI=cI, gD=gD, data=data,
              outMCMC=outMCMC) # this is biggish (OK with thinning)

  invisible(out)
}


set_boundaries <- function(ll){
  if (ll$spec=="unif" & !is.null(ll$min) & is.null(ll$a)) {ll$a <- ll$min}
  if (ll$spec=="unif" & !is.null(ll$max) & is.null(ll$b)) {ll$b <- ll$max}
  if(is.null(ll$a)){ll$a <- -Inf}
  if(is.null(ll$b)){ll$b <- Inf}
  invisible(ll)
}

density_generic = function(param, no.det=TRUE, allPriors){
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
  # if (no.det){ d_out <- sum(d_comp)} else {
  #   d_out <- sum(d_comp[-which(names(d_comp)=="Cstar")])} # not needed! (length of param is self-adapted)
  return(d_out)
}

sampler_generic = function(n=1, no.det=TRUE, allPriors){
  nPars <- length(allPriors)
  r_comp <- numeric(nPars)
  names(r_comp) <- names(allPriors)
  for (ind in 1:nPars){
    nam <- names(allPriors)[ind]
    pri <- allPriors[[nam]]
    pri[["n"]] <- n
    r_comp[nam] <- do.call(rtrunc, pri)
  }
  #if (!no.det){r_comp <- r_comp[-which(names(r_comp)=="Cstar")]}
  return(r_comp)
}

likelihood_generic <- function(param, river, ss, source.area, covariates, data, no.det=FALSE, ll.type="norm"){

  p <- 10^param["log_p0"]*exp(as.numeric(as.matrix(covariates) %*% as.matrix(param[grep("beta_",names(param))])))
  ConcMod <- evalConc2_cpp(river, ss, source.area, param["tau"]*3600, p, "AG")

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

  if (ll.type=="nbinom"){
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
                          beta.prior, sigma.prior, omega.prior, Cstar.prior){

  names.beta <- paste0("beta_",names(covariates))
  if (ll.type=="nbinom"){
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
  }
  invisible(list_density)
}
