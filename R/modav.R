## model averaging and unadjusted estimates

## function to get prevalences of all subgroups in any other subgroup
getPrev <- function(subMat){
  nSub <- ncol(subMat)
  subgrNams <- colnames(subMat)
  wS <- wC <- matrix(NA, nSub, nSub)
  for(i in 1:nSub){
    for(j in 1:nSub){
      ## prevalence of subgroup j in subgroup i
      wS[i,j] <- sum(subMat[,i] == 1 & subMat[,j] == 1, na.rm=TRUE)/
        sum(subMat[,i] == 1, na.rm=TRUE)
      ## prevalence of subgroup j in subgroup i complement
      wC[i,j] <- sum(subMat[,i] == 0 & subMat[,j] == 1, na.rm=TRUE)/
        sum(subMat[,i] == 0, na.rm=TRUE) 
    }
  }
  rownames(wS) <- colnames(wS) <- colnames(wC) <- subgrNams
  rownames(wC) <- paste(subgrNams, "compl", sep="-")
  list(wS=wS, wC=wC)
}

## get treatment effect estimates in the subgroups
predSub <- function(cf, vc, w){
  cvec <- c(1,w)
  delta <- sum(cvec*cf)
  sedelta <- as.numeric(sqrt(t(cvec)%*%vc%*%cvec))
  c(delta, sedelta)
}

## get treatment differences between subgr & complement
predIA <- function(cf, vc, w1, w2){
  ## beta+wS*delta-(beta+wC*delta)
  delta <- (w1-w2)*cf[2]
  sedelta <- abs(w1-w2)*sqrt(vc[2,2])
  c(delta, sedelta)
}

umfit <- function(resp, trt, subgr, covars = NULL, data, 
                  fitfunc = c("lm", "glm", "glm.nb", "survreg", "coxph", "rlm"),
                  type = c("modav", "unadj"),
                  event, exposure, 
                  level = 0.1, prior = 1, nullprior = 0, ...){
  ## argument checking
  stopifnot(is.numeric(level), level > 0 & level < 1, length(level) == 1)
  
  ## allow specification of non-character resp and trt
  ## (do not use this in non-interactive use of margMods)
  sizes <- checkTrtSub(data, trt, subgr)
  
  fitfunc <- match.arg(fitfunc)
  nSub <- length(subgr)
  
  
  if (!all(make.names(names(data)) == names(data)) & fitfunc == "rlm") {
    if(!is.null(covars)){
      vars.covars  = attr(terms(covars), "term.labels")
      vars.covars  = tidy.names(vars.covars)
      covars = reformulate(vars.covars)
    }
    names(data)  = tidy.names(names(data))
    subgr = tidy.names(subgr)
    message("rlm does not allow non-standard names for data.frame variables so they were transformed.")
  }
  
  if (!missing(event) & !(fitfunc %in% c("survreg", "coxph")))
    warning("You have specified an event variable with fitfunc other than survreg or coxph")
  if (!missing(exposure) & !(fitfunc %in% c("glm.nb")))
    warning("You have specified an exposure variable with fitfunc other than glm.nb")
  type <- match.arg(type)
  
  if(type == "modav"){
    ## additional checks for modav
    stopifnot(is.numeric(prior),
              length(prior) == 1 | length(prior) == nSub,
              is.numeric(nullprior),
              length(nullprior) == 1)
  }
 
  ## actual model fitting
  res <- fitMods(resp, trt, subgr, covars, data, fitfunc,
                 event, exposure, ...)

  ## calculate posterior model weigths (only needed for modav)
  bic <- sapply(res$fitmods, function(x){
    x$ests$BIC
  })
  bic <- c(bic, res$fitOverall$ests$BIC)
  bic <- bic-mean(bic)
  if(length(prior) == 1){
    pr <- c(rep(1, nSub), nullprior)
  } else {
    pr <- c(prior, nullprior)
  }
  postpr <- pr*exp(-0.5*bic)/sum(pr*exp(-0.5*bic))
  
  ## produce estimates for subgroups
  cis = getConfints(res, data, type, nSub, level, postpr)
  
  obj <- list()
  obj$fitMods <- res
  obj$trtEff <- cis$trtEff
  obj$trtEffDiff <- cis$trtEffDiff
  rownames(obj$trtEffDiff) <- NULL
  obj$n <- nrow(data)
  obj$subgroups <- res$subgroups
  obj$GroupSizes <- sizes
  obj$level <- level
  obj$overall <- cis$overall 
  names(postpr) <- c(subgr, "overall")
  obj$post.weights <- postpr
  obj$pvals <- cis$pvals
  obj$type <- type
  obj$fitfunc <- fitfunc
  obj$data <- data
  class(obj) <- "subtee"
  obj
}

unadj <- function(resp, trt, subgr, covars = NULL, data, 
                  fitfunc = c("lm", "glm", "glm.nb", "survreg", "coxph", "rlm"),
                  event, exposure, level = 0.1,...){
  call <- match.call()
  if(!is.character(call$resp))
    resp <- as.character(call$resp)
  if(!is.character(call$trt))
    trt <- as.character(call$trt)
  umfit(resp = resp, trt = trt, subgr = subgr, covars = covars,
        data = data, fitfunc = fitfunc,
        type = "unadj", 
        event = event, exposure = exposure,
        level = level, ...)
}

modav <- function(resp, trt, subgr, covars = NULL, data, 
                  fitfunc = c("lm", "glm", "glm.nb", "survreg", "coxph", "rlm"),
                  event, exposure, 
                  level = 0.1, prior = 1, nullprior = 0, ...){
  call <- match.call()
  if(!is.character(call$resp))
    resp <- as.character(call$resp)
  if(!is.character(call$trt))
    trt <- as.character(call$trt)
  umfit(resp = resp, trt = trt, subgr = subgr, covars = covars,
        data = data, fitfunc = fitfunc,
        type = "modav",
        event = event, exposure = exposure, 
        level = level, prior = prior, nullprior = nullprior, ...)
}

