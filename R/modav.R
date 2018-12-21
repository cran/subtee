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
  crit <- qnorm(1-level/2)
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
  ## extract overall estimates
  tmpO <- c(res$fitOverall$ests$ests,
            sqrt(res$fitOverall$ests$ests_vc))
  names(tmpO) <- NULL
  overall <- c(trtEff=tmpO[1],
               LB=tmpO[1]-crit*tmpO[2],
               UB=tmpO[1]+crit*tmpO[2])

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
  lst <- vector("list", nSub)
  
  if(type == "modav"){
    prevs <- getPrev(data[,subgr, drop=FALSE])
    tmpS <- tmpC <- tmpIA <- matrix(ncol=2, nrow=nSub+1)
    for(i in 1:nSub){
      for(j in 1:nSub){ # loop over predictions of all models
        resFit <- res$fitmods[[j]]$ests
        tmpS[j,] <- predSub(resFit$ests,
                            resFit$ests_vc,
                            w=prevs$wS[i,j])
        tmpC[j,] <- predSub(resFit$ests,
                            resFit$ests_vc,
                            w=prevs$wC[i,j])
        tmpIA[j,] <- predIA(resFit$ests, resFit$ests_vc,
                            prevs$wS[i,j], prevs$wC[i,j])
      }
      tmpS[nSub+1,] <- tmpC[nSub+1,] <- tmpO # add the overall estimate
      tmpIA[nSub+1,] <- c(0,0)
      est <- c(qmixnorm(0.5, tmpS[,1], tmpS[,2], postpr),
               qmixnorm(0.5, tmpC[,1], tmpC[,2], postpr),
               qmixnorm(0.5, tmpIA[,1], tmpIA[,2], postpr))
      ub <- c(qmixnorm(1-level/2, tmpS[,1], tmpS[,2], postpr),
              qmixnorm(1-level/2, tmpC[,1], tmpC[,2], postpr),
              qmixnorm(1-level/2, tmpIA[,1], tmpIA[,2], postpr))
      lb <- c(qmixnorm(level/2, tmpS[,1], tmpS[,2], postpr),
              qmixnorm(level/2, tmpC[,1], tmpC[,2], postpr),
              qmixnorm(level/2, tmpIA[,1], tmpIA[,2], postpr))
      lst[[i]] <- list(res$fitmods[[i]]$ests$pval[2],
                       est[1:2], ub[1:2], lb[1:2],
                       est[3], ub[3], lb[3])
    }
  }
  if(type == "unadj"){
    for(i in 1:nSub){
      est <- numeric(3)
      resFit <- res$fitmods[[i]]$ests
      tmp1 <- predSub(resFit$ests,
                      resFit$ests_vc,
                      w=1)
      tmp2 <- predSub(resFit$ests,
                      resFit$ests_vc,
                      w=0)
      tmp3 <- c(resFit$ests[2], sqrt(resFit$ests_vc[2,2]))
      est <- c(tmp1[1], tmp2[1], tmp3[1])
      ub <- est+crit*c(tmp1[2], tmp2[2], tmp3[2])
      lb <- est-crit*c(tmp1[2], tmp2[2], tmp3[2])
      lst[[i]] <- list(resFit$pval[2],
                       est[1:2], ub[1:2], lb[1:2],
                       est[3], ub[3], lb[3])
    }
  }

  pvals <- sapply(lst, function(x) x[[1]])
  trt_eff <- c(sapply(lst, function(x) x[[2]]))
  trt_ub <- c(sapply(lst, function(x) x[[3]]))
  trt_lb <- c(sapply(lst, function(x) x[[4]]))
  ia <- sapply(lst, function(x) x[[5]])
  ia_ub <- sapply(lst, function(x) x[[6]])
  ia_lb <- sapply(lst, function(x) x[[7]])
  fitmods <- lapply(res$fitmods, function(x) x$model)
  fitmods[["overall"]] <- res$fitOverall$model
  
  Group <- sapply(res$fitmods, function(x) x$subgrNam)
  cmpl <- rep(c("Subgroup", "Complement"), nSub)
  
  obj <- list()
  obj$fitmods <- fitmods
  obj$trtEff <- data.frame(Group = rep(Group, each=2),
                           Subset = cmpl,
                           LB = trt_lb, trtEff = trt_eff,
                           UB = trt_ub)
  obj$trtEffDiff <- data.frame(Group = Group,
                               LB = ia_lb, trtEffDiff = ia,
                               UB = ia_ub)
  rownames(obj$trtEffDiff) <- NULL
  obj$n <- nrow(data)
  obj$subgroups <- res$subgroups
  obj$GroupSizes <- sizes
  obj$level <- level
  obj$overall <- overall 
  names(postpr) <- c(subgr, "overall")
  obj$post.weights <- postpr
  names(pvals) <- subgr
  obj$pvals <- pvals
  obj$type <- type
  obj$fitfunc <- fitfunc
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

