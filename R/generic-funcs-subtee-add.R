getConfints = function(fitmods,
                       data,
                       type,
                       nSub,
                       level,
                       postpr,
                       boot_results){
  crit <- qnorm(1-level/2)
  ## extract overall estimates
  tmpO <- c(fitmods$fitOverall$ests$ests,
            sqrt(fitmods$fitOverall$ests$ests_vc))
  names(tmpO) <- NULL
  overall <- c(trtEff=tmpO[1],
               LB=tmpO[1]-crit*tmpO[2],
               UB=tmpO[1]+crit*tmpO[2])
  ## produce estimates for subgroups
  lst <- vector("list", nSub)
  if(type == "modav"){
    prevs <- getPrev(data[, fitmods$subgroups, drop = FALSE])
    tmpS <- tmpC <- tmpIA <- matrix(ncol=2, nrow=nSub+1)
    for(i in 1:nSub){
      for(j in 1:nSub){ # loop over predictions of all models
        resFit <- fitmods$fitmods[[j]]$ests
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
      lst[[i]] <- list(fitmods$fitmods[[i]]$ests$pval[2],
                       est[1:2], ub[1:2], lb[1:2],
                       est[3], ub[3], lb[3])
    }
  }
  if(type == "unadj"){
    for(i in 1:nSub){
      est <- numeric(3)
      resFit <- fitmods$fitmods[[i]]$ests
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
  if(type == "bagged"){
    SubgrEff = boot_results$SubgrEff
    ComplEff = boot_results$ComplEff
    InteractionEst = boot_results$InteractionEst
    for(i in 1:nSub){
      est <- numeric(3)
      resFit <- fitmods$fitmods[[i]]$ests
      tmp1 <- c(SubgrEff[i, "bagg_red"], SubgrEff[i, "se.bagg_red"])
      tmp2 <- c(ComplEff[i, "bagg_red"], ComplEff[i, "se.bagg_red"])
      tmp3 <- c(InteractionEst[i, "bagg_red"], InteractionEst[i, "se.bagg_red"])
      est <- c(tmp1[1], tmp2[1], tmp3[1])
      ub <- est+crit*c(tmp1[2], tmp2[2], tmp3[2])
      lb <- est-crit*c(tmp1[2], tmp2[2], tmp3[2])
      lst[[i]] <- list(resFit$pval[2],
                       est[1:2], ub[1:2], lb[1:2],
                       est[3], ub[3], lb[3])
    }
  }
  
  pvals <- sapply(lst, function(x) x[[1]])
  names(pvals) <- fitmods$subgroups
  trt_eff <- c(sapply(lst, function(x) x[[2]]))
  trt_ub <- c(sapply(lst, function(x) x[[3]]))
  trt_lb <- c(sapply(lst, function(x) x[[4]]))
  ia <- sapply(lst, function(x) x[[5]])
  ia_ub <- sapply(lst, function(x) x[[6]])
  ia_lb <- sapply(lst, function(x) x[[7]])
  
  Group <- sapply(fitmods$fitmods, function(x) x$subgrNam)
  cmpl <- rep(c("Subgroup", "Complement"), nSub)
  
  obj = list()
  obj$trtEff <- data.frame(Group = rep(Group, each=2),
                           Subset = cmpl,
                           LB = trt_lb, trtEff = trt_eff,
                           UB = trt_ub, stringsAsFactors = TRUE)
  obj$trtEffDiff <- data.frame(Group = Group,
                               LB = ia_lb, trtEffDiff = ia,
                               UB = ia_ub, stringsAsFactors = TRUE)
  rownames(obj$trtEffDiff) <- NULL
  obj$pvals = pvals
  obj$overall = overall
  obj
}


confint.subtee <- function(object, parm, level = 0.95, ...){
  if(level < 0.5) warning("confint uses 'level' for the confidence level required")
  level = 1 - level  # We use significance level in subtee
  fitmods <- object$fitMods
  data <- object$data
  type <- object$type
  nSub <- length(object$subgroups)
  crit <- qnorm(1-level/2)
  subgr <- object$subgroups
  ## calculate posterior model weigths (only needed for modav)
  postpr <- object$post.weights
  boot_results <- object$boot_results # will return NULL if not bagged object
  ## produce estimates for subgroups
  cis = getConfints(fitmods, data, type, nSub, level, postpr, boot_results)
  
  obj <- list()
  obj$fitmods <- fitmods
  obj$trtEff <- cis$trtEff
  obj$trtEffDiff <- cis$trtEffDiff
  obj$n <- object$n
  obj$subgroups <- object$subgroups
  obj$GroupSizes <- object$sizes
  obj$level <- level
  obj$overall <- cis$overall
  obj$post.weights <- object$postpr
  obj$pvals <- object$pvals
  obj$type <- object$type
  obj$fitfunc <- object$fitfunc
  if(type == "bagged"){
    obj$boot_results <- object$boot_results
  }
  class(obj) <- "subtee"
  obj
}