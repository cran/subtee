BIC_survreg = function(fit) AIC(fit, k = log(sum(fit$y[,2])))

## helper function to extract trt eff and interaction estimate and
## covariance matrix from all models and also the p-value of the test
## statistic
getEsts <- function(fit, trtNam, subgrNam){
  cf <- coef(fit)
  vc <- vcov(fit)
  cls <- class(fit)[1]
  sfit <- summary(fit)
  if(is.na(subgrNam)){
    ind1 <- match(trtNam, names(cf))
    bic <- switch(cls,
                  survreg = BIC_survreg(fit),
                  BIC(fit))
    return(list(ests = cf[ind1], ests_vc = vc[ind1,ind1],
                AIC=AIC(fit), BIC=bic))
  } else {
    p <- switch(cls,
                survreg = sfit$table[, 4],
                rlm = {
                  rat <- sfit$coefficients[,1]/sfit$coefficients[,2]
                  2*(1-pnorm(abs(rat)))
                }, 
                coxph = {
                  sfit$coefficients[, 5]
                }, {
                  sfit$coefficients[, 4]
                })
    bic <- switch(cls,
                  survreg = BIC_survreg(fit),
                  BIC(fit))
    ind1 <- match(trtNam, names(cf))
    ind2 <- match(sprintf("%s:%s", trtNam, subgrNam),
                  names(cf))
    ind <- c(ind1, ind2)
    return(list(ests = cf[ind], ests_vc = vc[ind,ind],
                pval = p[ind],
                AIC=AIC(fit), BIC=bic))
  }
}

## function to fit all subgroup models
## Inputs:
##        trt - binary trt variable
##        resp - response
##        subgr -  dataframe of candidate subgroups
##        covars - additional covariates included in the models
##        fitfunc - model fitting function;
##                  one of "lm", "glm", "survreg", "coxph" or "glm.nb"
##        event - event variable; need to be specified for fit functions survreg
##                and coxph
##        exposure - needs to be specified for fit function "glm.nb"
fitMods <- function(resp, trt, subgr, covars, data,
                    fitfunc = "lm", event, exposure, ...){

  ## create formula for model fitting
  ## look at special cases: survival data and overdispersed count data
  if(!is.element(fitfunc, c("survreg", "coxph", "glm.nb"))){
    form <- sprintf("%s ~ %s", resp, trt)
  }
  if(is.element(fitfunc, c("survreg", "coxph"))){
    if(missing(event))
      stop("need to specify event variable for survreg or coxph")
    if(!is.character(event))
      stop("event needs to be a character variable")
    form <- sprintf("Surv(%s, %s) ~ %s",
                    resp, event, trt)
  }
  if(fitfunc == "glm.nb"){
    if(missing(exposure))
      stop("need to specify exposure variable for glm.nb")
    if(!is.character(exposure))
      stop("exposure needs to be a character variable")
    form <- sprintf("%s ~ %s + offset(log(%s))",
                    resp, trt, exposure)
  }

  ## translate character to function object
  fitf <- get(fitfunc)

  ## fit all models
  nSub <- length(subgr)
  fitmods <- vector("list", nSub)
  if(is.null(covars)){
    progs <- NULL
  } else {
    progs <- attr(terms(covars), "term.labels")
  }
  for(i in 1:nSub){
    bsub <- sprintf(".subgroup__%s", i)
    assign(bsub, data[, subgr[i]])
    progDiff <- setdiff(setdiff(progs, subgr[i]), sprintf("`%s`", subgr[i]))
    progCovSub <- paste(c(progDiff, bsub), collapse = " + ")
    formSub <- sprintf("%s + %s + %s*%s", form, progCovSub,
                       trt, bsub)
    fit <- fitf(as.formula(formSub), data=data, ...)
    ests <- getEsts(fit, trt, bsub)
    if(any(is.na(ests$ests))){
      warning(sprintf("NA in treatment effect estimate in subgroup
model for variable \"%s\".", subgr[i]))
    }
    lst <- list(model=fit,
                ests=ests,
                subgrNam=subgr[i])
    fitmods[[i]] <- lst
  }
  names(fitmods) <- subgr

  ## fit overall model
  formOv <- form
  if (!is.null(progs)) {
    progCov <- paste(unique(progs),  collapse = " + ")
    formOv <- sprintf("%s + %s", form, progCov)
  }
  fit <- fitf(as.formula(formOv), data=data, ...)
  ests <- getEsts(fit, trt, NA)
  fitOverall <- list(model=fit, ests=ests)

  ## create output list
  out <- list()
  out$subgroups <- subgr
  out$fitmods <- fitmods
  out$fitOverall <- fitOverall
  out$trtNam <- trt
  out
}

