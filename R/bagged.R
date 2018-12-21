bagged <- function(resp, trt, subgr, covars = NULL, data, 
                   fitfunc = c("lm", "glm", "glm.nb", "survreg", "coxph", "rlm"),
                   event, exposure, 
                   level = 0.1,
                   B = 2000, mc.cores = 1, stratified = TRUE, 
                   select.by = c("BIC", "AIC"), ...) {
  ## argument checking
  ## allow specification of non-character resp and trt
  ## (do not use this in non-interactive use of margMods)
  select.by = match.arg(select.by)
  fitfunc  <- match.arg(fitfunc)
  
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
  call  <- match.call()
  dots. <- list(...)
  if(!is.character(call$resp))
    resp <- as.character(call$resp)
  if(!is.character(call$trt))
    trt <- as.character(call$trt)
  
  ## allow specification of non-character resp and trt
  ## (do not use this in non-interactive use of margMods)
  sizes <- checkTrtSub(data, trt, subgr)
  
  
  if(is.null(call$event)) event <- NULL
  if(is.null(call$exposure)) exposure <- NULL

  ## Analysis of original dataset ----------------------------------------------
  ## actual model fitting
  # Perform model fitting through the fitMods function. 
  args <- list(resp = resp, trt = trt,
               subgr = subgr, covars = covars, 
               data = data, fitfunc = fitfunc, 
               event = event, exposure = exposure)
  all.args <- c(args, dots.)
  res <- do.call(fitMods, args = all.args)
  crit <- qnorm(1 - level / 2)
  nSub <- length(subgr)
  ## extract overall estimates
  tmpO <- c(res$fitOverall$ests$ests,
            sqrt(res$fitOverall$ests$ests_vc))
  overallEff <- c(tmpO[1],
                  tmpO[1] - crit * tmpO[2],
                  tmpO[1] + crit * tmpO[2])
  
  # Obtain Treatment*subgroup interactions
  results <- do.call(rbind, lapply(res$subgroups, function(ii){
    est <- res$fitmods[[ii]]$ests
    c(est$ests[2], est$ests_vc[2,2], est$pval[2], est$AIC, est$BIC, est$ests[1])
  }))
  results = cbind(results, results[, 1] + results[, 2])
  colnames(results) <- c("estimate", "std.error", "p.value", "AIC", "BIC", "ComplEff", "SubgrEff")
  rownames(results) <- res$subgroups
  overall <- c(estimate = NA, 
               std.error = NA,
               p.value = NA,
               AIC = res$fitOverall$ests$AIC,
               BIC = res$fitOverall$ests$BIC, 
               ComplEff = res$fitOverall$ests$ests[1],
               SubgrEff = res$fitOverall$ests$ests[1])
  results <- rbind(results, overall)
  selected <- 1*(results[, select.by] == min(results[, select.by]))
  results <- data.frame(subgroup = c(res$subgroups, "none"), results, selected)
  
  ## Generate bootstrap datasets -----------------------------------------------
  # Order the dataset by treatment
  
  if (!any(names(data) == trt)){
    data[trt] = model.frame(as.formula(paste0("~ ", trt)), data = data)
  }
  trt_order = order(data[trt])
  ordered_data  <- data[trt_order, ]
  N <- nrow(ordered_data)
  if (!stratified){
    numberhits <- lapply(1:B, function(i){
      as.vector(rmultinom(1, N, rep(1, times = N)))
    })
  } else {
    numberhits <- lapply(1:B, function(i){
      N_t = sum(ordered_data[trt])
      N_c = N - N_t
      c(as.vector(rmultinom(1, N_c, rep(1, times = N_c))),
        as.vector(rmultinom(1, N_t, rep(1, times = N_t))))
    })
  }
  ## Analyze bootstrap samples -------------------------------------------------
  if (mc.cores > 1){
    boot_results <- parallel::mclapply(numberhits, function(x){
                                       fitMods_boot(x,
                                                    ordered_data = ordered_data, 
                                                    resp = resp, trt = trt,
                                                    subgr = subgr, covars = covars, 
                                                    fitfunc = fitfunc, 
                                                    event = event, exposure = exposure,
                                                    select.by = select.by,
                                                    dots. = dots.)
    }, mc.cores = mc.cores)
  } else {
    boot_results <- lapply(numberhits, function(x){
      fitMods_boot(x,
                   ordered_data = ordered_data, 
                   resp = resp, trt = trt,
                   subgr = subgr, covars = covars, 
                   fitfunc = fitfunc, 
                   event = event, exposure = exposure,
                   select.by = select.by,
                   dots. = dots.)
    })
  }
  
  ## Get the estimates of the interaction with the minimum SBC
  boot_results_ext <- t(sapply(boot_results, 
                               function(x) x[which.min(x[, select.by]), 
                                             c("subgroup", "estimate", "SubgrEff", "ComplEff")]))
  ## Transform it into a dataframe for better handling
  boot_results_ext <- data.frame(subgroup = unlist(boot_results_ext[, 1]), 
                                 estimate = unlist(boot_results_ext[, 2]), 
                                 SubgrEff = unlist(boot_results_ext[, 3]), 
                                 ComplEff = unlist(boot_results_ext[, 4]))
  
  InteractionEst = baggingSoupr("estimate", results, boot_results, boot_results_ext, subgr, numberhits, select.by, B, N, level)
  SubgrEff       = baggingSoupr("SubgrEff", results, boot_results, boot_results_ext, subgr, numberhits, select.by, B, N, level)
  ComplEff       = baggingSoupr("ComplEff", results, boot_results, boot_results_ext, subgr, numberhits, select.by, B, N, level)
  
  not_selected = paste(subgr[(which(InteractionEst$percent_selected[-(nSub+1)] == 0))])
  if (length(not_selected) > 0) message(paste0("The subgroup(s) ", paste(not_selected, collapse = ", "),
                                               " were not selected in the bootstrap samples and a corrected estimate is not available. Returning NA."))
  if(any(is.nan(InteractionEst$se.bagg_red)) & B < 2000) message("The variance for the bootstrap estimate of one or more subgroups could not be calculated because they were not selected in sufficient bootstrap samples. This may be due to a small number of bootstrap samples (B) or simply that the subgroup is not predictive.")
  
  postpr = unlist(InteractionEst["percent_selected"])
  model_fits <- results 
  rownames(model_fits) <- NULL
  ## produce estimates for subgroups
  lst <- vector("list", nSub)
  for(i in 1:nSub){
    est <- numeric(3)
    resFit <- res$fitmods[[i]]$ests
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
  obj$overall <- overallEff 
  names(postpr) <- c(subgr, "overall")
  obj$post.weights <- postpr
  names(pvals) <- subgr
  obj$pvals <- pvals
  obj$type <- "bagged"
  obj$fitfunc <- fitfunc
  obj$boot_results <- list(selected = names(which(selected == 1)),
                           unadjusted = results,
                           percSel = postpr[-(nSub+1)],
                           B = B,
                           InteractionEst = InteractionEst, 
                           SubgrEff = SubgrEff, 
                           ComplEff = ComplEff)
  class(obj) <- "subtee"
  obj
}



baggingSoupr = function(target, results, boot_results, boot_results_ext, 
                        subgr, numberhits, select.by, B, N, level){
  ## Calculate mean_boot_estimate pstar mean_u_boot_estimate
  boot_estimate <- sapply(boot_results, function(x) x[[target]])
  mean_boot_estimate <- rowMeans(boot_estimate)
  pstar <- as.vector(table(factor(boot_results_ext$subgroup, 
                                  levels = c(subgr, "none"))))/B
  mean_u_boot_estimate <- tapply(boot_results_ext[[target]], 
                                 factor(boot_results_ext$subgroup, 
                                        levels = c(subgr, "none")), mean)
  s_boot_estimate <- matrixStats::rowSds(sapply(boot_results, 
                                                function(x) x[[target]]))
  ## Put everything in a dataframe
  boot_mean <- data.frame(mean_boot_estimate, pstar, mean_u_boot_estimate, s_boot_estimate)
  
  # Combine results of the bootstrap samples with the estimates in the original
  # sample
  boot_mean <- cbind(boot_mean, results[, c(target, "std.error")])
  boot_mean$mean_u_estimate <- 2 * boot_mean[[target]] - boot_mean$mean_u_boot_estimate
  boot_mean$mean_u_bagg_estimate <- 2 * boot_mean$mean_boot_estimate - boot_mean$mean_u_boot_estimate

  ## Labels 
  ## mean_boot_estimate = 'Bagged estimate' 
  ## s_boot_estimate = 'Bootstrap Standard error' 
  ## mean_u_estimate = 'Bias reduced estimate' 
  ## mean_u_bagg_estimate = 'Bias reduced bagged estimate'
  
  res1 <- boot_mean[order(boot_mean$pstar, decreasing = T), ]
  res1
  
  ## Variance of aggregated estimates ------------------------------------------
  ## Create matrix with the Nbn-1 called boot_number The matrix has one row per
  ## bootstrap sample And one column per individual
  Nbn1 <- lapply(numberhits, function(x) x - 1)
  boot_number <- matrix(unlist(Nbn1), ncol = N, byrow = TRUE)
  
  ## Creat matrix with the delta*_b, that is, the bootstrap estimate per 
  ##  bootsrap sample. The estimate is 0 if the subgroup was not the one with 
  ##  minimum SBC Number of colums is the number of bootstrap samples B Number
  ##  of rows is the number of subgroups
  boot_estimate_mins <- boot_estimate
  minSBCRow <- sapply(boot_results, function(x) which.min(x[[select.by]]))
  for (i in 1:B) {
    boot_estimate_mins[-minSBCRow[i], i] <- 0
  }
  
  u_boot_estimate <- boot_estimate_mins
  
  ## Creat matrix with the delta*_b-delta* Number of colums is the number of
  ## bootstrap samples B Number of rows is the number of subgroups
  diff_u_estimate <- u_boot_estimate/pstar - as.vector(mean_u_boot_estimate)
  
  diff_u_bagg_estimate <- (2 * boot_estimate - u_boot_estimate/pstar) - 
    as.vector(boot_mean$mean_u_bagg_estimate)
  
  K <- length(subgr)+1
  diff_u_estimate_l <- list(0)
  
  Vk <- Ak <- V_baggk <- A_baggk <- se_k <- se_bagg_k <- vector()

  for (k in 1:K) {
    diff_u_estimate_l <- matrix(t(diff_u_estimate)[, k], ncol = N, nrow = B, 
                                byrow = F)
    diff_u_bagg_estimate_l <- matrix(t(diff_u_bagg_estimate)[, k], ncol = N, 
                                     nrow = B, byrow = F)
    cov_u <- vector(length = N)
    cov_u_bagg <- vector(length = N)
    for (n in 1:N) {
      cov_u[n] <- crossprod(boot_number[, n], diff_u_estimate_l[, n])
      cov_u_bagg[n] <- crossprod(boot_number[, n], diff_u_bagg_estimate_l[, n])
    }
    sum(cov_u^2)
    Ak[k] <- sum(cov_u^2)
    Vk[k] <- variance_u <- sum(cov_u^2)/(B * B)
    A_baggk[k] <- sum(cov_u_bagg^2)
    V_baggk[k] <- variance_u_bagg <- sum(cov_u_bagg^2)/(B * B)
    
    bias_u_diff <- colSums(diff_u_estimate_l^2)
    bias_u_diff_bagg <- colSums(diff_u_bagg_estimate_l^2)
    bias_u_diff_b <- bias_u_diff/B/B
    bias_u_diff_bagg_b <- bias_u_diff_bagg/B/B
    bias_u <- sum(bias_u_diff_b)
    bias_u_bagg <- sum(bias_u_diff_bagg_b)
    se_k[k] <- se_u <-   suppressWarnings(sqrt(variance_u - bias_u))
    se_bagg_k[k] <- se_u_bagg <- suppressWarnings(sqrt(variance_u_bagg - bias_u_bagg))
  }
  
  bagged_results <- data.frame(subgroup = results$subgroup, 
                               selected = as.numeric(pstar * 100), 
                               mean_boot_estimate = as.numeric(mean_boot_estimate), 
                               s_boot_estimate = as.numeric(s_boot_estimate), 
                               mean_u_estimate = as.numeric(boot_mean$mean_u_estimate), 
                               se_k = as.numeric(se_k), 
                               mean_u_bagg_estimate = as.numeric(boot_mean$mean_u_bagg_estimate),
                               se_bagg_k = as.numeric(se_bagg_k))
  
  rownames(bagged_results) <- NULL
  colnames(bagged_results) <- c("subgroup", "percent_selected", 
                                "bagg",     "se.bagg",
                                "boot_red", "se.boot_red",
                                "bagg_red", "se.bagg_red")
  bagged_results
}





fitMods_boot = function(numberhits, ordered_data, 
                        resp, trt,
                        subgr, covars, 
                        fitfunc, 
                        event, exposure,
                        select.by,
                        dots.) {
  data_boot <- ordered_data[which(numberhits > 0), ]
  ww <- numberhits[which(numberhits > 0)]
  args <- list(resp = resp, trt = trt,
               subgr = subgr, covars = covars, 
               data = data_boot, fitfunc = fitfunc, 
               event = event, exposure = exposure,
               weights = ww)
  all.args <- c(args, dots.)
  res <- do.call(fitMods, args = all.args)
  results <- do.call(rbind, lapply(res$subgroups, function(ii){
    x <- res$fitmods[[ii]]$model
    p <- length(x$coef)
    coef <- x$coefficients[p]
    complEff <- res$fitmods[[ii]]$ests$ests[1]
    se <- sqrt(diag(x$var)[p])
    tmp <- c(coef, se, 1 - pchisq((coef/se)^2, 1), AIC(x), BIC(x), complEff, complEff + coef)
    tmp
  }))
  rownames(results) <- res$subgroups
  colnames(results) <- c("estimate", "std.error", "p.value", "AIC", "BIC", "ComplEff", "SubgrEff")
  overall <- c(estimate = NA, 
               std.error = NA,
               p.value = NA,
               AIC = res$fitOverall$ests$AIC,
               BIC = res$fitOverall$ests$BIC,
               ComplEff = res$fitOverall$ests$ests,
               SubgrEff = res$fitOverall$ests$ests)
  results <- rbind(results, overall)
  selected <- 1*(results[, select.by] == min(results[, select.by]))
  results <- data.frame(subgroup = c(res$subgroups, "none"), results, selected)
  results
}
