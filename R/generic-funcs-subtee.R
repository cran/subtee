print.subtee <- function(x, digits = 4, ...){
  trtEff = x$trtEff
  trtEffDiff = x$trtEffDiff
  if(x$type == "bagged") {
    trtEff = trtEff[which(trtEff$Group == x$boot_results$selected), ]
    trtEffDiff = trtEffDiff[which(trtEffDiff$Group == x$boot_results$selected), ]
    rownames(trtEff) = rownames(trtEffDiff) = NULL
    if(x$boot_results$selected == "overall"){
      return(cat("No selected subgroup. In the bootstrap samples, the model with no interactions was selected",
             round(x$post.weights[x$boot_results$selected], 2),
             "% of the times.\n"))
    }
  }
  cat("Trt. Effect Estimates\n")
  print(trtEff, digits = digits)
  cat("\nDifference in Trt. Effect vs Complement\n")
  print(trtEffDiff, digits = digits)
  if(x$type == "bagged") {
    cat(sprintf("\n %s is the selected subgroup.\n",x$boot_results$selected),
        sprintf("It was selected in %s%% of %s bootstrap samples.\n",
                round(x$post.weights[x$boot_results$selected], 2),
                x$boot_results$B*(1-x$boot_results$boot_invalid/100)))
    if(x$boot_results$boot_invalid>0) 
      cat(sprintf(" %s%% the bootstrap samples could not be used\n",
                  x$boot_results$boot_invalid))
  }
  cat(sprintf("\nSubgroup Models fitted with \"%s\"\n", x$fitfunc))
  if (x$fitfunc=="coxph") cat("Effect estimates in terms of the log-hazard ratios\n")
}

summary.subtee <- function(object, ...){
  class(object) <- "summary.subtee"

  ests <- object$trtEff
  ## add p-value or posterior probability column
  if(object$type == "modav")
    pv <- object$post.weights
  if(object$type == "unadj")
    pv <- c(object$pvals, NA) ## no p-value for overall group
  if(object$type == "bagged")
    pv <- c(object$pvals, NA) ## no p-value for overall group
  
  n <- object$n
  GroupSize <- c(t(object$GroupSize[,1:2]))
  p <- rep(pv[-length(pv)], each=2)
  ests <- cbind(ests[,1:2], GroupSize = GroupSize, ests[,3:5], p=p)
  
  ests$Group <- factor(ests$Group, levels = c(levels(ests$Group), "overall"))
  overall <- object$overall
  ests$Subset <- factor(ests$Subset, levels=c(levels(ests$Subset), "overall"))
  over <- data.frame(Group="overall", Subset="overall",
                     GroupSize=n,
                     LB=overall[2],
                     trtEff=overall[1], UB=overall[3],
                     p=pv[length(pv)])
  ests <- rbind(ests, over)
  colnames(ests)[7] <- ifelse(object$type == "modav",
                              "post.mod.prob", "int.p.val")
  if (object$type == "bagged") ests$perc.selected = rep(object$post.weights, each=2)[-(nrow(ests)+1)]
  rownames(ests) <- 1:nrow(ests)
  ests
}

print.summary.subtee <- function(x, digits = 3, show.compl = FALSE, ...){
  print(x, digits = digits)
}

plot.subtee <- function(x, y = NULL, z = NULL, 
                       type = c("trtEff", "trtEffDiff"),
                       show.compl = FALSE,
                       xlab = "", ylab = "default", main = "default",
                       them,
                       point.size = 2.5, line.size = 1,
                       palette = "default", ...){
  type = match.arg(type)
  compare = !(is.null(y)&is.null(z))
  if(!compare){
    if(palette == "default") palette = "Set1"
    if(type == "trtEff"){
      p = plot_effect(x, show.compl,
                  xlab, ylab, main,
                  them,
                  point.size, line.size,
                  palette, ...)
    }
    if(type == "trtEffDiff"){
      p = plot_interaction(x, 
                       xlab, ylab, main,
                       them,
                       point.size, line.size,
                       palette, ...)
    }
  } else {
    if(palette == "default") palette = "Dark2"
    res.list = list(x)
    if(!is.null(y)) res.list[[2]] = y
    if(!is.null(z)) res.list[[length(res.list)+1]] = z
    
    if(type == "trtEff"){
      p = plot_compare_effect(res.list, show.compl,
                          xlab, ylab, main,
                          them,
                          point.size, line.size,
                          palette, ...)
    }
    if(type == "trtEffDiff"){
      p = plot_compare_interaction(res.list, 
                               xlab, ylab, main,
                               them,
                               point.size, line.size,
                               palette, ...)
    }
  }
  p
}
