## function used inside subbuild to generate subgroup column names
nameSubs <- function(x, varname, fac) {
  if (!fac) {
    n.groups <- length(x) - 1
    subnames <- numeric(n.groups)
    subnames[1] <- paste0(varname, "<=", x[2])
    subnames[n.groups] <- paste0(varname, ">", x[n.groups])
    if (n.groups > 2) {
      for (i in 2 : (n.groups - 1)) 
        subnames[i] <- paste0(x[i], "<", varname, "<=", x[i + 1])
    }
  } else {
    subnames <- paste0(varname, ".", substr(x, 1, 12)) ## use at most 12 chars
  }
  subnames
}

## function to perform pretty rounding (something similar is used in the cut function)
rndPretty <- function(x, digits){
  mid <- x[-c(1, length(x))]
  nb <- length(mid)
  
  for(dig in digits:max(12, digits)) {
    ## 0+ avoids printing signed zeros as "-0"
    ch.br <- formatC(0+mid, digits = dig, width = 1)
    if(all(ch.br[-1L] != ch.br[-nb]))
      break
  }
  c(x[1], ch.br, x[length(x)])
}


tidy.names <- function(x){
  x <- gsub("==", ".", x)
  x <- gsub(" ", "", x)
  x <- gsub("<=", ".leq.", x)
  x <- gsub(">=", ".geq.", x)
  x <- gsub("<", ".l.", x)
  x <- gsub(">", ".g.", x)
  x <- gsub("!=", ".neq.", x)
  x <- gsub("\"", "", x)
  x <- gsub("-", "m", x)
  x <- gsub("%in%", ".in.", x)
  x <- gsub("`", "", x)
  make.names(x, unique = TRUE)
}


## subbuild - building subgroups out of covariates
subbuild <- function(data, ..., n.cuts = 2, dig.lab = 3, dupl.rm = FALSE, 
                     make.valid.names = FALSE){
  
  mc <- match.call()
  tmp <- as.list(mc)[-(1:2)]
  subGr <- tmp[!(names(tmp) %in% c("n.cuts", "dig.lab", "dupl.rm", "make.valid.names"))]
  if(length(subGr) == 0){
    out <- buildAuto(data, n.cuts, dig.lab)
    
  } else{
    ll <- lapply(subGr, function(sub){ # find out which are logical expresssions and which covariates
      r <- as.numeric(eval(sub, data))
      if(all(r %in% c(0,1)))
        return(r)
      return(sub)
    })
    ind <- unlist(lapply(ll, is.numeric))
    ll1 <- ll[ind]
    out1 <- do.call(cbind, ll1)
    nams <- as.character(subGr[ind])
    if(length(dim(out1)) > 1)
      colnames(out1) <- nams
    # now build group definitions "by hand" for the rest
    ll2 <- ll[!ind]
    ll2.var <- lapply(ll2, function(x) eval(x, data))
    autoGr <- as.data.frame(ll2.var)
    out2 <- NULL
    if(ncol(autoGr) > 0){
      colnames(autoGr) <- as.character(ll2)
      out2 <- buildAuto(as.data.frame(autoGr), n.cuts, dig.lab)
    }
    out <- cbind(out1, out2)
  }
  
  if(dupl.rm){
    out.tmp1 <- out
    dupl.cols1 <- which(duplicated(out.tmp1, MARGIN = 2))
    out.tmp2 <- unique(out.tmp1, MARGIN = 2)
    dupl.cols2.mat <- apply(out.tmp2, 2, function(x){
      apply(out.tmp2, 2, function(y) all(y == (1-x)))
    })
    dupl.cols2 <- unlist(apply(which(dupl.cols2.mat, arr.ind = TRUE), 1, function(x) if(x[1] < x[2]) x[1]))
    if(length(dupl.cols2) > 0)
      out <- out.tmp2[, -dupl.cols2]
    rm.cols <- c(colnames(out.tmp1)[dupl.cols1], colnames(out.tmp2)[dupl.cols2])
    if(length(rm.cols) > 0)
      message(paste("removed duplicate columns:", paste(rm.cols, collapse = ", ")))
  }
  if(make.valid.names)
    colnames(out) <- tidy.names(colnames(out))
  as.data.frame(out)
}

# basically the old subbuild function, with some changes to catch too few values
# cases
buildAuto <- function(data, n.cuts, dig.lab){
  
  stopifnot(is.data.frame(data), is.numeric(n.cuts), is.numeric(dig.lab), 
            n.cuts > 0, dig.lab >= 0, length(n.cuts) == 1, length(dig.lab) == 1,
            n.cuts %% 1 == 0, dig.lab %% 1 == 0)
  
  k <- ncol(data)
  n <- nrow(data)
  
  types <- sapply(data, class)
  num <- types == "numeric" | types == "integer"
  fac <- !num
  data[, fac] <- lapply(data[, fac, drop = FALSE], factor)
  fac.lev <- sapply(data[, fac, drop = FALSE], levels, simplify = FALSE)
  var.names <- colnames(data)
  if(is.null(var.names))
    var.names <- paste0("x", 1:k)
  cut.names <- vector("list", k)
  names(cut.names) <- var.names
  cut.names[which(fac)] <- fac.lev
  
  vars.fac <- data 
  for(i in which(num)){
    n.cuts.tmp <- n.cuts
    cut.p <- seq(0, 1, 1/(n.cuts.tmp + 1))
    cutoffs.tmp <- quantile(data[, i], probs = cut.p, na.rm = TRUE, type = 1)+
      c(-Inf, rep(0, n.cuts.tmp), Inf)
    
    cutoffs <- rndPretty(cutoffs.tmp, dig.lab)
    while((length(unique(cutoffs)) < (n.cuts.tmp +2)) & (n.cuts.tmp > 0)){
      if(n.cuts.tmp >= 2){
        txt <- sprintf("insufficient number of unique values in %s to generate %i cutoffs,
                     reducing to %i.",
                       var.names[i], n.cuts.tmp, n.cuts.tmp - 1)
        warning(txt, call. = FALSE)
        n.cuts.tmp <- n.cuts.tmp - 1
        cut.p <- seq(0, 1, 1/(n.cuts.tmp + 1))
        cutoffs.tmp <- quantile(data[, i], probs = cut.p, na.rm = TRUE, type = 1) +
          c(-Inf, rep(0, n.cuts.tmp), Inf)
        cutoffs <- rndPretty(cutoffs.tmp, dig.lab)
      } else{
        cut.p <- seq(0, 1, 1/(n.cuts.tmp + 1)) 
        cutoffs.tmp <- quantile(data[, i], probs = cut.p, na.rm = TRUE, type = 1) +
                         c(-Inf, rep(0, n.cuts.tmp), Inf)
        cutoffs <- rndPretty(cutoffs.tmp, dig.lab)
        n.cuts.tmp <- 0
      }
    }
    vars.fac[, i] <- cut(data[, i], breaks = as.numeric(cutoffs))
    cut.names[[i]] <- cutoffs
  }
  naAct <- options()$na.action
  options(na.action="na.pass")
  subs.l <- lapply(as.list(vars.fac), function(x) model.matrix(~ x - 1))
  options(na.action=naAct)
  subs <- matrix(unlist(subs.l), nrow = n)
  subnames <- mapply(nameSubs, cut.names, var.names, fac)
  colnames(subs) <- unlist(subnames)
  subs
}

