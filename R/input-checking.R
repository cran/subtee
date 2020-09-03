checkTrtSub <- function(data, trt, subgr){
  ## input checking and calculating subgroup sizes (and NAs)
  ## check whether trt is a numeric 0,1 variable
  trttmp <- data[[trt]]
  chck <- function(tmp){
    mes <- "trt needs to be a numeric variable with values 0 and 1
  (indicating control and treatment)"
    if(!(class(tmp) %in% c("numeric", "integer"))){
      stop(mes)
    }
    srt <- sort(unique(tmp))
    if((srt[1] != 0 | srt[2] != 1) | length(srt) != 2)
      stop(mes)
  }
  if(!is.null(trttmp)){ # trt defined in data
    chck(trttmp)
  } else { # trt not defined in data
    trttmp <- get(trt)
    chck(trttmp)
  }

  ## check whether all variables in subgr are also in data
  if(!is.character(subgr))
    stop("\"subgr\" needs to be a character variable")
  ind <- match(subgr, colnames(data))
  if(any(is.na(ind))){
    stop(sprintf("Variable(s) \"%s\" specified in \"subgr\" not found in \"data\"",
                 paste(subgr[is.na(ind)], collapse ="\", \"")))
  }
  ## check for whether a subgroup defines the whole data-set and
  ## whether for all subgroups and complements we have patients in control and trt
  nSub <- length(subgr)
  sizes <- complSizes <- NAs <- numeric(nSub)
  for(i in 1:nSub){
    tmpvar <- data[,subgr[i]]
    len <- length(unique(tmpvar))
    if(len == 1){
      stop(sprintf("Subgroup \"%s\" contains all patients or no one.\n",
                   subgr[i]))
    }
    if(!is.numeric(tmpvar)){
      stop(sprintf("Subgroup \"%s\" needs to be numeric.\n",
                   subgr[i]))
    }
    if(!all(sort(unique(tmpvar)) == c(0,1))){
      stop(sprintf("Only values 0 or 1 are allowed (variable \"%s\" does not comply).\n",
                   subgr[i]))
    }
    tab <- table(tmpvar, trttmp)
    if(any(tab < 2)){
      stop(sprintf("Subgroup \"%s\" does not contain at least 1 patient on
  treatment and control in subgroup and complement.\n",
                   subgr[i]))
    }
    sizes[i] <- sum(tmpvar, na.rm=TRUE)
    complSizes[i] <- sum(tmpvar==0, na.rm=TRUE)
    NAs[i] <- sum(is.na(tmpvar))
  }
  out <- cbind(sizes, complSizes, NAs)
  rownames(out) <- subgr
  out
}