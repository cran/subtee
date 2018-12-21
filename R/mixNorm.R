## functions to deal with mixture of normals distributions
pmixnorm <- function(x, mean, sd, pi){
  out <- numeric(length(x))
  for(i in 1:length(mean)){
    out <- out + pi[i]*pnorm(x, mean[i], sd[i])
  }
  out
}

dmixnorm <- function(x, mean, sd, pi){
  out <- numeric(length(x))
  for(i in 1:length(mean)){
    out <- out + pi[i]*dnorm(x, mean[i], sd[i])
  }
  out
}

qmixnorm <- function(p, mean, sd, pi){
  out <- numeric(length(p))
  int <- c(min(mean)-3*max(sd), max(mean)+3*max(sd))
  sapply(p, function(x){
    uniroot(function(y)
            pmixnorm(y, mean, sd, pi)-x,
            int, extendInt = "yes")$root
  })
}
