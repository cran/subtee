\name{Simulated data-sets}
\alias{datbin}
\alias{datsurv}
\alias{datnorm}
\alias{datcount}
\docType{data}
\title{ Simulated example data-sets }
\description{
  Simulated test data-sets to illustrate and test methods in this
  package. The underlying simulation truth assumed that there is no
  treatment effect and no subgroup effect.
  
  The data set was simulated using the R code in the tests/data-sets.R
  file.  The data-sets (\code{datnorm}, \code{datbin}, \code{datsurv},
  \code{datcount}) only differ in the y variable (covariates are the
  same).
}
\usage{
data(datbin)

data(datnorm)

data(datcount)

data(datsurv)
}
\format{
  Data frames with 100 observations on the following variables.
  \describe{
    \item{\code{y}}{Response variable}
    \item{\code{treat}}{Categorical variable}
    \item{\code{height}}{Numeric variable}
    \item{\code{labvalue}}{Numeric variable}
    \item{\code{region}}{Categorical variable}
    \item{\code{smoker}}{Categorical variable}
    \item{\code{x1}}{Numeric variable}
    \item{\code{x2}}{Numeric variable}
  }
}
\keyword{datasets}
