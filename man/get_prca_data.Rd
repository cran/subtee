\name{get_prca_data}
\alias{get_prca_data}
\title{
  Downloads the prca dataset to use in the package's examples (Internet connection is required).
}
\description{
  Fetchs the Advanced prostate cancer data from Patrick Royston and Willi Sauerbrei' book Multivariable Model-building. Thw data comes from from a clinical trial with prostate carcinoma patients described in David P. Byar and Sylvan B. Green (1980).
}
\usage{
get_prca_data()
}

\value{
  A data frame with 475 rows and 15 variables:
  PATNR, AGE,  WT, SBP, DBP, SZ, AP, HG, SG, PF, HX, BM, STAGE, EKG, RX, SURVTIME, CENS, X_ST, X_D, X_T, X_T0
}
\references{
Royston, P., & Sauerbrei, W. (2008). Multivariable model-building: a pragmatic 
approach to regression anaylsis based on fractional polynomials for modelling 
continuous variables (Vol. 777). John Wiley & Sons.

Byar, D. P., & Green, S. B. (1980). The choice of treatment for cancer patients 
based on covariate information. Bulletin du cancer, 67(4), 477.
}

\seealso{
\code{\link{bagged}}, \code{\link{unadj}}, \code{\link{modav}}
}
\examples{
prca = get_prca_data()
head(prca)
}

