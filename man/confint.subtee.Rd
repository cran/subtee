\name{confint.subtee}
\alias{confint.subtee}
\title{
Confidence intervals for treatment effect estimates
}
\description{
Computes confidence intervals for subtee objects. This allows the recalculation of confidence intervals at a desired levels without fitting the models again, which is particularly useful for the results of the \code{\link{bagged}} function.
}
\usage{
\method{confint}{subtee}(object, parm, level = 0.95, ...)
}
\arguments{
  \item{object}{
An object of class \code{subtee}, usually a result of a call to \code{\link{modav}}, \code{\link{unadj}} or \code{\link{bagged}}.
}
\item{parm}{
Not used.
}
\item{level}{
the confidence level required. Note that this is 1 - the significance level used in \code{\link{modav}}, \code{\link{unadj}} or \code{\link{bagged}}.
}
\item{...}{
Not used.
}
}
\value{
An object of class \code{subtee}
}

\seealso{
\code{\link{plot.subtee}}, \code{\link{summary.subtee}}
}
\examples{
data(datnorm)
cand.groups <- subbuild(datnorm, height, labvalue, region, smoker)
fitd <- cbind(datnorm, cand.groups)
subgr <- colnames(cand.groups)
res <- modav(resp = "y", trt = "treat", subgr = subgr, data = fitd, 
             covars = ~ x1 + x2, fitfunc = "lm")
confint(res, level = 0.80)
}

