\name{summary.subtee}
\alias{summary.subtee}
\title{
Summarizing subgroup analyses estimates
}
\description{
Summary function for subtee objects. Shows estimates and confidence interval boundaries for all candidate subgroups (and complements).
}
\usage{
\method{summary}{subtee}(object, ...)
}
\arguments{
  \item{object}{
An object of class \code{subtee}, usually a result of a call to \code{\link{modav}}, \code{\link{unadj}} or \code{\link{bagged}}.
}
\item{...}{
Not used.
}
}
\value{
A dataframe containing information about treatment effects and group sizes in candidate subgroups.
}

\seealso{
\code{\link{plot.subtee}}
}
\examples{
data(datnorm)
cand.groups <- subbuild(datnorm, height, labvalue, region, smoker)
fitd <- cbind(datnorm, cand.groups)
subgr <- colnames(cand.groups)
res <- modav(resp = "y", trt = "treat", subgr = subgr, data = fitd, 
             covars = ~ x1 + x2, fitfunc = "lm")
summary(res)
}

