\name{bagged}
\alias{bagged}
\title{
Bootstrap estimates for interaction terms in exploratory subgroup analyses
}
\description{
  Perform exploratory subgroup analysis using bootstrap bias adjustment
  as described in Rosenkranz (2016).  The function fits a GLM or a Cox
  model in the data and then performs bootstrap samples to correct for
  bias.
}
\usage{
bagged(resp, trt, subgr, covars = NULL, data, 
       fitfunc = c("lm", "glm", "glm.nb", "survreg", "coxph", "rlm"),
       event, exposure, 
       level = 0.1,
       B = 2000, mc.cores = 1, stratified = TRUE, 
       select.by = c("BIC", "AIC"), ...)
}
\arguments{
\item{resp}{
Character giving the name of the response variable. The variable can be
either defined in the global environment or in the data-set \code{data}
specified below. For interactive use it is also possible to use unquoted
names (i.e. \code{bagged(resp,...)} instead of \code{bagged("resp",...)}),
avoid this for non-interactive use of the function.
}
\item{trt}{
  Character giving the name of the treatment variable. The variable can
  be either defined in the global environment or in the data-set
  \code{data} specified below. Note that the treatment variable itself
  needs to be defined as a numeric variable, with control coded as 0, and
  treatment coded as 1. For interactive use it is also possible to use unquoted
  names (as for the resp argument. see above).
}
  \item{subgr}{
Character vector giving the variable names in \code{data} to use as
subgroup identifiers. Note that the subgroup variables in \code{data}
need to be numeric 0-1 variables.
}
  \item{covars}{
Formula, specifying additional covariates to be included in the models
(need to be available in \code{data}).
}
  \item{data}{
    Data frame containing the variables referenced in \code{resp},
    \code{trt}, \code{subgr} and \code{covars}
    (and possibly \code{event} and \code{exposure}).
}
  \item{fitfunc}{
Model fitting functions. Currrently one of \code{'lm'}, \code{'glm'},
\code{'glm.nb'}, \code{'survreg'}, \code{'coxph'} or \code{'rlm'}.
}
  \item{event}{
Character giving the name of the event variable. Has to be specified
when using fit functions \code{'survreg'} and \code{'coxph'}. The variable can be
either defined in the global environment or in the data-set \code{data}.
}
  \item{exposure}{
    Character giving the name of the exposure variable, needed for
    negative binomial regression, when using fit functions
\code{'glm.nb'}. This is typically the time each patient is exposed to the drug.
    The fitted model uses the call
    \code{glm.nb(.~.+offset(log(exposure)))}. The variable
needs to be defined either in the global environment or in the data-set \code{data}.
}
\item{level}{
Significance level for confidence intervals will be calculated for
treatment effect estimates.
}
\item{B}{
 A numeric input. The number of bootstrap samples to perform.
}
\item{mc.cores}{
 A numeric input. This argument is passed to the \code{\link{mclapply}} function to perform computations in parallel.
 If \code{mc.cores = 1}, then \code{lapply} is used.
}
\item{stratified}{
  Should the bootstrap resampling be done stratifying by treatment group? (default: TRUE).
}
\item{select.by}{
  Should the model selection be done using BIC or AIC? (default: BIC).
}
\item{\dots}{
 Other arguments passed to the model fitting function.
}
}
\details{
In the generalized linear model case, \eqn{P} generalized linear models are fitted such that 
\deqn{
M_p: h(\mu_{pi}) = \alpha_p + \beta_p z_i +
(\gamma_p + \delta_p z_i) s_{pi} +  \sum_{k = 1}^{K} \tau_k x_{ik}
}

where \eqn{h} is the link function, \eqn{\mu_{pi} = E_p[Y_i]} is the expectation
of the response \eqn{Y_i} under model \eqn{M_p} and \eqn{x_{ik}} are additional 
covariates we control for. For survival data, a proportional
hazards model can be used:
\deqn{
M_p: \lambda_{pi}(t)= \lambda_{p0}(t) \exp\left\{\beta_p z_i + 
(\gamma_p + \delta_p z_i)s_{pi} +  \sum_{k = 1}^{K} \tau_k x_{ik} \right\}
}

The focus of estimation is the difference in the treatment effect
between a subgroup and its complement, the treatment by subgroup interaction
\eqn{\delta_p}.

Consider now \eqn{B} bootstrap samples from the original data. 
Let \eqn{(Y_{b1}^*, ..., Y_{bN}^*)} be a bootstrap sample from the original data.
Let \eqn{(z_{b1}^*, ..., z_{bN}^*)}, \eqn{(s_{b1}^*, ..., s_{bN}^*)}, and 
\eqn{(x_{b1k}^*, ..., x_{bNk}^*)} be the corresponding treatment indicators, group indicators, 
and covariates in the bootstrap samples, respectively. 
For each \eqn{p=1, ..., P} and \eqn{b=1, ..., B}
we fit the model:
\deqn{
h(E_p[Y_{bi}^*]) = \alpha^{*}_{bp} + \beta^{*}_{bp} z_{bi}^* + 
  (\gamma^{*}_{bp} + \delta^{*}_{bp} z_{bi}^*) s^{*}_{bpi}
  +  \sum_{k = 1}^{K} \tau_k x^*_{ik}
}

An estimator of \eqn{\delta_{p}} given that subgroup \eqn{S_p} provided the best fit 
can be calculated as
\deqn{
\overline{\delta}^{*}_p =
	\frac{\sum_{b=1}^B u_{bp} \hat{\delta}^{*}_{bp} }{\sum_{b=1}^B u_{bp}}
}
where \eqn{\hat{\delta}^{*}_{bp}} is the usual estimator of \eqn{\delta_{bp}} and
\eqn{u_{bp} = 1} if the subgroup \eqn{p} provides the best fit for bootstrap
sample \eqn{b} and 0 otherwise.

An bias-reduced estimator of \eqn{\delta_{p}} can be obtained as:
\deqn{
\check{\delta}^{*}_p = 2 \hat{\delta}_p - \overline{\delta}^{*}_p 
}

A bias-reduced estimator with decreased variability is obtained by
replacing the maximum likelihood estimator by the bagging estimator
\eqn{\hat{\delta}^{*}_{bp}}:

\deqn{
 \hat{\delta}^{*}_p = \frac{1}{B}\sum_{b=1}^B \hat{\delta}^{*}_{bp}
}
so that the bias reduced estimator is
\deqn{
\tilde{\delta}^{*}_p = 2 \hat{\delta}^*_p - \overline{\delta}^{*}_p
}
}
\value{
An object of class subtee. A list containing a dataframe (model_fits) with
the estimates using the original data, and a dataframe (bagged_results) with 
the bootstrap estimates with their percent of selection. The latter contains
the following columns: `percent_selected`: the relative proportion for selection of each subgroup,
`bagg`: the (uncorrected) bagged estimate  \eqn{\hat{\delta}^{*}_p}
`boot_red`: the bias reduced bootstrap estimate  \eqn{\check{\delta}^{*}_p}
`bagg_red`: the bias reduced bootstrap estimate with decreased variability by bagging \eqn{\tilde{\delta}^{*}_p}
and the respective standard deviations of the estimates.

}
\references{
Rosenkranz, G.(2016) "Exploratory subgroup analysis in clinical trials by model
selection", Biometrical Journal, 58, 1007-1259. doi: 10.1002/bimj.201500147
}

\seealso{
\code{\link{glm}}, \code{\link{coxph}}
}
\examples{
\dontrun{
## toy example calls using the simulated datnorm data-set without
## treatment and subgroup effect, see ?datnorm for details
data(datnorm)
head(datnorm)

## first need to create candidate subgroups (if not already defined in data-set)
## here generate candidate subgroups manually (need to be numeric 0-1 variables)
groups <- data.frame(labvalL.5=as.numeric(datnorm$labvalue < 0.5),
                     regUS=as.numeric(datnorm$region == "US"),
                     hgtL175=as.numeric(datnorm$height < 175))
fitdat <- cbind(datnorm, groups) # bind subgroup variables to main data
## subgroups of interest
subgr <- c("labvalL.5", "regUS", "hgtL175")
res <- bagged(resp = "y", trt = "treat", subgr = subgr, data = fitdat, 
                       covars = ~ x1 + x2, fitfunc = "lm")
res


## generate candidate subgroups using the subbuild function
## semi-automatically i.e. some groups specified directly (height and
## smoker), for region and labvalue subbuild generates subgroups (see
## ?subbuild)
cand.groups <- subbuild(datnorm, height < 175, smoker == 1, region, labvalue)
head(cand.groups)
fitdat <- cbind(datnorm, cand.groups) 
subgr <- colnames(cand.groups)
res <- bagged(resp = "y", trt = "treat", subgr = subgr, data = fitdat, 
                       covars = ~ x1 + x2, fitfunc = "lm")
res

## toy example call for binary data on simulated datbin data-set
data(datbin)
cand.groups <- subbuild(datbin, height < 175, smoker == 1, region, labvalue)
fitdat <- cbind(datbin, cand.groups) 
subgr <- colnames(cand.groups)
res <- bagged(resp = "y", trt = "treat", subgr = subgr, data = fitdat, 
                       covars = ~ x1 + x2, fitfunc = "glm", 
                       family = binomial(link = "logit"))
## scale of the treatment effect estimate: difference on log-odds scale
res

## toy example call for parametric and semi-parametric survival data on
## datsurv data-set
data(datsurv)
cand.groups <- subbuild(datsurv, height < 175, smoker == 1, region, labvalue)
fitdat <- cbind(datsurv, cand.groups)
subgr <- colnames(cand.groups)
res.survreg <- bagged(resp = "y", trt = "treat", subgr = subgr, data = fitdat,
                               covars = ~ x1 + x2,
                               fitfunc = "survreg", event = "event", dist = "exponential")
## parametric survival model (here exponential distribution)
## scale of treatment effect estimate: log scale (see ?survreg for details)
res.survreg

# Decreased B for a reduction in computational time
res.cox <- bagged(resp = "y", trt = "treat", subgr = subgr, data = fitdat,
                  fitfunc = "coxph", event = "event", B = 20) # B=2000 should be used
## scale of treatment effect estimate: difference in log-hazard rate
res.cox
## toy example call overdispersed count data on datcount data-set
data(datcount)
cand.groups <- subbuild(datcount, height < 175, smoker == 1, region, labvalue)
fitdat <- cbind(datcount, cand.groups)
subgr <- colnames(cand.groups)
# Decreased B for a reduction in computational time
res <- bagged(resp = "y", trt = "treat", subgr = subgr, data = fitdat,
              fitfunc = "glm.nb", exposure = "exposure", B = 20) # B=2000 should be used
## scale of treatment effect estimate: difference on log scale
res}
}
\keyword{models}

