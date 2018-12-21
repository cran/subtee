\name{unadj}
\alias{unadj}
\title{
Treatment effect estimation based on marginal subgroup models.
}
\description{
  Unadjusted estimation of treatment effects in subgroups.
  Fits separate (marginal) models for each candidate subgroup, i.e. including
  the subgroup as a main effect and interaction with treatment for each model.
}
\usage{
unadj(resp, trt, subgr, covars = NULL, data, 
      fitfunc = c("lm", "glm", "glm.nb", "survreg", "coxph", "rlm"),
      event, exposure, level = 0.1, ...)
}
\arguments{
  \item{resp}{
Character giving the name of the response variable. The variable can be
either defined in the global environment or in the data-set \code{data}
specified below. For interactive use it is also possible to use unquoted
names (i.e. \code{unadj(resp,...)} instead of \code{unadj("resp",...)}),
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
  Formula, specifying additional (prognostic) covariates to be included
  in the models (need to be available in \code{data}). It is crucial for
  the model averaging approach to include the important prognostic
  covariates (in particular if the corresponding prognostic covariate
  also defines a subgroup; otherwise models/subgroup might get
  upweighted just because the variable has prognostic value, but not
  because the treatment effect is modified).
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
  \item{\dots}{
other arguments passed to the model fitting function.
}
}
\details{
In the simple linear case (e.g when using fitfunc \code{\link{lm}}) for each of the \eqn{P} candidate subgroups the fitted model is of the form

\deqn{M_p : y_i \sim N(\mu_i^{(p)}, \sigma_p^2), i=1,...,n}{M_p : y_i ~  N(\mu_i^(p), \sigma_p^2), i=1,...,n}

where \deqn{\mu_i^{(p)} = 
\alpha_p + \beta_p z_i +
(\gamma_p + \delta_p z_i) s_{pi} +  \sum_{k = 1}^{K} \tau_k x_{ik}
}

where \eqn{s_i} denotes the subgroup indicators (the column vectors of \code{subgr}), \eqn{z_i} is the treatment indicator (from \code{trt}) and \eqn{x{.1}, ..., x{.K}} are additional covariates as specified in \code{covars}.
For other fitting functions the models are of similar form, including prognostic and predictive effects of subgroups.

A treatment effect (on the scale determined by \code{fitfunc}) for the candidate subgroups is estimated as \eqn{\hat{\beta} +  \hat{\delta_p}} and a treatment effect estimate for the complement is given by \eqn{\hat{\beta}}. 
Note that choosing subgroups based on these unadjusted treatment effect estimates may lead to overoptimistic conclusions in regards to the treatment effect in that subgroup. Naive estimates do not consider model
selection uncertainty and will often suffer from selection bias. 
}
\value{
  A list (object of class \code{subtee}). The most important entries are (i) \code{fitmods} containing all
  fitted subgroup models and the overall model (ii) \code{trtEff}
  containing the treatment effect estimates and CI for subgroup and subgroup
  complements. (iii) \code{trtEffDiff} containing the differences in
  treatment effect estimates (subgroup vs complement) and CI. 
}
\references{
Thomas, M., and Bornkamp, B. (2017) "Comparing Approaches to Treatment
Effect Estimation for Subgroups in Early Phase Clinical Trials."
Statistics in Biopharmaceutical Research, 9, 160-171,
doi: 10.1080/19466315.2016.1251490

Bornkamp, B., Ohlssen, D., Magnusson, B. P., and Schmidli, H. (2017)
"Model averaging for treatment effect estimation in subgroups."
Pharmaceutical Statistics,  16, 133-142,
doi: 10.1002/pst.1796

Raftery, A. E. (1995) "Bayesian model selection in social research."
Sociological Methodology, 25, 111-163.
}

\seealso{
\code{\link{summary.subtee}}, \code{\link{plot.subtee}},
\code{\link{lm}}, \code{\link{glm}}, \code{\link{glm.nb}},
\code{\link{survreg}}, \code{\link{coxph}}
}
\examples{
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
res <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitdat, 
             covars = ~ x1 + x2, fitfunc = "lm")
summary(res)
plot(res)
                     
## generate candidate subgroups using the subbuild function
## semi-automatically i.e. some groups specified directly (height and
## smoker), for region and labvalue subbuild generates subgroups (see
## ?subbuild).
cand.groups <- subbuild(datnorm, height < 175, smoker == 1, region, labvalue)
head(cand.groups)
fitdat <- cbind(datnorm, cand.groups) 
subgr <- colnames(cand.groups)
res <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitdat, 
             covars = ~ x1 + x2, fitfunc = "lm")
summary(res)
plot(res)

## toy example call for binary data on simulated datbin data-set
data(datbin)
cand.groups <- subbuild(datbin, height < 175, smoker == 1, region, labvalue)
fitdat <- cbind(datbin, cand.groups) 
subgr <- colnames(cand.groups)
res <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitdat, 
             covars = ~ x1 + x2, fitfunc = "glm", 
             family = binomial(link = "logit"))
## scale of the treatment effect estimate: difference on log-odds scale
summary(res)
plot(res)

## toy example call for parametric and semi-parametric survival data on
## datsurv data-set
data(datsurv)
cand.groups <- subbuild(datsurv, height < 175, smoker == 1, region, labvalue)
fitdat <- cbind(datsurv, cand.groups)
subgr <- colnames(cand.groups)
res.survreg <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitdat,
                     covars = ~ x1 + x2,
                     fitfunc = "survreg", event = "event", dist = "exponential")
## parametric survival model (here exponential distribution)
## scale of treatment effect estimate: log scale (see ?survreg for details)
summary(res.survreg)
plot(res.survreg)
res.cox <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitdat,
                 covars = ~ x1 + x2, fitfunc = "coxph", event = "event")
## scale of treatment effect estimate: difference in log-hazard rate
summary(res.cox)
plot(res.cox)

## toy example call overdispersed count data on datcount data-set
data(datcount)
cand.groups <- subbuild(datcount, height < 175, smoker == 1, region, labvalue)
fitdat <- cbind(datcount, cand.groups)
subgr <- colnames(cand.groups)
res <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitdat,
             covars = ~ x1 + x2, fitfunc = "glm.nb", exposure = "exposure")
## scale of treatment effect estimate: difference on log scale
summary(res)
plot(res)
}
\keyword{ models}

