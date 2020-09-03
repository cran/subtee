## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- show='hold'-------------------------------------------------------------
library(ggplot2)
library(subtee)
################################################################################
# We use the dataset from Rosenkranz (2016) https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201500147
# to illustrate the methods proposed in this work.
# The data comes from a clinical trial of an prostate cancer 
# treatment
# Data is loaded from  Royston, Patrick, and Willi Sauerbrei. 
# Multivariable model-building: a pragmatic approach to 
# regression anaylsis based on fractional polynomials for 
# modelling continuous variables. Vol. 777. John Wiley & Sons, 2008. 
# https://www.imbi.uni-freiburg.de/Royston-Sauerbrei-book
prca = get_prca_data()

## first create candidate subgroups
cand.groups <- subtee::subbuild(prca, dupl.rm = TRUE,
                               BM == 1, PF == 1, HX == 1,
                               STAGE == 4, AGE > 65, WT > 100)
fitdat <- cbind(prca, cand.groups)
subgr.names = names(cand.groups)
prog = as.formula(paste(" ~ ", paste0("`", names(cand.groups),"`", collapse = " + ")))

### Unadjusted estimates
res_unadj = unadj(resp = "SURVTIME", trt = "RX", subgr = subgr.names,
                  data = fitdat, covars = prog, 
                  event = "CENS", fitfunc = "coxph")
### ModelAveraging estimates
res_modav = modav(resp = "SURVTIME", trt = "RX", subgr = subgr.names,
                  data = fitdat, covars = prog, 
                  event = "CENS", fitfunc = "coxph")
### Bagged estimates
set.seed(321231) # set seed for reproducible results in the bootstrap samples
res_bagged = bagged(resp = "SURVTIME", trt = "RX", subgr = subgr.names,
                    data = fitdat, covars = prog, 
                    event = "CENS", fitfunc = "coxph",
                    select.by = "BIC", B = 200) #B = 2000)

## ---- fig.width = 5, fig.show='hold'------------------------------------------
ggplot(aes(y = Subset, x = trtEff, xmin = LB, xmax = UB, colour = Subset), 
       data = res_unadj$trtEff) + 
  geom_point(size = 2) +
  geom_errorbarh(size = 1, show.legend = FALSE, height = 0) +
  facet_grid(Group ~ .)

## ---- fig.width = 5, fig.show='hold'------------------------------------------
plot(res_unadj)

## ---- fig.width = 6, fig.height = 5, fig.show='hold'--------------------------
plot(res_unadj, show.compl = TRUE)

## ---- fig.width = 6, fig.show='hold'------------------------------------------
plot(res_bagged, show.compl = TRUE)

## ---- fig.width = 6, fig.height = 5, fig.show='hold'--------------------------
plot(res_unadj, res_modav, palette = "Dark2")

## ---- fig.width = 6, fig.height=5, fig.show='hold'----------------------------
plot(res_unadj, res_modav, show.compl = TRUE)

## ---- fig.width = 6, fig.show='hold'------------------------------------------
plot(res_unadj, res_modav, res_bagged, show.compl = TRUE)

## ---- fig.width = 5, fig.show='hold'------------------------------------------
plot(res_unadj, type = "trtEffDiff")

## ---- fig.width = 6, fig.show='hold'------------------------------------------
plot(res_unadj, res_modav, type = "trtEffDiff")

## ---- fig.width = 6, fig.show='hold'------------------------------------------
plot(res_unadj, res_modav, res_bagged, type = "trtEffDiff")

