library(testthat)
library(subtee)

########################################################################
## checks on the treatment variable
## "trt" should be a numeric (or integer) 0,1 variable
## all checks below should fail!
set.seed(1)
exdat <- data.frame(y = rnorm(100), trt = rbinom(100, size = 1, prob = 0.5),
                     base = rnorm(100),
                     region = factor(sample(c("US", "EU", "Japan"), 100, replace=TRUE)),
                     respo = factor(rbinom(100, 1, 0.3)))
covdat <- exdat[c("base", "region", "respo")]
csubgr <- subbuild(covdat)

dat2 <- data.frame(y=exdat$y)
dat2 <- cbind(dat2, csubgr)
dat2a <- dat2b <- dat2c <- dat2
subgr <- colnames(csubgr)
treatment <- exdat$trt
dat2a$trt <- factor(treatment)
trtmnta <- dat2a$trt
dat2b$trt <- NULL
dat2c$trt <- (treatment == 1) + 1
trtmntb <- (treatment == 1) + 10

test_that("Error thrown when treatment variable not numeric or integer (0, 1) variable", {
  expect_error(modav(resp = "y", trt = "trt", subgr = subgr, data = dat2a))
  expect_error(modav(resp = "y", trt = "trtmnta", subgr = subgr, data = dat2b))
  expect_error(modav(resp = "y", trt = "trt", subgr = subgr, data = dat2c))
  expect_error(modav(resp = "y", trt = "trtmntb", subgr = subgr, data = dat2b))
})


########################################################################
## checks on correct specification of the subgroup covariable
## catch cases where the subgroup contains out only 0s or only 1s, and
## cases where there are not at least 1 patient on treatment and
## control in subgroup and complement
## all checks below should fail!

dat2a <- dat2b <- dat2c <- dat2d <- dat2e <- dat2

dat2a$NULLGRP <- rep(0, nrow(dat2a))
dat2b$NULLGRP <- rep(1, nrow(dat2b))
dat2c$TRT2 <- treatment ## subgroups complete confounded with treatment

## catch cases where the subgroup variable is not a numeric variable
dat2d$NGRP <- factor(rbinom(nrow(dat2d), 1, 0.5))
## catch cases where the subgroup variable is not a 0,1 variable
dat2e$NGRP <- rnorm(nrow(dat2e))


test_that("Error thrown when there are not at least 1 patient on treatment and control in subgroup and complement", {
  expect_error(modav(resp = "y", trt = "treatment", subgr = c(subgr, "NULLGRP"), data = dat2a))
  expect_error(modav(resp = "y", trt = "treatment", subgr = c(subgr, "NULLGRP"), data = dat2b))
  expect_error(modav(resp = "y", trt = "treatment", subgr = c(subgr, "TRT2"), data = dat2c))
  expect_error(modav(resp = "y", trt = "treatment", subgr = c(subgr, "NGRP"), data = dat2d))
  expect_error(modav(resp = "y", trt = "treatment", subgr = c(subgr, "NGRP"), data = dat2e))
})


########################################################################
## using a subgroup variable in covars and subgr
set.seed(1)
exdat <- data.frame(y = rnorm(100),
                    trt = rbinom(100, size = 1, prob = 0.5),
                    marker = rbinom(100, 1, 0.25),
                    region = factor(sample(c("US", "EU", "Japan"), 100,
                                           replace=TRUE)), stringsAsFactors = FALSE)
covdat <- exdat[c("region")]
csubgr <- subbuild(covdat)
dat2 <- cbind(exdat, csubgr)
subgr <- c(colnames(csubgr), "marker")
res <- modav(resp = "y", trt = "trt", subgr = subgr,
             covars=~marker, data = dat2)


test_that("Subgroup variable can be used as covariate", {
  expect_equal(class(res), "subtee")
  expect_equal("marker" %in% summary(res)$Group, TRUE)
})

