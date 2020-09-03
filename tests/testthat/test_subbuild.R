library(testthat)
library(subtee)
########################################################################
## Example 1

datFr <- structure(list(PatID = c(1L, 1L, 1L, 2L, 3L, 4L, 4L, 5L, 5L, 
                          6L, 6L, 7L, 7L, 8L, 9L, 10L, 10L, 11L, 12L),
                        time = c(5L, 12L, 
                          25L, 8L, 6L, 4L, 19L, 5L, 16L, 5L, 14L, 25L, 35L, 31L, 11L, 19L, 
                          33L, 8L, 26L),
                        event = c(0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 
                          0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L),
                        start = c(0L, 8L, 17L, 
                          0L, 0L, 0L, 7L, 0L, 8L, 0L, 9L, 0L, 30L, 0L, 0L, 0L, 22L, 0L, 
                          0L),
                        stop = c(5L, 12L, 25L, 8L, 6L, 4L, 19L, 5L, 16L, 5L, 14L, 
                          25L, 35L, 31L, 11L, 19L, 33L, 8L, 26L),
                        treat = c(1L, 1L, 1L, 
                          1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L
                          )), .Names = c("PatID", "time", "event", "start", "stop", "treat"
                                ), class = "data.frame", row.names = c(NA, -19L))
subs <- subbuild(datFr, time, start)


test_that("Result is a dataframe of subgroups", {
  expect_equal(is.data.frame(subs), TRUE)
  expect_named(subs)
  expect_equal(all(apply(subs, 2, is.numeric)), TRUE)
  expect_equal(nrow(datFr), nrow(subs))
 }
)

########################################################################
## Example 2

set.seed(1)
## generate treatment and control response values
mu0 <- 1
mu1 <- 0.5
kappa <- 2
N<-100
## control and treatment
lambda0 <- rgamma(N,shape=1/kappa,rate=1/(kappa*mu0))
Y0 <- rep(NA,N)
for (i in 1:N) {
  Y0[i] <- rpois(1,lambda0[i])
}
lambda1 <- rgamma(N,shape=1/kappa,rate=1/(kappa*mu1))
Y1 <- rep(NA,N)
for (i in 1:N) {
  Y1[i] <- rpois(1,lambda1[i])
}
Y <- c(Y0,Y1)
TRT <- c(rep(0,N),rep(1,N))
## generate region indicator
reg <- rep("RoW",2*N)
regu <- runif(2*N)
reg[regu<0.4] <- "US"
reg[regu>0.7] <- "EU"
REGION <- as.factor(reg)
## number of relapses in past 2 years
mu <- 2*mu0
lambda <- rgamma(2*N,shape=1/kappa,rate=1/(kappa*mu))
NREL2 <- rep(NA,2*N)
for (i in 1:(2*N)) {
  NREL2[i] <- rpois(1,lambda[i])
}
## age covariate
AGE <- floor(rnorm(2*N,40,10))
## exposure
exposure <- rexp(2*N, 1)
## put things together
dat <- data.frame(Y=Y,TRT=TRT,REGION=REGION,NREL2=NREL2,
                  AGE=AGE,EXPOS=exposure)

### tests for different subgroup 
cand.groups1 <- subbuild(dat, REGION, NREL2)
cand.groups2 <- subbuild(dat, REGION)
cand.groups3 <- subbuild(dat, AGE, n.cuts = 2)

test_that("Results are dataframes of subgroups", {
  expect_equal(is.data.frame(cand.groups1), TRUE)
  expect_named(cand.groups1)
  expect_equal(all(apply(cand.groups1, 2, is.numeric)), TRUE)
  expect_equal(is.data.frame(cand.groups2), TRUE)
  expect_named(cand.groups2)
  expect_equal(all(apply(cand.groups2, 2, is.numeric)), TRUE)
  expect_equal(is.data.frame(cand.groups3), TRUE)
  expect_named(cand.groups3)
  expect_equal(all(apply(cand.groups3, 2, is.numeric)), TRUE)
}
)

test_that("Number of subgroups is correct", {
  expect_equal(ncol(cand.groups1) > ncol(cand.groups2), TRUE)
  expect_equal(ncol(cand.groups2) == ncol(cand.groups3), TRUE)
  expect_equal(ncol(cand.groups3) == 3, TRUE)
 }
)

## now all together
cand.groups <- subbuild(dat, REGION, NREL2, AGE)
subgr <- colnames(cand.groups)
dat <- cbind(dat, cand.groups)

est <- unadj("Y", "TRT", subgr=subgr, covars=~NREL2, data=dat,
             fitfunc="glm.nb", exposure="EXPOS")
estMA <- modav("Y", "TRT", subgr=subgr, covars=~NREL2, data=dat,
               fitfunc="glm.nb", exposure="EXPOS")
estMA2 <- modav("Y", "TRT", subgr=subgr, covars=~NREL2, data=dat,
                prior=1, nullprior=9, fitfunc="glm.nb",
                exposure="EXPOS")

test_that("Estimation functions with subgroup variables produce subtee outputs", {
  expect_equal(class(est) == "subtee", TRUE)
  expect_equal(class(estMA) == "subtee", TRUE)
  expect_equal(class(estMA2) == "subtee", TRUE)
 }
)
