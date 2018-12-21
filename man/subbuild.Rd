\name{subbuild}
\alias{subbuild}
\title{
Generating candidate subgroups based on an input data set
}
\description{
Takes categorical or continuous baseline covariate vectors and builds a
matrix of binary candidate subgroups indicators.
}
\usage{
subbuild(data, ..., n.cuts = 2, dig.lab = 3, dupl.rm = FALSE,
         make.valid.names = FALSE)
}
\arguments{
  \item{data}{
Data frame. Contains baseline covariates from which candidate subgroups are generated. Categorical covariates should be of type factor.
}
\item{\dots}{
  Optional subgroup definitions 
  or variable names, which are used to generate candidate subgroups.  
}
  \item{n.cuts}{
Integer. Number of cutoffs for each covariate.
}
  \item{dig.lab}{
Integer. Digit to which subgroup cutoffs are rounded.
}
  \item{dupl.rm}{
Logical. Remove duplicate subgroups. Note that this applies also to two
subgroup vectors a and b that satisfy a=1-b (i.e. only labels 0 and 1
exchanged), because these will give the same model fit (only the label
of "subgroup" and "complement" are exchanged).
}
  \item{make.valid.names}{
Logical. If \code{TRUE} subgroup names in the final output are transformed to be
syntactically valid names (see \code{?make.names})
}
}
\details{
The \dots argument allows manual specification of subgroups that should be included.
Subgroup definitions should be passed as (typically logical) expressions, that can be evaluated on \code{data} and result in binary subgroup indicator variables.

If only a variable name is specified in \dots, subgroup definitions based on this 
covariate are automatically generated in the following way:
For covariates of type factor or character candidate subgroups are the patients in each category.
For covariates of type numeric or integer, cutpoints for candidate subgroups are generated based on covariate quantiles. For each continuous covariates  '\code{n.cuts}' + 1 non-overlapping subgroups of (roughly) the same size are generated.

If no information about subgroups is supplied in \dots, candidate subgroups are 
automatically generated for all variables in \code{data}. 
Subgroup names are taken from the column names of the '\code{data}' data frame (or set to x1, x2,... if no names are supplied)

If \code{dupl.rm} is \code{TRUE} any duplicate columns (either subgroup or
complement is equal to another column) are removed from the final output.
}
\value{
A data frame of candidate subgroups.
}

\seealso{
\code{\link{bagged}}, \code{\link{modav}}, \code{\link{unadj}}
}
\examples{
data(datnorm)
## data frame of covariates considered for subgroup analysis
cov.dat <- datnorm[,c("height", "labvalue", "region", "smoker")]
## by default generate all subgroups for each categorical variable and
## use cut-offs based on quantiles for numeric variables
cand.groups <- subbuild(cov.dat)
head(cand.groups)
## alternatively use
cand.groups <- subbuild(datnorm, height, labvalue, region, smoker)
head(cand.groups)

## use more cutpoints
cand.groups2 <- subbuild(cov.dat, n.cuts = 4)
ncol(cand.groups)
ncol(cand.groups2)

## remove duplicate columns for smoker
cand.groups3 <- subbuild(cov.dat, dupl.rm = TRUE)
head(cand.groups3)
ncol(cand.groups3)

## syntactically valid names
cand.groups4 <- subbuild(cov.dat, make.valid.names = TRUE)
head(cand.groups4)

## manually specify subgroup definitions and which covariates to consider
cand.groups5 <- subbuild(cov.dat, region == "EU", height > 172, labvalue)
## note that for labvalue cut-offs are generated automatically based on quantiles
head(cand.groups5)

## further examples for manual specification of subgroups
cand.groups6 <- subbuild(cov.dat, region \%in\% c("Japan","EU"), smoker != 0)
## note that for labvalue cut-offs are generated automatically based on quantiles
head(cand.groups6)

## missing values in data-set are propagated through
cov.dat$labvalue[sample(1:nrow(cov.dat),10)] <- NA
cov.dat$region[sample(1:nrow(cov.dat),20)] <- NA
cov.dat$smoker[sample(1:nrow(cov.dat),10)] <- NA
cand.groups7 <- subbuild(cov.dat)
head(cand.groups7)

## if covariates in the data frame contain missing values consider
## imputing them for example with the rfImpute function from the
## randomForest package
}

