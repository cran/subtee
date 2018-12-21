\name{plot.subtee}
\alias{plot.subtee}

\title{
Plotting subgroup treatment effect estimates
}
\description{
Plotting function for objects of class 'subtee'. Visualizes estimates and confidence intervals for all candidate subgroups (and complements) in a forest plot. 
}
\usage{
\method{plot}{subtee}(x, y = NULL, z = NULL, 
     type = c("trtEff", "trtEffDiff"),
     show.compl = FALSE,
     xlab = "", ylab = "default", main = "default",
     them,
     point.size = 2.5, line.size = 1,
     palette = "default", ...)
}
\arguments{
  \item{x}{
An object of class \code{subtee}, usually a result of a call to \code{\link{modav}}, \code{\link{unadj}} or \code{\link{bagged}}.
}
  \item{y}{
An object of class \code{subtee}, usually a result of a call to \code{\link{modav}}, \code{\link{unadj}} or \code{\link{bagged}}. If this is used, a comparison between the estimates will be provided. 
}
  \item{z}{
An object of class \code{subtee}, usually a result of a call to \code{\link{modav}}, \code{\link{unadj}} or \code{\link{bagged}}. If this is used, a comparison between the estimates will be provided.
}
  \item{type}{
A character specifyng if treatment effects should be plotted (\code{"trtEff"}) or the treatment-subgroup interactions. (\code{"trtEffDiff"})
}
  \item{show.compl}{
Logical. If true estimates for candidate subgroup complements should be plotted as well. Only available if \code{type = "trtEff"}.
}
  \item{xlab}{
Character. Label for x-axis.
}
  \item{ylab}{
    Character. Label for y-axis. 
}
  \item{main}{
    Character. Title. The default is to provide a string with the
    level of the uncertainty intervals.
}
  \item{them}{
ggplot2 theme. Use \code{\link[ggplot2]{theme}()} if you just need to tweak the display, or provide a complete ggplot2 theme (e.g \code{\link[ggplot2]{theme_bw}}()).
}
  \item{point.size}{
Size for points, which denote point estimates of treatment effects. Default to 2.5.
}
  \item{line.size}{
Size for points, which denote confidence interval of treatment effects. Default to 1.
}
  \item{palette}{
A string providing a ggplot2 colour palette to use. This will be passed to the palette option in a \code{\link[ggplot2]{scale_colour_brewer}} sentence.
}
  \item{...}{
Not used.
}
}
\value{
Forest plot visualizing treatment effect estimates (if \code{type = "trtEff"}) or treatment-subgroup interactions in candidate subgroups (if \code{type = "trtEffDiff"}).
}

\seealso{
\code{\link{summary.subtee}}
}
\examples{
data(datnorm)
cand.groups <- subbuild(datnorm, height, labvalue, region, smoker)
fitd <- cbind(datnorm, cand.groups)
subgr <- colnames(cand.groups)
### Plot unadjusted estimates
res_unadj <- unadj(resp = "y", trt = "treat", subgr = subgr, data = fitd, 
             covars = ~ x1 + x2, fitfunc = "lm")
summary(res_unadj)
plot(res_unadj)
plot(res_unadj, show.compl = TRUE)
plot(res_unadj, type = "trtEffDiff")

### Compare unadjusted with model averaging estimates
res_modav <- modav(resp = "y", trt = "treat", subgr = subgr, data = fitd, 
             covars = ~ x1 + x2, fitfunc = "lm")
plot(res_unadj, res_modav, show.compl = TRUE)
plot(res_unadj, res_modav, type = "trtEffDiff")
}

