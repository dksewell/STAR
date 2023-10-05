\name{STAR-package}
\alias{STAR-package}
\alias{STAR}
\docType{package}
\title{
Simultaneous and Temporal Autoregressive Network Models
}
\description{
Link prediction for temporal networks.  This package implements the simultaneous and temporal dependence autoregressive models for dynamic networks.
}
\details{
\tabular{ll}{
Package: \tab STAR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-10-18\cr
License: \tab GPL (>= 2)\cr
}
There are two primary functions, one to fit the data and the other to plot it.  \code{STARnet()} can fit the model either accounting for or ignoring simultaneous dependence.  Estimation can be done by either variational Bayes or by a Gibbs sampler. \code{plot.STAR()} can provide trace plots, comparisons of posterior distributions, and plots of evidence of simultaneous dependence.  The \code{extractPostMeans} function will extract the posterior means from a \code{STAR} object.  There also exists a \code{summary} function for the \code{STAR} class.
}
\author{
Daniel K. Sewell, Amy Hahn

Maintainer: Daniel K. Sewell <daniel-sewell@uiowa.edu>
}
\references{
Sewell, D.K. (2017). 'Simultaneous and Temporal Autoregressive Network Models.'  Under Review.
}
\keyword{ network analysis; covariance structures; dynamic networks; homophily. }
\examples{
\dontrun{
library(STAR)
?STAR
?STARnet
data(knecht)
attach(knecht)
varNames = c("Intercept","gender","age","academEval","primary","alcohol","delinquency")
mod1 = STARnet(A,X,Xvarnames = varNames)
summary(mod1)
mod2 = STARnet(A,X,Xvarnames = varNames,RandEff=FALSE)
summary(mod2)
mod3 = STARnet(A,X,method="Gibbs",Xvarnames = varNames)
mod3 = extractPostMeans(mod3,burnin=100)
summary(mod3)
mod4 = STARnet(A,X,method="Gibbs",Xvarnames = varNames,RandEff=FALSE)
mod4 = extractPostMeans(mod4,burnin=100)
summary(mod4)

plot(mod1,plotType='epsBall',propVar = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3))
plot(mod1,mod2,plotType='Comparison')
plot(mod4,plotType='trace')
plot(mod3,mod4,plotType='Comparison')
}
}
