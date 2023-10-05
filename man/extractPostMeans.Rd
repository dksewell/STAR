\name{extractPostMeans}
\alias{extractPostMeans}
\title{
Extract posterior means
}
\description{
Takes a \code{STAR} object fit with a Gibbs sampler and computes the posterior means.
}
\usage{
extractPostMeans(STARObj, burnin = 0, thinning = 1)
}

\arguments{
\item{STARObj}{Object of class \code{STAR}.  \code{STARObj$method} must be 'Gibbs'.}
\item{burnin}{Burn in period at the beginning of the chain to remove before computing the posterior means.}
\item{thinning}{Thinning interval to reduce autocorrelation within the samples.}
}

\value{
An object of class \code{STAR}.  \code{PostMean} elements are added to the object:
\describe{
\item{BetaMean}
\item{ThetaMean}
\item{OmegaMean (if \code{RE} is \code{TRUE})}
\item{tauS (if \code{RE} is \code{TRUE})}
\item{tauR (if \code{RE} is \code{TRUE})}
\item{s2R (if \code{RE} is \code{TRUE})}
}
}

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
