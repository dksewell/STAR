\name{summary.STAR}
\alias{summary.STAR,plot.STAR}
\title{
Summarizing model fit from \code{STARnet}
}
\description{
\code{summary} method for class 'STAR'.
}
\usage{
summary(STARObj, alpha = 0.05, burnin = 0, thinning = 1,...)

print(sObj, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
\item{STARObj}{Object of class \code{STAR} (from \code{STARnet()}).}
\item{sObj}{Object of class \code{summary.STAR}.}
\item{alpha}{Probability outside of credible intervals.}
\item{burnin}{(Passed to \code{extractPostMeans} if not already run). Burn in period at the beginning of the chain to remove before computing the posterior means.}
\item{thinning}{(Passed to \code{extractPostMeans} if not already run). Thinning interval to reduce autocorrelation within the samples.}
\item{digits}{Number of digits to print.}
}

\value{
An object of class \code{summary.STAR}.
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
