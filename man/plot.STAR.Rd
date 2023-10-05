\name{plot.STAR}
\alias{plot.STAR}
\title{
Plot various facets of STAR objects.
}
\description{
\code{plot.STAR()} can provide trace plots, comparisons of posterior distributions, and plots of evidence of simultaneous dependence.
}
\usage{
plot(STARObj1, STARObj2 = NULL, burnin = 0, 
     plotType = c("trace","epsBall", "Comparison"),
     Titles = TRUE, savePlots = FALSE, prefix = "", 
     N = 5000, propVar = 0.5)	
}

\arguments{
\item{STARObj1}{An object of class \code{STAR}.}
\item{STARObj2}{An object of class \code{STAR}.  Optional, but if included \code{STARObj2$RE} should be \code{FALSE}.  Curves given in dotted line.}
\item{burnin}{How many samples to ignore at the beginning of the MCMC chain.  Ignored if \code{STARObj1$method} is 'VB'.}
\item{plotType}{Vector of 1,2, or 3 characters.  If 'trace' is included and \code{STARObj1$method} is 'Gibbs', trace plots are provided for each of the mean coefficients, including the burn-in period.  

If 'epsBall' is included, a plot is given providing the posterior probability that the random effects \eqn{(s_{1t},r_{1t},s_{2t},r_{2t})} fall within a ball of radius \eqn{\epsilon}, providing evidence for the existence of simultaneous dependence in the network data.  See Sewell '17 for more details.  

If 'Comparison' is included and \code{!is.null(STARObj2)}, then the posterior distributions for each of the mean coefficients is given for both accounting for and ignoring simultaneous dependence.  Solid lines gives the posterior distribution for \code{STARObj1}, and dotted lines for \code{STARObj2}.}
\item{Titles}{Logical.  Should plot titles be included.  }
\item{savePlots}{Logical.  Should plots be displayed immediately (\code{FALSE}) or should they be saved to the working directory (\code{TRUE}).}
\item{prefix}{Optional character vector if \code{savePlots} is \code{TRUE}.  Adds \code{prefix} to the beginning of the names of all plots being saved.}
\item{N}{If \code{STARObj1$method} is 'VB', the number of random draws to compute the posterior probabilities when \code{plotType} is \code{'epsBall'}.}
\item{propVar}{The proportion of variance accounted for by the simultaneous dependence (less reciprocal effects) for the reference distribution given when \code{plotType} is \code{'epsBall'}.}

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
