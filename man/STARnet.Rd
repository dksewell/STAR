\name{STARnet}
\alias{STARnet}
\title{
Simultaneous and temporal autoregressive model for dynamic networks
}
\description{
\code{STARnet()} can fit the model either accounting for or ignoring simultaneous dependence.  Estimation can be done by either variational Bayes or by a Gibbs sampler. }
\usage{
STARnet(A,X,lagFunctions = list(functions = 1:8, decay = NULL),
    RandEff = TRUE, method = c("VB", "Gibbs"),
    maxIter = 1000, eps = 1e-05, nSims = 1000, hyperparms = NULL,   
    Xvarnames = NULL)
}

\arguments{
\item{A}{Dynamic network data.  This should be in the form of a list of \eqn{n\times n} 'dgCMatrix'.  The total number of elements in this list should be the number of time points \eqn{T} + number of lags \eqn{L}.}
\item{X}{Covariate information.  This should be a list of \eqn{T} lists.  Each of these \eqn{T} lists should be a list of \eqn{p} \eqn{n\times n} 'matrix' objects with 0's along the diagonal, where \eqn{p} is the number of dyadic covariates, including the intercept (i.e., typically \code{X[[timePoint]][[1]]} will be a matrix of 1's with 0's along the diagonal).
}
\item{lagFunctions}{A list of lag functions you want to include.  You can use either numbers or abbreviations.  \code{decay} corresponds to the geometrically weighted degree statistic.
\describe{
\item{1 (outDeg)}{\eqn{A {\bf 1}{\bf 1}'}}
\item{2 (inDeg)}{\eqn{{\bf 1}{\bf 1}'A}}
\item{3 (stability)}{\eqn{A}}
\item{4 (reciprocity)}{\eqn{A'}}
\item{5 (trans1)}{\eqn{AA}}
\item{6 (trans2)}{\eqn{AA'}}
\item{7 (trans3)}{\eqn{A'A}}
\item{8 (3cycle)}{\eqn{A'A'}}
\item{9 (gwdeg)}{Snijders et al. 2006 eq. (17), using symmetric version of \eqn{A}.}
\item{10 (isolatesSendSym)}{Indicator of whether the sender is an isolate w.r.t. sending or receiving links.}
\item{11 (isolatesRecSym)}{Indicator of whether the receiver is an isolate w.r.t. sending or receiving  links.}
\item{12 (isolatesSend)}{Indicator of whether the sender is an isolate w.r.t. sending links.}
\item{13 (isolatesRec)}{Indicator of whether the receiver is an isolate w.r.t. receiving links.}
}
}
\item{RandEff}{Logical.  Should simultaneous dependence be accounted for?}
\item{method}{Variational Bayes or Gibbs sampler?}
\item{maxIter}{Maximum number of iterations for the VB algorithm.  Ignored if \code{method = 'Gibbs'.}}
\item{eps}{Maximum relative change in mean coefficients before convergence is reached in the VB algorithm.  Ignored if \code{method = 'Gibbs'}.}
\item{nSims}{Number of MCMC iterations.  Ignored if \code{method = 'VB'}.}
\item{hyperparms}{A list of hyperparameters.  If a hyperparameter is \code{NULL}, it will be set automatically.
\describe{
\item{mu_s}{Scalar. The prior mean of \eqn{\tau_s}}
\item{mu_r}{Scalar. The prior mean of \eqn{\tau_r}}
\item{mu_sr}{\eqn{2\times 2} matrix.  The prior mean of \eqn{\Omega}.}
\item{delta}{Scalar.  This value affects the variance of the above parameters.  Namely, \eqn{a_{s0}=}2+delta, and \eqn{b_{s0}=}mu_s*(\eqn{a_{s0}}-1), and similarly for \eqn{a_{r0}} and \eqn{b_{r0}}.  \eqn{a_{\Omega 0}=}5+delta, and \eqn{b_{\Omega 0}=}mu_sr*(\eqn{a_{\Omega 0}}-3).}
\item{aRec0}{Scalar.  Shape parameter for the prior of \eqn{\sigma^2_R}.}
\item{bRec0}{Scalar.  Scale parameter for the prior of \eqn{\sigma^2_R}.}
\item{sigB}{Scalar.  Variance for the prior of \eqn{\beta}.}
\item{sigTh}{Scalar.  Variance for the prior of \eqn{\theta}.}
\item{s2R}{Scalar.  Prior variance for the elements of \eqn{R_t}.  Only used if \code{method = 'VB'}.}
}
}
\item{Xvarnames}{Variable names corresponding to each matrix in the lists \code{X[[timePoint]]}.}

}

\details{
This function estimates the homophilic effects using a generalized linear mixed model that captures both temporal dependencies and simultaneous dependencies within network data.  The model is
\deqn{logit(P(A_{ijt}=1)) = X_{ijt}\beta + f(A_1,A_2,\ldots,A_{t-1})\theta + \epsilon_{ijt}.}
The temporal dependencies are captured by the autoregressive terms involving functions of the lagged network, and the simultaneous dependencies are captured in the \eqn{n(n-1)\times n(n-1)} covariance matrix of the \eqn{\epsilon_{ijt}}'s.
}

\value{
An object of class \code{STAR}.  If \code{method='VB'},
\describe{
\item{BetaMean}{Posterior mean of \eqn{\beta}.}
\item{ThetaMean}{Posterior mean of \eqn{\theta}.}
\item{BetaThetaCov}{Posterior covariance matrix of \eqn{(\beta',\theta')'}.}
\item{sShape,sScale}{Posterior shape and scale parameters of \eqn{\tau_s} (if \code{RE} is \code{TRUE}).}
\item{rShape,rScale}{Posterior shape and scale parameters of \eqn{\tau_r} (if \code{RE} is \code{TRUE}).}
\item{OmegaDOF,OmegaScaleMat}{Posterior degrees of freedom and scale matrix for \eqn{\Omega} (if \code{RE} is \code{TRUE}).}
\item{s1Means,s2Means,r1Means,r2Means,srCov}{Posterior means and covariance matrix for \eqn{s_{1t},s_{2t},r_{1t},r_{2t}} (if \code{RE} is \code{TRUE}).}
\item{RecMean}{Posterior means of \eqn{R_t}.}
\item{numIter}{Number of iterations performed before convergence or \code{maxIter} was reached.}
\item{convergence}{A vector of length \code{numIter}, where the \eqn{k}-th value corresponds to the maximum relative change in any of the mean parameters.}
}
If \code{method='Gibbs'},
\describe{
\item{Beta}{A \eqn{p\times}\code{nSims} matrix of posterior samples of \eqn{\beta}. }
\item{Theta}{A \eqn{q\times}\code{nSims} matrix of posterior samples of \eqn{\theta}. }
\item{Omega}{A \eqn{2\times 2\times}\code{nSims} array of posterior samples of \eqn{\Omega} (if \code{RE} is \code{TRUE}). }
\item{tauS}{A \code{nSims}-dim vector of posterior samples of \eqn{\tau_s} (if \code{RE} is \code{TRUE}). }
\item{tauR}{A \code{nSims}-dim vector of posterior samples of \eqn{\tau_r} (if \code{RE} is \code{TRUE}). }
\item{s1,s2,r1,r2}{A \eqn{n\times T\times}\code{nSims} array of posterior samples of \eqn{s_{kt}} or \eqn{r_{kt}} (if \code{RE} is \code{TRUE}). }
\item{s2R}{A \code{nSims}-dim vector of posterior samples of \eqn{\sigma^2_R} (if \code{RE} is \code{TRUE}). }
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
