% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ods-MSELE.R
\name{odsmle}
\alias{odsmle}
\title{MSELE estimator for analyzing the primary outcome in ODS design}
\usage{
odsmle(Y, X, beta, sig, pis, a, rs.size, size, strat)
}
\arguments{
\item{Y}{vector for the primary response}

\item{X}{the design matrix with a column of 1's for the intercept}

\item{beta}{starting parameter values for the regression coefficients that
relate Y to X.}

\item{sig}{starting parameter values for the error variance of the regression.}

\item{pis}{starting parameter values for the stratum probabilities (the
probability that Y belongs to certain stratum) e.g. pis = c(0.1, 0.8, 0.1).}

\item{a}{vector of cutpoints for the primary response (e.g., a = c(-2.5,2))}

\item{rs.size}{size of the SRS (simple random sample)}

\item{size}{vector of the stratum sizes of the supplemental samples (e.g. size
= c(50,0,50) represents that two supplemental samples each of size 50 are
taken from the upper and lower tail of Y.)}

\item{strat}{vector that indicates the stratum numbers (e.g. strat = c(1,2,3)
represents that there are three stratums).}
}
\value{
A list which contains the parameter estimates for the primary response
 model: \deqn{Y=\beta_{0}+\beta_{1}X+\epsilon,}{Y = beta0 + beta1*X +
 epsilon,} where epsilon has variance sig. The list contains the following
 components: \item{beta}{parameter estimates for beta} \item{sig}{estimates
 for sig} \item{pis}{estimates for the stratum probabilities}
 \item{grad}{gradient} \item{hess}{hessian} \item{converge}{whether the
 algorithm converges: True or False} \item{i}{Number of iterations}
}
\description{
\code{odsmle} provides a maximum semiparametric empirical likelihood estimator
(MSELE) for analyzing the primary outcome Y with respect to expensive exposure
and other covariates in ODS design (Zhou et al. 2002).
}
\details{
We assume that in the population, the primary outcome variable Y follows the
following model: \deqn{Y=\beta_{0}+\beta_{1}X+\epsilon,}{Y = beta0 + beta1*X +
epsilon,} where X are the covariates, and epsilon has variance sig. In ODS
design, a simple random sample is taken from the full cohort, then two
supplemental samples are taken from two tails of Y, i.e. (-Infty, mu_Y -
a*sig_Y) and (mu_Y + a*sig_Y, +Infty). Because ODS data has biased sampling
nature, naive regression analysis will yield biased estimates of the
population parameters. Zhou et al. (2002) describes a semiparametric empirical
likelihood estimator for estimating the parameters in the primary outcome
model.
}
\examples{
library(ODS)
# take the example data from the ODS package
# please see the documentation for details about the data set ods_data

Y <- ods_data[,1]
X <- cbind(rep(1,length(Y)), ods_data[,2:5])

# use the simple random sample to get an initial estimate of beta, sig #
# perform an ordinary least squares #
SRS <- ods_data[1:200,]
OLS.srs <- lm(SRS[,1] ~ SRS[,2:5])
OLS.srs.summary <- summary(OLS.srs)

beta <- coefficients(OLS.srs)
sig <- OLS.srs.summary$sigma^2
pis <- c(0.1,0.8,0.1)

# the cut points for this data is Y < 0.162, Y > 2.59.
a <- c(0.162,2.59)
rs.size <- 200
size <- c(100,0,100)
strat <- c(1,2,3)

odsmle(Y,X,beta,sig,pis,a,rs.size,size,strat)
}
