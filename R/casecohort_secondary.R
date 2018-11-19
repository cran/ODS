### code written by Yinghao Pan ###

### 2016-07-21 ###
### secondary analysis in case cohort study ###
### this code works for Z3 being discrete, and Z1, Z2, Z3 has more than 1 dimension ###
### an estimated likelihood approach which uses information in the full cohort ###
### use empirical distribution to estimate the conditional distribution of X given Z3 ###

### tildeT: event time
### C: Censoring time
### T: observation time     T=min(tildeT,C)
### Delta: event indicator    Delta=I(tildeT<C)
### Y2: secondary response
### X: expensive covariate
### Z1: other covariates that are in the linear regression model
### Z2: other covariates that are in the cox model
### Z3: other covariates that are related to the conditional distribution of X given Z

### define all functions needed ###
rep_row <- function(x,n) {
  matrix(rep(x, each=n), nrow=n)
}

#' Secondary analysis in case-cohort data
#'
#' \code{secondary_casecohort} performs the secondary analysis which describes
#' the association between a continuous secondary outcome and the expensive
#' exposure for case-cohort data.
#'
#' @export
#' @param SRS A data frame for subjects in the simple random sample. The first
#'   column is T: observation time for time-to-event outcome. The second column
#'   is Delta: the event indicator. The thid column is Y2: the continuous scale
#'   secondary outcome. The fourth column is X: the expensive exposure. Starting
#'   from the fifth column to the end are Z1, Z2 and Z3. Z1 is the set of
#'   covariates that are included in the linear regression model of the
#'   secondary outcome. Z2 is the set of covariates that are included in the Cox
#'   model (the proportional hazards model which relates the primary failure
#'   time to covariates). Z3 is the set of covariates that are related to the
#'   conditional distribution of X given other covariates.
#' @param CCH A data frame for subjects in the case-cohort sample. The
#'   case-cohort sample includes the simple random sample (SRS) and the
#'   supplemental cases. The data structure is the same as SRS.
#' @param NVsample A data frame for subjects in the non-validation sample.
#'   Subjects in the non-validation sample don't have the expensive exposure X
#'   measured. The data structure is the following: The first column is T. The
#'   second column is Delta. The thid column is Y2. Starting from the fourth
#'   column to the end are Z1, Z2 and Z3.
#' @param Z1.dim Dimension of Z1.
#' @param Z2.dim Dimension of Z2.
#' @param Z3.dim Dimension of Z3. Note here that the algorithm requires Z3 to be
#'   discrete and not high-dimensional, because we use the SRS sample to
#'   estimate the conditional distribution of X given other covariates.
#' @return A list which contains parameter estimates, estimated standard error
#'   for the primary outcome model: \deqn{\lambda
#'   (t)=\lambda_{0}(t)\exp{\gamma_{1}Y_{2}+\gamma_{2}X+\gamma_{3}Z2},}{lambda(t)=lambda0(t)exp{gamma1*Y2+gamma2*X+gamma3*Z2}}
#'    and the secondary outcome model:
#'   \deqn{Y_{2}=\beta_{0}+\beta_{1}X+\beta_{2}Z_{1}.}{Y2 = beta0 + beta1*X +
#'   beta2*Z1.} The list contains the following components:
#'   \item{gamma_paramEst}{parameter estimates for gamma in the primary outcome
#'   model} \item{gamma_stdErr}{estimated standard error for gamma in the
#'   primary outcome model} \item{beta_paramEst}{parameter estimates for beta in
#'   the secondary outcome model} \item{beta_stdErr}{estimated standard error
#'   for beta in the secondary outcome model}
#' @examples
#' \donttest{
#' library(ODS)
#' # take the example data from the ODS package
#' # please see the documentation for details about the data set casecohort_data_secondary
#' data <- casecohort_data_secondary
#'
#' # obtain SRS, CCH and NVsample from the original cohort data based on subj_ind
#' SRS <- data[data[,1]==1, 2:ncol(data)]
#' CCH <- data[data[,1]==1 | data[,1]==2, 2:ncol(data)]
#' NVsample <- data[data[,1]==0, 2:ncol(data)]
#'
#' # delete the fourth column (columns for X) from the non-validation sample
#' NVsample <- NVsample[,-4]
#'
#' Z1.dim <- 4
#' Z2.dim <- 3
#' Z3.dim <- 3
#' secondary_casecohort(SRS, CCH, NVsample, Z1.dim, Z2.dim, Z3.dim)
#' }

secondary_casecohort <- function(SRS, CCH, NVsample, Z1.dim, Z2.dim, Z3.dim) {

  p1 <- Z1.dim
  p2 <- Z2.dim
  p3 <- Z3.dim

  SRS.Z3 <- as.matrix(SRS[,(5+p1+p2):(4+p1+p2+p3)])

  n.SRS <- nrow(SRS)                       # size of the SRS sample
  Nc <- nrow(CCH)                          # size of the case-cohort sample
  N.NVsample <- nrow(NVsample)             # size of the non-validation sample
  N <- nrow(CCH) + nrow(NVsample)          # total number of subjects in the full cohort

  n1 <- sum(CCH$Delta==1)          # number of cases in the CCH sample
  n0 <- sum(CCH$Delta==0)          # number of controls in the CCH sample
  weight <- ifelse(CCH$Delta==1, 1, N/n.SRS)
  normweight <- weight/mean(weight)
  CCH <- cbind(CCH, normweight)
  CCH <- data.frame(CCH)

  ### get the initial parameter estimates of beta, gamma, sigma, baseline cumulative hazard and baseline hazard ###

  CCH.data <- data.matrix(CCH)         # a matrix form of the CCH data
  weighted.lm <- stats::lm(CCH.data[,3] ~ CCH.data[,4:(4+p1)], weights=CCH$normweight)
  beta0.initial <- as.numeric(stats::coef(weighted.lm)[1])
  beta1.initial <- as.numeric(stats::coef(weighted.lm)[2])
  beta2.initial <- as.numeric(stats::coef(weighted.lm)[3:length(stats::coef(weighted.lm))])
  sigma.initial <- summary(weighted.lm)$sigma

  weighted.cox <- survival::coxph(survival::Surv(CCH.data[,1], CCH.data[,2]) ~ CCH.data[,c(3,4,(5+p1):(4+p1+p2))], weights=CCH$normweight)
  gamma1.initial <- as.numeric(weighted.cox$coefficients[1])
  gamma2.initial <- as.numeric(weighted.cox$coefficients[2])
  gamma3.initial <- as.numeric(weighted.cox$coefficients[3:length(weighted.cox$coefficients)])

  base <- survival::basehaz(weighted.cox, centered=FALSE)
  cumulative.hazard <- base$hazard
  observation.times <- base$time
  #	plot(observation.times, cumulative.hazard)

  pair.cumulative.hazard <- function(x) {
    max(utils::tail(cumulative.hazard[observation.times<=x], 1), 0)
  }

  f2 <- function(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {

    SRS.X <- SRS$X[SRS.Z3[,1]==as.numeric(Z3[1]) & SRS.Z3[,2]==as.numeric(Z3[2]) & SRS.Z3[,3]==as.numeric(Z3[3])]
    den <- length(SRS.X)
    par <- c(gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)
    calmatrix <- cbind(rep(T,den), rep(Delta,den), rep(Y2,den), SRS.X, rep_row(Z1,den), rep_row(Z2,den), rep_row(Z3,den), rep(cumulative.hazard,den), rep_row(par, den))

    T.M <- calmatrix[,1]
    Delta.M <- calmatrix[,2]
    Y2.M <- calmatrix[,3]
    X.M <- calmatrix[,4]
    Z1.M <- calmatrix[,5:(4+p1)]
    Z2.M <- calmatrix[,(5+p1):(4+p1+p2)]
    Z3.M <- calmatrix[,(5+p1+p2):(4+p1+p2+p3)]
    cum.hazard.M <- calmatrix[,(5+p1+p2+p3)]
    gamma1.M <- calmatrix[,(6+p1+p2+p3)]
    gamma2.M <- calmatrix[,(7+p1+p2+p3)]
    gamma3.M <- calmatrix[1,(8+p1+p2+p3):(7+p1+2*p2+p3)]
    beta0.M <- calmatrix[,(8+p1+2*p2+p3)]
    beta1.M <- calmatrix[,(9+p1+2*p2+p3)]
    beta2.M <- calmatrix[1,(10+p1+2*p2+p3):(9+2*p1+2*p2+p3)]
    sigmasquare.M <- calmatrix[,(10+2*p1+2*p2+p3)]

    linearcomponent.cox <- gamma1.M*Y2.M+gamma2.M*X.M+Z2.M%*%gamma3.M
    part1 <- exp(-cum.hazard.M*exp(linearcomponent.cox))

    residual.linear <- Y2.M-beta0.M-beta1.M*X.M-Z1.M%*%beta2.M
    part2 <- 1/(sqrt(2*pi*sigmasquare.M))*exp(-residual.linear^2/(2*sigmasquare.M))

    return(log(sum(part1*part2)/den))
  }

  f2.vectorized <- function(x) {

    T <- x[1]
    Delta <- x[2]
    Y2 <- x[3]
    Z1 <- x[4:(3+p1)]
    Z2 <- x[(4+p1):(3+p1+p2)]
    Z3 <- x[(4+p1+p2):(3+p1+p2+p3)]
    cumulative.hazard <- x[(4+p1+p2+p3)]
    gamma1 <- x[(5+p1+p2+p3)]
    gamma2 <- x[(6+p1+p2+p3)]
    gamma3 <- x[(7+p1+p2+p3):(6+p1+2*p2+p3)]
    beta0 <- x[(7+p1+2*p2+p3)]
    beta1 <- x[(8+p1+2*p2+p3)]
    beta2 <- x[(9+p1+2*p2+p3):(8+2*p1+2*p2+p3)]
    sigmasquare <- x[(9+2*p1+2*p2+p3)]

    return(f2(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare))
  }

  LogLike <- function(x) {

    gamma1 <- x[1]
    gamma2 <- x[2]
    gamma3 <- x[3:(2+p2)]
    beta0 <- x[3+p2]
    beta1 <- x[4+p2]
    beta2 <- x[(5+p2):(4+p2+p1)]
    sigmasquare <- x[5+p2+p1]

    linearcomponent.cox <- gamma1*CCH.Y2+gamma2*CCH.X+CCH.Z2%*%gamma3
    residual.linear <- CCH.Y2-beta0-beta1*CCH.X-CCH.Z1%*%beta2

    part1 <- sum(CCH.Delta*linearcomponent.cox)
    part2 <- sum(CCH.cumulative.hazard.CCH*exp(linearcomponent.cox))
    part3 <- sum(residual.linear^2)/(2*sigmasquare)+Nc*log(2*pi*sigmasquare)/2

    LogLike.Vsample <- part1 - part2 - part3

    calmatrix.NVsample <- cbind(NVsample, rep_row(x, N.NVsample))
    LogLike.NVsample <- sum(apply(calmatrix.NVsample, 1, f2.vectorized))

    LogLike <- LogLike.Vsample + LogLike.NVsample
    return(-LogLike)
  }

  density.estimate <- function(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {
    return(exp(f2(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)))
  }

  partialdensity.estimate <- function(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {

    SRS.X <- SRS$X[SRS.Z3[,1]==as.numeric(Z3[1]) & SRS.Z3[,2]==as.numeric(Z3[2]) & SRS.Z3[,3]==as.numeric(Z3[3])]
    den <- length(SRS.X)
    par <- c(gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)
    calmatrix <- cbind(rep(T,den), rep(Delta,den), rep(Y2,den), SRS.X, rep_row(Z1,den), rep_row(Z2,den), rep_row(Z3,den), rep(cumulative.hazard,den), rep_row(par, den))

    T.M <- calmatrix[,1]
    Delta.M <- calmatrix[,2]
    Y2.M <- calmatrix[,3]
    X.M <- calmatrix[,4]
    Z1.M <- calmatrix[,5:(4+p1)]
    Z2.M <- calmatrix[,(5+p1):(4+p1+p2)]
    Z3.M <- calmatrix[,(5+p1+p2):(4+p1+p2+p3)]
    cum.hazard.M <- calmatrix[,(5+p1+p2+p3)]
    gamma1.M <- calmatrix[,(6+p1+p2+p3)]
    gamma2.M <- calmatrix[,(7+p1+p2+p3)]
    gamma3.M <- calmatrix[1,(8+p1+p2+p3):(7+p1+2*p2+p3)]
    beta0.M <- calmatrix[,(8+p1+2*p2+p3)]
    beta1.M <- calmatrix[,(9+p1+2*p2+p3)]
    beta2.M <- calmatrix[1,(10+p1+2*p2+p3):(9+2*p1+2*p2+p3)]
    sigmasquare.M <- calmatrix[,(10+2*p1+2*p2+p3)]

    linearcomponent.cox <- gamma1.M*Y2.M+gamma2.M*X.M+Z2.M%*%gamma3.M
    part1 <- exp(-cum.hazard.M*exp(linearcomponent.cox))

    residual.linear <- Y2.M-beta0.M-beta1.M*X.M-Z1.M%*%beta2.M
    part2 <- 1/(sqrt(2*pi*sigmasquare.M))*exp(-residual.linear^2/(2*sigmasquare.M))

    part3.gamma1 <- -cum.hazard.M*Y2.M*exp(linearcomponent.cox)
    part3.gamma2 <- -cum.hazard.M*X.M*exp(linearcomponent.cox)
    part3.gamma3 <- matrix(0, nrow=den, ncol=p2)
    for(k in 1:p2) {
      part3.gamma3[ ,k] <- -cum.hazard.M*Z2.M[, k]*exp(linearcomponent.cox)
    }

    part3.beta0 <- residual.linear/sigmasquare.M
    part3.beta1 <- residual.linear*X.M/sigmasquare.M
    part3.beta2 <- matrix(0, nrow=den, ncol=p1)
    for(k in 1:p1) {
      part3.beta2[,k] <- residual.linear*Z1.M[, k]/sigmasquare.M
    }
    part3.sigmasquare <- (-1/(2*sigmasquare.M)+residual.linear^2/(2*sigmasquare.M^2))*part2

    gradient.gamma1 <- sum(part1*part2*part3.gamma1)/den
    gradient.gamma2 <- sum(part1*part2*part3.gamma2)/den
    gradient.gamma3 <- rep(0, p2)
    for(k in 1:p2) {
      gradient.gamma3[k] <- sum(part1*part2*part3.gamma3[,k])/den
    }

    gradient.beta0 <- sum(part1*part2*part3.beta0)/den
    gradient.beta1 <- sum(part1*part2*part3.beta1)/den
    gradient.beta2 <- rep(0, p1)
    for(k in 1:p1) {
      gradient.beta2[k] <- sum(part1*part2*part3.beta2[,k])/den
    }
    gradient.sigmasquare <- sum(part1*part2*part3.sigmasquare)/den

    output <- c(gradient.gamma1, gradient.gamma2, gradient.gamma3, gradient.beta0, gradient.beta1, gradient.beta2, gradient.sigmasquare)
    return(output)
  }

  partialdensity <- function(T, Delta, Y2, X, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {

    linearcomponent.cox <- gamma1*Y2+gamma2*X+Z2%*%gamma3
    part1 <- exp(-cumulative.hazard*exp(linearcomponent.cox))

    residual.linear <- Y2-beta0-beta1*X-Z1%*%beta2
    part2 <- 1/(sqrt(2*pi*sigmasquare))*exp(-residual.linear^2/(2*sigmasquare))

    part3.gamma1 <- -cumulative.hazard*Y2*exp(linearcomponent.cox)
    part3.gamma2 <- -cumulative.hazard*X*exp(linearcomponent.cox)
    part3.gamma3 <- rep(0, p2)
    for(k in 1:p2) {
      part3.gamma3[k] <- -cumulative.hazard*Z2[k]*exp(linearcomponent.cox)
    }

    part3.beta0 <- residual.linear/sigmasquare
    part3.beta1 <- residual.linear*X/sigmasquare
    part3.beta2 <- rep(0, p1)
    for(k in 1:p1) {
      part3.beta2[k] <- residual.linear*Z1[k]/sigmasquare
    }
    part3.sigmasquare <- (-1/(2*sigmasquare)+residual.linear^2/(2*sigmasquare^2))*part2

    gradient.gamma1 <- part1*part2*part3.gamma1
    gradient.gamma2 <- part1*part2*part3.gamma2
    gradient.gamma3 <- rep(0, p2)
    for(k in 1:p2) {
      gradient.gamma3[k] <- part1*part2*part3.gamma3[k]
    }

    gradient.beta0 <- part1*part2*part3.beta0
    gradient.beta1 <- part1*part2*part3.beta1
    gradient.beta2 <- rep(0, p1)
    for(k in 1:p1) {
      gradient.beta2[k] <- part1*part2*part3.beta2[k]
    }
    gradient.sigmasquare <- part1*part2*part3.sigmasquare

    output <- c(gradient.gamma1, gradient.gamma2, gradient.gamma3, gradient.beta0, gradient.beta1, gradient.beta2, gradient.sigmasquare)
    return(output)

  }

  density <- function(T, Delta, Y2, X, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {

    linearcomponent.cox <- gamma1*Y2+gamma2*X+Z2%*%gamma3
    part1 <- exp(-cumulative.hazard*exp(linearcomponent.cox))

    residual.linear <- Y2-beta0-beta1*X-Z1%*%beta2
    part2 <- 1/(sqrt(2*pi*sigmasquare))*exp(-residual.linear^2/(2*sigmasquare))

    return(part1*part2)
  }

  score <- function(T, Delta, Y2, X, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {

    part1 <- partialdensity(T, Delta, Y2, X, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)/density.estimate(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)
    part2 <- density(T, Delta, Y2, X, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)*partialdensity.estimate(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)/(density.estimate(T, Delta, Y2, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)^2)
    denominator <- sum(SRS.Z3[,1]==as.numeric(Z3[1]) & SRS.Z3[,2]==as.numeric(Z3[2]) & SRS.Z3[,3]==as.numeric(Z3[3]))

    return((part1-part2)/denominator)
  }

  score.vectorized <- function(x) {

    T <- x[1]
    Delta <- x[2]
    Y2 <- x[3]
    X <- x[4]
    Z1 <- x[5:(4+p1)]
    Z2 <- x[(5+p1):(4+p1+p2)]
    Z3 <- x[(5+p1+p2):(4+p1+p2+p3)]
    cumulative.hazard <- x[(5+p1+p2+p3)]
    gamma1 <- x[(6+p1+p2+p3)]
    gamma2 <- x[(7+p1+p2+p3)]
    gamma3 <- x[(8+p1+p2+p3):(7+p1+2*p2+p3)]
    beta0 <- x[(8+p1+2*p2+p3)]
    beta1 <- x[(9+p1+2*p2+p3)]
    beta2 <- x[(10+p1+2*p2+p3):(9+2*p1+2*p2+p3)]
    sigmasquare <- x[(10+2*p1+2*p2+p3)]

    return(score(T, Delta, Y2, X, Z1, Z2, Z3, cumulative.hazard, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare))
  }


  Q <- function(X, Z1, Z2, Z3, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare) {

    NVsample.Z <- NVsample[NVsample.Z3[,1]==as.numeric(Z3[1]) & NVsample.Z3[,2]==as.numeric(Z3[2]) & NVsample.Z3[,3]==as.numeric(Z3[3]),]
    N.NVZ <- sum(NVsample.Z3[,1]==as.numeric(Z3[1]) & NVsample.Z3[,2]==as.numeric(Z3[2]) & NVsample.Z3[,3]==as.numeric(Z3[3]))
    par <- c(gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare)
    cal <- cbind(NVsample.Z[,1:3], rep(X, N.NVZ), NVsample.Z[,4:(p1+p2+p3+4)], rep_row(par, N.NVZ))
    result <- apply(cal, 1, score.vectorized)

    num.param <- p1+p2+5
    Qresult <- rep(0, num.param)
    for(k in 1:num.param) {
      Qresult[k] <- sum(result[c(seq(k, num.param*N.NVZ, by=num.param))])
    }

    return(Qresult)
  }

  Q.vectorized <- function(x) {

    X <- x[1]
    Z1 <- x[2:(1+p1)]
    Z2 <- x[(2+p1):(1+p1+p2)]
    Z3 <- x[(2+p1+p2):(1+p1+p2+p3)]
    gamma1 <- x[(2+p1+p2+p3)]
    gamma2 <- x[(3+p1+p2+p3)]
    gamma3 <- x[(4+p1+p2+p3):(3+p1+2*p2+p3)]
    beta0 <- x[(4+p1+2*p2+p3)]
    beta1 <- x[(5+p1+2*p2+p3)]
    beta2 <- x[(6+p1+2*p2+p3):(5+2*p1+2*p2+p3)]
    sigmasquare <- x[(6+2*p1+2*p2+p3)]

    return(Q(X, Z1, Z2, Z3, gamma1, gamma2, gamma3, beta0, beta1, beta2, sigmasquare))
  }


  ### sort the case cohort data and add the baseline cumulative hazard into the data frame ###
  CCH <- CCH[order(CCH$T) ,]
  cumulative.hazard.CCH <- sapply(CCH$T, pair.cumulative.hazard)
  CCH <- cbind(CCH, cumulative.hazard.CCH)

  CCH.T <- CCH$T
  CCH.Delta <- CCH$Delta
  CCH.Y2 <- CCH$Y2
  CCH.X <- CCH$X
  CCH.cumulative.hazard.CCH <- CCH$cumulative.hazard.CCH

  CCH.Z1 <- as.matrix(CCH[,5:(4+p1)])
  CCH.Z2 <- as.matrix(CCH[,(5+p1):(4+p1+p2)])
  CCH.Z3 <- as.matrix(CCH[,(5+p1+p2):(4+p1+p2+p3)])

  ### sort the non-validation sample data and add the baseline cumulative hazard into the data frame ###
  cumulative.hazard.NVsample <- rep(0, N.NVsample)
  NVsample <- NVsample[order(NVsample$T) ,]

  cumulative.hazard.NVsample <- sapply(NVsample$T, pair.cumulative.hazard)
  NVsample <- cbind(NVsample, cumulative.hazard.NVsample)

  NVsample.Z3 <- as.matrix(NVsample[,(4+p1+p2):(3+p1+p2+p3)])

  ### print the initial value ###
  initialvalue <- c(gamma1.initial, gamma2.initial, gamma3.initial, beta0.initial, beta1.initial, beta2.initial, sigma.initial^2)
#  print(initialvalue)

  CCH <- as.matrix(CCH)
  NVsample <- as.matrix(NVsample)

  ### starting maximizing the estimated log likelihood function ###
  optim <- stats::nlm(LogLike, initialvalue, hessian=TRUE)
#  print(optim)

  ### record the parameter estimates ###
  etahat <- optim$estimate
  gammahat <- etahat[1:(2+p2)]
  betahat <- etahat[(3+p2):(length(etahat)-1)]

  ### calculating the asymptotic covariance matrix ###
  I <- optim$hessian/N
  I.inverse <- solve(I)

  ### method II to get SIGMA ###
  cal <- cbind(SRS[,4:(4+p1+p2+p3)], rep_row(optim$estimate, n.SRS))
  score.matrix <- t(apply(cal,1,Q.vectorized))
  SIGMA <- stats::cov(score.matrix)

  SIGMA.MEL <- I.inverse+n.SRS/N*I.inverse%*%SIGMA%*%I.inverse

  estimatedstderr <- sqrt(diag(SIGMA.MEL)/N)
  estimatedstderrgamma <- estimatedstderr[1:(2+p2)]
  estimatedstderrbeta <- estimatedstderr[(3+p2):(length(etahat)-1)]

  result <- list(gamma_paramEst = gammahat, gamma_stdErr = estimatedstderrgamma, beta_paramEst = betahat, beta_stdErr = estimatedstderrbeta)

  return(result)
}

#' Data example for the secondary analysis in case-cohort design
#'
#' @format A data frame with 1000 rows and 15 columns: \describe{
#'   \item{subj_ind}{An indicator variable for each subject: 1 = SRS, 2 =
#'   supplemental cases, 0 = NVsample} \item{T}{observation time for failure
#'   outcome} \item{Delta}{event indicator} \item{Y2}{a continuous secondary
#'   outcome}\item{X}{expensive exposure} \item{Z11}{first covariate in the
#'   linear regression model} \item{Z12}{second covariate in the linear
#'   regression model}\item{Z13}{third covariate in the linear regression model}
#'   \item{Z14}{fourth covariate in the linear regression model}\item{Z21}{first
#'   covariate in the Cox model} \item{Z22}{second covariate in the Cox
#'   model}\item{Z23}{third covariate in the Cox model} \item{Z31}{first
#'   covariate that is related to the conditional distribution of X given other
#'   covariates}\item{Z32}{second covariate that is related to the conditional
#'   distribution} \item{Z33}{thid covariate that is related to the conditional
#'   distribution} }
#'
#' @source A simulated data set
"casecohort_data_secondary"
