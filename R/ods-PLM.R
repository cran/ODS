## the original matlab code is written by Guoyou Qin ##
## revised and rewritten by Yinghao Pan into R code ##
## implements the partial linear model for ODS data ##
## reference: Zhou, Qin and Longnecker (2011) ##

#' Generate logarithmically spaced vector
#'
#' \code{logspace} generates n logarithmically spaced points between 10^d1 and
#' 10^d2. The utility of this function is equivalent to logspace function in
#' matlab.
#'
#' @export
#' @param d1 first bound
#' @param d2 second bound
#' @param n number of points
#' @return a vector of n logarithmically spaced points between 10^d1 and 10^d2.
#' @examples
#' logspace(-6,7,30)

logspace <- function(d1, d2, n) {

  exp(log(10)*seq(d1, d2, length.out=n))

}

# function that rounds x to the nearest integer #
# different from the round function in package base, round(2.5)=2, round2(2.5)=3
round2 <- function(x) {trunc(x+sign(x)*0.5)}

#' Create knots at sample quantiles
#'
#' \code{quantileknots} creates knots at sample quantiles
#'
#' @export
#' @param x a vector. The knots are at sample quantiles of x.
#' @param nknots number of knots
#' @param boundstab parameter for boundary stability. The default is 0. If
#'   boundstab = 1, then nknots+2 knots are created and the first and last are
#'   deleted. This mitigates the extra variability of regression spline
#'   estimates near the boundaries.
#' @return a vector of knots at sample quantiles of x.
#' @examples
#' library(ODS)
#'
#' x <- c(1, 2, 3, 4, 5)
#' quantileknots(x, 3, 0)

quantileknots <- function(x, nknots, boundstab) {

  nargin <- length(as.list(match.call())) - 1

  if (nargin < 3) {boundstab = 0}
  x = unique(x)
  n = length(x)
  xsort = sort(x)
  loc = n*(1:(nknots+2*boundstab)) / (nknots+1+2*boundstab)
  knots=xsort[round2(loc)]
  knots=knots[(1 + boundstab) : (nknots + boundstab)]
  return (knots)
  #  REMOVE KNOTS NEAR BOUNDARIRES FOR
  #  STABILITY (= LOW VARIABILITY)
}

#' power basis functions of a spline of given degree
#'
#' \code{Bfct} returns the power basis functions of a spline of given degree.
#'
#' @export
#' @param x n by 1 matrix of the independent variable
#' @param degree the order of spline
#' @param knots the knots location
#' @param der the der-order derivative. The default is 0
#' @return n by (1+degree+nknots) matrix corresponding to the truncated power
#'   spline basis with knots and specified degree.
#' @examples
#' library(ODS)
#'
#' x <- matrix(c(1,2,3,4,5),ncol=1)
#' degree <- 2
#' knots <- c(1,3,4)
#'
#' Bfct(x, degree, knots)

Bfct <- function(x, degree, knots, der) {

  nargin <- length(as.list(match.call())) - 1
  if (nargin < 4) {der = 0}

  if (der > degree) {
    print('********************************************************')
    print('********************************************************')
    print('WARNING:  der > degree --- xm not returned by powerbasis')
    print('********************************************************')
    print('********************************************************')
    return (0);
  }

  n=dim(x)[1]
  nknots = length(knots)
  if (der == 0) {
    xm = matrix(1, n, 1)
  } else {
    xm = matrix(0, n, 1)
  }

  for (i in 1:degree) {
    if (i < der) {
      xm = cbind(xm, matrix(0, n, 1))
    } else {
      xm = cbind(xm, x^(i-der))
    }
  }

  if (nknots > 0) {
    for (i in 1:nknots) {
      xm = cbind(xm, (x-knots[i])^(degree-der)*(x > knots[i]))
    }
  }

  return (xm)

}

pI3_t2 <- function(y,D,theta0,sig0_sq0,Qh) {

  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  f_yx = stats::dnorm(y,mu_yx,var_yx^0.5)
  I3_t2 = 1/sig0_sq0 * t(D) %*% diag(as.vector(f_yx)/as.vector(Qh), nrow=length(Qh)) %*% (y-D %*% theta0)
  return (I3_t2)
}

pI3_t2s <- function(y,D,theta0,sig0_sq0,Qh) {

  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  f_yx = stats::dnorm(y,mu_yx,var_yx^0.5)
  I3_t2 = ((y-mu_yx)^2/sig0_sq0-1)*f_yx/Qh
  return (I3_t2)
}

pI3_t3s <- function(y,D,theta0,sig0_sq0,Qh) {

  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  f_yx = stats::dnorm(y,mu_yx,var_yx^0.5)
  I3_t2 = (-((y-mu_yx)^2/sig0_sq0^2)+((y-mu_yx)^2/sig0_sq0-1)^2/(2*sig0_sq0))*f_yx/Qh
  return (I3_t2)
}

pI3_t22 <- function(y,D,theta0,sig0_sq0,Qh) {

  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  f_yx = stats::dnorm(y,mu_yx,var_yx^0.5)
  I3_t2 = (-1/sig0_sq0*f_yx+1/(sig0_sq0^2)*(f_yx)*(y-D %*% theta0)^2)/Qh
  return (I3_t2)
}

pI3_t33 <- function(y,D,theta0,sig0_sq0,Qh) {

  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  f_yx = stats::dnorm(y,mu_yx,var_yx^0.5)
  I3_t2 = 1/sig0_sq0 * diag(as.vector(f_yx)/as.vector(Qh), nrow=length(Qh)) %*% (y-D %*% theta0)
  return (I3_t2)
}

#' Generalized cross-validation for ODS data
#'
#' \code{gcv_ODS} calculates the generalized cross-validation (GCV) for
#' selecting the smoothing parameter in the setting of outcome-dependent
#' sampling. The details can be seen in Zhou, Qin and Longnecker (2011) and its
#' supplementary materials.
#'
#' @export
#' @param X n by 1 matrix of the observed exposure variable
#' @param Y n by 1 matrix of the observed outcome variable
#' @param Z n by p matrix of the other covariates
#' @param n_f n_f = c(n0, n1, n2), where n0 is the SRS sample size, n1 is the
#'   size of the supplemental sample chosen from (-infty, mu_Y-a*sig_Y), n2 is
#'   the size of the supplemental sample chosen from (mu_Y+a*sig_Y, +infty).
#' @param eta a column matrix. eta = (theta^T pi^T v^T sig0_sq)^T where
#'   theta=(alpha^T, gamma^T)^T. We refer to Zhou, Qin and Longnecker (2011) for
#'   details of these notations.
#' @param q_s smoothing parameter
#' @param Cpt cut point a
#' @param mu_Y mean of Y in the population
#' @param sig_Y standard deviation of Y in the population
#' @param degree degree of the truncated power spline basis, default value is 2
#' @param nknots number of knots of the truncated power spline basis, default
#'   value is 10
#' @return the value of generalized cross-validation score
#' @examples
#' \donttest{
#' library(ODS)
#' # take the example data from the ODS package
#' # please see the documentation for details about the data set ods_data
#'
#' nknots = 10
#' degree = 2
#'
#' # get the initial value of the parameters from standard linear regression based on SRS data #
#' dataSRS = ods_data[1:200,]
#' YS = dataSRS[,1]
#' XS = dataSRS[,2]
#' ZS = dataSRS[,3:5]
#'
#' knots = quantileknots(XS, nknots, 0)
#' # the power basis spline function
#' MS = Bfct(as.matrix(XS), degree, knots)
#' DS = cbind(MS, ZS)
#' theta00 = as.numeric(lm(YS ~ DS -1)$coefficients)
#' sig0_sq00 = var(YS - DS %*% theta00)
#' pi00 = c(0.15, 0.15)
#' v00 = c(0, 0)
#' eta00 = matrix(c(theta00, pi00, v00, sig0_sq00), ncol=1)
#' mu_Y = mean(YS)
#' sig_Y = sd(YS)
#'
#' Y = matrix(ods_data[,1])
#' X = matrix(ods_data[,2])
#' Z = matrix(ods_data[,3:5], nrow=400)
#'
#' # In this ODS data, the supplemental samples are taken from (-Infty, mu_Y-a*sig_Y) #
#' # and (mu_Y+a*sig_Y, +Infty), where a=1 #
#' n_f = c(200, 100, 100)
#' Cpt = 1
#'
#' # GCV selection to find the optimal smoothing parameter #
#' q_s1 = logspace(-6, 7, 10)
#' gcv1 = rep(0, 10)
#'
#' for (j in 1:10) {
#'
#'   result = Estimate_PLMODS(X,Y,Z,n_f,eta00,q_s1[j],Cpt,mu_Y,sig_Y)
#'   etajj = matrix(c(result$alpha, result$gam, result$pi0, result$v0, result$sig0_sq0), ncol=1)
#'   gcv1[j] = gcv_ODS(X,Y,Z,n_f,etajj,q_s1[j],Cpt,mu_Y,sig_Y)
#' }
#'
#' b = which(gcv1 == min(gcv1))
#' q_s = q_s1[b]
#'
#' q_s
#'
#' # Estimation of the partial linear model in the setting of outcome-dependent sampling #
#' result = Estimate_PLMODS(X, Y, Z, n_f, eta00, q_s, Cpt, mu_Y, sig_Y)
#' result
#' }

gcv_ODS <- function(X,Y,Z,n_f,eta,q_s,Cpt,mu_Y,sig_Y,degree,nknots) {

  # check input
  nargin <- length(as.list(match.call())) - 1
  if (nargin <= 10) {degree = 2}
  if (nargin <= 11) {nknots = 10}

  ## Initializing ##

  n0 = n_f[1]
  n1 = n_f[2]
  n2 = n_f[3]
  n = sum(n_f)

  knots = quantileknots(X, nknots, 0)
  n3 = dim(Z)[1]
  p = dim(Z)[2]
  Psi = diag(c(rep(0, 1+degree), rep(1, nknots), rep(0, p)))
  M = Bfct(X, degree, knots)
  D = cbind(M, Z)

  low1 = mu_Y - 5*Cpt*sig_Y
  up1 = mu_Y - Cpt*sig_Y
  low2 = mu_Y + Cpt*sig_Y
  up2 = mu_Y + 5*Cpt*sig_Y

  Td = 1 + degree + nknots + p
  eta00 = eta
  theta0 = eta00[1:Td,]
  pi0 = eta00[(Td+1):(Td+2),]
  v0 = eta00[(Td+3):(Td+4),]
  sig0_sq0 = eta00[Td+5]

  I11 = -1/sig0_sq0*t(D) %*% D
  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  psi1 = stats::pnorm(up1,mu_yx,var_yx^0.5)-stats::pnorm(low1,mu_yx,var_yx^0.5)
  psi2 = stats::pnorm(up2,mu_yx,var_yx^0.5)-stats::pnorm(low2,mu_yx,var_yx^0.5)
  Qh = n0/n+ (n1/(n*pi0[1]))*psi1+(n2/(n*pi0[2]))*psi2
  psi1_f = cubature::adaptIntegrate(pI3_t2, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=dim(D)[2])$integral
  psi2_f = cubature::adaptIntegrate(pI3_t2, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=dim(D)[2])$integral
  psi1_s = t(D) %*% diag(cubature::adaptIntegrate(pI3_t22, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral) %*% D
  psi2_s = t(D) %*% diag(cubature::adaptIntegrate(pI3_t22, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral) %*% D

  Qh_f = -((n1/(n*pi0[1])+v0[1])*psi1_f+(n2/(n*pi0[2])+v0[2])*psi2_f)
  Qh_s1 = cubature::adaptIntegrate(pI3_t33, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
  Qh_s2 = cubature::adaptIntegrate(pI3_t33, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral

  Qh_s = -((n1/(n*pi0[1]))*psi1_s+(n2/(n*pi0[2]))*psi2_s)+t(D) %*% diag(((n1/(n*pi0[1]))*Qh_s1+(n2/(n*pi0[2]))*Qh_s2)^2) %*% D

  I22 = Qh_s
  V = I11+I22-n*q_s*Psi

  ## MLE ##
  ###############################

  mu_yx = D %*% theta0
  var_yx = sig0_sq0
  psi1 = stats::pnorm(up1,mu_yx,var_yx^0.5)-stats::pnorm(low1,mu_yx,var_yx^0.5)
  psi2 = stats::pnorm(up2,mu_yx,var_yx^0.5)-stats::pnorm(low2,mu_yx,var_yx^0.5)
  Qh = n0/n+ (n1/(n*pi0[1]))*psi1+(n2/(n*pi0[2]))*psi2+v0[1]*(psi1-pi0[1])+v0[2]*(psi2-pi0[2])
  like_mle = sum(log(stats::dnorm(Y,mu_yx,var_yx^0.5)))-sum(log(Qh))-(n1*log(pi0[1])+n2*log(pi0[2]))

  gcv1 = -1/n*like_mle/(1-1/n*sum(diag(solve(V,(I11+I22)))))^2

  return (gcv1)
}

#' Partial linear model for ODS data
#'
#' \code{Estimate_PLMODS} computes the estimate of parameters in a partial
#' linear model in the setting of outcome-dependent sampling. See details in
#' Zhou, Qin and Longnecker (2011).
#'
#' We assume that in the population, the primary outcome variable Y follows the
#' following partial linear model:
#' \deqn{E(Y|X,Z)=g(X)+Z^{T}\gamma}{E(Y|X,Z)=g(X)+Z^{T}*gamma} where X is the
#' expensive exposure, Z are other covariates. In ODS design, a simple random
#' sample is taken from the full cohort, then two supplemental samples are taken
#' from two tails of Y, i.e. (-Infty, mu_Y - a*sig_Y) and (mu_Y + a*sig_Y,
#' +Infty). Because ODS data has biased sampling nature, naive regression
#' analysis will yield biased estimates of the population parameters. Zhou, Qin
#' and Longnecker (2011) describes a semiparametric empirical likelihood
#' estimator for estimating the parameters in the partial linear model.
#'
#' @export
#' @param X n by 1 matrix of the observed exposure variable
#' @param Y n by 1 matrix of the observed outcome variable
#' @param Z n by p matrix of the other covariates
#' @param n_f n_f = c(n0, n1, n2), where n0 is the SRS sample size, n1 is the
#'   size of the supplemental sample chosen from (-infty, mu_Y-a*sig_Y), n2 is
#'   the size of the supplemental sample chosen from (mu_Y+a*sig_Y, +infty).
#' @param eta00 a column matrix. eta00 = (theta^T pi^T v^T sig0_sq)^T where
#'   theta=(alpha^T, gamma^T)^T. We refer to Zhou, Qin and Longnecker (2011) for
#'   details of these notations.
#' @param q_s smoothing parameter
#' @param Cpt cut point a
#' @param mu_Y mean of Y in the population
#' @param sig_Y standard deviation of Y in the population
#' @param degree degree of the truncated power spline basis, default value is 2
#' @param nknots number of knots of the truncated power spline basis, default
#'   value is 10
#' @param tol convergence criteria, the default value is 1e-6
#' @param iter maximum iteration number, the default value is 30
#' @return Parameter estimates and standard errors for the partial linear model:
#'   \deqn{E(Y|X,Z)=g(X)+Z^{T}\gamma}{E(Y|X,Z)=g(X)+Z^{T}*gamma} where the
#'   unknown smooth function g is approximated by a spline function with fixed
#'   knots. The results contain the following components: \item{alpha}{spline
#'   coefficient} \item{gam}{other linear regression coefficients}
#'   \item{std_gam}{standard error of gam} \item{cov_mtxa}{covariance matrix of
#'   alpha} \item{step}{numbers of iteration requied to acheive convergence}
#'   \item{pi0}{estimated notation pi} \item{v0}{estimated notation vtheta}
#'   \item{sig0_sq0}{estimated variance of error}
#' @examples
#' \donttest{
#' library(ODS)
#' # take the example data from the ODS package
#' # please see the documentation for details about the data set ods_data
#'
#' nknots = 10
#' degree = 2
#'
#' # get the initial value of the parameters from standard linear regression based on SRS data #
#' dataSRS = ods_data[1:200,]
#' YS = dataSRS[,1]
#' XS = dataSRS[,2]
#' ZS = dataSRS[,3:5]
#'
#' knots = quantileknots(XS, nknots, 0)
#' # the power basis spline function
#' MS = Bfct(as.matrix(XS), degree, knots)
#' DS = cbind(MS, ZS)
#' theta00 = as.numeric(lm(YS ~ DS -1)$coefficients)
#' sig0_sq00 = var(YS - DS %*% theta00)
#' pi00 = c(0.15, 0.15)
#' v00 = c(0, 0)
#' eta00 = matrix(c(theta00, pi00, v00, sig0_sq00), ncol=1)
#' mu_Y = mean(YS)
#' sig_Y = sd(YS)
#'
#' Y = matrix(ods_data[,1])
#' X = matrix(ods_data[,2])
#' Z = matrix(ods_data[,3:5], nrow=400)
#'
#' # In this ODS data, the supplemental samples are taken from (-Infty, mu_Y-a*sig_Y) #
#' # and (mu_Y+a*sig_Y, +Infty), where a=1 #
#' n_f = c(200, 100, 100)
#' Cpt = 1
#'
#' # GCV selection to find the optimal smoothing parameter #
#' q_s1 = logspace(-6, 7, 10)
#' gcv1 = rep(0, 10)
#'
#' for (j in 1:10) {
#'
#'   result = Estimate_PLMODS(X,Y,Z,n_f,eta00,q_s1[j],Cpt,mu_Y,sig_Y)
#'   etajj = matrix(c(result$alpha, result$gam, result$pi0, result$v0, result$sig0_sq0), ncol=1)
#'   gcv1[j] = gcv_ODS(X,Y,Z,n_f,etajj,q_s1[j],Cpt,mu_Y,sig_Y)
#' }
#'
#' b = which(gcv1 == min(gcv1))
#' q_s = q_s1[b]
#' q_s
#'
#' # Estimation of the partial linear model in the setting of outcome-dependent sampling #
#' result = Estimate_PLMODS(X, Y, Z, n_f, eta00, q_s, Cpt, mu_Y, sig_Y)
#' result
#' }

Estimate_PLMODS <- function(X,Y,Z,n_f,eta00,q_s,Cpt,mu_Y,sig_Y,degree,nknots,tol,iter) {

  ## check input, nargin is the number of paramaters ##

  nargin <- length(as.list(match.call())) - 1
  if (nargin <= 10) {degree = 2}
  if (nargin <= 11) {nknots = 10}
  if (nargin <= 12) {tol = 1e-6}
  if (nargin <= 13) {iter = 30}

  ## Initializing ##

  n0 = n_f[1]
  n1 = n_f[2]
  n2 = n_f[3]
  n = sum(n_f)

  knots = quantileknots(X, nknots, 0)
  n3 = dim(Z)[1]
  p = dim(Z)[2]
  Psi = diag(c(rep(0, 1+degree), rep(1, nknots), rep(0, p)))
  M = Bfct(X, degree, knots)
  D = cbind(M, Z)

  ## Newton-Raphson iteration ##

  eta0 = 0
  eta1 = eta00
  Td = 1 + degree + nknots + p
  theta1 = eta1[1:Td,]
  pi1 = eta1[(Td+1):(Td+2),]
  pi100 = eta00[(Td+1):(Td+2),]
  v1 = eta1[(Td+3):(Td+4),]
  sig0_sq1 = eta1[Td+5]
  sig0_sq100 = eta00[Td+5]

  low1 = mu_Y - 5*Cpt*sig_Y
  up1 = mu_Y - Cpt*sig_Y
  low2 = mu_Y + Cpt*sig_Y
  up2 = mu_Y + 5*Cpt*sig_Y

  step = 0
  while (sqrt(sum((eta1-eta0)^2)) > tol && step < iter) {
    step = step + 1
    eta1 = eta0
    theta0 = theta1
    sig0_sq0 = sig0_sq1
    I1 = 1/sig0_sq0*t(D) %*% (Y-D %*% theta0)
    I11 = -1/sig0_sq0*t(D) %*% D
    mu_yx = D %*% theta0
    var_yx = sig0_sq0
    psi1 = stats::pnorm(up1,mu_yx,var_yx^0.5)-stats::pnorm(low1,mu_yx,var_yx^0.5)
    psi2 = stats::pnorm(up2,mu_yx,var_yx^0.5)-stats::pnorm(low2,mu_yx,var_yx^0.5)
    Qh = n0/n+ (n1/(n*pi1[1])+v1[1])*psi1+(n2/(n*pi1[2])+v1[2])*psi2
    psi1_f = cubature::adaptIntegrate(pI3_t2, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=dim(D)[2])$integral
    psi2_f = cubature::adaptIntegrate(pI3_t2, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=dim(D)[2])$integral
    psi1_s = t(D) %*% diag(cubature::adaptIntegrate(pI3_t22, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral) %*% D
    psi2_s = t(D) %*% diag(cubature::adaptIntegrate(pI3_t22, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral) %*% D
    Qh_f = -((n1/(n*pi1[1])+v1[1])*psi1_f+(n2/(n*pi1[2])+v1[2])*psi2_f)
    Qh_s1 = cubature::adaptIntegrate(pI3_t33, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
    Qh_s2 = cubature::adaptIntegrate(pI3_t33, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
    Qh_s = -((n1/(n*pi1[1]))*psi1_s+(n2/(n*pi1[2]))*psi2_s)+t(D) %*% diag(((n1/(n*pi1[1]))*Qh_s1+(n2/(n*pi1[2]))*Qh_s2)^2) %*% D
    I2 = Qh_f
    I22 = Qh_s
    M_f = n*q_s*Psi %*% theta0
    M_s = n*q_s*Psi
    I_theta_f = I1 + I2 - M_f
    I_theta_s = I11 + I22 - M_s
    SC1_t = I_theta_f
    SC2_t = I_theta_s
    theta1 = theta0 - solve(SC2_t, SC1_t)
    psi_s1 = cubature::adaptIntegrate(pI3_t2s, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
    psi_s2 = cubature::adaptIntegrate(pI3_t2s, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
    SC1_s = 1/(2*sig0_sq0)*(sum((Y-mu_yx)^2/sig0_sq0-1)-n1/(n*pi1[1])*sum(psi_s1)-n2/(n*pi1[2])*sum(psi_s2))
    psi_ss1 = cubature::adaptIntegrate(pI3_t3s, low1, up1, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
    psi_ss2 = cubature::adaptIntegrate(pI3_t3s, low2, up2, D=D, theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh, fDim=n)$integral
    Q_ss = 1/(2*sig0_sq0)*(n1/(n*pi1[1])*(psi_s1)+n2/(n*pi1[2])*(psi_s2))

    SC2_s = -0.5/sig0_sq0^2*(sum((Y-mu_yx)^2/sig0_sq0-1)-n1/(n*pi1[1])*sum(psi_s1)-n2/(n*pi1[2])*sum(psi_s2))+0.5/sig0_sq0*(-sum((Y-mu_yx)^2/sig0_sq0^2)-n1/(n*pi1[1])*sum(psi_ss1)-n2/(n*pi1[2])*sum(psi_ss2)+sum(Q_ss^2))
    sig0_sq1 = sig0_sq0 - solve(SC2_s, SC1_s)
    pi0 = pi1
    if ((pi0[1])<0 | (pi0[1])>1) {pi0 = pi100}
    mu_yx = D %*% theta1
    var_yx = sig0_sq1
    psi1 = stats::pnorm(up1,mu_yx,var_yx^0.5)-stats::pnorm(low1,mu_yx,var_yx^0.5)
    psi2 = stats::pnorm(up2,mu_yx,var_yx^0.5)-stats::pnorm(low2,mu_yx,var_yx^0.5)
    K1_2 = n0/n + (n1/(n*pi0[1]))*psi1 + (n2/(n*pi0[2]))*psi2
    K_f = rbind(-n1/pi0[1]+sum((n1/(n*pi0[1]^2)*psi1)/K1_2), -n2/pi0[2]+sum((n2/(n*pi0[2]^2)*psi2)/K1_2))
    K_ff = rbind(cbind(t((n1/(n*pi0[1]^2)*psi1)/K1_2) %*% ((n1/(n*pi0[1]^2)*psi1)/K1_2), t((n1/(n*pi0[1]^2)*psi1)/K1_2) %*% ((n2/(n*pi0[2]^2)*psi2)/K1_2)), cbind(t((n1/(n*pi0[1]^2)*psi1)/K1_2) %*% ((n2/(n*pi0[2]^2)*psi2)/K1_2), t((n2/(n*pi0[2]^2)*psi2)/K1_2) %*% ((n2/(n*pi0[2]^2)*psi2)/K1_2)))
    K_s = rbind(cbind(n1/pi0[1]^2+(-2)*sum((n1/(n*pi0[1]^3)*psi1)/K1_2), 0), cbind(0, n2/pi0[2]^2+(-2)*sum((n2/(n*pi0[2]^3)*psi2)/K1_2))) + K_ff
    SC1_pi = K_f
    SC2_pi = K_s
    pi1 = pi0 - solve(SC2_pi, SC1_pi)
    v0 = v1
    if (abs(v0[1])>0.5 | abs(v0[2])>0.5) {v0 = matrix(0,2,1)}
    psi1 = stats::pnorm(up1,mu_yx,var_yx^0.5)-stats::pnorm(low1,mu_yx,var_yx^0.5)
    psi2 = stats::pnorm(up2,mu_yx,var_yx^0.5)-stats::pnorm(low2,mu_yx,var_yx^0.5)
    L1_2 = n0/n+ (n1/(n*pi1[1])+v0[1])*psi1+(n2/(n*pi1[2])+v0[2])*psi2
    L_f = -rbind(sum(((psi1-pi1[1]))/L1_2), sum(((psi2-pi1[2]))/L1_2))
    L_s = rbind(cbind(sum(((psi1-pi1[1])/L1_2)^2), sum(((psi1-pi1[1])/L1_2)*(((psi2-pi1[2]))/L1_2))), cbind(sum(((psi1-pi1[1])/L1_2)*(((psi2-pi1[2]))/L1_2)), sum(((psi2-pi1[2])/L1_2)^2)))
    SC1_v = L_f
    SC2_v = L_s
    v1 = v0-solve(SC2_v, SC1_v)

    eta0 = t(cbind(t(theta1), t(pi1), t(v1), t(sig0_sq1)))
  }

  #################################################################

  ############### Following is the part for covariance ############

  #################################################################


  theta0 = eta0[1:Td,]
  pi0 = eta0[(Td+1):(Td+2),]
  v0 = eta0[(Td+3):(Td+4),]
  sig0_sq0 = eta0[Td+5]
  V = I_theta_s

  U = 0
  D_i = t(D)
  for (ii in 1:n) {

    psi1_fi = cubature::adaptIntegrate(pI3_t2, low1, up1, D=t(D_i[,ii]), theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh[ii], fDim=dim(D)[2])$integral
    psi2_fi = cubature::adaptIntegrate(pI3_t2, low2, up2, D=t(D_i[,ii]), theta0=theta0, sig0_sq0=sig0_sq0, Qh=Qh[ii], fDim=dim(D)[2])$integral
    Qh_fi = ((n1/(n*pi0[1])+v0[1])*psi1_fi+(n2/(n*pi0[2])+v0[2])*psi2_fi)
    U = U+(1/sig0_sq0*D_i[,ii] %*% (Y[ii,]-t(D_i[,ii]) %*% theta0)-Qh_fi-q_s*Psi%*%theta0) %*% t(1/sig0_sq0*D_i[,ii] %*% (Y[ii,]-t(D_i[,ii]) %*% theta0)-Qh_fi-q_s*Psi %*% theta0)

  }
  cov_theta = diag(solve(V, U) %*% solve(V)) ^ 0.5
  cov_mtx = solve(V, U) %*% solve(V)

  alpha = eta0[1:(1+degree+nknots),]
  gam = eta0[(2+degree+nknots):(1+degree+nknots+p),]
  std_gam = cov_theta[(2+degree+nknots):(1+degree+nknots+p)]
  cov_mtxa = cov_mtx[1:(1+degree+nknots),1:(1+degree+nknots)]

  return (list(alpha=alpha, gam=gam, std_gam=std_gam, cov_mtxa=cov_mtxa, step=step, pi0=pi0, v0=v0, sig0_sq0=sig0_sq0))

}

