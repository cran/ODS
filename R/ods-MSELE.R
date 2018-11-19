## code written by Mark Weaver ##
## modified by Yinghao Pan ##
## implements the empirical likelihood method described in Zhou et al. (2002) ##

A.matrix <- function(alpha,dl.dpi)
{
  tot <- 1 + sum(exp(alpha))

  A <- matrix(0,ncol=length(alpha),nrow=length(alpha))
  for (i in 1:nrow(A)){
    for (j in 1:ncol(A)){

      d2pi <- rep(0,length(dl.dpi))

      if (i == j){
        for (l in 1:(length(d2pi)-1)){
          if(l== i){
            d2pi[l] <- (exp(alpha[i])*tot-exp(2*alpha[i]))*(tot-2*exp(alpha[i]))/tot^3}
          else{
            d2pi[l] <- (-exp(alpha[l]+alpha[i])*tot+2*exp(2*alpha[i]+alpha[l]))/tot^3}
        }

        d2pi[length(d2pi)] <-
          (-exp(alpha[i]) * tot + 2*exp(2*alpha[i]))/tot^3

      }

      else{
        for (l in 1:(length(d2pi)-1)){
          if(l ==i){
            d2pi[l] <- (2*exp(2*alpha[i]+alpha[j])-exp(alpha[i]+alpha[j])*tot)/tot^3}
          else if(l == j){
            d2pi[l] <- (2*exp(2*alpha[j]+alpha[i])-exp(alpha[i]+alpha[j])*tot)/tot^3}
          else {
            d2pi[l] <- 2*exp(alpha[l] + alpha[i] + alpha[j])/tot^3 }

        }
        d2pi[length(d2pi)] <- 2*exp(alpha[i]+alpha[j])/tot^3
      }

      A[i,j] <- t(as.matrix(d2pi)) %*% as.matrix(dl.dpi)
    }
  }

  A
}


cond.cdf <- function(x, beta, sig, a)
  #Calculates the conditional CDF for all cut points, including +- infinity
  #This is general for any number of strata
{
  a <- c(-Inf, a, Inf)
  nstrata <- length(a) -1
  Z <- matrix(0, nrow=nrow(x), ncol=length(a))

  for (j in 1:length(a))
  {
    Z[,j] <- ((a[j] - x%*%beta)/sqrt(sig))
  }

  probs <- matrix(0,nrow=nrow(x),ncol=nstrata)

  for (j in 1:nstrata)
  {
    probs[,j] <- stats::pnorm(Z[,j+1])-stats::pnorm(Z[,j])
  }
  probs
}


cpdf <- function(y,x,beta,sig)
  #  THE CONDITIONAL PDF
  #  Input:  y response vector, x the covariate matrix (with constant),
  #          current estimate for beta and sigma squared.
  #  Output:  Vector of conditional pdf values.
{
  func <- 1/(sqrt(2*pi*sig)) * exp(-1/(2*sig) * ( y - x%*%beta)^2)
  func
}


d2Pk.db <- function(a,x,beta,sig,strat)
  #This function is required for some that follow.  It gives the differences
  #required for the second derivatives of the conditional strata probs
  #Must multiply by -xtx/sig
{
  a <- c(-Inf,a,Inf)
  Z <- matrix(0,nrow=nrow(x),ncol=length(a))
  for(k in 1:ncol(Z)){
    if (abs(a[k]) != Inf){
      Z[,k] <- (a[k]-x%*%beta)*cpdf(y=a[k],x,beta,sig)}
    else Z[,k] <- 0
  }
  diff <- matrix(0, nrow=nrow(x),ncol=length(a)-1)
  for (k in 1:ncol(diff)){
    diff[,k] <- Z[,k+1]-Z[,k]
  }

  diff.use <- diff[,strat]
  diff.use
}


d2Pk.dbdsig <- function(a,x,beta,sig,strat)
  #  Derivative of the conditional stratum probabilities wrt sigma^2
  #  Actually, only gives the differences.  Must multiply by -x
{
  a <- c(-Inf,a,Inf)
  Z <- matrix(0,nrow=nrow(x),ncol=length(a))
  for(k in 1:ncol(Z)){
    if(abs(a[k]) != Inf){
      Z[,k] <- c(dfdsig(y=a[k],x,beta,sig))}
    else Z[,k] <- 0
  }
  diff <- matrix(0,nrow=nrow(x),ncol=(length(a) - 1))
  for(k in 1:ncol(diff)){
    diff[,k] <- Z[,k+1] - Z[,k]}

  diff.use <- diff[,strat]
  diff.use
}


d2Pk.dsig <- function(a,x,beta,sig,strat)
  #  Derivative of the conditional stratum probabilities wrt sigma^2
{
  a <- c(-Inf,a,Inf)
  Z <- matrix(0,nrow=nrow(x),ncol=length(a))
  for(k in 1:ncol(Z)){
    if(abs(a[k]) != Inf){
      Z[,k] <- -(a[k]-x%*%beta)*
        (cpdf(y=a[k],x,beta,sig)*(-1)/(2*sig^2)+dfdsig(y=a[k],x,beta,sig)*1/(2*sig))
    }
    else Z[,k] <- 0
  }
  diff <- matrix(0,nrow=nrow(x),ncol=(length(a) - 1))
  for(k in 1:ncol(diff)){
    diff[,k] <- Z[,k+1] - Z[,k]}

  diff.use <- diff[,strat]
  diff.use
}


d2f.dbdsig <- function(y,x,beta,sig)
  #  Second deriv of cpdf wrt beta and sigma^2
{
  pdf <- cpdf(y,x,beta,sig)
  df <- dfdsig(y,x,beta,sig)
  temp <- (y-x%*% beta)
  func <- -x * c(temp*pdf/(sig^2)) + x*c(temp*df/sig)

  func
}


d2f.dsig <- function(y,x,beta,sig)
  #  Second derivative of cpdf wrt sigma ^2
{
  pdf <- cpdf(y,x,beta,sig)
  temp <- (y - x%*% beta)^2
  func <- c(pdf) * ((1/(2*sig^2)-temp/(sig^3)) +
                      (-1/(2*sig) + temp/(2*sig^2))^2)

  func
}


d2lc.dbdb <- function(x,beta,sig,a,N.edf,rhos,strat)
  #  The second derivative matrix with respect to beta
{
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  temp3 <- dpdb(x,beta,sig,a,N.edf,rhos,strat)
  func <- -crossprod(x)/sig + d2pdb(x,beta,sig,a,N.edf,rhos,strat) -
    crossprod(temp3/c(p.i^2),temp3)

  func
}

d2lc.dbdpis <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
{
  temp1 <- d2p.dbdpi(x,beta,sig,pis,a,N.edf,rhos,strat)
  temp2 <- dpdpis(x,beta,sig,pis,a,N.edf,rhos,strat)
  temp3 <- dpdb(x,beta,sig,a,N.edf,rhos,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  temp4 <- temp2/c(p.i^2)
  temp5 <- crossprod(temp3,temp4)
  func <- temp1 - temp5
  func
}


d2lc.dbdsig <- function(y, x,beta,sig,a,N.edf,rhos,strat)
  #  Second derivative of ods likelihood wrt beta and sigma^2
{
  d2f <- d2f.dbdsig(y,x,beta,sig)
  df <- dfdsig(y,x,beta,sig)
  pdf <- cpdf(y,x,beta,sig)
  df.db <- dfdb(y,x,beta,sig)
  d2p <- d2p.dbdsig(x,beta,sig,a,N.edf,rhos,strat)
  dp <- dpdsig(x,beta,sig,a,N.edf,rhos,strat)
  dp.db <- dpdb(x,beta,sig,a,N.edf,rhos,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  first <- (d2f*c(pdf) - df.db*c(df))/c(pdf^2)
  first <- apply(first,2,sum)
  second <- (d2p*c(p.i) - dp.db*c(dp))/c(p.i^2)
  second <- apply(second,2,sum)
  func <- first + second

  func
}


d2lc.dpidpi<- function(x,beta,sig,pis,a,N.edf,rhos,strat)
{
  temp1 <- d2p.dpis(x,beta,sig,pis,a,N.edf,rhos,strat)
  temp1 <- as.matrix(temp1)
  temp2 <- dpdpis(x,beta,sig,pis,a,N.edf,rhos,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  temp3 <- temp2/c(p.i^2)
  temp4 <- crossprod(temp2,temp3)

  rhos <- rhos[strat]/pis[strat]
  func <- temp1 - temp4  + diag(rhos, nrow=length(rhos))
  func
}


d2lc.dpisdsig <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
  #  Second deriv of ods likelihood wrt sigma^2
{
  d2p <- d2p.dpisdsig(x,beta,sig,pis,a,N.edf,rhos,strat)
  dp <- dpdpis(x,beta,sig,pis,a,N.edf,rhos,strat)
  dp.ds <- dpdsig(x,beta,sig,a,N.edf,rhos,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  func <- (d2p * c( p.i) - dp * c(dp.ds))/c(p.i^2)
  func <- as.matrix(func)
  func <- apply(func,2,sum)

  func
}


d2lc.dsig <- function(y, x,beta,sig,a,N.edf,rhos,strat)
  #  Second derivative of ods likelihood wrt sigma ^2
{
  d2f <- d2f.dsig(y,x,beta,sig)
  df <- dfdsig(y,x,beta,sig)
  pdf <- cpdf(y,x,beta,sig)
  dp <- dpdsig(x,beta,sig,a,N.edf,rhos,strat)
  d2p <- d2p.dsig(x,beta,sig,a,N.edf,rhos,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  first <- (c(d2f) * c(pdf)- c(df^2))/c(pdf^2)
  second <- (c(d2p)*c(p.i) - c(dp^2))/c(p.i^2)
  func <- sum(first) + sum(second)

  func
}


d2p.dbdpi <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
  #The second derivative matrix of p.i wrt to beta and pis
{
  dp.db <- dpdb(x,beta,sig,a,N.edf,rhos,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat]/pis[strat])
  dPk <- as.matrix(dPk.db(a,x,beta,sig,strat))
  probs <- cond.cdf(x,beta,sig,a)
  probs <- as.matrix(probs[,strat])

  d2p <- matrix(0,nrow=length(beta),ncol=length(strat))
  for(k in 1:length(rhos)){
    func <- rhos[k]*c(p.i)*dPk[,k]*(-x) + 2*rhos[k]*c(probs[,k])*dp.db
    func <- as.matrix(func)
    d2p[,k] <- apply(func,2,sum)
  }
  d2p
}


d2p.dbdsig <- function(x,beta,sig,a,N.edf,rhos,strat)
  #  second deriv of NPMLE wrt to beta and sigma^2
{
  dPkdb <- dPk.db(a,x,beta,sig,strat)
  dPkds <- dPk.dsig(a,x,beta,sig,strat)
  d2Pk <- d2Pk.dbdsig(a,x,beta,sig,strat)
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat])

  func <- -x*
    (2*c(p.i^3)*c(dPkds%*%(-rhos))*c(dPkdb%*%(-rhos))+c(p.i^2)*c(d2Pk%*%(-rhos)))

  func
}


d2p.dpis <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
  #This function computes the second derivative matrix for p.is wrt pis
  #It is general and can be used for any number of strata
{
  probs <- cond.cdf(x,beta,sig,a)
  probs <- as.matrix(probs[,strat])
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat]/pis[strat])
  d2p <- matrix(0,nrow=length(rhos),ncol=length(rhos))
  for(i in 1:nrow(d2p)){
    for(j in 1:ncol(d2p)){
      if(i==j){
        d2p[i,j] <-
          sum(probs[,i]*(-2*rhos[i]/pis[strat[i]]*p.i+2*(rhos[i]^2)*probs[,i]*(p.i^2)))}
      else d2p[i,j] <- sum(2*rhos[i]*rhos[j]*c(probs[,i])*c(probs[,j])*c(p.i^2))
    }
  }
  d2p
}


d2p.dpisdsig <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
  #  Second deriv of NPMLE wrt pis and sigma^2
{
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat])
  probs <- cond.cdf(x,beta,sig,a)
  temp <- as.matrix(probs[,strat])
  dPk <- as.matrix(dPk.dsig(a,x,beta,sig,strat))
  for (k in 1:ncol(temp)){
    temp[,k] <- (rhos[k]/pis[strat[k]])*dPk[,k]*c(p.i^2) +
      2*(rhos[k]/pis[strat[k]])*c(temp[,k])*c(p.i^3)*c(dPk%*%(-rhos))
  }

  temp
}


d2p.dsig <- function(x,beta,sig,a,N.edf,rhos,strat)
  #  Second derivative of NPMLE wrt sigma^2
{
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat])
  dPk <- dPk.dsig(a,x,beta,sig,strat)
  d2Pk <- d2Pk.dsig(a,x,beta,sig,strat)
  func <- 2*c(p.i^3)*(dPk%*%(-rhos))^2 + c(p.i^2)*(d2Pk%*%(-rhos))

  func
}


d2pdb <- function(x,beta,sig,a,N.edf,rhos,strat)
  #This function gives the second derivative matrix for the p.i's wrt beta
  #Actually gives second derivatives with p.is divided out as required
  #This function is general
{

  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)

  x1.prime <- x*c(p.i) *
    c((d2Pk.db(a,x,beta,sig,strat)%*%as.matrix(rhos[strat]))/sig)
  x2.prime <- x*c(p.i^2) *
    c(2*(dPk.db(a,x,beta,sig,strat)%*%as.matrix(rhos[strat]))^2)

  func <- crossprod(x1.prime,x) + crossprod(x2.prime,x)

  func
}


dPk.db <- function(a,x,beta,sig,strat)
  #This function provides the matrix of differences between the cpdf in
  #the intervals.  It is necessary for functions which follow
  #Note that this is not the entire dPk.db.  Must multiply by -x!!
{
  a <- c(-Inf,a,Inf)
  Z <- matrix(0,nrow=nrow(x),ncol=length(a))
  for(k in 1:ncol(Z)){
    Z[,k] <- cpdf(y=a[k],x,beta,sig)
  }
  diff <- matrix(0, nrow=nrow(x),ncol=length(a)-1)
  for (k in 1:ncol(diff)){
    diff[,k] <- Z[,k+1]-Z[,k]
  }

  diff.use <- diff[,strat]

  diff.use
}


dPk.dsig <- function(a,x,beta,sig,strat)
  #  Derivative of the conditional stratum probabilities wrt sigma^2
{
  a <- c(-Inf,a,Inf)
  Z <- matrix(0,nrow=nrow(x),ncol=length(a))
  for(k in 1:ncol(Z)){
    if(abs(a[k]) != Inf){
      Z[,k] <- -(a[k]-x%*%beta)*cpdf(y=a[k],x,beta,sig)*1/(2*sig)
    }
    else Z[,k] <- 0
  }
  diff <- matrix(0,nrow=nrow(x),ncol=(length(a) - 1))
  for(k in 1:ncol(diff)){
    diff[,k] <- Z[,k+1] -Z[,k]}

  diff.use <- diff[,strat]
  diff.use
}


dfdb <- function(y,x,beta,sig)
  #The derivative of the conditional pdf wrt beta
  #General for any number of strata (does not depend on strata)
{

  func <- c(cpdf(y,x,beta,sig))* x* c((y-x%*%beta)/sig)

  func
}


dfdsig <- function(y,x,beta,sig)
  # The derivative of the cpdf wrt sigma2
{
  pdf <- cpdf(y,x,beta,sig)
  temp <- (y - x%*% beta)^2
  func <- c(-1/(2*sig) + temp/(2*sig^2)) * pdf
  func
}


dlc.db <- function(y,x,beta,sig,a,N.edf,rhos,strat)
  #This function gets the part of the gradient coresponding to the betas
  #It is general, as long as the component functions are general!

{
  temp2 <- dfdb(y,x,beta,sig)/c(cpdf(y,x,beta,sig))
  func1 <- apply(temp2,2,sum)

  func2 <- dpdb(x,beta,sig,a,N.edf,rhos,strat)/
    c( npmle(x,beta,sig,a,N.edf,rhos,strat))
  func2 <- apply(func2,2,sum)

  func <- func1 + func2
  func
}


dlc.dpis <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
  # Computes the part of the gradient corresponding to the pis
  # This is general, as long as the components are all general!
{
  func <- dpdpis(x,beta,sig,pis,a,N.edf,rhos,strat)/
    c(npmle(x,beta,sig,a,N.edf,rhos,strat))
  func <- as.matrix(func)
  func <- apply(func,2,sum)
  func <- func - rhos[strat]
  func
}


dlc.dsig <- function(y, x,beta,sig,a,N.edf,rhos,strat)
  # Gradient coresponding to sigma^2
{
  first <- dfdsig(y,x,beta,sig)/c(cpdf(y,x,beta,sig))
  first <- sum(first)
  second <- dpdsig(x,beta,sig,a,N.edf,rhos,strat)/
    c(npmle(x,beta,sig,a,N.edf,rhos,strat))
  second <- sum(second)
  func <- first + second

  func
}


dpdb <- function(x,beta,sig,a,N.edf,rhos,strat)
  #This function gives the first derivative of p.i wr to beta
  #This is general for any number of strata
{
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat])
  func <- -(-x * c(dPk.db(a,x,beta,sig,strat)%*%(rhos))*c(p.i^2))

  func
}

dpdpis <- function(x,beta,sig,pis,a,N.edf,rhos,strat)
  #This function gives the first derivative of p.i wr to the pis
  #This is now general for any number of strata
{
  temp2 <- rhos[strat] / pis[strat]
  temp3 <- cond.cdf(x,beta,sig,a)
  temp3 <- as.matrix(temp3[,strat])

  for(k in 1:ncol(temp3)){
    temp3[,k] <- temp3[,k]*temp2[k]
  }
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  func <- temp3 * c(p.i)^2

  func
}


dpdsig <- function(x,beta,sig,a,N.edf,rhos,strat)
  #  Derivative of NPMLE wrt sigma^2
{
  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)
  rhos <- as.matrix(rhos[strat])
  func <- c(p.i^2) * (dPk.dsig(a,x,beta,sig,strat) %*% (-rhos))

  func
}


dpis.dalpha <- function (alpha)
  # First, the derivatives of the pis wrt the alphas, the new parms.
{
  den <- (1 + sum(exp(alpha)))^2

  func <- matrix(0,nrow=(length(alpha)+1),ncol=length(alpha))

  for(i in 1:ncol(func)){
    for (j in 1:ncol(func)){
      if (i == j){
        func[i,j] <- (exp(alpha[i])*(1+sum(exp(alpha))) - exp(2*alpha[i]))/ den
      }
      else{
        func[i,j] <- -exp(alpha[i] + alpha[j]) / den }
    }
  }

  for(i in 1:ncol(func)){
    func[nrow(func),i] <- -exp(alpha[i])/den
  }

  func
}


grad.lc <- function(y,x,beta,sig,pis,a,N.edf,rhos,strat)
  #  The Gradient function
{
  temp1 <- dlc.db(y,x,beta,sig,a,N.edf,rhos,strat)
  temp2 <- dlc.dpis(x,beta,sig,pis,a,N.edf,rhos,strat)
  temp3 <- dlc.dsig(y,x,beta,sig,a,N.edf,rhos,strat)
  func <- as.matrix(c(temp1,temp3,temp2))
  func
}


hessian.lc <- function(y,x,beta,sig,pis,a,N.edf,rhos,strat)
  #The Hessian matrix.  this is a general function
{
  cols <- length(beta)+length(pis[strat])+1
  b <- length(beta)
  b1 <- length(beta)+1
  b2 <- b1+1
  hess <- matrix(0,ncol=cols,nrow=cols)
  hess[1:b,1:b]<- d2lc.dbdb(x,beta,sig,a,N.edf,rhos,strat)
  hess[1:b,b1] <- d2lc.dbdsig(y,x,beta,sig,a,N.edf,rhos,strat)
  hess[1:b,b2:cols]<- d2lc.dbdpis(x,beta,sig,pis,a,N.edf,rhos,strat)
  hess[b1,1:b] <- t(hess[1:b,b1])
  hess[b1,b1] <- d2lc.dsig(y,x,beta,sig,a,N.edf,rhos,strat)
  hess[b1,b2:cols] <- d2lc.dpisdsig(x,beta,sig,pis,a,N.edf,rhos,strat)
  hess[b2:cols,1:b]<- t(hess[1:b,b2:cols])
  hess[b2:cols,b1] <- t(hess[b1,b2:cols])
  hess[b2:cols,b2:cols]<- d2lc.dpidpi(x,beta,sig,pis,a,N.edf,rhos,strat)

  hess
}


npmle <- function(x,beta,sig,a,N.edf,rhos,strat)
  #The npmle of G_x.  General for any number of strata
{
  cdf <- cond.cdf(x,beta,sig,a)     #Matrix of conditional interval probs
  cdf <- as.matrix(cdf[,strat])

  rhos <- as.matrix(rhos[strat])
  p.i <- 1/(N.edf + cdf%*%rhos)

  p.i
}


#Main call routines
#Main call routines

#'MSELE estimator for analyzing the primary outcome in ODS design
#'
#'\code{odsmle} provides a maximum semiparametric empirical likelihood estimator
#'(MSELE) for analyzing the primary outcome Y with respect to expensive exposure
#'and other covariates in ODS design (Zhou et al. 2002).
#'
#'We assume that in the population, the primary outcome variable Y follows the
#'following model: \deqn{Y=\beta_{0}+\beta_{1}X+\epsilon,}{Y = beta0 + beta1*X +
#'epsilon,} where X are the covariates, and epsilon has variance sig. In ODS
#'design, a simple random sample is taken from the full cohort, then two
#'supplemental samples are taken from two tails of Y, i.e. (-Infty, mu_Y -
#'a*sig_Y) and (mu_Y + a*sig_Y, +Infty). Because ODS data has biased sampling
#'nature, naive regression analysis will yield biased estimates of the
#'population parameters. Zhou et al. (2002) describes a semiparametric empirical
#'likelihood estimator for estimating the parameters in the primary outcome
#'model.
#'
#'@export
#'@param Y vector for the primary response
#'@param X the design matrix with a column of 1's for the intercept
#'@param beta starting parameter values for the regression coefficients that
#'  relate Y to X.
#'@param sig starting parameter values for the error variance of the regression.
#'@param pis starting parameter values for the stratum probabilities (the
#'  probability that Y belongs to certain stratum) e.g. pis = c(0.1, 0.8, 0.1).
#'@param a vector of cutpoints for the primary response (e.g., a = c(-2.5,2))
#'@param rs.size size of the SRS (simple random sample)
#'@param size vector of the stratum sizes of the supplemental samples (e.g. size
#'  = c(50,0,50) represents that two supplemental samples each of size 50 are
#'  taken from the upper and lower tail of Y.)
#'@param strat vector that indicates the stratum numbers (e.g. strat = c(1,2,3)
#'  represents that there are three stratums).
#'@return A list which contains the parameter estimates for the primary response
#'  model: \deqn{Y=\beta_{0}+\beta_{1}X+\epsilon,}{Y = beta0 + beta1*X +
#'  epsilon,} where epsilon has variance sig. The list contains the following
#'  components: \item{beta}{parameter estimates for beta} \item{sig}{estimates
#'  for sig} \item{pis}{estimates for the stratum probabilities}
#'  \item{grad}{gradient} \item{hess}{hessian} \item{converge}{whether the
#'  algorithm converges: True or False} \item{i}{Number of iterations}
#' @examples
#'library(ODS)
#' # take the example data from the ODS package
#' # please see the documentation for details about the data set ods_data
#'
#' Y <- ods_data[,1]
#' X <- cbind(rep(1,length(Y)), ods_data[,2:5])
#'
#' # use the simple random sample to get an initial estimate of beta, sig #
#' # perform an ordinary least squares #
#' SRS <- ods_data[1:200,]
#' OLS.srs <- lm(SRS[,1] ~ SRS[,2:5])
#' OLS.srs.summary <- summary(OLS.srs)
#'
#' beta <- coefficients(OLS.srs)
#' sig <- OLS.srs.summary$sigma^2
#' pis <- c(0.1,0.8,0.1)
#'
#' # the cut points for this data is Y < 0.162, Y > 2.59.
#' a <- c(0.162,2.59)
#' rs.size <- 200
#' size <- c(100,0,100)
#' strat <- c(1,2,3)
#'
#' odsmle(Y,X,beta,sig,pis,a,rs.size,size,strat)

odsmle <- function(Y,X,beta,sig,pis,a,rs.size,size,strat)
  # This function get the semiparametric MLE given an ODS sample has already
  # been obtained.

{

  #  strat <- c(1,3)

  b <- length(beta)
  b1 <- b+1
  b2 <- b + 2


  k1 <- length(pis) - 1
  alpha <- rep(0,k1)
  den.alpha <- 1 - sum(pis[1:k1])
  for(i in 1:k1)
  {
    alpha[i] <- log(pis[i]/den.alpha)
  }

  #  alpha <- rep(0,2)
  #  alpha[1] <- log(pis[1]/(1-pis[1]-pis[3]))
  #  alpha[2] <- log(pis[3]/(1-pis[1]-pis[3]))

  grad <- grad.lc(Y,X,beta,sig,pis,a,
                  N.edf=rs.size,rhos=size/pis,strat)
  hess <- hessian.lc(Y,X,beta,sig,pis,a,
                     N.edf=rs.size,rhos=size/pis,strat)

  dpi.da <- dpis.dalpha(alpha)
  temp.grad <- t(dpi.da) %*% as.matrix(grad[b2:length(grad)])
  A <- A.matrix(alpha,grad[b2:length(grad)])
  temp.hess <-t(dpi.da)%*%hess[b2:length(grad),b2:length(grad)]%*%dpi.da
  temp2.hess <- hess[1:b1,b2:length(grad)]%*%dpi.da
  grad <- c(grad[1:b1],temp.grad)
  hess <- rbind(cbind(hess[1:b1,1:b1],temp2.hess),
                cbind(t(temp2.hess),(temp.hess+A)))




  eta <- c(beta,sig,alpha)

  #eta <- c(beta,sig,pis[strat])

  converge <- F

  det <- prod(svd(hess,nu=0,nv=0)$d)
  if(det <= 1.0e-06){
    return(list(beta=beta,sig=sig,pis=pis,grad=grad,hess=hess,converge=converge))
  }

  delta <- max(abs(grad))
  if ( delta < 1.0e-04)
  {
    converge <- T
    return(list(beta=beta,sig=sig,pis=pis,grad=grad,hess=hess,converge=converge))
  }

  for(i in 1:20)
  {
    eta.new <- eta - solve(hess)%*%grad
    gamma <- eta - eta.new
    eta <- eta.new
    beta <- eta[1:length(beta)]
    sig <- eta[(length(beta)+1)]
    #  pis[strat] <- eta[(length(beta)+2):length(eta)]
    #  pis[2] <- 1-pis[1]-pis[3]
    alpha <- eta[(b1+1):length(eta)]
    for(j in 1:k1){
      pis[j] <- exp(alpha[j])/(1 + sum(exp(alpha)))
    }
    pis[length(pis)] <- 1 - sum(pis[1:k1])

    # pis[1] <- exp(alpha[1])/(1+exp(alpha[1])+exp(alpha[2]))
    #  pis[3] <- exp(alpha[2])/(1+exp(alpha[1])+exp(alpha[2]))
    #  pis[2] <- 1-pis[1]-pis[3]
    grad <- grad.lc(Y,X,beta,sig,pis,a,
                    N.edf=rs.size,rhos=size/pis,strat)
    hess <- hessian.lc(Y,X,beta,sig,pis,a,
                       N.edf=rs.size,rhos=size/pis,strat)

    dpi.da <- dpis.dalpha(alpha)
    temp.grad <- t(dpi.da) %*% as.matrix(grad[b2:length(grad)])
    A <- A.matrix(alpha,grad[b2:length(grad)])
    temp.hess <-t(dpi.da)%*%hess[b2:length(grad),b2:length(grad)]%*%dpi.da
    temp2.hess <- hess[1:b1,b2:length(grad)]%*%dpi.da
    grad <- c(grad[1:b1],temp.grad)
    hess <- rbind(cbind(hess[1:b1,1:b1],temp2.hess),
                  cbind(t(temp2.hess),(temp.hess+A)))


    det <- prod(svd(hess,nu=0,nv=0)$d)
    if(det <= 1.0e-06){
      return(list(beta=beta,sig=sig,pis=pis,grad=grad,hess=hess,converge=converge))
    }

    delta <- max(abs(grad))
    gamma <- max(abs(gamma))
    if ( delta < 1.0e-04 || gamma < 1.0e-04)
    {
      converge <- T
      return(list(beta=beta,sig=sig,pis=pis,grad=grad,hess=hess,converge=converge,i=i))
    }
  }


  list(beta=beta,sig=sig,pis=pis,grad=grad,hess=hess,converge=converge)
}

#' standard error for MSELE estimator
#'
#' \code{se.spmle} calculates the standard error for MSELE estimator in Zhou et
#' al. 2002
#'
#' @export
#' @param y vector of the primary response
#' @param x the design matrix with a column of 1's for the intercept
#' @param beta final estimates of the regression coefficients obtained from
#'   odsmle
#' @param sig final estimate of the error variance obtained from odsmle
#' @param pis final estimates of the stratum probabilities obtained from odsmle
#' @param a vector of cutpoints for the primary response (e.g., a = c(-2.5,2))
#' @param N.edf should be the size of the SRS (simple random sample)
#' @param rhos which is size/pis, where size is a vector representing the
#'   stratum sizes of supplemental samples. e.g. size = c(100, 0, 100), and pis
#'   are the final estimates obtained from odsmle.
#' @param strat vector that indicates the stratum numbers of supplemental
#'   samples, except that you should only list stratum with size > 0. (e.g. if
#'   the supplemental size is c(100, 0, 100), then the strat vector should be
#'   c(1,3))
#' @param size.nc total size of the validation sample (SRS plus supplemental
#'   samples)
#' @return A list which contains the standard error estimates for betas in the
#'   model : \deqn{Y=\beta_{0}+\beta_{1}X+\epsilon,}{Y = beta0 + beta1*X +
#'   epsilon,} where epsilon has variance sig.
#' @examples
#' library(ODS)
#' # take the example data from the ODS package
#' # please see the documentation for details about the data set ods_data
#'
#' Y <- ods_data[,1]
#' X <- cbind(rep(1,length(Y)), ods_data[,2:5])
#'
#' # use the simple random sample to get an initial estimate of beta, sig #
#' # perform an ordinary least squares #
#' SRS <- ods_data[1:200,]
#' OLS.srs <- lm(SRS[,1] ~ SRS[,2:5])
#' OLS.srs.summary <- summary(OLS.srs)
#'
#' beta <- coefficients(OLS.srs)
#' sig <- OLS.srs.summary$sigma^2
#' pis <- c(0.1,0.8,0.1)
#'
#' # the cut points for this data is Y < 0.162, Y > 2.59.
#' a <- c(0.162,2.59)
#' rs.size <- 200
#' size <- c(100,0,100)
#' strat <- c(1,2,3)
#'
#' # obtain the parameter estimates
#' ODS.model = odsmle(Y,X,beta,sig,pis,a,rs.size,size,strat)
#'
#' # calculate the standard error estimate
#' y <- Y
#' x <- X
#' beta <- ODS.model$beta
#' sig <- ODS.model$sig
#' pis <- ODS.model$pis
#' a <- c(0.162,2.59)
#' N.edf <- rs.size
#' rhos <- size/pis
#' strat <- c(1,3)
#' size.nc <- length(y)
#'
#' se = se.spmle(y, x, beta, sig, pis, a, N.edf, rhos, strat, size.nc)
#'
#' # summarize the result
#' ODS.tvalue <- ODS.model$beta / se
#' ODS.pvalue <- 2 * pt( - abs(ODS.tvalue), sum(rs.size, size)-2)
#'
#' ODS.results <- cbind(ODS.model$beta, se, ODS.tvalue, ODS.pvalue)
#' dimnames(ODS.results)[[2]] <- c("Beta","SEbeta","tvalue","Pr(>|t|)")
#' row.names(ODS.results) <- c("(Intercept)","X","Z1","Z2","Z3")
#'
#' ODS.results

se.spmle <- function(y,x,beta,sig,pis,a,N.edf,rhos,strat,size.nc)
{
  #  Computes the standard error estimator

  p.i <- npmle(x,beta,sig,a,N.edf,rhos,strat)

  temp1 <- dfdb(y,x,beta,sig)/c(cpdf(y,x,beta,sig))
  temp1a <- dfdsig(y,x,beta,sig)/c(cpdf(y,x,beta,sig))
  temp2 <- dpdb(x,beta,sig,a,N.edf,rhos,strat)/c(p.i)
  temp2a <- dpdsig(x,beta,sig,a,N.edf,rhos,strat)/c(p.i)
  dlcdb <- temp1 + temp2
  dlcds <- temp1a + temp2a

  temp3 <- dpdpis(x,beta,sig,pis,a,N.edf,rhos,strat)/c(p.i)
  dlcdpis <- temp3 - rhos[strat]/size.nc

  dlc <- cbind(dlcdb,dlcds,dlcdpis)

  U <- crossprod(dlc)/size.nc

  A <- -(hessian.lc(y,x,beta,sig,pis,a,N.edf,rhos,strat))/size.nc

  A.inv <- solve(A)

  var <- (A.inv %*% U %*% A.inv)/size.nc

  se <- sqrt(diag(var))

  se.beta <- se[1:ncol(x)]

  se.beta
}

#' Data example for analyzing the primary response in ODS design
#'
#' Data example for analyzing the primary response in ODS design (zhou et al.
#' 2002)
#'
#' @format A matrix with 400 rows and 5 columns. The first 200 observations are
#'   from the simple random sample, while 2 supplemental samples each with size
#'   100 are taken from one standard deviation above the mean and below the
#'   mean, i.e. (Y1 < 0.162) and (Y1 > 2.59). \describe{ \item{Y1}{primary
#'   outcome for which the ODS sampling scheme is based on} \item{X}{expensive
#'   exposure} \item{Z1}{a simulated covariate} \item{Z2}{a simulated covariate}
#'   \item{Z3}{a simulated covariate} }
#'
#' @source A simulated data set
"ods_data"
