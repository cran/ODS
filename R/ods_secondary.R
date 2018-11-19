# code written by Yinghao Pan

# July 7, 2016
# Project2: secondary outcome analysis in continuous ODS design
# use the dimension of Z as a parameter in all functions
# Method2: Using data from all observations
# augmented inverse probability weighted estimator
# Y1: primary outcome
# Y2: secondary outcome
# X: expensive exposure
# Z: other covariates, has more than 1 dimensions

# repeat a vector by n times --------------------------------------
# suppose a vector x has length p, the resulting matrix has dim n*p
# each row is x

rep_row <- function(x,n) {
  matrix(rep(x, each=n), nrow=n)
}

# score function --------------------------------------------------

secondary_ODS_score <- function(y1,y2,x,z,beta0,beta1,beta2,gamma0,gamma1,gamma2,Q11,Q12,Q21,Q22,weight,r,phi0,phi1,phi2,phi3,varestimate,Z.dim) {

  D <- matrix(c(1,x,z,rep(0, 2*(Z.dim+2)),1,x,z),nrow=2,byrow=TRUE)
  Q <- matrix(c(Q11,Q12,Q21,Q22),nrow=2)
  Qinverse <- solve(Q)
  Qinverse11 <- Qinverse[1,1]
  Qinverse12 <- Qinverse[1,2]
  Qinverse21 <- Qinverse[2,1]
  Qinverse22 <- Qinverse[2,2]
  residual1 <- as.vector(r/weight*(y1-beta0-beta1*x-z%*%beta2))
  residual2 <- as.vector(r/weight*(y2-gamma0-gamma1*x-z%*%gamma2))
  residual <- matrix(c(residual1,residual2),nrow=2)
  if (r==0) {
    matrix1 <- matrix(0, nrow=2*(2+Z.dim))
  } else {
    matrix1 <- t(D)%*%solve(Q)%*%residual
  }

  expectationX <- as.vector(phi0+phi1*y1+phi2*y2+z%*%phi3)
  expectationXsquare <- varestimate + expectationX^2
  expectationU1 <- Qinverse11*(y1-beta0-beta1*expectationX-z%*%beta2)+Qinverse12*(y2-gamma0-gamma1*expectationX-z%*%gamma2)
  expectationU2 <- Qinverse21*(y1-beta0-beta1*expectationX-z%*%beta2)+Qinverse22*(y2-gamma0-gamma1*expectationX-z%*%gamma2)
  expectationxU1 <- (Qinverse11*(y1-beta0-z%*%beta2)+Qinverse12*(y2-gamma0-z%*%gamma2))*expectationX-(beta1*Qinverse11+gamma1*Qinverse12)*expectationXsquare
  expectationxU2 <- (Qinverse21*(y1-beta0-z%*%beta2)+Qinverse22*(y2-gamma0-z%*%gamma2))*expectationX-(beta1*Qinverse21+gamma1*Qinverse22)*expectationXsquare

  vector <- (1-r/weight)*c(expectationU1, expectationxU1, z*as.vector(expectationU1), expectationU2, expectationxU2, z*as.vector(expectationU2))
  matrix2 <- matrix(vector,nrow=4+2*Z.dim)

  answer <- matrix1 + matrix2
  return(answer)
}


# vectorized score function ---------------------------------------

secondary_ODS_score_vectorized <- function(x) {

  Z.dim <- x[length(x)]
  Y1 <- x[1]
  Y2 <- x[2]
  X <- x[3]
  Z <- x[4:(3+Z.dim)]
  beta0 <- x[(4+Z.dim)]
  beta1 <- x[(5+Z.dim)]
  beta2 <- x[(6+Z.dim):(5+2*Z.dim)]
  gamma0 <- x[(6+2*Z.dim)]
  gamma1 <- x[(7+2*Z.dim)]
  gamma2 <- x[(8+2*Z.dim):(7+3*Z.dim)]
  Q11 <- x[(8+3*Z.dim)]
  Q12 <- x[(9+3*Z.dim)]
  Q21 <- x[(10+3*Z.dim)]
  Q22 <- x[(11+3*Z.dim)]
  weight <- x[(12+3*Z.dim)]
  r <- x[(13+3*Z.dim)]
  phi0 <- x[(14+3*Z.dim)]
  phi1 <- x[(15+3*Z.dim)]
  phi2 <- x[(16+3*Z.dim)]
  phi3 <- x[(17+3*Z.dim):(16+4*Z.dim)]
  varestimate <- x[(17+4*Z.dim)]

  return(secondary_ODS_score(Y1,Y2,X,Z,beta0,beta1,beta2,gamma0,gamma1,gamma2,Q11,Q12,Q21,Q22,weight,r,phi0,phi1,phi2,phi3,varestimate,Z.dim))
}

# hessian matrix --------------------------------------------------

secondary_ODS_hessian <- function(y1,y2,x,z,beta0,beta1,beta2,gamma0,gamma1,gamma2,Q11,Q12,Q21,Q22,weight,r,phi0,phi1,phi2,phi3,varestimate,Z.dim) {

  D <- matrix(c(1,x,z,rep(0, 2*(Z.dim+2)),1,x,z),nrow=2,byrow=TRUE)
  Q <- matrix(c(Q11,Q12,Q21,Q22),nrow=2)
  Qinverse <- solve(Q)
  Qinverse11 <- Qinverse[1,1]
  Qinverse12 <- Qinverse[1,2]
  Qinverse21 <- Qinverse[2,1]
  Qinverse22 <- Qinverse[2,2]
  calmatrix <- -Qinverse*r/weight
  if (r==0) {
    matrix1 <- matrix(0, nrow=4+2*Z.dim, ncol=4+2*Z.dim)
  } else {
    matrix1 <- t(D)%*%calmatrix%*%D
  }

  expectationX <- phi0+phi1*y1+phi2*y2+z%*%phi3
  expectationXsquare <- varestimate + expectationX^2

  kroneckerA <- -matrix(c(Qinverse11, Qinverse12, Qinverse21, Qinverse22), nrow=2, byrow=TRUE)
  vector <- matrix(c(1, expectationX, z), nrow=2+Z.dim)
  kroneckerB <- vector %*% t(vector)
  kroneckerB[2,2] <- expectationXsquare

  matrix2 <- (1-r/weight)*kronecker(kroneckerA, kroneckerB)

  answer <- matrix1 + matrix2
  return(answer)
}

# vectorized hessian function ------------------------------------

secondary_ODS_hessian_vectorized <- function(x) {

  Z.dim <- x[length(x)]
  Y1 <- x[1]
  Y2 <- x[2]
  X <- x[3]
  Z <- x[4:(3+Z.dim)]
  beta0 <- x[(4+Z.dim)]
  beta1 <- x[(5+Z.dim)]
  beta2 <- x[(6+Z.dim):(5+2*Z.dim)]
  gamma0 <- x[(6+2*Z.dim)]
  gamma1 <- x[(7+2*Z.dim)]
  gamma2 <- x[(8+2*Z.dim):(7+3*Z.dim)]
  Q11 <- x[(8+3*Z.dim)]
  Q12 <- x[(9+3*Z.dim)]
  Q21 <- x[(10+3*Z.dim)]
  Q22 <- x[(11+3*Z.dim)]
  weight <- x[(12+3*Z.dim)]
  r <- x[(13+3*Z.dim)]
  phi0 <- x[(14+3*Z.dim)]
  phi1 <- x[(15+3*Z.dim)]
  phi2 <- x[(16+3*Z.dim)]
  phi3 <- x[(17+3*Z.dim):(16+4*Z.dim)]
  varestimate <- x[(17+4*Z.dim)]

  return(secondary_ODS_hessian(Y1,Y2,X,Z,beta0,beta1,beta2,gamma0,gamma1,gamma2,Q11,Q12,Q21,Q22,weight,r,phi0,phi1,phi2,phi3,varestimate,Z.dim))
}

#' Secondary analysis in ODS design
#'
#' \code{secondary_ODS} performs the secondary analysis which describes the
#' association between a continuous scale secondary outcome and the expensive
#' exposure for data obtained with ODS (outcome dependent sampling) design.
#'
#' @export
#' @param SRS A data matrix for subjects in the simple random sample. The first
#'   column is Y1: the primary outcome for which the ODS scheme is based on. The
#'   second column is Y2: a secondary outcome. The third column is X: the
#'   expensive exposure. Starting from the fourth column to the end is Z: other
#'   covariates.
#' @param lowerODS A data matrix for supplemental samples taken from the lower
#'   tail of Y1 (eg. Y1 < a). The data structure is the same as SRS.
#' @param upperODS A data matrix for supplemental samples taken from the upper
#'   tail of Y1 (eg. Y1 > b). The data structure is the same as SRS.
#' @param NVsample A data matrix for subjects in the non-validation sample.
#'   Subjects in the non-validation sample don't have the expensive exposure X
#'   measured. The data structure is the same as SRS, but the third column
#'   (which represents X) has values NA.
#' @param cutpoint A vector of length two that represents the cut off points for
#'   the ODS design. eg. cutpoint <- c(a,b). In the ODS design, a simple random
#'   sample is taken from the full cohort, then two supplemental samples are
#'   taken from \{Y1 < a\} and \{Y1 > b\}, respectively.
#' @param Z.dim Dimension of the covariates Z.
#' @return A list which contains parameter estimates, estimated standard error
#'   for the primary outcome model:
#'   \deqn{Y_{1}=\beta_{0}+\beta_{1}X+\beta_{2}Z,}{Y1 = beta0 + beta1*X + beta2*Z,}
#'   and the secondary outcome model:
#'   \deqn{Y_{2}=\gamma_{0}+\gamma_{1}X+\gamma_{2}Z.}{Y2 = gamma0 + gamma1*X + gamma2*Z.}
#'   The list contains the following components: \item{beta_paramEst}{parameter estimates
#'   for beta in the primary outcome model} \item{beta_stdErr}{estimated
#'   standard error for beta in the primary outcome model}
#'   \item{gamma_paramEst}{parameter estimates for gamma in the secondary
#'   outcome model} \item{gamma_stdErr}{estimated standard error for gamma in
#'   the secondary outcome model}
#' @examples
#' library(ODS)
#' # take the example data from the ODS package
#' # please see the documentation for details about the data set ods_data_secondary
#' data <- ods_data_secondary
#'
#' # divide the original cohort data into SRS, lowerODS, upperODS and NVsample
#' SRS <- data[data[,1]==1,2:ncol(data)]
#' lowerODS <- data[data[,1]==2,2:ncol(data)]
#' upperODS <- data[data[,1]==3,2:ncol(data)]
#' NVsample <- data[data[,1]==0,2:ncol(data)]
#'
#' # obtain the cut off points for ODS design. For this data, the ODS design
#' # uses mean plus and minus one standard deviation of Y1 as cut off points.
#' meanY1 <- mean(data[,2])
#' sdY1 <- sd(data[,2])
#' cutpoint <- c(meanY1-sdY1, meanY1+sdY1)
#'
#' # the data matrix SRS has Y1, Y2, X and Z. Hence the dimension of Z is ncol(SRS)-3.
#' Z.dim <- ncol(SRS)-3
#'
#' secondary_ODS(SRS, lowerODS, upperODS, NVsample, cutpoint, Z.dim)

secondary_ODS <- function(SRS, lowerODS, upperODS, NVsample, cutpoint, Z.dim) {

  Vsample <- rbind(SRS, lowerODS, upperODS)                                    # validation sample
  orgdata <- rbind(Vsample,NVsample)                                           # reorganize the data set

  N <- nrow(orgdata)                                                           # full cohort size

  # calculate the selection probability
  n0 <- nrow(SRS)                                                              # SRS sample size
  n1 <- nrow(lowerODS)                                                         # sample size taken from low tail
  n3 <- nrow(upperODS)                                                         # sample size taken from upper tail
  N1 <- sum(orgdata[,1]<=cutpoint[1])                                          # total number of subjects in low tail
  N2 <- sum(orgdata[,1]<cutpoint[2] & orgdata[,1]>cutpoint[1])                 # total number of subjects in the middle
  N3 <- sum(orgdata[,1]>=cutpoint[2])                                          # total number of subjects in upper tail
  n01 <- sum(SRS[,1]<=cutpoint[1])                                             # total number of SRS subjects in low tail
  n02 <- sum(SRS[,1]<cutpoint[2] & SRS[,1]>cutpoint[1])                        # total number of SRS subjects in the middle
  n03 <- sum(SRS[,1]>=cutpoint[2])                                             # total number of SRS subjects in the upper tail

  # observed selection probability
  p1 <- (n01+n1)/N1
  p2 <- n02/N2
  p3 <- (n03+n3)/N3

  observedweight <- rep(0,nrow(orgdata))
  for (j in 1:nrow(orgdata)) {
    p <- 0
    if (orgdata[j,1]<=cutpoint[1]) {
      p <- p1
    } else if (orgdata[j,1]>=cutpoint[2]) {
      p <- p3
    } else {
      p <- p2
    }
    observedweight[j] <- p
  }

  # find a starting value for the newton-raphson algorithm
  # regress Y1 on X and Z
  # regress Y2 on X and Z

  fit1 <- stats::lm(SRS[,1] ~ SRS[,3:(Z.dim+3)])
  fit2 <- stats::lm(SRS[,2] ~ SRS[,3:(Z.dim+3)])

  beta0estimate.initial <- as.numeric(stats::coef(fit1)[1])
  beta1estimate.initial <- as.numeric(stats::coef(fit1)[2])
  beta2estimate.initial <- as.numeric(stats::coef(fit1)[3:length(stats::coef(fit1))])

  gamma0estimate.initial <- as.numeric(stats::coef(fit2)[1])
  gamma1estimate.initial <- as.numeric(stats::coef(fit2)[2])
  gamma2estimate.initial <- as.numeric(stats::coef(fit2)[3:length(stats::coef(fit2))])

  etahat <- c(beta0estimate.initial,beta1estimate.initial,beta2estimate.initial,gamma0estimate.initial,gamma1estimate.initial,gamma2estimate.initial)

  # the covariance matrix of Y1 and Y2
  Q <- stats::cov(orgdata[,1:2])
  Q11 <- Q[1,1]
  Q12 <- Q[1,2]
  Q21 <- Q[2,1]
  Q22 <- Q[2,2]

  # regression of X on Y1, Y2, Z

  fit3 <- stats::lm(SRS[,3] ~ SRS[,c(1,2,4:(Z.dim+3))])
  phi0 <- as.numeric(stats::coef(fit3)[1])
  phi1 <- as.numeric(stats::coef(fit3)[2])
  phi2 <- as.numeric(stats::coef(fit3)[3])
  phi3 <- as.numeric(stats::coef(fit3)[4:length(stats::coef(fit3))])

  varestimate <- sum(fit3$residual^2)/(n0-Z.dim-1)

  missingindicator <- c(rep(1,n0+n1+n3),rep(0,N-n0-n1-n3))

  # newton-Raphson algorithm -------------------------------------

  while (TRUE) {

    cal <- cbind(orgdata, rep_row(etahat, nrow(orgdata)), rep(Q11,nrow(orgdata)),rep(Q12,nrow(orgdata)),rep(Q21,nrow(orgdata)),rep(Q22,nrow(orgdata)),observedweight,missingindicator,rep(phi0,nrow(orgdata)),rep(phi1,nrow(orgdata)),rep(phi2,nrow(orgdata)),rep_row(phi3,nrow(orgdata)),rep(varestimate,nrow(orgdata)),rep(Z.dim,nrow(orgdata)))

    allscorematrix <- matrix(0, 2*(Z.dim+2), nrow(orgdata))

    hessianmatrix <- matrix(0, 2*(Z.dim+2), 2*(Z.dim+2))

    for (j in 1:nrow(orgdata)) {

      allscorematrix[,j] <- secondary_ODS_score_vectorized(cal[j,])

      hessianmatrix <- hessianmatrix + secondary_ODS_hessian_vectorized(cal[j,])
    }

    scorematrix <- rowSums(allscorematrix)

    oldetahat <- etahat
    etahat <- etahat - t(solve(hessianmatrix)%*%scorematrix)
    hatdiff <- etahat - oldetahat
    convergemargin <- 10^(-10)

    if(sqrt(sum(hatdiff^2))<convergemargin) break
  }

#  print(etahat)
  betahat <- etahat[1:(length(etahat)/2)]
  gammahat <- etahat[(length(etahat)/2+1):length(etahat)]

  # calculate the consistent estimator for the asymptotic covariance matrix

  I <- -1/N*hessianmatrix

  SIGMA <- matrix(0,2*(Z.dim+2),2*(Z.dim+2))

  for (j in 1:nrow(orgdata)) {

    calvector <- allscorematrix[,j]

    SIGMA <- SIGMA + calvector%*%t(calvector)
  }

  SIGMA <- SIGMA/N
  SIGMAGEE2 <- solve(I)%*%SIGMA%*%t(solve(I))
#  print(SIGMAGEE2)
  estimatedstderr <- sqrt(diag(SIGMAGEE2)/N)

  estimatedstderrbeta <- estimatedstderr[1:(length(estimatedstderr)/2)]
  estimatedstderrgamma <- estimatedstderr[(length(estimatedstderr)/2+1):length(estimatedstderr)]

  result <- list(beta_paramEst = betahat, beta_stdErr = estimatedstderrbeta, gamma_paramEst = gammahat, gamma_stdErr = estimatedstderrgamma)

  return(result)
}

#' Data example for the secondary analysis in ODS design
#'
#' @format A matrix with 3000 rows and 7 columns: \describe{ \item{subj_ind}{An
#'   indicator variable for each subject: 1 = SRS, 2 = lowerODS, 3 = upperODS, 0
#'   = NVsample} \item{Y1}{primary outcome for which the ODS sampling scheme is
#'   based on} \item{Y2}{a secondary outcome} \item{X}{expensive exposure}
#'   \item{Z1}{a simulated covariate} \item{Z2}{a simulated covariate}
#'   \item{Z3}{a simulated covariate} }
#'
#' @source A simulated data set
"ods_data_secondary"
