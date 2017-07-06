#simulate_grpLM
#
# This file contains function to simulate data according to a linear model with grouped features
# -----------------------

#' Simulate Data According to a linear model with grouped feautres
#'
#' The predictors are split into group carrying different amount of information (different coeffcient sizes)
#' The simulation is based on multivariate normal distributed predictors.
#'
#' @param n number of samples
#' @param p number of features
#' @param beta_best coefficient amplitude in most-informative group
#' @param blockSize block size of correlated blocks
#' @param block_cor overall correlation between all features
#' @param equiCor coefficient amplitude in most-informative group
#' @param G number of groups of features
#' @param diffFac differential signal per group \eqn{\beta_{g+1} = \beta_g *diffFac}
#' @param SLovG (overall) sparsity level, i.e. percentage of zero groups
#' @param  SLinG sparsity level within each group, i.e. percentage of zero feautres within group
#' @param model model for response (e.g. linear, based on logs)
#' @param  family family of response, e.g. gaussian, binomial
#' @param  sigma2 noise variance
#' @param onlypos boolean determining wether within a group all non-zero coeffcients should have the same postive value (if TRUE) or wether half get a negative sign (if FALSE)
#'
#' @details Special cases would be orthogonal design (blockSize=1,equiCor=0) or equicorrelated design (blockSize=1, equiCor=rho).
#'
#' @return A list containing
#' \item{X}{nxp matrix of feautres with given correlation structure}
#' \item{y}{vector of length p of simulated response}
#' \item{annot}{factor containing group assignemnt of featrues}
#' \item{beta0}{true coefficients in the model}
#' \item{family}{family used for response}
# -----------------------

simulate_grpLM<-function(n=100,p=1000, beta_best=3, G=10, block_cor=0.3, blockSize=10, equiCor=0.1,
                           SLinG=0.7, SLovG=0.5, model="linear", diffFac=0.7, seed=167, sigma2=0.1, family="gaussian", onlypos=T){

  set.seed(seed)

  stopifnot(p%%G==0) #number of features should be a multiple of the number of groups
  stopifnot(p%%blockSize==0) #number of features should be a multiple of the blockSize
  stopifnot((SLovG*G)%%1==0) #number of groups should be a multiple of number of zero_groups
  stopifnot((SLinG*p/G)%%1==0) #number of features per group should be a multiple of number of zero feautres per group

  #construct design
  if(block_cor!=0 | equiCor!=0){
    block<-matrix(block_cor,nrow=blockSize, ncol=blockSize)
    diag(block)<-1
    Sigma<-Matrix::bdiag(rep(list(block),p/blockSize))
    if(equiCor!=0) Sigma[Sigma==0]<-equiCor
    X<- mvtnorm::rmvnorm(n,rep(0,p),as.matrix(Sigma))
  } else X<-matrix(rnorm(p*n,0,1), ncol=p, nrow=n)

  colnames(X)<-paste("Feature",1:p, sep="_")

  #choose non-zero features and amplitude
  #base-level signal per group
  betabase<-sapply(1:G, function(i) beta_best*diffFac^(i-1))
  #set weakest  proportion of SLovG of group to zero
  betabase[(G-SLovG*G+1):G]<-0
  #within each group set a proportion of SLinG of feature to zero
  beta0<-rep(betabase, each=p/G)
  zeros<-as.numeric(sapply(1:G, function(i) (i-1)*(p/G)+sample(1:(p/G), SLinG*(p/G))))
  beta0[zeros]<-0
  names(beta0)<-colnames(X)
  #use postive and dengative coefficents (50-50)
  if(!onlypos) {
    neg<-sample(which(beta0!=0), round(sum(beta0!=0)/2))
    beta0[neg]<- -beta0[neg]
  }

  #simulate response according to model and family
  if(model=="linear" & family=="gaussian")
    y<-X%*%beta0+rnorm(n,0,sigma2)
  else if(model=="log" & family=="gaussian")
    y<-log(X+1000)%*%beta0+rnorm(n,0,sigma2)
  else if(model=="linear" & family=="binomial"){
    exponent<-X%*%beta0
    prob<-exp(exponent)/(1+exp(exponent))
    y<-rbinom(length(prob),1,prob)
  }
  else if(model=="linear" & family=="binomial"){
    exponent<-log(X+1000)%*%beta0
    prob<-exp(exponent)/(1+exp(exponent))
    y<-rbinom(length(prob),1,prob)
  }
  else stop("model needs to be either linear or log")

  #annotation by groups
  annot<-as.factor(rep(1:G, each=p/G))

  return(list(X=X, y=y, annot=annot, beta0=beta0, family=family))
}


simulateExplicit <- function(n,p, beta, sigma2, seed, block_cor=0, equiCor=0, blockSize=0){

  stopifnot(p==length(beta)) #number of features should be a multiple of the number of groups
  stopifnot(block_cor==0 | p%%blockSize==0) #number of features should be a multiple of the blockSize

  set.seed(seed)

  #construct design
  if(block_cor!=0 | equiCor!=0){
    block<-matrix(block_cor,nrow=blockSize, ncol=blockSize)
    diag(block)<-1
    Sigma<-Matrix::bdiag(rep(list(block),p/blockSize))
    if(equiCor!=0) Sigma[Sigma==0]<-equiCor
    X<- mvtnorm::rmvnorm(n,rep(0,p),as.matrix(Sigma))
  } else X<-matrix(rnorm(p*n,0,1), ncol=p, nrow=n)

  #simulate response
  y<-X%*%beta+rnorm(n,0,sigma2)

  return(list(y=y, X=X))
}
