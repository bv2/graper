# This file contains function to simulate data according to a linear model with grouped features
# - simulateExplicit uses fixed coefficients
# - makeExampleData simulates coeffients from the Bayesian model

#' @title Simulate Data According to a linear model with given coefficients
#' @name simulateExplicit
#' @param n number of samples
#' @param p number of features
#' @param beta model coefficients
#' @param sigma2 noise variance
#' @param seed random seed
#' @param exp_decay_cor covariance matrix with exponential decay of covariance
#' @param equiCor coefficient amplitude in most-informative group
#' @param block_cor correlation within feature blocks
#' @param blockSize size of correlated blocks (if block_cor>0)
#' @return list with simulated response y and design matrix X
#' @export

simulateExplicit <- function(n, p, beta, sigma2, seed, exp_decay_cor =0, block_cor = 0, equiCor = 0, blockSize = 0) {

    stopifnot(p == length(beta))  #number of features should be a multiple of the number of groups
    stopifnot(block_cor == 0 | p%%blockSize == 0)  #number of features should be a multiple of the blockSize

    set.seed(seed)

    stopifnot(exp_decay_cor == 0 | ( block_cor == 0 & equiCor == 0))
    # construct design
    if(exp_decay_cor>0){
      # option A: construct covariance matrix with exponential decay of covariance
    pp <- expand.grid(1:p, 1:p)
    Sigma <- matrix(exp(-1/exp_decay_cor*abs(pp[,1]-pp[,2])), nrow=p)
    X <- mvtnorm::rmvnorm(n, rep(0, p), Sigma)
    } else if (block_cor != 0 | equiCor != 0) {
        # option B: construct covariance matrix with block correlation
        block <- matrix(block_cor, nrow = blockSize, ncol = blockSize)
        diag(block) <- 1
        Sigma <- Matrix::bdiag(rep(list(block), p/blockSize))
        if (equiCor != 0)
            Sigma[Sigma == 0] <- equiCor
        X <- mvtnorm::rmvnorm(n, rep(0, p), as.matrix(Sigma))
    } else {
      #option C: indep.
      X <- matrix(rnorm(p * n, 0, 1), ncol = p, nrow = n)
    }

    # simulate response
    y <- X %*% beta + rnorm(n, 0, sqrt(sigma2))

    return(list(y = y, X = X))
}

#' Simulate Data According to a linear model with given coefficients and toeplitz design matrix
#' @param n number of samples
#' @param p number of features
#' @param beta model coefficients
#' @param sigma2 noise variance
#' @param seed random seed
#' @param rho Toeplitz parameter
#' @return list with simulated response y and design matrix X
#' @export
simulateData_toeplitz <- function(n, p, beta, sigma2, seed, rho) {

    stopifnot(p == length(beta))  #number of features should be a multiple of the number of groups
    set.seed(seed)

    # construct design
    Sigma = toeplitz(rho^(0:(p-1)))
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)

    # simulate response
    y <- X %*% beta + rnorm(n, 0, sqrt(sigma2))

    return(list(y = y, X = X))
}


#' @title Simulate example data with groups of equal size
#' @name makeExampleData
#' @description Simulate data from the Bayesian model with specified parameters $\gamma$ and $\tau$.
#' @param n number of samples
#' @param p number of features
#' @param g number of groups
#' @param gammas vector of length g, sepcifying the slab precision of the prior on beta per group
#' @param pis vector of length g, sepcifying the probability of s to be 1 (slab)
#' @param tau noise precision
#' @param rho correlation of design matrix (Toeplitz structure)
#' @return list containin the design matrix X, the response y, the feautre annotation to
#'  groups annot as well as the different parameters in the Bayesian model and the correlation strength rho
#' @export
makeExampleData <- function(n=100, p=200, g=4, gammas=c(0.1,1,10,100), pis=c(0.5,0.5,0.5,0.5), tau=1, rho=0) {

    #checks
    stopifnot(p%%g==0)
    stopifnot(length(gammas)==g)
    stopifnot(length(pis)==g)

    makeExampleDataWithUnequalGroups(n=n, pg=rep(p/g,g), gammas=gammas, pis=pis, tau=tau, rho=rho)
}

#' @title Simulate example data with groups of unequal size
#' @name makeExampleDataWithUnequalGroups
#' @description Simulate data from the Bayesian model with specified parameters $\gamma$ and $\tau$.
#' @param n number of samples
#' @param pg vector of length g (desired number of groups) with number of features per group
#' @param gammas vector of length g, sepcifying the slab precision of the prior on beta per group
#' @param pis vector of length g, sepcifying the probability of s to be 1 (slab)
#' @param tau noise precision
#' @param rho correlation of design matrix (Toeplitz structure)
#' @return list containin the design matrix X, the response y, the feautre annotation to
#'  groups annot as well as the different parameters in the Bayesian model and the correlation strength rho
#' @export
makeExampleDataWithUnequalGroups <- function(n=100, pg=c(100,100,10,10), gammas=c(0.1,10,0.1,10), pis=c(0.5,0.5,0.5,0.5), tau=1, rho=0) {

    #checks
    g <- length(pg)
    p <- sum(pg)
    stopifnot(length(gammas)==g)
    stopifnot(length(pis)==g)

    # construct design
    Sigma <- toeplitz(rho^(0:(p-1)))
    X <- matrix(rnorm(n*p),n,p) %*% chol(Sigma)
    X <- scale(X)

    # simulate coefficients
    beta_tilde <- Reduce(c,lapply(1:g, function(k) rnorm(pg[k],0,sqrt(1/gammas[k]))))
    s <- Reduce(c,lapply(1:g, function(k) rbinom(pg[k],1,pis[k])))
    beta <- s* beta_tilde

    #simulate response
    y <- rnorm(n, X%*%beta, 1/sqrt(tau))
    annot <- rep(1:length(pg), times=pg)
    list(X=X, y=y, annot=annot, gammas=gammas, pis=pis, tau=tau, rho=rho, g=g, p=p, n=n, beta=beta, s=s, beta_tilde=beta_tilde)
}
