#' @title Simulate example data from the grpRR model
#' @name makeExampleData
#' @description Simulate data from the grpRR model with groups of equal size and pre-specified parameters gamma, pi and tau.
#' @param n number of samples
#' @param p number of features
#' @param g number of groups
#' @param gammas vector of length g, specifying the slab precision
#'  of the prior on beta per group
#' @param pis vector of length g, specifying the probability of s to be 1 (slab)
#' @param tau noise precision
#' @param rho correlation of design matrix (Toeplitz structure)
#' @param response "gaussian" for continuous response from a linear regression model,
#'  "bernoulli"  for a binary response from a logistic regression model.
#' @param intercept model intercept (default: 0)
#' @return list containing the design matrix \code{X},
#'  the response \code{y}, the feature annotation to
#'  groups \code{annot} as well as the different parameters
#'  in the Bayesian model and the correlation strength rho
#' @export
#' @examples
#' dat <- makeExampleData()

makeExampleData <- function(n=100, p=200, g=4,
                            gammas=c(0.1,1,10,100), pis=c(0.5,0.5,0.5,0.5),
                            tau=1, rho=0, response = "gaussian", intercept = 0) {

    #checks
    stopifnot(p%%g==0)
    stopifnot(length(gammas)==g)
    stopifnot(length(pis)==g)

    makeExampleDataWithUnequalGroups(n=n, pg=rep(p/g,g), gammas=gammas, pis=pis,
                                     tau=tau, rho=rho, response=response, intercept = intercept)
}

#' @title Simulate example data from the grpRR model with groups of unequal size
#' @name makeExampleDataWithUnequalGroups
#' @description Simulate data from the grpRR model with groups of unequal size and pre-specified parameters gamma, pi and tau.
#' @param n number of samples
#' @param pg vector of length g (desired number of groups) with number of features per group
#' @param gammas vector of length g, specifying the slab precision of the prior on beta per group
#' @param pis vector of length g, specifying the probability of s to be 1 (slab)
#' @param tau noise precision (only relevant for gaussian response)
#' @param rho correlation of design matrix (Toeplitz structure)
#' @param response "gaussian" for continuous response from a linear regression model,
#'  "bernoulli"  for a binary response from a logistic regression model.
#' @param intercept model intercept (default: 0)
#' @return list containin the design matrix \code{X}, the response \code{y}, the feature annotation to
#'  groups \code{annot} as well as the different parameters in the Bayesian model
#'   and the correlation strength rho
#' @export
#' @importFrom stats toeplitz rnorm rbinom
#' @examples
#' dat <- makeExampleDataWithUnequalGroups()

makeExampleDataWithUnequalGroups <- function(n=100, pg=c(100,100,10,10),
                                             gammas=c(0.1,10,0.1,10),
                                             pis=c(0.5,0.5,0.5,0.5),
                                             tau=1, rho=0, response = "gaussian",
                                             intercept = 0) {

    #checks
    g <- length(pg)
    p <- sum(pg)
    stopifnot(length(gammas)==g)
    stopifnot(length(pis)==g)
    if(!response %in% c("gaussian", "bernoulli")){
      stop("Response needs to be 'gaussian' or 'bernoulli'.")
    }

    # construct design
    Sigma <- stats::toeplitz(rho^(0:(p-1)))
    X <- matrix(stats::rnorm(n*p),n,p) %*% chol(Sigma)
    X <- scale(X)

    # simulate coefficients
    beta_tilde <- Reduce(c,lapply(seq_len(g),
                                  function(k) stats::rnorm(pg[k],0,sqrt(1/gammas[k]))))
    s <- Reduce(c,lapply(seq_len(g),
                         function(k) stats::rbinom(pg[k],1,pis[k])))
    beta <- s* beta_tilde

    #simulate response
    if(response == "gaussian"){
      y <- stats::rnorm(n, X%*%beta + intercept, 1/sqrt(tau))
    } else y <- stats::rbinom(n, 1, 1/(1+exp(- (X%*%beta +intercept)) ))

    #group annotations
    annot <- rep(seq_along(pg), times=pg)

    list(X=X, y=y, annot=annot, gammas=gammas, pis=pis, tau=tau, rho=rho,
         g=g, p=p, n=n, beta=beta, s=s, beta_tilde=beta_tilde, intercept=intercept)
}
