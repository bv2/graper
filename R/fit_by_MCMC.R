#' fit_by_MCMC
#'
#' Function to fit the hier. Bayes model used for grpRR using MCMC for inference. Only implemented for the normal prior and gaussian response at the moment.
#' @param X Design matrix of size n x p
#' @param y Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param iter maximum number of iterations
#' @param chains number of chains to run
#' @param warmup number of iterations for burn-in period
#' @param ... other parameters that can be passed to stan
#' @return SummaryList containing model coefficients, penalty factors, run time, fit model
#' @import rstan
#' @export 


fit_by_MCMC <- function(X,y, annot, iter=12000, warmup=2000, chains=3, verbose=T,thin=1,...){

  y <- as.vector(y)
  y <- scale(y, center = T, scale = F)
  X <- scale(X, center = T, scale = F)
  stopifnot(length(y) == nrow(X))
  stopifnot(length(annot) == ncol(X))
  
  n <- nrow(X)
  p <- ncol(X)
  G <- length(unique(annot))
    
    
  stanmodelcode = "
  data {                      // ***Data block
  int<lower=1> N;           // Sample size
  int<lower=1> P;           // Dimension of model matrix  
  int<lower=1> G;           // Number of groups
  
  matrix[N, P] X;           // Model Matrix
  vector[N] y;              // Target variable
  int<lower=1> annot[P];          // coding group membership per feature
  vector[P] zeros;          // zero vector
  }
  
  parameters {                // ***Parameters block
  vector[P] beta;           // Coefficient vector
  real<lower=0> tau;        // Noise precision
  real<lower=0> gamma[G];           //Parameters controlling shrinkage
  }
  
  model {                     // ***Model block
  vector[N] mu;
  vector[P] long_gamma;     //transfomed value of gamma
  int gr;
  
  mu <- X * beta;           // Creation of linear predictor
  for (i in 1:P) {
  gr = annot[i];
  long_gamma[i] = gamma[gr];  //transform gamma using the normal cdf
  }
  
  // priors
  beta ~  normal(zeros, long_gamma);
  tau ~ gamma(0.001,0.001);
  gamma ~ gamma(0.001,0.001);
  
  
  // likelihood
  y ~ normal(mu, 1/tau);
  }
  "
  
  dat4MCMC <- list(N = n, P = p, G = G, y = as.vector(y) , X = X, annot = annot, zeros= rep(0,p))


  ### Run the model and examine results ###
  tmp <- Sys.time()
  fit_MCMC = stan(model_code=stanmodelcode, data=dat4MCMC,  iter=iter,
               warmup=warmup, thin=thin, chains=chains, verbose=verbose, ...)
  time_MCMC <- Sys.time() - tmp
  
  # Get all relevant quantities in same format as for other methods
  MCMCbetas <- get_posterior_mean(fit_MCMC, pars="beta")[,chains+1]
  MCMCintercept <- attr(y, "scaled:center")-attr(X, "scaled:center")%*%MCMCbetas
  MCMCgammas <- get_posterior_mean(fit_MCMC, pars="gamma")[,chains+1]
  MCMCsummary <- list()
  MCMCsummary$runtime <- as.numeric(time_MCMC)
  MCMCsummary$pf <- as.numeric(MCMCgammas)
  MCMCsummary$beta <- MCMCbetas
  MCMCsummary$intercept <- MCMCintercept
  MCMCsummary$sparsity <- NULL
  MCMCsummary$out <- fit_MCMC
  
  return(MCMCsummary)
}