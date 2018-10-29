#' @title Fit a regression model with grpRR
#' @name grpRR
#' @description Fit a regression model with grpRR given a matrix of predictors (\code{X}), a response vector (\code{y}) and
#' a vector of group memberships for each predictor in \code{X} (\code{annot}). For each group a different strength of penalization is determined adaptively.
#' @param X design matrix of size n (samples) x p (features)
#' @param y response vector of size n
#' @param annot factor of length p indicating group membership of each feature (column) in X
#' @param family Likelihood model for the response,
#'  either "gaussian" for linear regression or "binomial" for logistic regression
#' @param factoriseQ if set to TRUE, the variational distribution is assumed
#'  to fully factorize across features (faster, default). FALSE uses a multivariate variational distribution.
#' @param spikeslab if set to TRUE, a spike and slab prior on the coefficients (default).
#' @param d_tau hyper-parameters for prior of tau (noise precision)
#' @param r_tau hyper-parameters for prior of tau (noise precision)
#' @param d_gamma hyper-parameters for prior of gamma (coefficients' prior precision)
#' @param r_gamma hyper-parameters for prior of gamma (coefficients' prior precision)
#' @param r_pi hyper-parameters for Beta prior of the mixture probabilities in the spike and slab prior
#' @param d_pi hyper-parameters for Beta prior of the mixture probabilities in the spike and slab prior
#' @param max_iter maximum number of iterations
#' @param th convergence threshold for the evidence lower bound (ELB)
#' @param intercept whether to include an intercept into the model
#' @param calcELB whether to calculate the evidence lower bound (ELB)
#' @param verbose  whether to print out intermediate messages during fitting
#' @param freqELB frequency at which the evidence lower bound (ELB) is to be calculated,
#'  i.e. each freqELB-th iteration
#' @param n_rep number of repetitions with different random initializations  to be fit
#' @param standardize whether to standardize the predictors to unit variance
#' @param init_psi initial value for the spike variables
#' @param nogamma if TRUE, the normal prior will have same variance for all groups
#' (only relevant for spikeslab = TRUE)
#' @details The function trains the grpRR model given a matrix of predictors (\code{X}), a response vector (\code{y}) and
#' a vector of group memberships for each predictor in \code{X} (\code{annot}).
#' For each feature group as specified in \code{annot} a penalty factor and sparsity level is learnt.
#'
#'  By default it uses a Spike-and-Slab prior on the coefficients and uses a
#'  fully factorized variational distribution in the inference.
#'  This provides a fast way to train the model. Using \code{spikeslab=FALSE} a
#'  ridge regression like model can be fitted using a normal instead of the spike and slab prior.
#'  Setting \code{factoriseQ = FALSE} gives a more exact inference
#'  scheme based on a multivariate variational distribution, but can be much slower.
#'
#'  As the optimization is non-convex is can
#'  be helpful to use multiple random initializations by setting \code{n_rep} to a value larger 1. The returned model is then chosen
#'  as the optimal fit with respect to the evidence lower bound (ELB).
#'
#'  Depending on the response vector a linear regression model (\code{family = "gaussian"}) or a logistic regression model
#'  (\code{family = "binomial"}) is fitted. Note, that the implementation of logistic regression is still experimental.
#'
#' @return A grpRR object containing
#' \describe{
#' \item{EW_beta}{estimated model coefficients in liner/logistic regression}
#' \item{EW_s}{estimated posterior-inclusion probabilities for each feature}
#' \item{intercept}{estimated intercept term}
#' \item{annot}{annotation vector of features to the groups as specified when calling \code{\link{grpRR}}}
#' \item{EW_gamma}{estimated penalty factor per group}
#' \item{EW_pi}{estimated sparsity level per group (from 1 (dense) to 0 (sparse))}
#' \item{EW_tau}{estimated noise precision}
#' \item{sigma2_tildebeta_s1, EW_tildebeta_s1, alpha_gamma, alpha_tau, beta_tau, Sigma_beta, alpha_pi, beta_pi}{parameters of
#'  the variational distributions of beta, gamma, tau and pi}
#' \item{ELB}{final value of the evidence lower bound}
#' \item{ELB_trace}{values of the  evidence lower bound for all iterations}
#' \item{Options}{other options used when calling \code{\link{grpRR}}}
#' }
#' @useDynLib grpRR
#' @import Rcpp
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#'
#' # fit a sparse model with spike and slab prior
#' fit <- grpRR(dat$X, dat$y, dat$annot)
#' fit # print fitted object
#' beta <- coef(fit, include_intercept = FALSE) # model coeffients
#' pips <- getPIPs(fit) # posterior inclusion probabilities
#' pf <- fit$EW_gamma # penalty factors per group
#' sparsities <- fit$EW_pi # sparsity levels per group
#'
#' # fit a dense model without spike and slab prior
#' fit <- grpRR(dat$X, dat$y, dat$annot, spikeslab = FALSE)
#'
#' # fit a dense model using a multivariate variational distribution
#' fit <- grpRR(dat$X, dat$y, dat$annot, factoriseQ = TRUE, spikeslab = FALSE)


grpRR <- function(X, y, annot, factoriseQ = TRUE, spikeslab = TRUE,
                      intercept = TRUE, family = "gaussian",
                      standardize=TRUE, n_rep=1,
                      max_iter = 3000, th = 0.01,
                      d_tau = 0.001, r_tau = 0.001,
                      d_gamma = 0.001, r_gamma = 0.001,
                      r_pi = 1, d_pi = 1,
                      calcELB = TRUE, verbose = TRUE,
                      freqELB = 1, nogamma=FALSE,
                      init_psi=1) {

    stopifnot(ncol(X) == length(annot))

    # nogamma only of use when spikeslab
    if(!spikeslab & !nogamma) nogamma <- FALSE

    annot <- factor(annot, levels = unique(annot))

    # get data dimension
    p <- ncol(X)  #no of features
    n <- nrow(X)  #no of samples

    # get group structure
    g <- length(unique(annot))
    NoPerGroup <- vapply(unique(annot), function(x){
      sum(annot == x)
      }, numeric(1))
    names(NoPerGroup) <- unique(annot)
    message(paste("Fitting a model with", g, "groups,", n, "samples and", p , "features."))

    if(standardize){
        X <- scale(X, center = FALSE, scale=TRUE)
        sf <- attr(X, "scaled:scale")
    } else sf <- rep(1,p)

    if(family == "binomial"){
      if(calcELB){
        calcELB <- FALSE
        warning("The implementation of logistic regression is still experimental.
                ELB calculations are not yet implemented for logistic regression.")
      }
    }

    if(!calcELB & n_rep >1) {
            warning("For model selection with multiple trials calcELB needs to be set to TRUE.
                    Only using a single trial now.")
            n_rep <-1
    }

    reslist <- lapply(seq_len(n_rep), function(rep){
        message("Fitting with random init ", rep)

        if (family == "gaussian") {
            # remove intercept effect by centering X and y
            if (intercept) {
                X <- scale(X, center = TRUE, scale = FALSE)
                y <- scale(y, center = TRUE, scale = FALSE)
            }

            # call C function depending on FacType and spikeslap argument
            if (spikeslab) {
                if (!factoriseQ)
                    warning("Using fully factorized approach with a spike and slab prior")
                # initialize slab mean and spike prob.
                mu_init <- rnorm(p)
                # psi_init <- runif(p)
                psi_init <- rep(init_psi,p)
                if(!nogamma)
                res <- grRRCpp_sparse_ff(X, y, annot, g, NoPerGroup, d_tau, r_tau,
                                         d_gamma, r_gamma, r_pi, d_pi, max_iter,
                                         th, calcELB, verbose, freqELB, mu_init,
                                         psi_init)
                else
                res <- grRRCpp_sparse_ff_nogamma(X, y, annot, g, NoPerGroup,
                                                 d_tau, r_tau, d_gamma, r_gamma,
                                                 r_pi, d_pi, max_iter, th, calcELB,
                                                 verbose, freqELB, mu_init, psi_init)
            } else {
                if (factoriseQ) {
                    # initialize coefficients mean randomly
                    mu_init <- rnorm(p)
                    res <- grRRCpp_dense_ff(X, y, annot, g, NoPerGroup, d_tau, r_tau,
                                            d_gamma, r_gamma, max_iter, th, calcELB,
                                            verbose,freqELB, mu_init)
                    } else {
                      message("You are using no factorization of the variational distribution.
                              This might take some time to compute.
                              Set factoriseQ = TRUE for fast solution.")
                      res <- grRRCpp_dense_nf(X, y, annot, g, NoPerGroup, d_tau, r_tau, d_gamma,
                                              r_gamma, max_iter, th, calcELB, verbose, freqELB)
                    }
            }

        # calculate intercept
            if (intercept) {
                res$intercept <- attr(y, "scaled:center") - sum(attr(X, "scaled:center")*res$EW_beta)
                } else {
                  res$intercept <- NULL
                }

            # give proper names
            if(!nogamma) {
              rownames(res$EW_gamma) <- unique(annot)
            }
            res

        } else if (family == "binomial") {
            if (spikeslab){
                if (!factoriseQ)
                    warning("Using fully factorized approach with a spike and slab prior")
                # initialize slab mean and spike prob. randomly
                mu_init <- rnorm(p)
                # psi_init <- runif(p)
                psi_init <- rep(init_psi,p)
                res <- grpRRCpp_sparse_logistic_ff(X, y, annot, g, NoPerGroup, d_gamma, r_gamma,
                                                   r_pi, d_pi, max_iter, th, calcELB,verbose,
                                                   freqELB, mu_init, psi_init, intercept)
            } else {
                if (factoriseQ){
                    # initialize coefficients mean randomly
                    mu_init <- rnorm(p)
                    res <- grpRRCpp_logistic_ff(X, y, annot, g, NoPerGroup, d_gamma, r_gamma,
                                                max_iter, th, calcELB, verbose, freqELB,
                                                mu_init, intercept)
                    } else {
                      warning("factoriseQ=FALSE is not maintained currently for the logistic model.
                              No intercept option and ELB available.")
                      res <- grpRRCpp_logistic_nf(X, y, annot, g, NoPerGroup,
                                                  d_gamma, r_gamma, max_iter,
                                                  th, calcELB, verbose, freqELB)
                    }
            }
            if(!intercept) {
              res$intercept <- NULL
            }
            res
        }
        else stop("Family not implemented. Needs to be either binomial or gaussian.")
})

    if(n_rep==1) {
        res <- reslist[[1]]
    }  else {
            best_idx <- which.max(vapply(reslist, function(l) l$ELB, numeric(1)))
            if(is.na(best_idx) | is.null(best_idx)) {
                warning("Model selection based on ELB encountered errors.
                        Returned model is picked arbitrarily!")
                best_idx <- 1
            }
            res <- reslist[[best_idx]]
    }

    # revert coefficients to original scale
        if(standardize){
          res$EW_beta <- res$EW_beta/sf
          if(!factoriseQ) {
            res$Sigma_beta <- diag(1/sf) %*% res$Sigma_beta %*% diag(1/sf)
          } else {
            res$Sigma_beta <-diag(1/(sf^2) * diag(as.matrix(res$Sigma_beta)))
          }
        }
        res$annot <- annot
        res$Options <- list(factoriseQ = factoriseQ,
                            spikeslab = spikeslab,
                            d_tau = d_tau, r_tau = r_tau,
                            d_gamma = d_gamma, r_gamma =r_gamma,
                            r_pi = r_pi, d_pi = d_pi,
                            max_iter = max_iter, th = th,
                            intercept = intercept,
                            calcELB = calcELB,
                            verbose = verbose,
                            freqELB = freqELB,
                            family = family,
                            nogamma=nogamma,
                            standardize=standardize,
                            featurenames = colnames(X))

    #remove ELB slot if not calculated
    if(all(is.na(res$ELB_trace) | !is.finite(res$ELB_trace))){
      res$ELB_trace <- NULL
    }
    class(res) <- "grpRR"
    return(res)
}
