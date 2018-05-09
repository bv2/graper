# R function to call the different C functions defined in fit_grpRR.cpp fit_grpRR This file contains the main function to fit a
# grpRR model.

# ---------------------------
#'  Fit grpRR model
#'
#'  Main function to fit a grpRR model.
#' @param X Design matrix of size n x p
#' @param y Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param factoriseQ If true, the variational distribution is assumed to fully factorize across features (rougher approx., but faster). If spikeslab=F, this is always done.
#' @param spikeslab If true, a spike and slab prior is used instead of a normal prior
#' @param nogamma If true, the normal prior will have same variance for all groups (only relevant for SS)
#' @param d_tau hyper-parameters for prior of tau (noise precision)
#' @param r_tau hyper-parameters for prior of tau (noise precision)
#' @param d_gamma hyper-parameters for prior of gamma (coeffient's prior precision)
#' @param r_gamma hyper-parameters for prior of gamma (coeffient's prior precision)
#' @param r_pi hyper-parameters for prior of pi (Spike-Slab-Bernoulli variable prior probabiliy of being 1)
#' @param d_pi hyper-parameters for prior of pi (Spike-Slab-Bernoulli variable prior probabiliy of being 1)
#' @param max_iter maximum number of iterations
#' @param th convergence threshold for ELB
#' @param intercept boolean, indicating wether to fit an intercept
#' @param calcELB boolean, indicating wether to calculate ELB
#' @param verbose boolean, indicating wether to print out intermediate messages during fitting
#' @param freqELB determines frequency at which ELB is to be calculated, i.e. each feqELB-th iteration
#' @return List of fitted parameters .....


#' @useDynLib grpRR
#' @import Rcpp
#' @export


fit_grpRR <- function(X, y, annot, factoriseQ = TRUE, spikeslab = TRUE, d_tau = 0.001, r_tau = 0.001, d_gamma = 0.001, r_gamma = 0.001,
    r_pi = 1, d_pi = 1, max_iter = 1000, th = 1e-05, intercept = TRUE, calcELB = TRUE, verbose = TRUE, freqELB = 10, family = "gaussian",
    nogamma=F, standardize=TRUE) {

    stopifnot(ncol(X) == length(annot))

    #nogamma only of use when spikeslab
    if(!spikeslab & !nogamma) nogamma <- FALSE

    # check structure of annot: needs to be 1:g with 1 <-> frist group etc
    annot <- factor(annot, levels = unique(annot))

    # get data dimension
    p <- ncol(X)  #no of features
    n <- nrow(X)  #no of samples

    # get group structure
    g <- length(unique(annot))
    NoPerGroup <- sapply(unique(annot), function(x) sum(annot == x))
    names(NoPerGroup) <- unique(annot)
    message(paste("Fitting a model with", g, "groups", n, "samples and", p , "features."))

    if(standardize){
        X <- scale(X, center = F, scale=T)
        sf <- attr(X, "scaled:scale")
    } else sf <- rep(1,p)

    if (family == "gaussian") {
        # remove intercept effect by centering X and y
        if (intercept) {
            X <- scale(X, center = T, scale = F)
            y <- scale(y, center = T, scale = F)
        }

        # call C function depending on FacType and spikeslap arguments
        # eventually this should be repeated x-times and model with lowest ELBO chosen, for now keep it for robustness anaylsis
        if (spikeslab) {
            if (!factoriseQ)
                warning("Using fully factorized approach with a spike and slab prior")
            # initialize slab mean and spike prob.
            mu_init <- rnorm(p)
            # psi_init <- runif(p)
            psi_init <- rep(0.5,p)
            if(!nogamma)
            res <- grRRCpp_sparse_ff(X, y, annot, g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB,
                verbose, freqELB, mu_init, psi_init)
            else
            res <- grRRCpp_sparse_ff_nogamma(X, y, annot, g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB,
                verbose, freqELB, mu_init, psi_init)
        } else {
            if (factoriseQ) {
                # initialize coefficients mean randomly
                mu_init <- rnorm(p)
                res <- grRRCpp_dense_ff(X, y, annot, g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose,
                  freqELB, mu_init) 
                } else res <- grRRCpp_dense_nf(X, y, annot, g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose,
                freqELB)
        }

        # calculate intercept
        if (intercept)
            intercept <- attr(y, "scaled:center") - sum(attr(X, "scaled:center")*res$EW_beta) else intercept <- NULL

        # give proper names
        if(!nogamma) rownames(res$EW_gamma) <- unique(annot)
        # return mean of approximate posterior (other quantities of interes: tau, lower bound on model evidence etc)
        res <- append(res, list(intercept = intercept))
    } else if (family == "binomial") {
        # in case intercept =TRUE  this is removed iteratively during training in terms of the variational approximation
        if (spikeslab){
            if (!factoriseQ) 
                warning("Using fully factorized approach with a spike and slab prior")
            # initialize slab mean and spike prob. randomly
            mu_init <- rnorm(p)
            psi_init <- runif(p)
            res <- grpRRCpp_sparse_logistic_ff(X, y, annot, g, NoPerGroup, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB,
                verbose, freqELB, mu_init, psi_init, intercept)
        } else {
            if (factoriseQ){
                # initialize coefficients mean randomly
                mu_init <- rnorm(p)
                res <- grpRRCpp_logistic_ff(X, y, annot, g, NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB, mu_init, intercept)
                } else {
                	warning("factoriseQ=FALSE is not maintained currently for the logistic model. No intercept option and ELBO available.")
                	res <- grpRRCpp_logistic_nf(X, y, annot, g, NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB)
                }
        }
        if(!intercept) res$intercept <- NULL
    }
    else stop("Family not implemented. Needs to be either binomial or gaussian.")
    
    # revert coefficients to original scale
        res$EW_beta <- res$EW_beta/sf
        res$Sigma_beta <- diag(sf) %*% res$Sigma_beta %*% diag(sf)
    
    return(res)
}
