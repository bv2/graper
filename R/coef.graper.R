#' @title Get estimated coefficients from a graper object
#' @name coef.graper
#' @description Function to obtain estimated coefficients
#' from a fitted graper model.
#' @param object fitted graper model as obtained from  \code{\link{graper}}
#' @param include_intercept whether to include the estimated
#' intercept value in the output
#' @param ... other arguments
#' @return 1-Column matrix of estimated coefficients.
#' @export
#' @importFrom methods is
#' @examples
#' # create data
#' dat <- makeExampleData()
#' fit <- graper(dat$X, dat$y, dat$annot)
#' coef(fit)

coef.graper <- function(object, include_intercept = TRUE, ...){
    # sanity check
    if(!is(object, "graper")) {
        stop("object needs to be a graper object.")
    }

    if(is.null(object$intercept)){
        object$intercept <- 0
    }
    coefs  <- rbind("(intercept)" = object$intercept, object$EW_beta)
    if(!is.null( object$Options$featurenames)) {
        rownames(coefs) <- c("(intercept)", object$Options$featurenames)
    }
    if(include_intercept) {
        coefs
    } else {
        coefs[-1,,drop=FALSE]
    }
}

