#' @title Get estimated coefficients from a grpRR object
#' @name coef.grpRR
#' @description Function to obtain estimated coefficients from a fitted grpRR model.
#' @param object fitted grpRR model as obtained from  \code{\link{grpRR}}
#' @param include_intercept whether to include the estimated intercept value in the output
#' @param ... other arguments
#' @return 1-Column matrix of estimated coefficients.
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#' fit <- grpRR(dat$X, dat$y, dat$annot)
#' coef(fit)

coef.grpRR <- function(object, include_intercept = TRUE, ...){
  # sanity check
  if(class(object) != "grpRR") {
    stop("object needs to be a grpRR object.")
  }

  if(is.null(object$intercept)){
    object$intercept <- 0
  }
  coefs  <- rbind("(intercept)" = object$intercept, object$EW_beta)
  if(!is.null( object$Options$featurenames)){
    rownames(coefs) <- c("(intercept)", object$Options$featurenames)
  }
  if(include_intercept) {
    return(coefs)} else return(coefs[-1,])

}

