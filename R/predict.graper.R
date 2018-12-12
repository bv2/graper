#' @title Predict response on new data
#' @name predict.graper
#' @description Function to predict the response on a
#' new data set using a fitted graper model.
#' @param object fitted graper model as obtained from  \code{\link{graper}}
#' @param newX Predictor matrix of size n_test
#'(number of new test samples) x p (number of predictors)
#' (same feature structure as used in \code{\link{graper}})
#' @param type type of prediction returned, either:
#' \itemize{
#'  \item{\strong{response}:}{returns the linear predictions
#'   for linear regression and class probabilities
#'   for logistic regression}
#'  \item{\strong{link}:}{returns the linear predictions}
#'  \item{\strong{inRange}:}{returns linear predictions for linear
#'   and class memberships for logistic regression}
#' }
#' @param ... other arguments
#' @importFrom methods is
#' @return A vector with predictions.
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#' ntrain <- dat$n / 2
#' fit <- graper(dat$X[seq_len(ntrain), ],
#'               dat$y[seq_len(ntrain)], dat$annot)
#' ypred <- predict(fit, dat$X[seq_len(ntrain) + dat$n / 2, ])
#'
#' dat <- makeExampleData(response="bernoulli")
#' ntrain <- dat$n / 2
#' fit <- graper(dat$X[seq_len(ntrain), ],
#'               dat$y[seq_len(ntrain)],
#'               dat$annot, family="binomial")
#' ypred <- predict(fit, dat$X[seq_len(ntrain) + dat$n / 2, ])

predict.graper <- function(object, newX,
                        type = c("inRange","response", "link"), ...){

    # sanity check
    if(ncol(newX) != nrow(object$EW_beta)) {
        stop("Number of columns in newX need to agree with number
            of predictors in the graper object.")
    }
    # sanity check
    if(!is(object, "graper")) {
        stop("object needs to be a graper object.")
    }
    # Get type of predictions wanted
    type = match.arg(type)

    if(is.null(object$intercept)) {
        object$intercept <- 0
    }
    if(object$Options$family == "gaussian"){
        pred <- object$intercept + newX %*% object$EW_beta
    } else {
        predexp <- object$intercept + newX %*%  object$EW_beta
        probs <- exp(predexp)/(1 + exp(predexp))
        if(type == "link") {
            pred <- predexp
        } else if(type == "response"){
            pred <- probs
        } else {
            pred <- round(probs)
        }
    }
    return(pred)
}
