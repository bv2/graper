# utils This file contains simple helper functions
# - MarginalCoefficient: Calculate Marginal coefficeints in a univariate GLM -
# ImputeByMean: Impute missing values of a vector by its mean

#' @title Get marginal coefficient from a (generalized) linear model
#' @name MarginalCoefficient
#' @description Function to compute marignal regression coefficients
#' @param response response vector
#' @param data design matrix
#' @param family liklihood model for the response
#' @export
#' @return a matrix containing estimates, standard errors, statistics
#'  and p-values (rows) for each feature (columns) in the desing matrix
#' @importFrom stats glm
#' @examples
#' dat <- makeExampleData(response="bernoulli")
#' MarginalCoefficient(dat$y, dat$X, family= "binomial")
MarginalCoefficient <- function(response, data, family = "gaussian") {
    apply(data, 2, function(c) {
        lm.fit <- stats::glm(response ~ c, family = family)
        s <- summary(lm.fit)
        s$coefficients[2, ]
    })
}

#' @title Impute by mean
#' @name ImputeByMean
#' @description Function to impute by mean
#' @param x vector to with value to impute
#' @export
#' @return vector with imputed values
#' @examples
#' x <- c(1,0.4,1.3, NA)
#' ImputeByMean(x)
ImputeByMean <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
}

#' @title Evaluate regression model on test data
#' @name EvaluateModel
#' @description Function to evaluate model fits on test data
#' @param beta_est Estimated model coefficients
#' @param intercept Estiamted intercept
#' @param X_test design matrix test data of size
#'  number of test samples times features (used for fitting)
#' @param y_test response vectorof length number of  test samples
#' @param beta0 true model coefficients (if known)
#' @param family liklihood model for the response, either
#'  "gaussian" for linear regression or "binomial" for logisitc regression
#' @importFrom GRridge auc roc
#' @export
#' @return A list containting various performance measures such as RMSE, FNR, FPR etc
#' @examples
#' dat <- makeExampleData()
#' # take a null model with all coefficeints set to zero for evaluation
#' beta <- rep(0, dat$p)
#' eval.out <- EvaluateModel(beta, intercept = 0, X_test = dat$X,
#'  y_test = dat$y, beta0 = dat$beta, family = "gaussian")

EvaluateModel <- function(beta_est, intercept, X_test,
                          y_test, beta0 = NULL, family) {
    if (is.null(intercept))
        intercept <- 0
    intercept <- as.numeric(intercept)
    p <- ncol(X_test)
    n_test <- nrow(X_test)
    stopifnot(family %in% c("gaussian", "binomial"))
    RMSE_test <- FNR <- sensitivity <- specificity <- FPR <- BrierScore <- NULL
    ROC <- AUC <- predprob <- pred_gauss <- precision <- recall <- F1score <- NULL

    # if 'true features' known evaluate sensitivity and specificity
    if (!is.null(beta0)) {
        stopifnot(length(beta0) == length(beta_est))
        # Feature selection performance
        FPR <- sum(beta0 == 0 & beta_est != 0)/sum(beta0 == 0)
        specificity <- 1 - FPR
        FNR <- sum(beta0 != 0 & beta_est == 0)/sum(beta0 != 0)
        sensitivity <- 1 - FNR
        precision <- sum(beta0 != 0 & beta_est != 0)/sum(beta_est != 0)
        recall <- sum(beta0 != 0 & beta_est != 0)/sum(beta0 != 0)
        F1score <- 2*precision*recall/(precision + recall)
    } else NULL

    # Prediction performance on test set
    if (family == "gaussian") {
        # RMSE
        pred_gauss <- intercept + X_test %*% beta_est
        RMSE_test <- sqrt(sum((y_test - pred_gauss)^2)/length(y_test))
    } else if (family == "binomial") {
        # Brier score
      # to avoid NaN for too large numbers
        predexp <- pmin(intercept + X_test %*% beta_est, 500)
        predprob <- exp(predexp)/(1 + exp(predexp))
        BrierScore <- sum((y_test - predprob)^2)/length(y_test)

        # ROC (function from grridge used)
        cutoffs <- rev(seq(0, 1, by = 0.01))
        ROC <- GRridge::roc(predprob, y_test, cutoffs)

        # AUC (function from grridge used)
        AUC <- GRridge::auc(ROC)
    }


    return(list(beta_est = beta_est, RMSE_test = RMSE_test, BrierScore = BrierScore,
                ROC = ROC, AUC = AUC, specificity = specificity,
                sensitivity = sensitivity, FNR = FNR, FPR = FPR, precision=precision,
                recall=recall, F1score=F1score, beta0 = beta0, y_test = y_test,
                predprob = predprob, pred_gauss = pred_gauss, intercept = intercept))
}
