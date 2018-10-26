# Functions to extract specific statistics from the results list produced by \code{\link{runMethods}}

#' @title get penalty factors
#' @name getPenaltyFactors
#' @description get penalty factors from various methods
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export

getPenaltyFactors <- function(AllFits) {
    G <- AllFits$G
    pfmat <- vapply(AllFits$summaryList, function(l) {
        pfs <- l$pf
        if (!is.null(pfs))
            as.numeric(pfs) else rep(NA, G)
    }, numeric(G))
    rownames(pfmat) <- unique(AllFits$groupnames)
    return(pfmat)
}

#' @title get sparsity levels
#' @name getSparsityLevel
#' @description get sparsity levels (1=dense, 0=sparse) from various methods
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
getSparsityLevel <- function(AllFits) {
    G <- AllFits$G
    sparsity_mat <- vapply(AllFits$summaryList, function(l) {
        sparsity <- l$sparsity
        if (!is.null(sparsity))
            as.numeric(sparsity) else rep(NA, G)
    }, numeric(G))
    rownames(sparsity_mat) <- unique(AllFits$groupnames)
    return(sparsity_mat)
}

#' @title get coefficients
#' @name getCoefficients
#' @description get coefficients estimated from various methods
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export

getCoefficients <- function(AllFits) {
    p <- AllFits$p
    coefmat <- vapply(AllFits$summaryList, function(l) {
        coef <- l$beta
        if (!is.null(coef))
            as.numeric(coef) else rep(NA, p)
    }, numeric(p))
    rownames(coefmat) <- AllFits$varnames
    return(coefmat)
}

#' @title get intercept
#' @name getIntercept
#' @description get intercept estimated from various methods
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export
getIntercept <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        intercept <- l$intercept
        if (is.null(intercept))
            intercept <- NA
        intercept
    }, numeric(1))
}

#' @title get run times
#' @name getRunTime
#' @description get run times of various methods
#' @param AllFits object as produced by \code{\link{runMethods}}
#' @export

getRunTime <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        runtime <- l$runtime
        if (is.null(runtime))
            runtime <- NA
        runtime
    }, numeric(1))
}

#' @title get RMSE
#' @name getRMSE
#' @description get RMSE of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export

getRMSE <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        RMSE <- l$RMSE
        if (is.null(RMSE))
            RMSE <- NA
        RMSE
    }, numeric(1))
}

#' @title get Brier Score
#' @name getBS
#' @description get Brier Score of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export

getBS <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        BS <- l$BS
        if (is.null(BS))
            BS <- NA
        BS
    }, numeric(1))
}



#' @title get AUC
#' @name getAUC
#' @description get AUC of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export

getAUC <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        AUC <- l$AUC
        if (is.null(AUC))
            AUC <- NA
        AUC
    }, numeric(1))
}


#' @title get MSE
#' @name getMSE
#' @description get MSE of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export

getMSE <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        MSE <- l$RMSE^2
        if (is.null(MSE))
            MSE <- NA
        MSE
    }, numeric(1))
}

#' @export

#' @title get FNR
#' @name getFNR
#' @description get false negative rates of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#'
getFNR <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        FNR <- l$FNR
        if (is.null(FNR))
            FNR <- NA
        FNR
    }, numeric(1))
}


#' @title get FPR
#' @name getFPR
#' @description get false positive rates of various methods
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export

getFPR <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        fpr <- l$FPR
        if (is.null(fpr))
            fpr <- NA
        fpr
    }, numeric(1))
}

#' @title get l1 error on beta
#' @name getl1error_beta
#' @description get absolute error on beta
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export
#'
getl1error_beta <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
        l1error_beta <- l$l1error_beta
        if (is.null(l1error_beta))
            l1error_beta <- NA
        l1error_beta
    }, numeric(1))
}

#' @title get l1 error on intercept
#' @name getl1error_intercept
#' @description get absolute error on intercept
#' @param AllFits object as produced by \code{\link{evaluateFits}}
#' @export

getl1error_intercept <- function(AllFits) {
    vapply(AllFits$summaryList, function(l) {
      l1error_intercept <- l$l1error_intercept
        if (is.null(l1error_intercept))
          l1error_intercept <- NA
        l1error_intercept
    }, numeric(1))
}
