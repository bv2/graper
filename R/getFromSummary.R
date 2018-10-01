# Functions to extract specific statistics from the results list produced by RunMethods

#' @title get penalty factors
#' @name getPenaltyFactors
#' @description get penalty factors from various methods
#' @param AllFits object as produced by RunMethods
#' @export

getPenaltyFactors <- function(AllFits) {
    G <- AllFits$G
    pfmat <- sapply(AllFits$summaryList, function(l) {
        pfs <- l$pf
        if (!is.null(pfs))
            pfs else rep(NA, G)
    })
    rownames(pfmat) <- unique(AllFits$groupnames)
    return(pfmat)
}

#' @title get sparsity levels
#' @name getSparsityLevel
#' @description get sparsity levels (1=dense, 0=sparse) from various methods
#' @param AllFits object as produced by RunMethods
#' @export
getSparsityLevel <- function(AllFits) {
    G <- AllFits$G
    sparsity_mat <- sapply(AllFits$summaryList, function(l) {
        sparsity <- l$sparsity
        if (!is.null(sparsity))
            sparsity else rep(NA, G)
    })
    rownames(sparsity_mat) <- unique(AllFits$groupnames)
    return(sparsity_mat)
}

#' @title get coefficients
#' @name getCoefficients
#' @description get coefficients estimated from various methods
#' @param AllFits object as produced by RunMethods
#' @export

getCoefficients <- function(AllFits) {
    p <- AllFits$p
    coefmat <- sapply(AllFits$summaryList, function(l) {
        coef <- l$beta
        if (!is.null(coef))
            as.numeric(coef) else rep(NA, p)
    })
    rownames(coefmat) <- AllFits$varnames
    return(coefmat)
}

#' @title get intercept
#' @name getIntercept
#' @description get intercept estimated from various methods
#' @param AllFits object as produced by RunMethods
#' @export
getIntercept <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        intercept <- l$intercept
        if (is.null(intercept))
            intercept <- NA
        intercept
    })
}

#' @title get run times
#' @name getRunTime
#' @description get run times of various methods
#' @param AllFits object as produced by RunMethods
#' @export

getRunTime <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        runtime <- l$runtime
        if (is.null(runtime))
            runtime <- NA
        runtime
    })
}

#' @title get RMSE
#' @name getRMSE
#' @description get RMSE of various methods
#' @param AllFits object as produced by RunMethods
#' @export

getRMSE <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        RMSE <- l$RMSE
        if (is.null(RMSE))
            RMSE <- NA
        RMSE
    })
}

#' @title get Brier Score
#' @name getBS
#' @description get Brier Score of various methods
#' @param AllFits object as produced by RunMethods
#' @export

getBS <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        BS <- l$BS
        if (is.null(BS))
            BS <- NA
        BS
    })
}



#' @title get AUC
#' @name getAUC
#' @description get AUC of various methods
#' @param AllFits object as produced by RunMethods
#' @export

getAUC <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        AUC <- l$AUC
        if (is.null(AUC))
            AUC <- NA
        AUC
    })
}


#' @title get MSE
#' @name getMSE
#' @description get MSE of various methods
#' @param AllFits object as produced by RunMethods
#' @export

getMSE <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        MSE <- l$RMSE^2
        if (is.null(MSE))
            MSE <- NA
        MSE
    })
}

#' @export

#' @title get FNR
#' @name getFNR
#' @description get false negative rates of various methods
#' @param AllFits object as produced by RunMethods
#' @export
#'
getFNR <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        FNR <- l$FNR
        if (is.null(FNR))
            FNR <- NA
        FNR
    })
}


#' @title get FPR
#' @name getFPR
#' @description get false positive rates of various methods
#' @param AllFits object as produced by RunMethods
#' @export

getFPR <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        fpr <- l$FPR
        if (is.null(fpr))
            fpr <- NA
        fpr
    })
}

#' @title get l1 error on beta
#' @name getl1error_beta
#' @description get absolute error on beta
#' @param AllFits object as produced by RunMethods
#' @export
#'
getl1error_beta <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        l1error_beta <- l$l1error_beta
        if (is.null(l1error_beta))
            l1error_beta <- NA
        l1error_beta
    })
}

#' @title get l1 error on intercept
#' @name getl1error_intercept
#' @description get absolute error on intercept
#' @param AllFits object as produced by RunMethods
#' @export

getl1error_intercept <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        intercept <- l$intercept
        if (is.null(intercept))
            intercept <- NA
        intercept
    })
}
