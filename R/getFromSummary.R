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

#' @export

getCoefficients <- function(AllFits) {
    p <- AllFits$p
    coefmat <- sapply(AllFits$summaryList, function(l) {
        coef <- l$beta
        if (!is.null(coef)) 
            coef else rep(NA, p)
    })
    rownames(coefmat) <- AllFits$varnames
    return(coefmat)
}

#' @export

getIntercept <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        intercept <- l$intercept
        if (is.null(intercept)) 
            intercept <- NA
        intercept
    })
}

#' @export

getRunTime <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        runtime <- l$runtime
        if (is.null(runtime)) 
            runtime <- NA
        runtime
    })
}

#' @export

getRMSE <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        RMSE <- l$RMSE
        if (is.null(RMSE)) 
            RMSE <- NA
        RMSE
    })
}

#' @export

getBS <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        BS <- l$BS
        if (is.null(BS)) 
            BS <- NA
        BS
    })
}

#' @export

getROC <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        ROC <- l$ROC
        if (is.null(ROC)) 
            ROC <- NA
        ROC
    })
}

#' @export

getAUC <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        AUC <- l$AUC
        if (is.null(AUC)) 
            AUC <- NA
        AUC
    })
}


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

getFNR <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        FNR <- l$FNR
        if (is.null(FNR)) 
            FNR <- NA
        FNR
    })
}

#' @export

getFPR <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        fpr <- l$FPR
        if (is.null(fpr)) 
            fpr <- NA
        fpr
    })
}

#' @export

getl1error_beta <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        l1error_beta <- l$l1error_beta
        if (is.null(l1error_beta)) 
            l1error_beta <- NA
        l1error_beta
    })
}

#' @export

getl1error_intercept <- function(AllFits) {
    sapply(AllFits$summaryList, function(l) {
        intercept <- l$intercept
        if (is.null(intercept)) 
            intercept <- NA
        intercept
    })
}
