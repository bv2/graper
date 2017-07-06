
getPenaltyFactors <- function(AllFits) {
  G <- AllFits$G
  pfmat <- sapply(AllFits$summaryList, function(l) {
        pfs <- l$pf
        if(!is.null(pfs)) pfs else rep(NA,G)
      })
  rownames(pfmat) <- unique(AllFits$groupnames)
  return(pfmat)
}

getSparsityLevel <- function(AllFits) {
  G <- AllFits$G
  sparsity_mat <- sapply(AllFits$summaryList, function(l) {
    sparsity <- l$sparsity
    if(!is.null(sparsity)) sparsity else rep(NA,G)
  })
  rownames(sparsity_mat) <- unique(AllFits$groupnames)
  return(sparsity_mat)
}


getCoefficients <- function(AllFits) {
  p <- AllFits$p
  coefmat <- sapply(AllFits$summaryList, function(l) {
    coef <- l$beta
    if(!is.null(coef)) coef else rep(NA,p)
  })
  rownames(coefmat) <- AllFits$varnames
  return(coefmat)
}

getIntercept <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    intercept <- l$intercept
    if(is.null(intercept)) intercept <- NA
    intercept
  })
}

getRunTime <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    runtime <- l$runtime
    if(is.null(runtime)) runtime <- NA
    runtime
  })
}

getRMSE <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    RMSE <- l$RMSE
    if(is.null(RMSE)) RMSE <- NA
    RMSE
  })
}

getMSE <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    MSE <- l$RMSE^2
    if(is.null(MSE)) MSE <- NA
    MSE
  })
}

getFNR <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    FNR <- l$FNR
    if(is.null(FNR)) FNR <- NA
    FNR
  })
}

getFPR <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    fpr <- l$FPR
    if(is.null(fpr)) fpr <- NA
    fpr
  })
}

getl1error_beta <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    l1error_beta <- l$l1error_beta
    if(is.null(l1error_beta)) l1error_beta <- NA
    l1error_beta
  })
}

getl1error_intercept <- function(AllFits) {
  sapply(AllFits$summaryList, function(l) {
    intercept <- l$intercept
    if(is.null(intercept)) intercept <- NA
    intercept
  })
}
