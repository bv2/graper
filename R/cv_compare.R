#' cv_compare
#'
#' Function to run serveral different methods for high-dimensional regression and evaluate them in a cross-validated fashion
#' @param X Design matrix of size n x p
#' @param y Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param family Liklihood model to use for response, either binomial or gaussian
#' @param nfolds Number of fold for evaluation
#' @param ncores Number of cores to use
#' @param plot_cv boolean whether to plot summary from evaluation
#' @param seed optional seed for the choice of folds
#' @param ... Other parameters that can be passed to RunMethods
#'
#' @return List of fitted models and two data frames with coeffcients and penalty factors
#' @import ggplot2
#' @import parallel
#' @export

cv_compare <- function(X, y, annot, family="gaussian",
                       ncores=1, nfolds=10, plot_cv=TRUE,
                       seed=NULL, parallel=FALSE, saveFits=FALSE,...){

  # split observations into folds
  if(!is.null(seed)) set.seed(seed)
  foldid <- sample(rep(seq(nfolds), length=nrow(X)))

  # function for each fold in cross-validation: train method on all excpet one fold, evaluate on the hold-out fold
  runPerFold <- function(foldidx){
    # split in train and test data
    use4test <- foldid==foldidx
    ytrain <- y[ !use4test]
    ytest <- y[use4test]
    Xtrain <- X[ !use4test,]
    Xtest <- X[use4test,]

    #fit models
    AllFits <- RunMethods(Xtrain = as.matrix(Xtrain),
                          ytrain =as.vector(ytrain),
                          annot = annot, ...)

    # save the fit
    if(saveFits) save(AllFits, file=paste0("AllFits",foldidx,".RData"))

    # extract relevant parameter from the model
    pf_mat <- getPenaltyFactors(AllFits)
    sparsity_mat <- getSparsityLevel(AllFits)
    beta_mat <- getCoefficients(AllFits)
    intercepts <- getIntercept(AllFits)

    # evaluate prediciton performance
    AllFits <- evaluateFits(AllFits, Xtest=as.matrix(Xtest), ytest=ytest)
    runtime <- getRunTime(AllFits)
    FNR <- getFNR(AllFits)
    FPR <- getFPR(AllFits)
    l1error_intercept <- getl1error_intercept(AllFits)
    l1error_beta <- getl1error_beta(AllFits)

    if (family=="gaussian"){
      RMSE <- getRMSE(AllFits)
      l <- list(FPR=FPR, FNR=FNR, RMSE=RMSE, pf_mat=pf_mat, beta_mat=beta_mat,
                intercepts=intercepts, sparsity_mat=sparsity_mat, annot=AllFits$annot, runtime=runtime,
                l1error_intercept, l1error_beta)
    } else if(family=="binomial"){
      BS <- getBS(AllFits)
      AUC <- getAUC(AllFits)
      ROC <- getROC(AllFits)
      l <- list(FPR=FPR, FNR=FNR, BS=BS, AUC=AUC, ROC=ROC, pf_mat=pf_mat, beta_mat=beta_mat,
                intercepts=intercepts, sparsity_mat=sparsity_mat, annot=AllFits$annot, runtime=runtime,
                l1error_intercept, l1error_beta)
    }
    else stop("Family not implemented")
    return(l)
  }


  if(parallel){
  resultList <- parallel::mclapply(1:nfolds, runPerFold, mc.cores = ncores)
  } else resultList <- lapply(1:nfolds, runPerFold)

  # plot results
  if(plot_cv) plotMethodComparison(resultList, family = family)

  return(resultList)
}


