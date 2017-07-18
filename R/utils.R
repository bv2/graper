# utils
#
# This file contains simple helper functions
#         - MarginalCoefficient: Calculate Marginal coefficeints in a univariate GLM
#         - ImputeByMean: Impute missing values of a vector by its mean

#' MarginalCoefficient
#'
#' Function to compute marignal regression coefficients
#' @export 
MarginalCoefficient<-function(response,data, family="gaussian"){
  apply(data,2, function (c) {
    lm.fit<-glm(response~c, family=family)
    s<-summary(lm.fit)
    #print(s)
    s$coefficients[2,]
  })
}

#' ImputeByMean
#'
#' Function to impute by mean
#' @export 
ImputeByMean <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}

#' EvaluateModel
#'
#' Function to evaluate model fits
#' @export 
EvaluateModel<-function(beta_est, intercept, X_test,y_test, beta0=NULL,family) {
  if(is.null(intercept)) intercept <- 0
  intercept <- as.numeric(intercept)
  p<-ncol(X_test)
  n_test<-nrow(X_test)
  stopifnot(family %in% c("gaussian", "binomial"))
  RMSE_test<-FNR <- sensitivity <- specificity <- FPR <-BrierScore<-ROC<-AUC<-predprob<-pred_gauss<-NULL

  #if "true features" known evaluate sensitivity and specificity
  if(!is.null(beta0)){
    stopifnot(length(beta0)==length(beta_est))
    #Feature selection performance
    FPR <- sum(beta0==0 & beta_est!=0)/sum(beta0==0)
    specificity<-1-FPR
    FNR<- sum(beta0!=0 & beta_est==0)/sum(beta0!=0)
    sensitivity<-1-FNR
  } else  NULL

  #Prediction performance on test set
  if(family=="gaussian")  {
    #RMSE
    pred_gauss<-intercept+X_test%*%beta_est
    RMSE_test<-sqrt(sum((y_test-pred_gauss)^2)/length(y_test))
  } else if(family=="binomial") {
    #Brier score
    predexp<-pmin(intercept+X_test%*%beta_est, 500) # to avoid NaN for too large numbers
    predprob<-exp(predexp)/(1+exp(predexp))
    BrierScore<- sum((y_test- predprob)^2)/length(y_test)

    #ROC (function from grridge used)
    cutoffs <- rev(seq(0,1,by=0.01))
    ROC<-GRridge::roc(predprob, y_test, cutoffs)

    #AUC (function from grridge used)
    AUC<-GRridge::auc(ROC)
  }


  return(list(beta_est=beta_est,
              RMSE_test=RMSE_test,
              BrierScore=BrierScore,
              ROC=ROC,
              AUC=AUC,
              specificity=specificity,
              sensitivity=sensitivity,
              FNR=FNR,
              FPR=FPR,
              beta0=beta0,
              y_test=y_test,
              predprob=predprob,
              pred_gauss=pred_gauss,
              intercept=intercept
  ))
}
