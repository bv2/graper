# evaluation_functions
#
# This file contains functions used for comparison to other methods and evaluation of the produced fit
#         - RunMethods: Run a number of different methods for estimation in the model including
#                       - Bayes Linear Regression with groupwise prior
#                       - Bayes Linear Regression with groupwise prior and fully-factorized inference
#                       - Ridge
#                       - Lasso
#                       - Group Lasso
#         - evalResult: Compare Estimates produced in RunMethods by MSE and estimation error
#
# ---------------------------
#' RunMethods
#'
#' Function to run serveral different methods for high-dimensional regression
#' @param Xtrain Design matrix of size n x p
#' @param ytrain Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param beta0 True coefficients in the linear model if known, NULL otherwise (default)
#' @param trueintercept True intercept in the linear model if known, NULL otherwise (default)
#' @param max_iter maximum number of iterations
#' @param intercept boolean, indicating wether to fit an intercept
#' @param plotit boolean, indicating wether to produce diagnositc plots
#' @param standardize boolean, indicating wether features for Ridge and Lasso fit should be standardized. Note this does not affect GRridge and grouplasso where standardization is default.
#' @param calcELB boolean, indicating wether to calculate ELB
#' @param verbose boolean, indicating wether to print out intermediate messages during fitting
#' @param compareGRridge boolean, indicating wether to fit a GRridge model (might cause errors, look for stable version!)
#' @param freqELB determines frequency at which ELB is to be calculated, i.e. each feqELB-th iteration
#' @return List of fitted models and data frame DFresult with coeffcients and penalty factors
#' @import ggplot2
#' @import glmnet
#' @import gridExtra
#' @import GRridge

#ToDo: Option to standardize the input and then re-adjust estimates coefficeints to not standardized
RunMethods<-function(Xtrain, ytrain, annot, beta0=NULL, trueintercept=NULL, max_iter=500, intercept=T, plotit=T, standardize=F, verbose=F, compareGRridge=F,
                     freqELB=10, calcELB=T, include_nonfacQ=T, family="gaussian", constantXcol=F, compareGroupLasso=T, includeRF=T, th=1e-7){

  if(!standardize) warning("Group Lasso and Grridge are standardized despite of standardized = F")

  if(intercept & family=="binomial") warning("Intercept not yet implemented for logistic regression")

  if(constantXcol) stopifnot(var(Xtrain[,1])==0)
  if(constantXcol) {
    if(intercept) {warning("Intercept not used, as constant X columns included") ;intercept<-F}
    penaltyFac<-c(0,rep(1,ncol(Xtrain)-1))
    } else penaltyFac<-rep(1,ncol(Xtrain))

  #turn annot to facotr
  annot<-as.factor(annot)

  #not fully factorized, normal prior
  if(include_nonfacQ){
  tmp<-Sys.time()
  result<-fit_grpRR(Xtrain, ytrain, annot=annot, factoriseQ=F, spikeslab= F, max_iter = max_iter, intercept=intercept,
                    verbose = verbose, freqELB = freqELB, calcELB = calcELB, family = family, th=th)
  timeNF<-Sys.time()-tmp
  result$runningTime<-timeNF
  #set intercept to zero as not yet implemented
  if(family=="binomial" & intercept) result$intercept<-0
  if(plotit) plotVBFit(result, whichParam = c("ELB", "tau", "gamma"))
  gamma_est<-result$EW_gamma[annot]
  beta_est<-result$EW_beta
  }

  #fully factorized, normal prior
  tmp<-Sys.time()
  resultFF<-fit_grpRR(Xtrain, ytrain, annot=annot, factoriseQ = T, spikeslab= F, max_iter = max_iter, intercept=intercept,
                      verbose = verbose, freqELB = freqELB, calcELB = calcELB, family = family, th=th)
  #set intercept to zero as not yet implemented
  if(family=="binomial" & intercept) resultFF$intercept<-0
  timeFF<-Sys.time()-tmp
  resultFF$runningTime<-timeFF
  if(plotit) plotVBFit(result, whichParam = c("ELB", "tau", "gamma"))
  gamma_estFF<-resultFF$EW_gamma[annot]
  beta_estFF<-resultFF$EW_beta

  #Only implemented for gaussian yet
  if(family=="gaussian"){
  #fully factorized, spike and slab
  tmp<-Sys.time()
  resultSS<-fit_grpRR(Xtrain, ytrain, annot=annot, factoriseQ = T, spikeslab= T, max_iter = max_iter, intercept=intercept,
                      verbose = verbose, freqELB = freqELB, calcELB = calcELB, th=th)
  timeSS<-Sys.time()-tmp
  resultSS$runningTime<-timeSS
  if(plotit) plotVBFit(result, whichParam = c("ELB", "tau", "gamma"))
  gamma_estSS<-resultSS$EW_gamma[annot]
  beta_estSS<-resultSS$EW_beta
  beta_estSScutoff<-ifelse(resultSS$EW_s<0.5,0,resultSS$EW_beta)
  pi_estSS<-resultSS$EW_pi[annot]
  includeSS<-T
  } else includeSS<-F


  #ridge regression
  RidgeFit<-cv.glmnet(Xtrain, ytrain, alpha=0, intercept=intercept, standardize=standardize, family=family, penalty.factor=penaltyFac)
  beta_ridge<-as.vector(coef(RidgeFit, RidgeFit$lambda.min))

  #Lasso
  LassoFit<-cv.glmnet(Xtrain, ytrain, alpha=1, intercept=intercept, standardize=standardize, family=family, penalty.factor=penaltyFac)
  beta_lasso<-as.vector(coef(LassoFit, LassoFit$lambda.min))

  #EN
  ENFit<-cv.glmnet(Xtrain, ytrain, alpha=0.2, intercept=intercept, standardize=standardize, family=family, penalty.factor=penaltyFac)
  beta_EN<-as.vector(coef(ENFit, ENFit$lambda.min))

  #Random Forest
  if(includeRF) rf.out<-randomForest(x = Xtrain, y = ytrain)

  # #Ridge with PF by average marginal coefficients...better use significance instead of effect size....
  # marg<-MarginalCoefficient(ytrain, scale(Xtrain), family = family)
  # avMargGroup<-sapply(annot, function(g) mean(abs(marg[1,annot==g])))
  # pf_margCoeff<-1/abs(avMargGroup)  #estimates
  # RidgeavMargFit<-cv.glmnet(Xtrain,ytrain, alpha=0, intercept=intercept, penalty.factor=pf_margCoeff,standardize=standardize, family=family)
  # beta_RidgeavMarg<-as.vector(coef(RidgeavMargFit$glmnet.fit, RidgeavMargFit$lambda.min))
  #
  # #Lasso with PF by average marginal coefficients
  # marg<-MarginalCoefficient(ytrain, Xtrain)
  # avMargGroup<-sapply(annot, function(g) mean(abs(marg[1,annot==g])))
  # pf_margCoeff<-1/abs(avMargGroup)  #estimates
  # LassoavMargFit<-cv.glmnet(Xtrain,ytrain, alpha=1, intercept=intercept, penalty.factor=pf_margCoeff,standardize=standardize, family=family)
  # beta_LassoavMarg<-as.vector(coef(LassoavMargFit$glmnet.fit, LassoavMargFit$lambda.min))

  #group lasso
  if(compareGroupLasso){
  GrpLassoFit<-grpreg::cv.grpreg(Xtrain, ytrain,
                                 group=as.factor(annot), penalty="grLasso",
                                 intercept=intercept, family=family)
  beta_grplasso<-coef(GrpLassoFit, GrpLassoFit$lambda.min)
  }

  #grridge
  if(compareGRridge){
    partition<-CreatePartition(as.factor(annot))
    if(intercept) MessagesGR<-capture.output(GRfit<-grridge(t(Xtrain), as.numeric(ytrain), list(partition), unpenal=~1))
    else MessagesGR<-capture.output(GRfit<-grridge(t(Xtrain), as.numeric(ytrain), list(partition), unpenal=~0))
    if(verbose) MessagesGR
    lambda_gr<-GRfit$lambdamultvec[,2]
    beta_gr<-GRfit$betas
  }

  #zero model
  beta_zeromodel<-rep(0,ncol(Xtrain))
  intercept_zeromodel<-mean(ytrain)

  groupnames<-as.character(annot)
  varnames<-colnames(Xtrain)
  if(length(varnames)==0) varnames<-factor(paste("Feautre",1:ncol(Xtrain), sep=""))

  if(intercept){
    beta_estFF<-c(resultFF$intercept,beta_estFF)
    gamma_estFF<-c(0,gamma_estFF)
    beta_zeromodel<-c(intercept_zeromodel,beta_zeromodel)
    if(include_nonfacQ){
      beta_est<-c(result$intercept,beta_est)
      gamma_est<-c(0,gamma_est)
    }
    if(compareGRridge){
      beta_gr<-c(GRfit$predobj$GroupRegul@unpenalized,beta_gr)
      lambda_gr<-c(0,lambda_gr)
    }
    if(includeSS){
      beta_estSS<-c(resultSS$intercept,beta_estSS)
      beta_estSScutoff<-c(resultSS$intercept,beta_estSScutoff)
      pi_estSS<-c(0,pi_estSS)
      gamma_estSS<-c(0,gamma_estSS)
    }
    # pf_margCoeff<-c(0,pf_margCoeff)
    groupnames<-c("intercept",groupnames)
    varnames<-c("intercept",varnames)
    if(!is.null(beta0)) beta0<-c(trueintercept,beta0)
  } else{
    # beta_RidgeavMarg<-beta_RidgeavMarg[-1]
    beta_ridge<-beta_ridge[-1]
    # beta_LassoavMarg<-beta_LassoavMarg[-1]
    beta_lasso<-beta_lasso[-1]
    beta_EN<-beta_EN[-1]
    if(compareGroupLasso) beta_grplasso<-beta_grplasso[-1]
  }

  #result dataframe

  DFresult<-data.frame(beta_estFF=beta_estFF,
                       beta_ridge=beta_ridge,
                       beta_lasso=beta_lasso,
                       beta_EN=beta_EN,
                       # beta_LassoavMarg=beta_LassoavMarg,
                       # beta_RidgeavMarg=beta_RidgeavMarg,
                       beta_zeromodel=beta_zeromodel,
                       gamma_estFF=gamma_estFF,
                       # avMarg=pf_margCoeff,
                       id= varnames,
                       groupMembership=groupnames)

  if(compareGroupLasso) DFresult$beta_grplasso=beta_grplasso
  if(compareGRridge){
    DFresult$beta_gr=beta_gr
    DFresult$lambda_gr=lambda_gr
  }
  if(include_nonfacQ){
    DFresult$beta_est=beta_est
    DFresult$gamma_est=gamma_est
  }
  if(includeSS){
    DFresult$beta_estSS=beta_estSS
    DFresult$beta_estSScutoff=beta_estSScutoff
    DFresult$pi_estSS=pi_estSS
    DFresult$gamma_estSS=gamma_estSS
  }

  if(!is.null(beta0))  {
    DFresult$beta_true<-beta0
    if(!intercept) DFresult$sparsity_true<-sapply(unique(annot), function(gr) sum(beta0[annot==gr]==0))[annot]
    else DFresult$sparsity_true<-c(0,sapply(unique(annot), function(gr) sum(beta0[-1][annot==gr]==0)/sum(annot==gr))[annot])
  }


  listOfFits<-list(DFresult=DFresult,
                    ENFit=ENFit,
                   # LassoavMargFit=LassoavMargFit,
                   LassoFit=LassoFit,
                   # RidgeavMargFit=RidgeavMargFit,
                   RidgeFit=RidgeFit,
                   resultFF=resultFF,
                   intercept=intercept,
                   beta0=beta0,
                   beta_zeromodel=beta_zeromodel)

  if(compareGroupLasso){
    listOfFits$GrpLassoFit=GrpLassoFit
  }
  if(compareGRridge){
    listOfFits$GRfit=GRfit
  }
  if(include_nonfacQ){
    listOfFits$result=result
  }
  if(includeSS){
    listOfFits$resultSS=resultSS
  }
  if(includeRF){
    listOfFits$RFout=rf.out
  }
  return(listOfFits)

}


# ---------------------------
#'  evalResult
#'
#' Function to evaluate results (MSE, Peanlty Factors...)
#' @param allFits List as produced by \code{\link{runMethods}}
#' @param Xtest Design matrix of size n' x p (same feature structure as used in runMethods)
#' @param ytest Response vector of size n'
#' @param plotit Boolean, indication wether to produce plots of MSE and penalty factors
#' @return Datframe containing error measures per method

evalResult<-function(allFits, Xtest, ytest, plotit=T, family="gaussian", plotbeta=F, saveit=F, filenm=""){

  DFresult<-allFits$DFresult
  UseIntercept<-allFits$intercept
  beta0<-allFits$beta0

  # Plot group-wise penalties
  DFGroupPenalties<-data.frame(gamma_dense_ff=DFresult$gamma_estFF,
                               # avMarg=DFresult$avMarg,
                               group=factor((DFresult$groupMembership), level=unique(DFresult$groupMembership)))

  if(!is.null(beta0)){
    DFGroupPenalties$beta0=beta0
    DFGroupPenalties$sparsity_true=DFresult$sparsity_true
  }
  if("GRfit" %in% names(allFits)) DFGroupPenalties$gamma_GR=DFresult$lambda_gr
  if("result" %in% names(allFits)) DFGroupPenalties$gamma_dense=DFresult$gamma_est
  if("resultSS" %in% names(allFits)) {
    DFGroupPenalties$gamma_sparse=DFresult$gamma_estSS
    DFGroupPenalties$pi_sparse=DFresult$pi_estSS
  }


  dfgp <- melt(DFGroupPenalties,id.vars = c("group"))
  ggrelPenal<-ggplot(dfgp, aes(x=group, y=value, group=variable, fill=variable))+geom_bar(stat="summary", position="dodge", fun.y="mean")+
    ggtitle("Penalty factor per group") +facet_wrap(~variable, scales = "free_y")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  if(plotit) print(ggrelPenal)

  #Extract estimatedmodel coefficinets
  if(is.null(beta0)) DFresult<-DFresult[,!grepl("beta_true", colnames(DFresult))]
  betaDF<-DFresult[,grepl("beta", colnames(DFresult))]


  #For gaussian family calculate MSE as measure of prediciton performance
  if(family=="gaussian"){
    if(UseIntercept) MSE<-apply(betaDF, 2, function(beta) 1/length(ytest)*sum((cbind(1,Xtest)%*%beta- as.vector(ytest))^2))
    else MSE<-apply(betaDF, 2, function(beta) 1/length(ytest)*sum((Xtest%*%beta- as.vector(ytest))^2))
    if("RFout" %in% names(allFits)) MSE<-c(MSE, RF=1/length(ytest)*sum((predict(allFits$RFout, Xtest)-ytest)^2))

    EvalDF<-data.frame(method= names(MSE), MSE=MSE)
    ggMSE<-ggplot(EvalDF, aes(x=method, y=MSE, fill=method))+geom_bar(stat="identity")+ggtitle("MSE")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    if(plotit) print(ggMSE)
  }
  #For binomial family calculate ROC, aUC and Brier Score
  else if(family=="binomial"){
    #ROC curve
    if(UseIntercept) {
      interceptDF = betaDF[1,]
      betaDF_withoutIntercept = betaDF[-1,]
    } else{
      interceptDF=rep(0,ncol(betaDF))
      betaDF_withoutIntercept = betaDF
    }

    ROC<-lapply(1:ncol(betaDF),function( idx) EvaluateModel(betaDF_withoutIntercept[,idx], intercept=interceptDF[idx],
                                                            Xtest,ytest,
                                                            beta0=beta0, family="binomial")$ROC)
    if(plotit){
    colors4method<-rainbow(ncol(betaDF))
    names(colors4method)<-colnames(betaDF)
    plot(NA,xlim=c(0,1), ylim=c(0,1), xlab="FPR", ylab="TPR")
    for(idx in 1:ncol(betaDF)) lines(ROC[[idx]][1,], ROC[[idx]][2,], col=colors4method[colnames(betaDF)[idx]])
    legend(x=0.4,y=0.3, legend=names(colors4method), fill=colors4method)
    }
    #AUC
    AUC<-sapply(1:ncol(betaDF),function( idx) EvaluateModel(betaDF_withoutIntercept[,idx], intercept=interceptDF[idx],
                                                            Xtest,ytest,
                                                           beta0=beta0, family="binomial")$AUC)

    #Test prediction performance (Brier Score)
    BS<-sapply(1:ncol(betaDF),function( idx) EvaluateModel(betaDF_withoutIntercept[,idx], intercept=interceptDF[idx],
                                                           Xtest,ytest,
                                                          beta0=beta0, family="binomial")$BrierScore)
    EvalDF<-data.frame(method= colnames(betaDF), AUC=AUC, BS=BS)
    ggAUC<-ggplot(EvalDF, aes(x=method, y=AUC, fill=method))+geom_bar(stat="identity")+ggtitle("AUC")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    ggBS<-ggplot(EvalDF, aes(x=method, y=BS, fill=method))+geom_bar(stat="identity")+ggtitle("Brier Score")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    if(plotit) grid.arrange(ggAUC, ggBS, ncol=2)
  }
  #No other families implemented
  else("Family not known")


  #evaluate estimates
  if(!is.null(beta0)){
    #l1-Error in estimation of coeffcients
    L1DiffBeta<-apply(betaDF, 2, function(beta) sum(abs(beta[-1]- betaDF$beta_true[-1])))
    #l1-Error in estimation of intercept
    if(UseIntercept) InterceptDiff<-apply(betaDF, 2, function(beta) sum(abs(beta[1]- betaDF$beta_true[1]))) else InterceptDiff<-NA

    ggL1<-ggplot(EvalDF, aes(x=method, y=L1DiffBeta, fill=method))+geom_bar(stat="identity")+ggtitle("L1 error in estimated coeffcients")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    if(UseIntercept)ggInterceptDiff<-ggplot(EvalDF, aes(x=method, y=InterceptDiff, fill=method))+geom_bar(stat="identity")+ggtitle("abs error in estimated intercept")+
      theme(axis.text.x = element_text(angle = 60, hjust = 1))

     if("RFout" %in% names(allFits)) {
       #at the moment only works if RF is last
       stopifnot(EvalDF$method[nrow(EvalDF)]=="RF")
       L1DiffBeta<-c(L1DiffBeta,NA)
       InterceptDiff<-c(InterceptDiff,NA)
       }

    EvalDF<-cbind(EvalDF, L1DiffBeta=L1DiffBeta,
                       InterceptDiff=InterceptDiff)

  if(plotit){
    if(UseIntercept) grid.arrange(ggL1,ggInterceptDiff, ncol=2)
    else print(ggL1)
  }
  }

  betaDF$group<-allFits$DFresult$groupMembership
  if(UseIntercept) betaDF$Feature<-allFits$DFresult$id
  if(UseIntercept) betaDF$FeatureNo<-0:ncol(Xtest) else betaDF$FeatureNo<-1:ncol(Xtest)
  moltenbetaDF<-melt(betaDF[betaDF$Feature!="intercept",], id.vars = c("Feature", "group", "FeatureNo"))
  ggbeta<-ggplot(moltenbetaDF, aes(x=as.numeric(FeatureNo), y=value, fill=group))+geom_bar(stat="identity", position = "dodge")+facet_wrap(~variable, scales = "free_y")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))


  if(plotbeta){
    print(ggbeta)
  }
  if(saveit) save(list(EvalDF=EvalDF, DFGroupPenalties=dfgp, betaDF=moltenbetaDF), file = paste(filenm,"_evalgrpRR.RData", sep=""))

  return(list(EvalDF=EvalDF, DFGroupPenalties=dfgp, betaDF=moltenbetaDF))
  }







