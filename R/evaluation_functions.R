#' @title Run various regression methods
#' @name runMethods
#' @description Function to run serveral different methods for high-dimensional regression
#' @param Xtrain Design matrix of size n x p
#' @param ytrain Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param beta0 True coefficients in the linear model if known, NULL otherwise (default)
#' @param trueintercept True intercept in the linear model if known, NULL otherwise (default)
#' @param max_iter maximum number of iterations
#' @param family liklihood model for the response, either "gaussian" for linear regression
#'  or "binomial" for logisitc regression
#' @param intercept boolean, indicating wether to fit an intercept
#' @param standardize boolean, indicating wether features for Ridge and Lasso fit should be standardized.
#'  Note this does not affect GRridge and grouplasso where standardization is default.
#' @param calcELB boolean, indicating wether to calculate ELB
#' @param freqELB determines frequency at which ELB is to be calculated, i.e. each feqELB-th iteration
#' @param th convergence threshold on the ELBO for the variational Bayes algorithm
#' @param n_rep number of repeated fits with variational Bayes using different random intilizations
#' @param verbose boolean, indicating wether to print out intermediate messages during fitting
#' @param compareGRridge boolean, indicating wether to fit a GRridge model
#' @param include_nonfacQ include a VB method with multivariate variational distributon
#'  (can be very timme consuming for large data sets)
#' @param compareIPF boolean, indicating whether to fit a IPFLasso
#' @param compareSparseGroupLasso boolean, indicating whether to fit a sparse group lasso
#' @param compareGroupLasso boolean, indicating whether to fit a group lasso
#' @param compareAdaLasso boolean, indicating whether to fit an adpative lasso
#' @param includeRF boolean, indicating whether to fit a random forest
#' @param include_varbvs boolean, indicating whether to fit varbvs
#' @param include_nogamma boolean, indicating whether to fit a grpRR model without different slab parameters
#' @param include_grpRR_SS_ungrouped boolean, indicating whether to fit a grpRR model without group annotations
#' @param verbose_progress boolean, indicating whether to print details on the progress
#' @return List of fitted models and two data frames with coeffcients and penalty factors
#' @importFrom glmnet cv.glmnet
#' @importFrom varbvs varbvs
#' @importFrom randomForest randomForest
#' @importFrom grpreg cv.grpreg
#' @importFrom SGL cvSGL
#' @importFrom GRridge CreatePartition grridge
#' @importFrom ipflasso cvr2.ipflasso
#' @importFrom stats coef
#' @return a list with
#'  - a summaryList containing the complete individual fits (out) as well as other
#' statistics (coefficients, runtime, sparsity, intercept, penalty factors)
#' - the details on the data (sample size n, predictor number p, covariate annotation annot,
#'  group number G, feature and annotation names)
#' @export

runMethods <- function(Xtrain, ytrain, annot, beta0 = NULL, trueintercept = NULL, max_iter = 5000,
                       family = "gaussian", intercept = TRUE, standardize = TRUE,
                       freqELB = 10, calcELB = TRUE, th = 0.01,
                       n_rep=1, verbose = FALSE, verbose_progress = TRUE,
                       compareGRridge = FALSE, include_nonfacQ = FALSE,
                       compareSparseGroupLasso =TRUE, compareIPF = TRUE,
                       compareGroupLasso = TRUE, includeRF = T,
                       compareAdaLasso = TRUE, include_varbvs=FALSE,
                       include_nogamma=FALSE, include_grpRR_SS_ungrouped= FALSE) {

  if (!standardize)
    warning("Group Lasso and GRridge are alwyas standardized, despite of standardized = FALSE.")

  # sanity checks
  stopifnot(nrow(Xtrain) == length(ytrain))
  if (!is.null(beta0)) stopifnot(ncol(Xtrain) == length(beta0))

  # turn annot to factor
  annot <- as.factor(annot)

  # extract important parameters and names
  p <- ncol(Xtrain)
  n <- nrow(Xtrain)
  G <- length(unique(annot))
  groupnames <- as.character(annot)
  varnames <- colnames(Xtrain)
  if (length(varnames) == 0)
    varnames <- factor(paste("Feature", 1:ncol(Xtrain), sep = ""))

  #### RUN DIFFERENT METHODS ####
  summaryList <- list()

  # grpRR: multivariate variational distribution, normal prior (dense)
  if (include_nonfacQ) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting grpRR model...")
    grpRR <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = F, spikeslab = F,
                       max_iter = max_iter, intercept = intercept,
                       verbose = verbose, freqELB = freqELB, calcELB = calcELB,
                       family = family, th = th, standardize=standardize, n_rep=n_rep)
    timeNF <- difftime(Sys.time(), tmp, units = "secs")
    grpRR_summary <- list()
    grpRR_summary$runtime <- as.numeric(timeNF)
    grpRR_summary$pf <- as.numeric(grpRR$EW_gamma)
    grpRR_summary$beta <- grpRR$EW_beta
    grpRR_summary$intercept <- grpRR$intercept
    grpRR_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
    grpRR_summary$out <- grpRR
    rm(grpRR)
    summaryList$grpRR <- grpRR_summary
  }

  # grpRR_FF : fully factorized variational distribution, normal prior (dense)
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting grpRR_FF model...")
  grpRR_FF <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = T, spikeslab = F,
                        max_iter = max_iter, intercept = intercept,
                        verbose = verbose, freqELB = freqELB, calcELB = calcELB,
                        family = family, th = th, standardize=standardize, n_rep=n_rep)
  timeFF <- difftime(Sys.time(), tmp, units = "secs")
  grpRR_FF_summary <- list()
  grpRR_FF_summary$runtime <- as.numeric(timeFF)
  grpRR_FF_summary$pf <- as.numeric(grpRR_FF$EW_gamma)
  grpRR_FF_summary$beta <- grpRR_FF$EW_beta
  grpRR_FF_summary$intercept <- grpRR_FF$intercept
  grpRR_FF_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
  grpRR_FF_summary$out <- grpRR_FF
  rm(grpRR_FF)
  summaryList$grpRR_FF <- grpRR_FF_summary

  # grpRR_SS: fully factorized variational distribution, spike and slab prior (sparse)
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting grpRR_SS model...")
  grpRR_SS <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = T, spikeslab = T,
                        max_iter = max_iter, intercept = intercept,
                        verbose = verbose, freqELB = freqELB, calcELB = calcELB,
                        th = th, family = family, standardize=standardize, n_rep=n_rep)
  timeSS <- difftime(Sys.time(), tmp, units = "secs")
  grpRR_SS_summary <- list()
  grpRR_SS_summary$runtime <- as.numeric(timeSS)
  grpRR_SS_summary$pf <- as.numeric(grpRR_SS$EW_gamma)
  grpRR_SS_summary$beta <- grpRR_SS$EW_beta
  grpRR_SS_summary$intercept <- grpRR_SS$intercept
  grpRR_SS_summary$sparsity <- grpRR_SS$EW_pi
  grpRR_SS_summary$out <- grpRR_SS
  summaryList$grpRR_SS <- grpRR_SS_summary

  # set factos with a low posteriori inclusion probability to zero
  grpRR_SScutoff_summary <- list()
  grpRR_SScutoff_summary$runtime <- as.numeric(timeSS)
  grpRR_SScutoff_summary$pf <- as.numeric(grpRR_SS$EW_gamma)
  grpRR_SScutoff_summary$beta <- ifelse(grpRR_SS$EW_s < 0.5, 0, grpRR_SS$EW_beta)
  grpRR_SScutoff_summary$intercept <- grpRR_SS$intercept
  grpRR_SScutoff_summary$sparsity <- grpRR_SS$EW_pi
  grpRR_SScutoff_summary$out <- NULL
  summaryList$grpRR_SScutoff <- grpRR_SScutoff_summary
  rm(grpRR_SS)

  # grpRR_SS_nogamma: fully factorized variational distribution, spike and slab prior (sparse) without different slab precisions
  if(include_nogamma){
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting grpRR_SS model without gamma...")
    grpRR_SS_nogamma <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = TRUE,
                                  spikeslab = TRUE, max_iter = max_iter, intercept = intercept,
                                  verbose = verbose, freqELB = freqELB, calcELB = calcELB, th = th,
                                  family = family,  nogamma = TRUE, standardize=standardize, n_rep=n_rep)
    timeSS_nogamma <- difftime(Sys.time(), tmp, units = "secs")
    grpRR_SS_nogamma_summary <- list()
    grpRR_SS_nogamma_summary$runtime <- as.numeric(timeSS_nogamma)
    grpRR_SS_nogamma_summary$pf <- rep(as.numeric(grpRR_SS_nogamma$EW_gamma), G)
    grpRR_SS_nogamma_summary$beta <- grpRR_SS_nogamma$EW_beta
    grpRR_SS_nogamma_summary$intercept <- grpRR_SS_nogamma$intercept
    grpRR_SS_nogamma_summary$sparsity <- grpRR_SS_nogamma$EW_pi
    grpRR_SS_nogamma_summary$out <- grpRR_SS_nogamma
    summaryList$grpRR_SS_nogamma <- grpRR_SS_nogamma_summary
    rm(grpRR_SS_nogamma)
  }

  # grpRR_SS_ungrouped: fully factorized variational distribution, spike and slab prior
  # (sparse) without group annotations
  if(include_grpRR_SS_ungrouped){
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting grpRR_SS model without group annotations...")
    grpRR_SS_ungrouped <- fit_grpRR(Xtrain, ytrain, annot = rep(1,ncol(Xtrain)), factoriseQ = TRUE,
                                    spikeslab = TRUE, max_iter = max_iter, intercept = intercept,
                                    verbose = verbose, freqELB = freqELB, calcELB = calcELB, th = th,
                                    family = family,  nogamma = TRUE, standardize=standardize, n_rep=n_rep)
    timeSS_ungrouped <- difftime(Sys.time(), tmp, units = "secs")
    grpRR_SS_ungrouped_summary <- list()
    grpRR_SS_ungrouped_summary$runtime <- as.numeric(timeSS_ungrouped)
    grpRR_SS_ungrouped_summary$pf <- rep(as.numeric(grpRR_SS_ungrouped$EW_gamma), G)
    grpRR_SS_ungrouped_summary$beta <- grpRR_SS_ungrouped$EW_beta
    grpRR_SS_ungrouped_summary$intercept <- grpRR_SS_ungrouped$intercept
    grpRR_SS_ungrouped_summary$sparsity <- rep(grpRR_SS_ungrouped$EW_pi, G)
    grpRR_SS_ungrouped_summary$out <- grpRR_SS_ungrouped
    summaryList$grpRR_SS_ungrouped <- grpRR_SS_ungrouped_summary
    rm(grpRR_SS_ungrouped)
  }

  # ridge regression
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting ridge regression...")
  RidgeFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0, intercept = intercept,
                                standardize = standardize, family = family)
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  beta_ridge <- as.vector(stats::coef(RidgeFit, RidgeFit$lambda.min))[-1]
  Ridge_summary <- list()
  Ridge_summary$runtime <- as.numeric(tmp)
  Ridge_summary$pf <- rep(RidgeFit$lambda.min, G)
  Ridge_summary$beta <- beta_ridge
  if(intercept) Ridge_summary$intercept <- as.vector(stats::coef(RidgeFit, RidgeFit$lambda.min))[1]
  Ridge_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
  Ridge_summary$out <- RidgeFit
  rm(RidgeFit, beta_ridge)
  summaryList$Ridge <- Ridge_summary

  # Lasso
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting Lasso...")
  LassoFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 1, intercept = intercept, standardize = standardize, family = family)
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  beta_lasso <- as.vector(stats::coef(LassoFit, LassoFit$lambda.min))[-1]
  Lasso_summary <- list()
  Lasso_summary$runtime <- as.numeric(tmp)
  Lasso_summary$pf <- rep(LassoFit$lambda.min, G)
  Lasso_summary$beta <- beta_lasso
  if(intercept) Lasso_summary$intercept <- as.vector(stats::coef(LassoFit, LassoFit$lambda.min))[1]
  Lasso_summary$sparsity <- vapply(unique(annot),
                                   function(gr) sum(beta_lasso[annot == gr] != 0)/sum(annot == gr),
                                   numeric(1))
  Lasso_summary$out <- LassoFit
  rm(LassoFit, beta_lasso)
  summaryList$Lasso <- Lasso_summary

  # elastic net
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting elastic net...")
  ENFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0.2, intercept = intercept, standardize = standardize, family = family)
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  beta_EN <- as.vector(stats::coef(ENFit, ENFit$lambda.min))[-1]
  ElasticNet_summary <- list()
  ElasticNet_summary$runtime <- as.numeric(tmp)
  ElasticNet_summary$pf <- rep(ENFit$lambda.min, G)
  ElasticNet_summary$beta <- beta_EN
  if(intercept) ElasticNet_summary$intercept <- as.vector(stats::coef(ENFit, ENFit$lambda.min))[1]
  ElasticNet_summary$sparsity <- vapply(unique(annot),
                                        function(gr) sum(beta_EN[annot == gr] != 0)/sum(annot == gr),
                                        numeric(1))
  ElasticNet_summary$out <- ENFit
  rm(ENFit, beta_EN)
  summaryList$ElasticNet <- ElasticNet_summary

  # varbvs
  if(include_varbvs){
    if(!intercept) warning("varbvs always fits an intercept")
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting varbvs...")
    varbvsFit <- varbvs::varbvs(X=Xtrain, Z=NULL, y=ytrain, family = family)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    #if (intercept)
    beta_varbvs <- apply(varbvsFit$alpha * varbvsFit$mu,1, function(b) sum(b * varbvsFit$w))
    varbvs_summary <- list()
    varbvs_summary$runtime <- as.numeric(tmp)
    varbvs_summary$pf <- rep(sum(varbvsFit$sa * varbvsFit$w),G)
    varbvs_summary$beta <- beta_varbvs
    if(intercept) varbvs_summary$identercept <- sum(varbvsFit$mu.cov * varbvsFit$w)
    varbvs_summary$sparsity <- vapply(unique(annot),
                                      function(gr) mean(varbvsFit$pip[annot == gr]),
                                      numeric(1))
    varbvs_summary$out <- varbvsFit
    rm(varbvsFit, beta_varbvs)
    summaryList$varbvs <- varbvs_summary
  }

  # Random Forest
  if (includeRF) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting random forest...")
    if(family=="gaussian") rf.out <- randomForest::randomForest(x = Xtrain, y = ytrain)
    else if(family=="binomial") rf.out <- randomForest::randomForest(x = Xtrain, y = as.factor(ytrain))
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    RandomForest_summary <- list()
    RandomForest_summary$runtime <- as.numeric(tmp)
    RandomForest_summary$pf <- NULL
    RandomForest_summary$beta <- NULL
    RandomForest_summary$intercept <- NULL
    RandomForest_summary$sparsity <- NULL
    RandomForest_summary$out <- rf.out
    rm(rf.out)
    summaryList$RandomForest <- RandomForest_summary
  }

  # group lasso
  if (compareGroupLasso) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting group Lasso...")
    GrpLassoFit <- try(grpreg::cv.grpreg(Xtrain, ytrain, group = as.factor(annot),
                                         penalty = "grLasso", intercept = intercept,
                                         family = family))
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (class(GrpLassoFit) == "try-error") {
      warning("Group Lasso encountered errors, not included in the comparison!")
    } else {
      beta_GrpLasso <- as.vector(stats::coef(GrpLassoFit, GrpLassoFit$lambda.min))[-1]
      GroupLasso_summary <- list()
      GroupLasso_summary$runtime <- as.numeric(tmp)
      GroupLasso_summary$pf <- rep(GrpLassoFit$lambda.min, G)
      GroupLasso_summary$beta <- beta_GrpLasso
      if(intercept) GroupLasso_summary$intercept <- as.vector(stats::coef(GrpLassoFit, GrpLassoFit$lambda.min))[1]
      GroupLasso_summary$sparsity <- vapply(unique(annot),
                                            function(gr) sum(beta_GrpLasso[annot == gr] != 0)/sum(annot == gr),
                                            numeric(1))
      GroupLasso_summary$out <- GrpLassoFit
      rm(GrpLassoFit, beta_GrpLasso)
      summaryList$GroupLasso <- GroupLasso_summary
    }
  }

  # sparse group lasso
  if (compareSparseGroupLasso) {
    if(intercept) warning("Sparse group lasso does not fit an intercept, need to center beforehand")
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting sparse group Lasso...")
    SpGrpLassoFit <- try(SGL::cvSGL(list(x=Xtrain, y=ytrain),
                                    index = as.factor(annot), standardize = standardize,
                                    type = ifelse(family=="gaussian", "linear", "logit")))
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (class(SpGrpLassoFit) == "try-error") {
      warning("Sparse Group Lasso encountered errors, not included in the comparison!")
    } else {
      beta_SpGrpLasso <- SpGrpLassoFit$fit$beta[, which.min(SpGrpLassoFit$lldiff)]
      SpGroupLasso_summary <- list()
      SpGroupLasso_summary$runtime <- as.numeric(tmp)
      SpGroupLasso_summary$pf <- rep(SpGrpLassoFit$lambdas[which.min(SpGrpLassoFit$lldiff)], G)
      SpGroupLasso_summary$beta <- beta_SpGrpLasso
      SpGroupLasso_summary$intercept <- NULL
      SpGroupLasso_summary$sparsity <- vapply(unique(annot),
                                              function(gr) sum(beta_SpGrpLasso[annot == gr] != 0)/sum(annot == gr),
                                              numeric(1))
      SpGroupLasso_summary$out <- SpGrpLassoFit
      rm(SpGrpLassoFit, beta_SpGrpLasso)
      summaryList$SparseGroupLasso <- SpGroupLasso_summary
    }
  }

  # GRridge
  if (compareGRridge) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting GRridge...")
    partition <- GRridge::CreatePartition(as.factor(annot))
    if (intercept) {
      #notes itself which type of respone (new version!)
      MessagesGR <- utils::capture.output(GRfit <- try(GRridge::grridge(t(Xtrain), as.numeric(ytrain), list(partition), unpenal = ~1)))
    } else {
      #notes itself which type of respone (new version!)
      MessagesGR <- utils::capture.output(GRfit <- try(GRridge::grridge(t(Xtrain), as.numeric(ytrain), list(partition), unpenal = ~0)))
    }
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (verbose)
      MessagesGR
    if (class(GRfit) == "try-error") {
      warning("GRridge encountered errors, not included in the comparison!")
    } else {
      GRridge_summary <- list()
      GRridge_summary$runtime <- as.numeric(tmp)
      GRridge_summary$pf <- as.numeric(GRfit$lambdamults[[1]])
      GRridge_summary$beta <- GRfit$betas
      if(intercept) {
        if(family=="gaussian") GRridge_summary$intercept <- GRfit$predobj$GroupRegul@unpenalized
        else GRridge_summary$intercept <- GRfit$predobj$GroupRegul@unpenalized
      }
      GRridge_summary$sparsity <- rep(1, G)
      GRridge_summary$out <- GRfit
      rm(GRfit)
      summaryList$GRridge <- GRridge_summary
    }
  }

  # zero model
  tmp <- Sys.time()
  if(verbose_progress) message(" ### Fitting zero model...")
  if (intercept) intercept_zeromodel <- mean(ytrain) else intercept_zeromodel <- 0
  tmp <- difftime(Sys.time(), tmp, units = "secs")
  NullModel_summary <- list()
  NullModel_summary$runtime <- as.numeric(tmp)
  NullModel_summary$pf <- rep(Inf, G)
  NullModel_summary$beta <- rep(0, p)
  NullModel_summary$intercept <- intercept_zeromodel
  NullModel_summary$sparsity <- rep(0, G)
  NullModel_summary$out <- NULL
  summaryList$NullModel <- NullModel_summary

  # True Model
  if (!is.null(beta0)) {
    if(verbose_progress) message(" ### Including true model...")
    TrueModel_summary <- list()
    TrueModel_summary$runtime <- 0
    TrueModel_summary$pf <- 1/vapply(unique(annot),
                                     function(gr) mean((beta0[annot == gr])^2),
                                     numeric(1))
    TrueModel_summary$beta <- beta0
    TrueModel_summary$intercept <- trueintercept
    TrueModel_summary$sparsity <- vapply(unique(annot),
                                         function(gr) sum(beta0[annot == gr] != 0)/sum(annot == gr),
                                         numeric(1))
    TrueModel_summary$out <- NULL
    summaryList$TrueModel <- TrueModel_summary
  }

  # IPF -Lasso
  if (compareIPF) {
    # penalty factors to consider for cross-validation
    # (unclear how to choose, take a very rough grid here to make it applicable to larger number of groups)
    lambda_1d <- c(0.1, 0.5, 1, 2, 10)
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting IPF-Lasso...")
    pfgrid <- expand.grid(rep(list(lambda_1d), G))
    pflist <- lapply(seq_len(nrow(pfgrid)), function(i) pfgrid[i, ])
    type.measure <- ifelse(family == "gaussian", "mse", "class")
    ipf.out <- try(ipflasso::cvr2.ipflasso(Xtrain, ytrain, alpha = 1, standardize = standardize,
                                           family = family, type.measure = type.measure,
                                           blocks = lapply(unique(annot), function(gr) which(annot == gr)),
                                           pflist = pflist, nfolds = 10, ncv = 3))
    #using same cv parameter as standard glmnet leads to errors, needs to be ncv>1
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (class(ipf.out) == "try-error") {
      warning("ipf-lasso encountered errors, not included in the comparison!")
    } else {
      IPFLasso_summary <- list()
      IPFLasso_summary$runtime <- as.numeric(tmp)
      IPFLasso_summary$pf <- pflist[[ipf.out$ind.bestpf]]
      IPFLasso_summary$beta <- ipf.out$coeff[-1,ipf.out$ind.bestlambda]
      IPFLasso_summary$intercept <- ipf.out$coeff[1,ipf.out$ind.bestlambda]
      IPFLasso_summary$sparsity <-  vapply(unique(annot),
                                           function(gr) sum(IPFLasso_summary$beta[annot == gr] != 0)/sum(annot == gr),
                                           numeric(1))
      IPFLasso_summary$out <- ipf.out
      rm(ipf.out)
      summaryList$IPFLasso <- IPFLasso_summary
    }

  }

  # Adaptive Lasso
  if (compareAdaLasso) {
    tmp <- Sys.time()
    if(verbose_progress) message(" ### Fitting adaptive Lasso...")
    ## Ridge Regression to create the Adaptive Weights Vector
    RidgeFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0, intercept = intercept, standardize = standardize, family = family)
    wRidge <- pmin(1/abs((stats::coef(RidgeFit, s = RidgeFit$lambda.min))), 1e+300)
    wRidge <- wRidge[-1]
    adaLassoFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 1, intercept = intercept,
                                     standardize = standardize, family = family, penalty.factor = wRidge)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    beta_adalasso <- as.vector(stats::coef(adaLassoFit, adaLassoFit$lambda.min))[-1]
    adaLasso_summary <- list()
    adaLasso_summary$runtime <- as.numeric(tmp)
    adaLasso_summary$pf <- vapply(unique(annot),
                                  function(gr) mean(adaLassoFit$lambda.min * wRidge[annot == gr]),
                                  numeric(1))
    adaLasso_summary$beta <- beta_adalasso
    if(intercept) adaLasso_summary$intercept <- as.vector(stats::coef(adaLassoFit, adaLassoFit$lambda.min))[1]
    adaLasso_summary$sparsity <- vapply(unique(annot),
                                        function(gr) sum(beta_adalasso[annot == gr] != 0)/sum(annot == gr),
                                        numeric(1))
    adaLasso_summary$out <- adaLassoFit
    rm(adaLassoFit, beta_adalasso)
    summaryList$adaptiveLasso <- adaLasso_summary
  }

  return(list(summaryList = summaryList, groupnames = groupnames, varnames = varnames,
              family = family, n = n, p = p, G = G, annot = annot))
}



# ---------------------------

#' @title Evaluate fits from various regression methods
#' @name evaluateFits
#' @description Function to evaluate results on test data
#' @param allFits List as produced by \code{\link{runMethods}}
#' @param Xtest Design matrix of size n' x p (same feature structure as used in \code{\link{runMethods}})
#' @param ytest Response vector of size n'
#' @return List as prodcused by \code{\link{runMethods}} with additional predicition performance slots
#' @export
#' @importFrom stats predict


evaluateFits <- function(allFits, Xtest, ytest) {

  stopifnot(nrow(Xtest) == length(ytest))
  stopifnot(ncol(Xtest) == length(allFits$summaryList[[1]]$beta))

  family <- allFits$family
  summaryList <- allFits$summaryList
  ytest <- as.vector(ytest)
  ntest <- length(ytest)
  beta0 <- summaryList$TrueModel$beta
  intercept0 <- summaryList$TrueModel$intercept

  ###### Prediction Performance

  # For gaussian family calculate RMSE as measure of prediciton performance
  if (family == "gaussian") {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      intercept <- summary$intercept
      if (!is.null(beta)) {
        # for cases with linear coefficients
        RMSE <- EvaluateModel(beta, intercept = intercept, Xtest, ytest, beta0 = beta0, family = "gaussian")$RMSE_test
        summary$RMSE <- RMSE
      }
      summary
    })
    if ("RandomForest" %in% names(summaryList))
      summaryList$RandomForest$RMSE <- sqrt(1/length(ytest) * sum((stats::predict(summaryList$RandomForest$out, Xtest) - ytest)^2))
    if ("varbvs" %in% names(summaryList))
      summaryList$varbvs$RMSE <- sqrt(1/length(ytest) * sum((stats::predict(summaryList$varbvs$out, Xtest) - ytest)^2))
  }

  # For binomial family calculate AUC and Brier Score
  if (family == "binomial") {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      intercept <- summary$intercept
      if (!is.null(beta)) {
        # for cases with linear coeeficients
        eval.out <- EvaluateModel(beta, intercept = intercept, Xtest, ytest, beta0 = beta0, family = "binomial")
        summary$AUC <- eval.out$AUC
        summary$BS <- eval.out$BrierScore
      }
      # TODO add varbvs and RandomForest
      summary
    })
  }

  ###### Feature Recovery
  if (!is.null(beta0)) {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      intercept <- summary$intercept
      if (!is.null(beta)) {
        # for cases without linear coeeficients e.g. Random Forest
        eval.out <- EvaluateModel(beta, intercept = intercept, Xtest, ytest, beta0 = beta0, family = family)
        summary$FPR <- eval.out$FPR
        summary$FNR <- eval.out$FNR
        summary$precision <- eval.out$precision
        summary$recall <- eval.out$recall
        summary$F1score <- eval.out$F1score
      }
      summary
    })
  }


  ###### Error on estimate

  # l1-Error in estimation of coeffcients
  if (!is.null(beta0)) {
    summaryList <- lapply(summaryList, function(summary) {
      beta <- summary$beta
      if (!is.null(beta))
        summary$l1error_beta <- sum(abs(beta - beta0))
      summary
    })
  }

  # l1-Error in estimation of intercept
  if (!is.null(intercept0)) {
    summaryList <- lapply(summaryList, function(summary) {
      intercept <- summary$intercept
      if (!is.null(intercept))
        summary$l1error_intercept <- sum(abs(intercept - intercept0))
      summary
    })
  }

  allFits$summaryList <- summaryList

  return(allFits)
}

#' @title Compare various regression method via cross-validation
#' @name cv_compare
#' @description  Function to run serveral different methods for high-dimensional regression and evaluate them in a cross-validated fashion
#' @param X Design matrix of size n x p
#' @param y Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param family Liklihood model to use for response, either binomial or gaussian
#' @param nfolds Number of fold for evaluation
#' @param ncores Number of cores to use
#' @param plot_cv boolean whether to plot summary from evaluation
#' @param seed optional seed for the choice of folds
#' @param parallel boolean: Run cross-validation in parallel?
#' @param saveFits boolean: Save the fit of each fold?
#' @param ... Other parameters that can be passed to \code{\link{runMethods}}
#'
#' @return List of fitted models and two data frames with coeffcients and penalty factors
#' @import ggplot2
#' @import parallel
#' @export

cv_compare <- function(X, y, annot, family="gaussian",
                       ncores=1, nfolds=10, plot_cv=TRUE,
                       seed=NULL, parallel=FALSE, saveFits=FALSE, ...){

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
    AllFits <- runMethods(Xtrain = as.matrix(Xtrain),
                          ytrain =as.vector(ytrain),
                          annot = annot, family = family,
                          ...)

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
                l1error_intercept=l1error_intercept, l1error_beta=l1error_beta)
    } else if(family=="binomial"){
      BS <- getBS(AllFits)
      AUC <- getAUC(AllFits)
      l <- list(FPR=FPR, FNR=FNR, BS=BS, AUC=AUC, pf_mat=pf_mat, beta_mat=beta_mat,
                intercepts=intercepts, sparsity_mat=sparsity_mat, annot=AllFits$annot, runtime=runtime,
                l1error_intercept=l1error_intercept, l1error_beta=l1error_beta)
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



