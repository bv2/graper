#  evaluation_functions
# This file contains functions used for comparison to other methods and evaluation of the produced fit
#   - RunMethods: Run a number of different methods for estimation in the model including
#        - Bayes Model with groupwise normal prior
#         - Bayes Model with groupwise normal prior and fully-factorized inference
#         - Bayes Modeln with groupwise spike and slab prior and fully-factorized inference
#         - Bayes Model with groupwise spike and slab prior and fully-factorized inference and feature selection
#         - Ridge
#         - Lasso
#         - Group Lasso
#         - Random forest
#         - Adaptive Lasso
#         - Elastic Net
#         - IPF-Lasso
#         - GRridge
#   - evalResult: Compare Estimates produced in RunMethods by MSE and estimation error
#   - plotMethodComparison: Plot comparison results from evalResult
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
#' @param compareIPF boolean, indicating wether to fit a IPFLasso (might cause errors, look for stable version!)
#' @param compareAdaLasso boolean, indicating wether to fit an adpative lasso
#' @param freqELB determines frequency at which ELB is to be calculated, i.e. each feqELB-th iteration
#' @return List of fitted models and two data frames with coeffcients and penalty factors
#' @import ggplot2
#' @export

# ToDo: Option to standardize the input and then re-adjust estimates coefficeints to not standardized
RunMethods <- function(Xtrain, ytrain, annot, beta0 = NULL, trueintercept = NULL, max_iter = 2000, intercept = T, plotit = F, standardize = T,
    verbose = F, compareGRridge = F, freqELB = 10, calcELB = T, include_nonfacQ = T, family = "gaussian", constantXcol = F, compareGroupLasso = T,
    includeRF = T, th = 1e-07, compareIPF = T, compareAdaLasso = T) {

    if (!standardize)
        warning("Group Lasso and Grridge are standardized despite of standardized = F")

    stopifnot(nrow(Xtrain) == length(ytrain))

    if (!is.null(beta0))
        stopifnot(ncol(Xtrain) == length(beta0))


    if (constantXcol)
        stopifnot(var(Xtrain[, 1]) == 0)
    if (constantXcol) {
        if (intercept) {
            warning("Intercept not used, as constant X columns included")
            intercept <- F
        }
        penaltyFac <- c(0, rep(1, ncol(Xtrain) - 1))
    } else penaltyFac <- rep(1, ncol(Xtrain))


    # turn annot to facotr
    annot <- as.factor(annot)
    # need to strucute annot - TO DO stopifnot(all(order(annot) == 1:p))

    # extract important parameters and names
    p <- ncol(Xtrain)
    n <- nrow(Xtrain)
    G <- length(unique(annot))
    groupnames <- as.character(annot)
    varnames <- colnames(Xtrain)
    if (length(varnames) == 0)
        varnames <- factor(paste("Feautre", 1:ncol(Xtrain), sep = ""))

    #### RUN DIFFERENT METHODS
    summaryList <- list()

    # grpRR: not fully factorized, normal prior
    if (include_nonfacQ) {
        tmp <- Sys.time()
        grpRR <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = F, spikeslab = F, max_iter = max_iter, intercept = intercept,
            verbose = verbose, freqELB = freqELB, calcELB = calcELB, family = family, th = th)
        timeNF <- difftime(Sys.time(), tmp, units = "secs")
        if (plotit)
            plotVBFit(grpRR, whichParam = c("ELB", "tau", "gamma"))

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

    # grpRR_FF : fully factorized, normal prior
    tmp <- Sys.time()
    grpRR_FF <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = T, spikeslab = F, max_iter = max_iter, intercept = intercept,
        verbose = verbose, freqELB = freqELB, calcELB = calcELB, family = family, th = th)
    timeFF <- difftime(Sys.time(), tmp, units = "secs")
    if (plotit)
        plotVBFit(grpRR_FF, whichParam = c("ELB", "tau", "gamma"))

    grpRR_FF_summary <- list()
    grpRR_FF_summary$runtime <- as.numeric(timeFF)
    grpRR_FF_summary$pf <- as.numeric(grpRR_FF$EW_gamma)
    grpRR_FF_summary$beta <- grpRR_FF$EW_beta
    grpRR_FF_summary$intercept <- grpRR_FF$intercept
    grpRR_FF_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
    grpRR_FF_summary$out <- grpRR_FF
    rm(grpRR_FF)
    summaryList$grpRR_FF <- grpRR_FF_summary



    # grpRR_SS: fully factorized, spike and slab This part is only implemented for gaussian so far
        includeSS <- T

        tmp <- Sys.time()
        grpRR_SS <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = T, spikeslab = T, max_iter = max_iter, intercept = intercept,
            verbose = verbose, freqELB = freqELB, calcELB = calcELB, th = th, family = family)
        timeSS <- difftime(Sys.time(), tmp, units = "secs")
        if (plotit)
            plotVBFit(grpRR_SS, whichParam = c("ELB", "tau", "gamma"))

        grpRR_SS_summary <- list()
        grpRR_SS_summary$runtime <- as.numeric(timeSS)
        grpRR_SS_summary$pf <- as.numeric(grpRR_SS$EW_gamma)
        grpRR_SS_summary$beta <- grpRR_SS$EW_beta
        grpRR_SS_summary$intercept <- grpRR_SS$intercept
        grpRR_SS_summary$sparsity <- grpRR_SS$EW_pi
        grpRR_SS_summary$out <- grpRR_SS
        summaryList$grpRR_SS <- grpRR_SS_summary

        # set factos with low a posteriori probability to zero
        grpRR_SScutoff_summary <- list()
        grpRR_SScutoff_summary$runtime <- as.numeric(timeSS)
        grpRR_SScutoff_summary$pf <- as.numeric(grpRR_SS$EW_gamma)
        grpRR_SScutoff_summary$beta <- ifelse(grpRR_SS$EW_s < 0.5, 0, grpRR_SS$EW_beta)
        grpRR_SScutoff_summary$intercept <- grpRR_SS$intercept
        grpRR_SScutoff_summary$sparsity <- grpRR_SS$EW_pi
        grpRR_SScutoff_summary$out <- NULL
        summaryList$grpRR_SScutoff <- grpRR_SScutoff_summary

        rm(grpRR_SS)

    # grpRR_SS: fully factorized, spike and slab This part is only implemented for gaussian so far
        tmp <- Sys.time()
        grpRR_SS_nogamma <- fit_grpRR(Xtrain, ytrain, annot = annot, factoriseQ = T, spikeslab = T, max_iter = max_iter, intercept = intercept,
            verbose = verbose, freqELB = freqELB, calcELB = calcELB, th = th, family = family,  nogamma = TRUE)
        timeSS_nogamma <- difftime(Sys.time(), tmp, units = "secs")
        if (plotit)
            plotVBFit(grpRR_SS, whichParam = c("ELB", "tau", "gamma"))

        grpRR_SS_nogamma_summary <- list()
        grpRR_SS_nogamma_summary$runtime <- as.numeric(timeSS_nogamma)
        grpRR_SS_nogamma_summary$pf <- rep(as.numeric(grpRR_SS_nogamma$EW_gamma), G)
        grpRR_SS_nogamma_summary$beta <- grpRR_SS_nogamma$EW_beta
        grpRR_SS_nogamma_summary$intercept <- grpRR_SS_nogamma$intercept
        grpRR_SS_nogamma_summary$sparsity <- grpRR_SS_nogamma$EW_pi
        grpRR_SS_nogamma_summary$out <- grpRR_SS_nogamma
        summaryList$grpRR_SS_nogamma <- grpRR_SS_nogamma_summary

    # ridge regression
    tmp <- Sys.time()
    RidgeFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0, intercept = intercept, standardize = standardize, family = family, penalty.factor = penaltyFac)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (intercept)
        beta_ridge <- as.vector(coef(RidgeFit, RidgeFit$lambda.min))[-1] else beta_ridge <- as.vector(coef(RidgeFit, RidgeFit$lambda.min))

    Ridge_summary <- list()
    Ridge_summary$runtime <- as.numeric(tmp)
    Ridge_summary$pf <- rep(RidgeFit$lambda.min, G)
    Ridge_summary$beta <- beta_ridge
    Ridge_summary$intercept <- ifelse(intercept, as.vector(coef(RidgeFit, RidgeFit$lambda.min))[1], NULL)
    Ridge_summary$sparsity <- rep(1, G)  #dense - no sparsity per groups
    Ridge_summary$out <- RidgeFit
    rm(RidgeFit, beta_ridge)
    summaryList$Ridge <- Ridge_summary



    # Lasso
    tmp <- Sys.time()
    LassoFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 1, intercept = intercept, standardize = standardize, family = family, penalty.factor = penaltyFac)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (intercept)
        beta_lasso <- as.vector(coef(LassoFit, LassoFit$lambda.min))[-1] else beta_lasso <- as.vector(coef(LassoFit, LassoFit$lambda.min))

    Lasso_summary <- list()
    Lasso_summary$runtime <- as.numeric(tmp)
    Lasso_summary$pf <- rep(LassoFit$lambda.min, G)
    Lasso_summary$beta <- beta_lasso
    Lasso_summary$intercept <- ifelse(intercept, as.vector(coef(LassoFit, LassoFit$lambda.min))[1], NULL)
    Lasso_summary$sparsity <- sapply(unique(annot), function(gr) sum(beta_lasso[annot == gr] != 0)/sum(annot == gr))
    Lasso_summary$out <- LassoFit
    rm(LassoFit, beta_lasso)
    summaryList$Lasso <- Lasso_summary


    # EN
    tmp <- Sys.time()
    ENFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0.2, intercept = intercept, standardize = standardize, family = family, penalty.factor = penaltyFac)
    tmp <- difftime(Sys.time(), tmp, units = "secs")
    if (intercept)
        beta_EN <- as.vector(coef(ENFit, ENFit$lambda.min))[-1] else beta_EN <- as.vector(coef(ENFit, ENFit$lambda.min))

    ElasticNet_summary <- list()
    ElasticNet_summary$runtime <- as.numeric(tmp)
    ElasticNet_summary$pf <- rep(ENFit$lambda.min, G)
    ElasticNet_summary$beta <- beta_EN
    ElasticNet_summary$intercept <- ifelse(intercept, as.vector(coef(ENFit, ENFit$lambda.min))[1], NULL)
    ElasticNet_summary$sparsity <- sapply(unique(annot), function(gr) sum(beta_EN[annot == gr] != 0)/sum(annot == gr))
    ElasticNet_summary$out <- ENFit
    rm(ENFit, beta_EN)
    summaryList$ElasticNet <- ElasticNet_summary

    # Random Forest
    if (includeRF) {
        tmp <- Sys.time()
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
        GrpLassoFit <- try(grpreg::cv.grpreg(Xtrain, ytrain, group = as.factor(annot), penalty = "grLasso", intercept = intercept,
            family = family))
        tmp <- difftime(Sys.time(), tmp, units = "secs")
        if (class(GrpLassoFit) == "try-error") {
            warning("Group Lasso encountered errors, not included in the comparison!")
        } else {
            if (intercept)
                beta_GrpLasso <- as.vector(coef(GrpLassoFit, GrpLassoFit$lambda.min))[-1] else beta_GrpLasso <- as.vector(coef(GrpLassoFit, GrpLassoFit$lambda.min))

            GroupLasso_summary <- list()
            GroupLasso_summary$runtime <- as.numeric(tmp)
            GroupLasso_summary$pf <- rep(GrpLassoFit$lambda.min, G)
            GroupLasso_summary$beta <- beta_GrpLasso
            GroupLasso_summary$intercept <- ifelse(intercept, as.vector(coef(GrpLassoFit, GrpLassoFit$lambda.min))[1], NULL)
            GroupLasso_summary$sparsity <- sapply(unique(annot), function(gr) sum(beta_GrpLasso[annot == gr] != 0)/sum(annot == gr))
            GroupLasso_summary$out <- GrpLassoFit
            rm(GrpLassoFit, beta_GrpLasso)
            summaryList$GroupLasso <- GroupLasso_summary

        }
    }

    # grridge
    if (compareGRridge) {
        tmp <- Sys.time()
        partition <- GRridge::CreatePartition(as.factor(annot))
        if (intercept) {
          if(family=="gaussian") MessagesGR <- capture.output(GRfit <- try(GRridge::grridgelin(t(Xtrain), as.numeric(ytrain), list(partition), unpenal = ~1)))
            else MessagesGR <- capture.output(GRfit <- try(GRridge::grridge(t(Xtrain), as.numeric(ytrain), list(partition), unpenal = ~1)))
        } else {
          if(family=="gaussian") MessagesGR <- capture.output(GRfit <- try(GRridge::grridgelin(t(Xtrain), as.numeric(ytrain), list(partition), unpenal = ~0)))
            else  MessagesGR <- capture.output(GRfit <- try(GRridge::grridge(t(Xtrain), as.numeric(ytrain), list(partition), unpenal = ~0)))
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
            if(family=="gaussian") GRridge_summary$intercept <- ifelse(intercept, coef(GRfit$predobj$GroupRegul)[1], NULL)
              else GRridge_summary$intercept <- ifelse(intercept, GRfit$predobj$GroupRegul@unpenalized, NULL) 
            GRridge_summary$sparsity <- rep(1, G)
            GRridge_summary$out <- GRfit
            rm(GRfit)
            summaryList$GRridge <- GRridge_summary

        }
    }

    # zero model
    if (intercept) {
        tmp <- Sys.time()
        intercept_zeromodel <- mean(ytrain)
        tmp <- difftime(Sys.time(), tmp, units = "secs")

        NullModel_summary <- list()
        NullModel_summary$runtime <- as.numeric(tmp)
        NullModel_summary$pf <- rep(Inf, G)
        NullModel_summary$beta <- rep(0, p)
        NullModel_summary$intercept <- intercept_zeromodel
        NullModel_summary$sparsity <- rep(0, G)
        NullModel_summary$out <- NULL
        summaryList$NullModel <- NullModel_summary

    }

    if (!is.null(beta0)) {
        TrueModel_summary <- list()
        TrueModel_summary$runtime <- 0
        # TO DO What the best correspondence?
        TrueModel_summary$pf <- 1/sapply(unique(annot), function(gr) mean(abs(beta0[annot == gr])))
        TrueModel_summary$beta <- beta0
        TrueModel_summary$intercept <- trueintercept
        TrueModel_summary$sparsity <- sapply(unique(annot), function(gr) sum(beta0[annot == gr] != 0)/sum(annot == gr))
        TrueModel_summary$out <- NULL
        summaryList$TrueModel <- TrueModel_summary
    }

    # IPF -Lasso NOTE alwyas fit an intercept error because cv.glmnet and its element cv.glmnet$glmnet.fit haev different lambda
    # sequences? Can happen if some cvsd are NA nas = is.na(cvsd)???
    if (compareIPF) {
        # penalty factors to consider for cross-validation (unclear how to choose)
        lambda_1d <- seq(1, 10, 2)
        tmp <- Sys.time()
        pfgrid <- expand.grid(rep(list(lambda_1d), G))
        pflist <- lapply(seq_len(nrow(pfgrid)), function(i) pfgrid[i, ])
        type.measure <- ifelse(family == "gaussian", "mse", "class")
        ipf.out <- try(ipflasso::cvr2.ipflasso(Xtrain, ytrain, alpha = 1, standardize = standardize, family = family, type.measure = type.measure,
            blocks = lapply(unique(annot), function(gr) which(annot == gr)), pflist = pflist, nfolds = 10, ncv = 1))  #using same cv parameter as standard glmnet
        tmp <- difftime(Sys.time(), tmp, units = "secs")

        if (class(ipf.out) == "try-error") {
            warning("ipf-lasso encountered errors, not included in the comparison!")
        } else {
            IPFLasso_summary <- list()
            IPFLasso_summary$runtime <- as.numeric(tmp)
            IPFLasso_summary$pf <- ind.bestpf
            IPFLasso_summary$beta <- NULL
            IPFLasso_summary$intercept <- NULL
            IPFLasso_summary$sparsity <- NULL
            IPFLasso_summary$out <- rf.out
            rm(ipf.out)
            summaryList$IPFLasso <- IPFLasso_summary
        }

    }

    # Adaptive Lasso
    if (compareAdaLasso) {
        tmp <- Sys.time()
        ## Ridge Regression to create the Adaptive Weights Vector
        RidgeFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 0, intercept = intercept, standardize = standardize, family = family, penalty.factor = penaltyFac)
        wRidge <- pmin(1/abs((coef(RidgeFit, s = RidgeFit$lambda.min))), 1e+300)
        if (intercept)
            wRidge <- wRidge[-1]

        ## Adaptive Lasso
        adaLassoFit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = 1, intercept = intercept, standardize = standardize, family = family, penalty.factor = penaltyFac *
            wRidge)
        tmp <- difftime(Sys.time(), tmp, units = "secs")

        if (intercept)
            beta_adalasso <- as.vector(coef(adaLassoFit, adaLassoFit$lambda.min))[-1] else beta_adalasso <- as.vector(coef(adaLassoFit, adaLassoFit$lambda.min))

        adaLasso_summary <- list()
        adaLasso_summary$runtime <- as.numeric(tmp)
        adaLasso_summary$pf <- sapply(unique(annot), function(gr) mean(adaLassoFit$lambda.min * penaltyFac * wRidge[annot == gr]))
        adaLasso_summary$beta <- beta_adalasso
        adaLasso_summary$intercept <- ifelse(intercept, as.vector(coef(adaLassoFit, adaLassoFit$lambda.min))[1], NULL)
        adaLasso_summary$sparsity <- sapply(unique(annot), function(gr) sum(beta_adalasso[annot == gr] != 0)/sum(annot == gr))
        adaLasso_summary$out <- adaLassoFit
        rm(adaLassoFit, beta_adalasso)
        summaryList$adaptiveLasso <- adaLasso_summary
    }
    # #Ridge with PF by average marginal coefficients...better use significance instead of effect size....
    # marg<-MarginalCoefficient(ytrain, scale(Xtrain), family = family) avMargGroup<-sapply(annot, function(g)
    # mean(abs(marg[1,annot==g]))) pf_margCoeff<-1/abs(avMargGroup) #estimates RidgeavMargFit<-cv.glmnet(Xtrain,ytrain, alpha=0,
    # intercept=intercept, penalty.factor=pf_margCoeff,standardize=standardize, family=family)
    # beta_RidgeavMarg<-as.vector(coef(RidgeavMargFit$glmnet.fit, RidgeavMargFit$lambda.min)) #Lasso with PF by average marginal
    # coefficients marg<-MarginalCoefficient(ytrain, Xtrain) avMargGroup<-sapply(annot, function(g) mean(abs(marg[1,annot==g])))
    # pf_margCoeff<-1/abs(avMargGroup) #estimates LassoavMargFit<-cv.glmnet(Xtrain,ytrain, alpha=1, intercept=intercept,
    # penalty.factor=pf_margCoeff,standardize=standardize, family=family)
    # beta_LassoavMarg<-as.vector(coef(LassoavMargFit$glmnet.fit, LassoavMargFit$lambda.min))

    return(list(summaryList = summaryList, groupnames = groupnames, varnames = varnames, family = family, n = n, p = p, G = G, annot = annot))
}

# ---------------------------
#'  evaluateFits
#'
#' Function to evaluate results on test data
#' @param allFits List as produced by \code{\link{runMethods}}
#' @param Xtest Design matrix of size n' x p (same feature structure as used in runMethods)
#' @param ytest Response vector of size n'
#' @return List as prodcused by \code{\link{runMethods}} with additional predicition performance slots
#' @export


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
                # for cases without linear coeeficients e.g. Random Forest
                RMSE <- EvaluateModel(beta, intercept = intercept, Xtest, ytest, beta0 = beta0, family = "gaussian")$RMSE_test
                summary$RMSE <- RMSE
            }
            summary
        })
        if ("RandomForest" %in% names(summaryList))
            summaryList$RandomForest$RMSE <- sqrt(1/length(ytest) * sum((predict(summaryList$RandomForest$out, Xtest) - ytest)^2))
    }

    # For binomial family calculate ROC, AUC and Brier Score
    if (family == "binomial") {
        summaryList <- lapply(summaryList, function(summary) {
            beta <- summary$beta
            intercept <- summary$intercept
            if (!is.null(beta)) {
                # for cases without linear coeeficients e.g. Random Forest
                eval.out <- EvaluateModel(beta, intercept = intercept, Xtest, ytest, beta0 = beta0, family = "binomial")
                summary$AUC <- eval.out$AUC
                summary$ROC <- eval.out$ROC
                summary$BS <- eval.out$BrierScore
            }
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


#' #'  plotMethodComparison
#' Function to plot method comparison across several runs
#' @param resultList List as in simulation_setting1.Rmd
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @export

plotMethodComparison <- function(resultList, plotbeta = F, family = "gaussian") {
    # get results in dataframe format
    if(family=="gaussian"){
    eval_summary <- melt(lapply(resultList, function(l) rbind(FPR = l$FPR, FNR = l$FNR, RMSE = l$RMSE, l1error_beta = l$l1error_beta)),
        varnames = c("measure", "method"), level = "run")
    } else {
      eval_summary <- melt(lapply(resultList, function(l) rbind(FPR = l$FPR, FNR = l$FNR, BS = l$BS, AUC = l$AUC, l1error_beta = l$l1error_beta)),
                           varnames = c("measure", "method"), level = "run")
    }

    pf_summary <- lapply(seq_along(resultList), function(i) cbind(melt(resultList[[i]]$pf_mat, varnames = c("group", "method"),
        value.name = "penalty_factor"), Lrun = i)) %>% bind_rows()

    sparsity_summary <- lapply(seq_along(resultList), function(i) cbind(melt(resultList[[i]]$sparsity_mat, varnames = c("group",
        "method"), value.name = "sparsity_level"), Lrun = i)) %>% bind_rows()

    beta_summary <- lapply(seq_along(resultList), function(i) cbind(melt(resultList[[i]]$beta_mat, varnames = c("feature", "method"),
        value.name = "beta"), Lrun = i)) %>% bind_rows()
    # the folowing only works if annot is names properly, needs to be fixes
    beta_summary$group <- sapply(1:nrow(beta_summary), function(i) resultList[[beta_summary$Lrun[i]]]$annot[beta_summary$feature[i]])

    intercepts_summary <- melt(lapply(resultList, function(l) t(l$intercepts)), varnames = c("const", "method"), value.name = "intercept",
        level = "run")[, 2:4]

    runtime_summary <- melt(lapply(resultList, function(l) t(l$runtime)), varnames = c("const", "method"), value.name = "runtime",
        level = "run")[, 2:4]

    gg_pf <- ggplot(pf_summary, aes(x = as.factor(group), y = penalty_factor, fill = as.factor(group), group = as.factor(group))) + geom_boxplot() + facet_wrap(~method,
        scales = "free_y") + ggtitle("Penalty Factors per group") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    print(gg_pf)

    gg_sparse <- ggplot(sparsity_summary, aes(x = as.factor(group), y = sparsity_level, fill = as.factor(group), group = as.factor(group))) + geom_boxplot() + facet_wrap(~method,
        scales = "free_y") + ggtitle("Sparsity Level per group (1=dense)") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    print(gg_sparse)

    gg_perf <- ggplot(filter(eval_summary, method != "TrueModel"), aes(x = method, y = value, fill = method)) + geom_boxplot() +
        ggtitle("Method comparison") + facet_wrap(~measure, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        ggtitle("Performance measures")
    print(gg_perf)

    if (plotbeta) {
        gg_beta <- ggplot(beta_summary, aes(x = beta, fill = group, group = group)) + geom_histogram(alpha = 0.6, position = "identity") +
            facet_wrap(~method, scales = "free") + ggtitle("Distribution of estimated coefficients per group")
        print(gg_beta)
    }

    gg_runtime <- ggplot(runtime_summary, aes(x = method, y = runtime, group = method, fill = method)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 60,
        hjust = 1)) + ggtitle("Runtime") + ylab("secs")
    print(gg_runtime)

    if(family=="binomial"){
        fprMat <- lapply(resultList, function(res) sapply(res$ROC, function(r) {
            if(!is.na(r)) r["FPR",] else rep(NA, 101)}))
        tprMat <- lapply(resultList, function(res) sapply(res$ROC, function(r) {
            if(!is.na(r)) r["TPR",] else rep(NA, 101)}))

        fprDF <- melt(fprMat, varnames=c("cut", "method"), value.name = "FPR")
        tprDF <- melt(tprMat, varnames=c("cut", "method"), value.name = "TPR")
        rocDF <- merge.data.frame(fprDF, tprDF, by=c("method", "cut", "L1"))

        ggROC <- ggplot(rocDF, aes(x=FPR, y=TPR, col=method)) + geom_line() +facet_wrap(~L1)
        print(ggROC)

    }

}
