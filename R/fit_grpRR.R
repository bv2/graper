# R function to call the different C functions defined in fit_grpRR.cpp
#
# fit_grpRR
#
# This file contains the main function to fit a grpRR model.
# fit_gprRR
# ---------------------------

#'  Fit grpRR model
#'
#'  Main function to fit a grpRR model.
#' @param X Design matrix of size n x p
#' @param y Response vector of size n
#' @param annot Factor of length p indicating group membership of each feature
#' @param factoriseQ If true, the variational distribution is assumed to fully factorize across features (rougher approx., but faster). If spikeslab=F, this is always done.
#' @param spikeslab If true, a spike and slab prior is used instead of a normal prior
#' @param d_tau hyper-parameters for prior of tau (noise precision)
#' @param r_tau hyper-parameters for prior of tau (noise precision)
#' @param d_gamma hyper-parameters for prior of gamma (coeffient's prior precision)
#' @param r_gamma hyper-parameters for prior of gamma (coeffient's prior precision)
#' @param r_pi hyper-parameters for prior of pi (Spike-Slab-Bernoulli variable prior probabiliy of being 1)
#' @param d_pi hyper-parameters for prior of pi (Spike-Slab-Bernoulli variable prior probabiliy of being 1)
#' @param max_iter maximum number of iterations
#' @param th convergence threshold for ELB
#' @param intercept boolean, indicating wether to fit an intercept
#' @param calcELB boolean, indicating wether to calculate ELB
#' @param verbose boolean, indicating wether to print out intermediate messages during fitting
#' @param freqELB determines frequency at which ELB is to be calculated, i.e. each feqELB-th iteration
#' @return List of fitted parameters .....


#' @useDynLib grpRR
#' @import Rcpp
#' @export


fit_grpRR<-function(X,y,annot, factoriseQ=T, spikeslab= F, d_tau=0.001, r_tau=0.001, d_gamma=0.001, r_gamma=0.001,
                    r_pi=1, d_pi=1, max_iter=1000, th=1e-7, intercept=T, calcELB=T, verbose=F, freqELB=10, family="gaussian"){

  stopifnot(ncol(X)==length(annot))

  #check structure of annot: needs to be 1:g with 1 <-> frist group etc
  annot<-as.factor(annot)

  #get data dimension
  p<-ncol(X) #no of features
  n<-nrow(X) #no of samples

  #get group structure
  g<-length(unique(annot))
  NoPerGroup<-sapply(unique(annot), function(x) sum(annot==x))
  names(NoPerGroup)<-unique(annot)

  if(family=="gaussian"){
  if(intercept){
    X<-scale(X, center=T, scale=F)
    y<-scale(y, center=T, scale =F)
  }

    #call C function depending on FacType and spikeslap arguments
  if(spikeslab){
    if(!factoriseQ) warning("Using fully factorized approach with a spike and slab prior")
    res<-grRRCpp_sparse_ff(X,y, annot,g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB, verbose, freqELB)
  } else{
    if(factoriseQ) res<-grRRCpp_dense_ff(X,y, annot,g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB)
    else res<-grRRCpp_dense_nf(X,y, annot,g, NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB)
  }

  #calculate intercept
  if(intercept) intercept<-attr(y, "scaled:center")-attr(X, "scaled:center")%*%res$EW_beta else intercept<-NULL

  #return mean of approximate posterior (other quantities of interes: tau, lower bound on model evidence etc)
  return(append(res, list(intercept=intercept)))
  }

  else if (family=="binomial"){
    if(spikeslab) stop("spikeslab not yet implemented")
    if(intercept) {warning("inercept not yet implemented"); intercept<-F}


    if(factoriseQ) res<-grpRRCpp_logistic_ff(X,y, annot,g, NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB)
    else res<-grpRRCpp_logistic_nf(X,y, annot,g, NoPerGroup,  d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB)

  }
  else stop("Family not implemented")
}
