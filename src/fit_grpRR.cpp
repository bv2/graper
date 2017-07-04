//C Functions calling the main function from classes

#include <RcppArmadillo.h>
#include "grpRRclass.hpp"
#include "grpRRclass_fullyfac.hpp"
#include "grpRR_sparseclass_fullyfac.hpp"
#include "grpRRclass_logistic.hpp"
#include "grpRRclass_logistic_fullyfac.hpp"


// [[Rcpp::depends(RcppArmadillo)]]

//Fitting a normal prior model with only partially factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grRRCpp_dense_nf(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau =0.001, double r_tau =0.001,
   double d_gamma =0.001, double r_gamma =0.001, int max_iter=1000, double th=1e-7, bool calcELB=true, bool verbose=true, int freqELB=10){

    grpRR MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
    List result =MyModel.fitModel();

    return(result); // transforms an arbitrary object into a SEXP.
}


//Fitting a normal prior model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grRRCpp_dense_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau =0.001, double r_tau =0.001,
                            double d_gamma =0.001, double r_gamma =0.001, int max_iter=1000, double th=1e-7, bool calcELB=true, bool verbose=true, int freqELB=10){

  grpRR_dense_ff MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
  List result =MyModel.fitModel();

  return(result);
}


//Fitting a spike and slab prior model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grRRCpp_sparse_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau =0.001, double r_tau =0.001,
                          double d_gamma =0.001, double r_gamma =0.001, double r_pi=1, double d_pi=1, int max_iter=1000, double th=1e-7, bool calcELB=true, bool verbose=true, int freqELB=10){

  grpRR_sparse_ff MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB, verbose, freqELB);
  List result =MyModel.fitModel();

  return(result);
}

//Fitting a normal prior logistic model with only partially factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grpRRCpp_logistic_nf(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                            double d_gamma =0.001, double r_gamma =0.001, int max_iter=1000, double th=1e-7,
                            bool calcELB=true, bool verbose=true, int freqELB=10){

  grpRR_logistic_nf MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
  List result =MyModel.fitModel();

  return(result);
}


//Fitting a normal prior logistic model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grpRRCpp_logistic_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                            double d_gamma =0.001, double r_gamma =0.001, int max_iter=1000, double th=1e-7,
                            bool calcELB=true, bool verbose=true, int freqELB=10){

  grpRR_logistic_ff MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
  List result =MyModel.fitModel();

  return(result);
}

