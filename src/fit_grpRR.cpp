//C Functions calling the main function from classes

#include <RcppArmadillo.h>
#include "grpRRclass.hpp"
#include "grpRRclass_fullyfac.hpp"
#include "grpRR_sparseclass_fullyfac.hpp"
#include "grpRRclass_logistic.hpp"
#include "grpRRclass_logistic_fullyfac.hpp"
#include "grpRR_sparseclass_logistic_fullyfac.hpp"


// [[Rcpp::depends(RcppArmadillo)]]

//Fitting a normal prior model with only partially factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grRRCpp_dense_nf(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
   double d_gamma, double r_gamma, int max_iter, double th, bool calcELB, bool verbose, int freqELB){

    grpRR MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
    List result =MyModel.fitModel();

    return(result); // transforms an arbitrary object into a SEXP.
}


//Fitting a normal prior model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grRRCpp_dense_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
                            double d_gamma, double r_gamma, int max_iter, double th, bool calcELB, bool verbose, int freqELB, arma::vec mu_init){

  grpRR_dense_ff MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB, mu_init);
  List result =MyModel.fitModel();

  return(result);
}


//Fitting a spike and slab prior model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grRRCpp_sparse_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
                          double d_gamma, double r_gamma, double r_pi, double d_pi, int max_iter, double th, bool calcELB, bool verbose,
                             int freqELB, arma::vec mu_init, arma::vec psi_init){

  grpRR_sparse_ff MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB, verbose, freqELB, mu_init, psi_init);
  List result =MyModel.fitModel();

  return(result);
}

//Fitting a normal prior logistic model with only partially factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grpRRCpp_logistic_nf(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                            double d_gamma, double r_gamma, int max_iter, double th,
                            bool calcELB, bool verbose, int freqELB){

  grpRR_logistic_nf MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
  List result =MyModel.fitModel();

  return(result);
}


//Fitting a normal prior logistic model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grpRRCpp_logistic_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                            double d_gamma, double r_gamma, int max_iter, double th,
                            bool calcELB, bool verbose, int freqELB, arma::vec mu_init, bool intercept ){

  grpRR_logistic_ff MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB, mu_init, intercept);
  List result =MyModel.fitModel();

  return(result);
}

//Fitting a spike and slab prior logistic model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List grpRRCpp_sparse_logistic_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                                double d_gamma, double r_gamma, int max_iter, double th,
                                bool calcELB, bool verbose, int freqELB, arma::vec mu_init, arma::vec psi_init, bool intercept){
    
    grpRR_sparse_logistic_ff MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB, mu_init, psi_init,intercept);
    List result =MyModel.fitModel();
    
    return(result);
}
