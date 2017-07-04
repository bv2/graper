#include <armadillo>
#include  <boost/math/special_functions/digamma.hpp>

# include <iostream>
# include <chrono>
using get_time = std::chrono::steady_clock ;

using  namespace  arma;
using namespace Rcpp;


class  grpRR_logistic_ff {

private:
  // values remaining constant
  mat  X;
  vec y;
  Row<int> annot;
  double yty;
  int p,n,g;
  vec NoPerGroup;
  double d_gamma, r_gamma;
  int max_iter;
  double th;
  bool calcELB, verbose;
  int freqELB;
  List ListOfOuterX;
  List ListXrowSquared;
  vec term4betamu;


  // changing values
  double ELB;
  vec  alpha_gamma, beta_gamma;
  vec xi;
  sp_mat Sigma_beta;
  vec sigma2_beta;
  vec mu_beta;
  vec EW_gamma;
  double diff;
  int n_iter;
  vec EW_betasq;
  double EW_leastSquares;
  vec ELB_trace;
public:

  //initaliser list
  grpRR_logistic_ff(mat X, vec y, Row<int> annot, int g, vec NoPerGroup,
       double d_gamma =0.001, double r_gamma =0.001, int max_iter=1000, double th=1e-7, bool calcELB=true, bool verbose=true,
       int freqELB =10):
  X(X)                                // design matrix
  , y(y)                                // response vector
  , annot(annot)                        // assignement of each feautre to a group
  , d_gamma(d_gamma)                    // hyperparameters of gamma distribution for gamma
  , r_gamma(r_gamma)                    // hyperparameters of gamma distribution for gamma
  , p(X.n_cols)                         //number of samples
  , n(X.n_rows)                         //number of samples
  , g(g)                                // number of groups
  , NoPerGroup(NoPerGroup)               //number of features per group
  , max_iter(max_iter)                  // maximal number of iterations
  , th(th)                              //threshold for ELBO to stop iterations
  , calcELB(calcELB)                    //whether to calculate ELBO
  , verbose(verbose)                    //whether to print intermediate messages
  , xi(n)                               // variational parameter, initialised to 0, should better be close to yX\beta
  , ELB(-std::numeric_limits<double>::infinity())                           //evidence lower bound
  , sigma2_beta(p)                      //variance parameter of normal distribution for beta
  , Sigma_beta(speye(p,p))              //diagonal covariance matrix
  , mu_beta(p)                          //initialise by 0
  , EW_gamma(g)                         //initialise by expected value of a gamma distribution, one value per group
  , alpha_gamma(g)                      //parameter of gamma distribution for tau (stays same in each iteration)
  , beta_gamma(g)
  , diff(th+1)                          // to ensure it is larger than th at the beginning
  , n_iter(0)                           // counter of iterations
  , freqELB(freqELB)                    // freuqency of ELB calculation: each freqELB-th iteration ELBO is calculated
  , ELB_trace(max_iter)
  , ListOfOuterX(n)
  , term4betamu(p)
  , ListXrowSquared(n)
  {
    EW_gamma.fill(r_gamma/d_gamma);
    alpha_gamma=r_gamma+NoPerGroup/2;
    mu_beta.fill(0);
    term4betamu.fill(0);
    xi.fill(0);
    //calculate often used quantities -slow
    for(int i = 0; i< n; i++) {
      ListOfOuterX(i) = X.row(i).t()*X.row(i);
      term4betamu = term4betamu + (y(i)-0.5)*X.row(i).t();
      ListXrowSquared(i)=X.row(i).t()%X.row(i).t();
      }
  }

  //helper functions
  double sigmoid (double x){
    return(1/(1+exp(-x)));
  }
  double lambda (double x) {
    if(x!=0) return (1/(2*x)*(sigmoid(x)-0.5));
      else return(1/8);
  }

  //main function: fit  model
  List fitModel() {
    while(n_iter<max_iter && (std::abs(diff)>th || isinf(diff) || isnan(diff))){
      iterate();
    }

    if(diff<th){
      Rcout << "ELB converged" << endl;
      ELB_trace=ELB_trace(span(0,n_iter-1));
    }
    else{
      Rcout << "Maximum numbers of iterations reached - no convergence or ELB not calculated" << endl;
    }

    List results=List::create(Named("EW_beta")=mu_beta, Named("EW_gamma")=EW_gamma, Named("ELB")=ELB,
                                    Named("alpha_gamma")=alpha_gamma,
                                    Named("beta_gamma")=beta_gamma, Named("Sigma_beta")=Sigma_beta, Named("ELB_trace")=ELB_trace);


    return(results);
  }

  //function to do one iteration
  void iterate(){
    n_iter=n_iter+1;                          //increasing counter by 1
    if(verbose) Rcout << "iteration " << n_iter << endl;

    update_param_beta();
    update_exp_beta();
    update_param_gamma();
    update_exp_gamma();
    update_param_xi();


    //optional: calculate ELB every freqELB-th step
    if(calcELB & n_iter%freqELB==0) calculate_ELBO();
    ELB_trace(n_iter-1)=ELB;

  }


  //function to calculate updated parameters for beta variational distirbution
  void update_param_beta(){
    if(verbose) Rcout << "Updating beta.." << endl;
    auto start_beta=get_time::now();

    vec gamma_annot(p);
    for(int i = 0; i< p; i++) {
      gamma_annot(i)=EW_gamma(annot(i)-1);      // minus one as annot starts counting at 1 instead of 0
    }

    vec term1(p);
    term1.fill(0);
    vec lambda_xi(n);
    for(int k = 0; k< n; k++) {
      lambda_xi(k) = lambda(xi(k));
      vec XrowSquared = ListXrowSquared(k);
      term1 = term1+lambda_xi(k) * XrowSquared;
    }
    sigma2_beta = 1/(gamma_annot+2*term1);
    Sigma_beta.diag() = sigma2_beta;


    sp_mat XiD = speye(n,n);
    XiD.diag() = lambda_xi;
    mat XTxiX = X.t() * XiD * X;

    for(int i = 0; i< p; i++){
      mu_beta(i) = sigma2_beta(i) * as_scalar(0.5 * X.col(i).t() * (2*y-1) - 2* ((XTxiX.row(i)*mu_beta - XTxiX(i,i)*mu_beta(i)) ));
    }

    auto time_beta = get_time::now() - start_beta;
    if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_beta).count()<<" ms "<<endl;
  }

  //function to calculate updated parameters for tau variational distribution
  void update_param_xi(){

    for(int i = 0; i< n; i++) {
      xi(i) = sqrt(as_scalar(X.row(i)*(Sigma_beta + mu_beta*mu_beta.t())*X.row(i).t()));
    }
  }

  //function to calculate updated parameters for gamma variational distribution
  void update_param_gamma(){

    beta_gamma.fill(d_gamma);

    for(int i = 0; i< p; i++){
      int k = annot[i]-1;                          // minus one as annot stars counting at 1 instead of 0
      beta_gamma[k]=beta_gamma[k]+0.5*EW_betasq[i];
    }
  }

  //function to update expected values of beta
  void update_exp_beta(){
    EW_betasq=square(mu_beta)+sigma2_beta;
  }


  //function to update expected values of gamma
  void update_exp_gamma(){
    EW_gamma=alpha_gamma/beta_gamma;
  }


  //function to calculate ELBO
  void calculate_ELBO(){
    if(n_iter==1) Rcout<<"ELB not implemented"<<endl;
  //   if(verbose) Rcout<<"Calculating ELB.."<<endl;
  //   auto start_ELB=get_time::now();
  //
  //   double ELB_old = ELB;
  //
  //   vec lgamma_alpha_gamma(g);
  //   vec digamma_alpha_gamma(g);
  //   for(int i =0; i<g;i++){
  //     lgamma_alpha_gamma(i)=lgamma((alpha_gamma(i)));                     // important to directly use log gamma to avoid numerical overflow
  //     digamma_alpha_gamma(i)=boost::math::digamma(alpha_gamma(i));
  //   }
  //
  //   //expected values required in addition (log of Gamma r.v.)
  //   double EW_logtau = boost::math::digamma(alpha_tau)-log(beta_tau);
  //   vec EW_loggamma = digamma_alpha_gamma-log(beta_gamma);
  //
  //   //to get EW_loggamma[annot] and EW_gamma[annot]
  //   vec EW_loggamma_annot(p);
  //   vec EW_gamma_annot(p);
  //   for(int i = 0; i< p; i++){
  //     int k = annot(i)-1;                          // minus one as annot stars counting at 1 instead of 0
  //     EW_loggamma_annot(i) = EW_loggamma(k);
  //     EW_gamma_annot(i) = EW_gamma(k);
  //   }
  //
  //   //expectation under variational density of log joint distribution
  //   double exp_logcondDy=n/2*EW_logtau -0.5*EW_tau*EW_leastSquares-n/2*log(2*M_PI);
  //   double exp_logcondDbeta=accu(0.5*EW_loggamma_annot-0.5*EW_gamma_annot%EW_betasq)-p/2*log(2*M_PI);
  //   double exp_logDgamma=accu((r_gamma-1)*EW_loggamma-d_gamma * EW_gamma)-g*lgamma(r_gamma)+g*r_gamma*log(d_gamma);
  //   double exp_logDtau=(r_tau-1)*EW_logtau-d_tau* EW_tau-lgamma(r_tau)+r_tau*log(d_tau);
  //   double exp_Djoint=exp_logcondDy+exp_logcondDbeta+exp_logDgamma+exp_logDtau;
  //
  //   //entropy of variational distribution
  //   double logdet_Sigma = real(log_det(Sigma_beta));      //replace log(det) by log_det to avoid numeric issues of Inf
  //                                                         // Are there faster ways? Better use inverse from above? Reuse Cholesky?
  //   double entropy_beta=p/2*(log(2*M_PI)+1)+0.5*logdet_Sigma;
  //   double entropy_gamma=accu(alpha_gamma-log(beta_gamma)+log(lgamma_alpha_gamma)+(1-alpha_gamma)%digamma_alpha_gamma); //replace log(tgamma) by lgamma to avoid numeric issues of Inf
  //   double entropy_tau=alpha_tau-log(beta_tau)+lgamma(alpha_tau)+(1-alpha_tau)*boost::math::digamma(alpha_tau);
  //
  //   //evidence lower bound
  //   ELB=exp_Djoint+entropy_beta +entropy_gamma+entropy_tau;
  //   diff=ELB-ELB_old;
  //
  //   auto time_ELB = get_time::now() - start_ELB;
  //   if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_ELB).count()<<" ms "<<endl;
  //
  //   if(verbose){
  //     Rcout<<"ELB="<<ELB<<endl;
  //     Rcout<<"ELB improved by "<<diff<<endl;
  //     Rcout<<endl;
  //   }
  //
   }

};
