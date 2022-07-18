#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace RcppArmadillo;


float normal_pdf(double x, double m, double s){
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

// [[Rcpp::export]]
NumericVector jpr(NumericVector eps, NumericVector htt, double h0m, double h0V, double h_mu, double h_rho, double h_sig, double theta, double tau2, NumericVector z, int T, double c) {
  double h_rho2 = pow(h_rho,2);
  NumericVector rnum1 = rnorm(1);
  NumericVector rnum2 = rnorm(T);
  NumericVector rnum3 = log(runif(T));
  
  for (int i=0;i<T;++i){
    double ht_mu;
    double ht_sig;
    
    if(i==0){
      double h0t_sig = 1/(h0m/h0V + h_rho2/h_sig);
      double h0t_mu = h0t_sig * (h0m/h0V + h_rho * (htt[0]-h_mu)/h_sig);
      double h0d = h0t_mu + (sqrt(h0t_sig) * rnum1[0]);
      
      ht_mu = ((1-h_rho)*h_mu + h_rho * (h0d + htt[1]))/(1+h_rho2);
      ht_sig = h_sig/(1+h_rho2);
    }else if(i==(T-1)){
      ht_mu = h_mu + (h_rho * htt[i]);
      ht_sig = h_sig;
    }else{
      ht_mu = ((1-h_rho)*h_mu + h_rho*(htt[i-1] + htt[i+1]))/(1+h_rho2);
      ht_sig = h_sig/(1+h_rho2);
    }
    
    double ht_prop = htt[i] + (sqrt(c) * rnum2[i]);
    double liki_new = log(normal_pdf(eps[i],theta*z[i]*exp(ht_prop),sqrt(tau2*z[i])*exp(ht_prop))) + log(normal_pdf(ht_prop,ht_mu,sqrt(ht_sig)));
    double liki_old = log(normal_pdf(eps[i],theta*z[i]*exp(htt[i]),sqrt(tau2*z[i])*exp(htt[i]))) + log(normal_pdf(htt[i],ht_mu,sqrt(ht_sig)));
    
    double alpha = liki_new-liki_old;
    if(rnum3[i]<alpha){
      htt[i] = ht_prop;
    }
  }
  return wrap(htt);
}
