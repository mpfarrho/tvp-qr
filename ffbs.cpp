#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace RcppArmadillo;
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::export]]
NumericVector ffbs(NumericMatrix y, NumericMatrix Z,NumericMatrix Ht, NumericMatrix Qtt,int m, int p, int t, NumericVector B0, NumericMatrix V0) {
  static Rcpp::Function asVector("as.vector");
  arma::vec bp = Rcpp::as<arma::vec>(B0);
  arma::mat Vp = Rcpp::as<arma::mat>(V0);
  arma::mat Qt =Rcpp::as<arma::mat>(Qtt);
  arma::mat bt(t,m);
  arma::mat Vt(pow(m,2),t);
  arma::mat R(p,m);
  arma::mat H(t*m,p);
  arma::mat cfe(p,1);
  arma::mat yt(p,1);
  arma::mat f(p,p);
  arma::mat inv_f(p,p);
  arma::mat btt(m,1);
  arma::mat Vtt(m,m);
  arma::cube test(m,m,t);
  
  for (int i=1;i<(t+1);++i){
    arma::mat Qt =diag(Qtt(i-1,_));
    R= Ht(Range((i-1)*p,(i*p)-1),_);
    H = Z(Range((i-1)*p,(i*p)-1),_);
    yt=y(_,i-1);
    cfe= yt-H*bp;
    f = H*Vp*trans(H)+R;
    inv_f =  inv(f);
    btt = bp+Vp*trans(H)*inv_f*cfe; 
    Vtt = Vp-Vp*trans(H)*inv_f*H*Vp;
    if ((i-1)<t){
      bp=btt;
      Vp=Vtt+Qt;
    }
    bt.row(i-1) = trans(btt);
    test.slice(i-1)=Vtt;
  }
  arma::mat bdraw(t,m);
  arma::mat Y = arma::randn(m, 1);
  arma::mat bmean(m,1);
  arma::mat bvar(m,m);
  bdraw.row(t-1)=trans(btt+ trans(arma::chol(Vtt))*Y);
  
  //backward recurssions
  arma::mat bf(1,m);
  for (int i=1;i<t;++i){
    arma::mat Qt =diag(Qtt(t-i,_));
    bf=trans(bdraw.row(t-i));
    btt=trans(bt.row(t-i-1));
    Vtt=test.slice(t-i-1);
    f=Vtt+Qt;
    inv_f=inv(f);
    
    cfe=bf-btt;
    bmean=btt+Vtt*inv_f*cfe;
    bvar=Vtt-Vtt*inv_f*Vtt;
    arma::mat Y = arma::randn(m, 1);
    bdraw.row(t-i-1)=trans(bmean+trans(arma::chol(bvar))*Y);
  }
  bdraw=trans(bdraw);
  return wrap(bdraw);
}

