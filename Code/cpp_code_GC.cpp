# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define pi 3.14159265358979323846

// [[Rcpp::export]]
double matern(double d, double sigma2, double phi, double nu){
  double hphi = d / phi;
  if(d==0){
    return sigma2;
  }else{
    return sigma2 / (Rf_gammafn(nu) * pow(2,(nu-1)) ) * pow(hphi,nu) * Rf_bessel_k(hphi,nu,1);
  }
}

// [[Rcpp::export]]
arma::mat matern_stx_cpp(double sigma2, double phi, double rho, double psi, double nu, arma::vec t1, arma::vec t2, arma::mat X1, arma::mat X2, arma::mat locs1, arma::mat locs2, bool sq=false){
  double r = locs1.n_rows;
  double c = locs2.n_rows;
  arma::mat Sigma(r,c);
  double d;
  
  if(sq){
    for(int i=0; i < r; i++){
      for(int j=0; j <= i; j++){
        if(i==j){
          Sigma(i,j) = sigma2;
        }else{
          d = pow(pow(locs1(i,0)/phi - locs2(j,0)/phi, 2) + pow(locs1(i,1)/phi - locs2(j,1)/phi, 2) + pow(t1(i)/rho - t2(j)/rho,2) + sum(pow(X1.row(i)/psi - X2.row(j)/psi,2)),0.5);
          Sigma(i,j) = Sigma(j,i) = matern(d, sigma2, 1, nu);
        }
      }
    }
  }else{
    for(int i=0; i < r; i++){
      for(int j=0; j < c; j++){
        d = pow(pow(locs1(i,0)/phi - locs2(j,0)/phi, 2) + pow(locs1(i,1)/phi - locs2(j,1)/phi, 2) + pow(t1(i)/rho - t2(j)/rho,2) + sum(pow(X1.row(i)/psi - X2.row(i)/psi,2)),0.5);
        Sigma(i,j) = matern(d, sigma2, 1, nu);
      }
    }
  }
  return Sigma;
}  

// [[Rcpp::export]]
arma::mat matern_st_cpp(double sigma2, double phi, double rho, double nu, arma::vec t1, arma::vec t2, arma::mat locs1, arma::mat locs2, bool sq=false){
  double r = locs1.n_rows;
  double c = locs2.n_rows;
  arma::mat Sigma(r,c);
  double d;
  
  if(sq){
    for(int i=0; i < r; i++){
      for(int j=0; j <= i; j++){
        if(i==j){
          Sigma(i,j) = sigma2;
        }else{
          d = pow(pow(locs1(i,0)/phi - locs2(j,0)/phi, 2) + pow(locs1(i,1)/phi - locs2(j,1)/phi, 2) + pow(t1(i)/rho - t2(j)/rho,2),0.5);
          Sigma(i,j) = Sigma(j,i) = matern(d, sigma2, 1, nu);
        }
      }
    }
  }else{
    for(int i=0; i < r; i++){
      for(int j=0; j < c; j++){
        d = pow(pow(locs1(i,0)/phi - locs2(j,0)/phi, 2) + pow(locs1(i,1)/phi - locs2(j,1)/phi, 2) + pow(t1(i)/rho - t2(j)/rho,2),0.5);
        Sigma(i,j) = matern(d, sigma2, 1, nu);
      }
    }
  }
  return Sigma;
}  

// [[Rcpp::export]]
arma::mat matern_s_cpp(double sigma2, double phi, double nu, arma::mat locs1, arma::mat locs2, bool sq=false){
  double r = locs1.n_rows;
  double c = locs2.n_rows;
  arma::mat Sigma(r,c);
  double d;
  
  if(sq){
    for(int i=0; i < r; i++){
      for(int j=0; j <= i; j++){
        if(i==j){
          Sigma(i,j) = sigma2;
        }else{
          d = pow(pow(locs1(i,0)/phi - locs2(j,0)/phi, 2) + pow(locs1(i,1)/phi - locs2(j,1)/phi, 2),0.5);
          Sigma(i,j) = Sigma(j,i) = matern(d, sigma2, 1, nu);
        }
      }
    }
  }else{
    for(int i=0; i < r; i++){
      for(int j=0; j < c; j++){
        d = pow(pow(locs1(i,0)/phi - locs2(j,0)/phi, 2) + pow(locs1(i,1)/phi - locs2(j,1)/phi, 2),0.5);
        Sigma(i,j) = matern(d, sigma2, 1, nu);
      }
    }
  }
  return Sigma;
}  

// [[Rcpp::export]]  
List krig_STX_cpp(arma::vec par,
                 arma::mat S1, arma::vec t1, arma::mat X1,
                 arma::mat S2, arma::vec t2, arma::mat X2, 
                 arma::colvec Y2, arma::colvec Xb1, arma::colvec Xb2){
  // arma::mat d11 = as<arma::mat>(d["d11"]);
  int N = Y2.n_elem;
  int n = S1.n_rows;
  
  double sigma2 = par(0); // space-time variance
  double phi = par(1); // spatial range
  double rho = par(2); // temporal range
  double psi = par(3); // covariate range
  double tau2 = par(4); // nugget
  
  arma::mat S22(N, N);
  arma::mat S12(1, N);
  arma::mat err(N,N);  
  err.eye();
  err *= tau2;
  
  S22 = matern_stx_cpp(sigma2, phi, rho, psi, 0.5, t2, t2, X2, X2, S2, S2, true) + err;
  
  arma::mat Q = inv(S22);
  arma::vec var(n); 
  arma::vec pred(n);
  arma::mat QS12(N, 1); 
  

  
  for(int i=0; i < n; i++){
    S12 = matern_stx_cpp(sigma2, phi, rho, psi, 0.5, t1.row(i), t2, X1.row(i), X2, S1.row(i), S2, false);

    QS12 = Q*S12.t();
  
    pred(i) = as<double>(wrap(S12*Q*(Y2-Xb2) + Xb1(i)));
    
    var(i) =  as<double>(wrap(sigma2 + tau2 - S12*QS12));
  }
  return List::create(
    Named("pred") = pred,
    Named("var") = var);
}

  