#include <RcppArmadillo.h>
#include <cassert>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat center(arma::mat M) {
  return M - arma::repmat(arma::mean(M), M.n_rows, 1);
}

// [[Rcpp::export]]
arma::mat mycov(arma::mat M1, arma::mat M2) {
// Armadillo's covariance function requires that both matrices be of the same
// dimension. Mine, like the one in base R, does not.
  assert(M1.n_rows == M2.n_rows);
  int n = M1.n_rows;
  return center(M1).t() * center(M2) / (n - 1);
}

// [[Rcpp::export]]
arma::vec foo(arma::vec z) {
  return arma::ones(z.n_elem) - z;
}


class BCS {
public:
  BCS(const arma::vec&, const arma::vec&, const arma::vec&, const arma::mat&);
//private:
  int n;
  double kappa_n, q;
  arma::mat W, W_tilde, boot_raw, Sigma_II, M_EE, M_IE;
  arma::vec W_tilde_bar, w1_z0, w1_z1, w_z_bar, s_II;
};
//Class constructor
BCS::BCS(const arma::vec& y, const arma::vec& Tobs, const arma::vec& z,
         const arma::mat& normal_draws) {
  q = arma::mean(z);
  kappa_n = sqrt(log(n));
  arma::vec w1 = Tobs;
  arma::vec w2 = y;
  arma::vec w3 = 2.0 * y % Tobs;
  arma::vec w4 = pow(y, 2.0);
  W = arma::join_horiz(arma::join_horiz(w1, w2), arma::join_horiz(w3, w4));
  W_tilde = arma::repmat(z, 1, W.n_cols) % center(W);
  W_tilde_bar = arma::mean(W_tilde).t();
  M_EE = arma::cov(W_tilde - q * W);
  arma::vec z0 = arma::ones(z.n_elem) - z;
  w1_z0 = w1 % z0 / (1 - q);
  w1_z1 = w1 % z / q;
  w_z_bar << arma::mean(w1_z0) << arma::endr
          << arma::mean(w1_z1);
  Sigma_II = arma::cov(arma::join_horiz(w1_z0, w1_z1));
  s_II = arma::sqrt(Sigma_II.diag());
  M_IE = arma::join_vert(mycov(w1_z0, q * W) - mycov(w1_z0, W_tilde),
                         mycov(w1_z1, q * W) - mycov(w1_z1, W_tilde));
}

// [[Rcpp::export]]
List testy(arma::vec y, arma::vec Tobs, arma::vec z){
  BCS myBCS(y, Tobs, z, z);
  return List::create(Named("W_tilde_bar") = myBCS.W_tilde_bar,
                      Named("M_EE") = myBCS.M_EE,
                      Named("w1_z0") = myBCS.w1_z0,
                      Named("w1_z1") = myBCS.w1_z1,
                      Named("w_z_bar") = myBCS.w_z_bar,
                      Named("Sigma_II") = myBCS.Sigma_II,
                      Named("s_II") = myBCS.s_II,
                      Named("M_IE") = myBCS.M_IE);
}
