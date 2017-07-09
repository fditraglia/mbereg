#include <RcppArmadillo.h>
#include <stdexcept>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat center(arma::mat M) {
  return M - arma::repmat(arma::mean(M), M.n_rows, 1);
}

// [[Rcpp::export]]
arma::mat mycov(arma::mat M1, arma::mat M2) {
// Armadillo's covariance function requires that both matrices be of the same
// dimension. Mine, like the one in base R, does not.
  if(M1.n_rows != M2.n_rows) {
    throw std::invalid_argument("received non-square matrix");
  }
  int n = M1.n_rows;
  return center(M1).t() * center(M2) / (n - 1);
}

// [[Rcpp::export]]
double SS(arma::vec v) {
// Returns the sum of squares of a vector v
  arma::vec v_sq = arma::pow(v, 2.0);
  return arma::as_scalar(arma::sum(v_sq));
}

// [[Rcpp::export]]
double SS_neg(arma::vec v) {
// Returns the sum of squares of the subset of v whose elements are negative
  arma::vec v_neg = v.elem(arma::find(v < 0.0));
  return SS(v_neg);
}

// [[Rcpp::export]]
arma::mat sqrtm_cpp(arma::mat M){
// Square root of symmetric, positive semi-definite matrix via eigen-decomposition
  double tol = 1e-06;
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, M);
  if(arma::any(eigval < -tol * fabs(eigval.max()))) {
    throw std::invalid_argument("Matrix is not positive semi-definite!");
  }
  arma::vec eigval_trunc = arma::max(eigval, arma::zeros(eigval.n_elem));
  return eigvec * arma::diagmat(arma::sqrt(eigval_trunc));
}

// [[Rcpp::export]]
arma::mat cov2cor_cpp(arma::mat V){
  arma::vec Is = arma::sqrt(1 / V.diag());
  if(!arma::is_finite(Is)){
    throw std::invalid_argument("One or more variances is zero or NaN!");
  }
  return arma::diagmat(Is) * V * arma::diagmat(Is);
}

// [[Rcpp::export]]
arma::mat foo(arma::mat V) {
  arma::vec Is = arma::sqrt(1 / V.diag());
  if(!arma::is_finite(Is)){
    throw std::invalid_argument("One or more variances is zero or NaN!");
  }
  return arma::diagmat(Is) * V * arma::diagmat(Is);
}


class BCS {
public:
  BCS(const arma::vec&, const arma::vec&, const arma::vec&, const arma::mat&);
  arma::mat get_Nu(double a1, double beta) {
    double theta1 = beta / (1 - a1);
    double theta2 = beta * theta1;
    arma::vec nu1, nu2;
    nu1 << -theta1 << arma::endr
        << 1.0 << arma::endr
        << 0.0 << arma::endr
        << 0.0 << arma::endr;
    nu2 << theta2 << arma::endr
        << 0.0 << arma::endr
        << -theta1 << arma::endr
        << 1.0 << arma::endr;
    return arma::join_horiz(nu1, nu2);
  }
  arma::mat get_Sigma_EE(double a1, double beta) {
    arma::mat Nu = get_Nu(a1, beta);
    return Nu.t() * M_EE * Nu;
  }
  arma::mat get_s_EE(double a1, double beta) {
    arma::mat Sigma_EE = get_Sigma_EE(a1, beta);
    return arma::sqrt(Sigma_EE.diag());
  }
  arma::mat get_Sigma(double a1, double beta) {
    arma::mat Nu = get_Nu(a1, beta);
    arma::mat Sigma_EE = get_Sigma_EE(a1, beta);
    arma::mat Sigma_IE = M_IE * Nu;
    return arma::join_vert(arma::join_horiz(Sigma_II, Sigma_IE),
                           arma::join_horiz(Sigma_IE.t(), Sigma_EE));
  }
  arma::vec get_mbar_I(double a1) {
    arma::vec ones_vec = arma::ones(w_z_bar.n_elem);
    return ones_vec - a1 * ones_vec - w_z_bar;
  }
  arma::vec get_mbar_E(double a1, double beta) {
    arma::mat Nu = get_Nu(a1, beta);
    return Nu.t() * W_tilde_bar;
  }
  double get_Qn(double a1, double beta) {
    arma::vec mbar_I = get_mbar_I(a1);
    double Tn_I = SS_neg(sqrt(n) * mbar_I / s_II);
    arma::vec mbar_E = get_mbar_E(a1, beta);
    arma::vec s_EE = get_s_EE(a1, beta);
    double Tn_E = SS(sqrt(n) * mbar_E / s_EE);
    return Tn_I + Tn_E;
  }
  arma::vec get_phi(double a1) {
    arma::vec mbar_I = get_mbar_I(a1);
    arma::vec test_stats = sqrt(n) * mbar_I / s_II;
    arma::vec phi(mbar_I.n_elem, arma::fill::zeros);
    phi.elem(find(test_stats > kappa_n)).fill(arma::datum::inf);
    return phi;
  }
  double get_Qn_DR(double a1, double beta) {
    arma::vec phi = get_phi(a1);
    return 0;
  }
  double get_Qn_PR_i(double a1, double beta, int i) {
    return 0;
  }
//private:
  int n;
  double kappa_n, q;
  arma::mat W, W_tilde, boot_draws, Sigma_II, M_EE, M_IE;
  arma::vec W_tilde_bar, w1_z0, w1_z1, w_z_bar, s_II;
};
//Class constructor
BCS::BCS(const arma::vec& y, const arma::vec& Tobs, const arma::vec& z,
         const arma::mat& normal_draws) {
  n = y.n_elem;
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
  double a1 = 0.5;
  double beta = 1.0;
  arma::mat Nu = myBCS.get_Nu(a1, beta);
  arma::mat Sigma_EE = myBCS.get_Sigma_EE(a1, beta);
  arma::mat Sigma = myBCS.get_Sigma(a1, beta);
  double Qn = myBCS.get_Qn(a1, beta);
  arma::vec phi = myBCS.get_phi(a1);
  return List::create(Named("W_tilde_bar") = myBCS.W_tilde_bar,
                      Named("M_EE") = myBCS.M_EE,
                      Named("w1_z0") = myBCS.w1_z0,
                      Named("w1_z1") = myBCS.w1_z1,
                      Named("w_z_bar") = myBCS.w_z_bar,
                      Named("Sigma_II") = myBCS.Sigma_II,
                      Named("s_II") = myBCS.s_II,
                      Named("M_IE") = myBCS.M_IE,
                      Named("Nu") = Nu,
                      Named("Sigma_EE") = Sigma_EE,
                      Named("Sigma") = Sigma,
                      Named("Qn") = Qn,
                      Named("phi") = phi);
}
