#include <RcppArmadillo.h>
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
