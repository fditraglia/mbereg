#include <RcppArmadillo.h>
#include <stdexcept>
#include "brent.hpp"
#include "helper.h"

using namespace Rcpp;

class BCSclass {
public:
  BCSclass(const arma::vec&, const arma::vec&, const arma::vec&);
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
  arma::vec get_Qn_DR(double a1, double beta, arma::mat normal_draws) {
    arma::mat Sigma = get_Sigma(a1, beta);
    arma::mat Omega_sqrt = sqrtm_cpp(cov2cor_cpp(Sigma));
    arma::mat boot_draws = Omega_sqrt * normal_draws.t();
    arma::vec phi = get_phi(a1);
    int B = normal_draws.n_rows;
    arma::vec Tn_DR_I(B, arma::fill::zeros);
    arma::vec Tn_DR_E(B, arma::fill::zeros);
    for(int b = 0; b < B; ++b){
      arma::vec draw_b = boot_draws.col(b);
      // Zero indexing in C++!
      Tn_DR_I(b) = SS_neg(draw_b.rows(0, 1) + phi);
      Tn_DR_E(b) = SS(draw_b.rows(2, 3));
    }
    return Tn_DR_I + Tn_DR_E;
  }
  double get_Qn_PR(double a1, double beta, arma::vec normal_draw) {
    arma::mat Sigma = get_Sigma(a1, beta);
    arma::mat Omega_sqrt = sqrtm_cpp(cov2cor_cpp(Sigma));
    arma::vec boot_draw = Omega_sqrt * normal_draw;
    arma::vec mbar_I = get_mbar_I(a1);
    arma::vec ell_I = sqrt(n) * mbar_I / (kappa_n * s_II);
    arma::vec mbar_E = get_mbar_E(a1, beta);
    double Tn_PR_I = SS_neg(boot_draw.rows(0, 1) + ell_I); // zero-indexing!
    arma::vec s_EE = get_s_EE(a1, beta);
    arma::vec ell_E = sqrt(n) * mbar_E / (kappa_n * s_EE);
    double Tn_PR_E = SS(boot_draw.rows(2, 3) + ell_E); // zero-indexing!
    return Tn_PR_I + Tn_PR_E;
  }
//private:
  int n;
  double kappa_n, q;
  arma::mat W, W_tilde, Sigma_II, M_EE, M_IE;
  arma::vec W_tilde_bar, w1_z0, w1_z1, w_z_bar, s_II;
};
//Class constructor
BCSclass::BCSclass(const arma::vec& y, const arma::vec& Tobs, const arma::vec& z) {
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
List testy(arma::vec y, arma::vec Tobs, arma::vec z, arma::mat normal_draws){
  BCSclass myBCS(y, Tobs, z);
  double a1 = 0.1;
  double beta = 1.0;
  arma::mat Nu = myBCS.get_Nu(a1, beta);
  arma::mat Sigma_EE = myBCS.get_Sigma_EE(a1, beta);
  arma::mat Sigma = myBCS.get_Sigma(a1, beta);
  double Qn = myBCS.get_Qn(a1, beta);
  arma::vec phi = myBCS.get_phi(a1);
  arma::vec Tn_DR = myBCS.get_Qn_DR(a1, beta, normal_draws);
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
                      Named("phi") = phi,
                      Named("Tn_DR") = Tn_DR);
}

// Functor for the optimization of Qn from which we calculate Tn
class Qn_Functor: public brent::func_base
{
public:
  Qn_Functor (BCSclass myBCS_, double beta_null_) : myBCS(myBCS_), beta_null(beta_null_) {}
  double operator() (double a1) {
    return myBCS.get_Qn(a1, beta_null);
    }
private:
  BCSclass myBCS;
  double beta_null;
};

// Functor for the optimization of Qn_DR from which we calculate Tn_DR
class Qn_PR_Functor: public brent::func_base
{
public:
  Qn_PR_Functor (BCSclass myBCS_, double beta_null_,
                 arma::vec normal_draw_) : myBCS(myBCS_),
                                          beta_null(beta_null_),
                                          normal_draw(normal_draw_){}
  double operator() (double a1) {
    return myBCS.get_Qn_PR(a1, beta_null, normal_draw);
    }
private:
  BCSclass myBCS;
  double beta_null;
  arma::vec normal_draw;
};

// [[Rcpp::export]]
double test_Qn(double a1, arma::vec y, arma::vec Tobs, arma::vec z){
  BCSclass myBCS(y, Tobs, z);
  Qn_Functor myQn(myBCS, 1.0);
  return myQn(a1);
}

// [[Rcpp::export]]
List test_Qn_opt(arma::vec y, arma::vec Tobs, arma::vec z){
  BCSclass myBCS(y, Tobs, z);
  Qn_Functor myQn(myBCS, 1.0);
  double a1_star;
  double Qn_star = brent::local_min(0.0, 0.99, 0.0001, myQn, a1_star);
  return List::create(Named("minimum") = a1_star, Named("objective") = Qn_star);
}

// [[Rcpp::export]]
List BCS_test_cpp(double beta_null, arma::vec y, arma::vec Tobs, arma::vec z,
                  arma::mat normal_draws) {

  // Preliminary calculations and intitializations
  BCSclass BCS(y, Tobs, z);
  double opt_tol = 0.0001;
  double a1_lower = 0.0;
  double a1_upper = 0.99;

  // Add beta_null == 0 case later

  // Calculate Tn
  Qn_Functor Qn(BCS, beta_null);
  double a1_star;
  double Tn = brent::local_min(a1_lower, a1_upper, opt_tol, Qn, a1_star);

  // Calculate bootstrap draws of Tn_DR
  arma::vec Tn_DR = BCS.get_Qn_DR(a1_star, beta_null, normal_draws);

  // Calculate bootstrap draws of Tn_PR
  int B = normal_draws.n_rows;
  arma::vec Tn_PR(B, arma::fill::zeros);
  for(int b = 0; b < B; ++b){
    arma::vec normal_draw_b = normal_draws.row(b).t();
    Qn_PR_Functor Qn_DR_b(BCS, beta_null, normal_draw_b);
    double a1_star_b;
    Tn_PR(b) =  brent::local_min(a1_lower, a1_upper, opt_tol, Qn_DR_b, a1_star_b);
  }

  return List::create(Named("a1_star") = a1_star,
                      Named("Tn") = Tn,
                      Named("Tn_DR") = Tn_PR,
                      Named("Tn_PR") = Tn_DR);
}

