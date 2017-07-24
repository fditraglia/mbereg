#ifndef HELPER_H
#define HELPER_H

arma::mat center(arma::mat M);

arma::mat mycov(arma::mat M1, arma::mat M2);

double SS(arma::vec v);

double SS_neg(arma::vec v);

arma::mat sqrtm_cpp(arma::mat M);

arma::mat cov2cor_cpp(arma::mat V);

#endif
