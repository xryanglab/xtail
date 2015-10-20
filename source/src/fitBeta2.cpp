//copied from DESeq2
//alpha_hatSEXP are matrix
// include RcppArmadillo and Rcpp
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// user includes
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

// fit the Negative Binomial GLM.
// note: the betas are on the natural log scale
// change: alpha_hatSEXP are n * 2 matrix
// [[Rcpp::export]]
Rcpp::List fitBeta2(SEXP ySEXP, SEXP xSEXP, SEXP nfSEXP, SEXP alpha_hatSEXP, SEXP contrastSEXP, SEXP beta_matSEXP, SEXP lambdaSEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP useQRSEXP) {
  arma::mat y = Rcpp::as<arma::mat>(ySEXP);
  arma::mat nf = Rcpp::as<arma::mat>(nfSEXP);
  arma::mat x = Rcpp::as<arma::mat>(xSEXP);
  int y_n = y.n_rows;
  int y_m = y.n_cols;
  int x_p = x.n_cols;
  arma::mat alpha_hat = Rcpp::as<arma::mat>(alpha_hatSEXP);
  arma::mat beta_mat = Rcpp::as<arma::mat>(beta_matSEXP);
  arma::mat beta_var_mat = arma::zeros(beta_mat.n_rows, beta_mat.n_cols);
  arma::mat contrast_num = arma::zeros(beta_mat.n_rows, 1);
  arma::mat contrast_denom = arma::zeros(beta_mat.n_rows, 1);
  arma::mat hat_matrix = arma::zeros(x.n_rows, x.n_rows);
  arma::mat hat_diagonals = arma::zeros(y.n_rows, y.n_cols);
  arma::colvec lambda = Rcpp::as<arma::colvec>(lambdaSEXP);
  arma::colvec contrast = Rcpp::as<arma::colvec>(contrastSEXP);
  int maxit = Rcpp::as<int>(maxitSEXP);
  arma::colvec yrow, nfrow, beta_hat, mu_hat, z;
  arma::mat w, ridge, sigma;
  arma::colvec alpha_mu;
  // vars for QR
  bool useQR = Rcpp::as<bool>(useQRSEXP);
  arma::colvec gamma_hat, big_z;
  arma::vec big_w_diag;
  arma::mat weighted_x_ridge, q, r, big_w;
  // deviance, convergence and tolerance
  double dev, dev_old, conv_test;
  double tol = Rcpp::as<double>(tolSEXP);
  double large = 30.0;
  Rcpp::NumericVector iter(y_n);
  Rcpp::NumericVector deviance(y_n);
  // bound the estimated count, as weights include 1/mu
  double minmu = 0.1;
  for (int i = 0; i < y_n; i++) {
    Rcpp::checkUserInterrupt();
    nfrow = nf.row(i).t();
    yrow = y.row(i).t();
    beta_hat = beta_mat.row(i).t();
    mu_hat = nfrow % exp(x * beta_hat);
    for (int j = 0; j < y_m; j++) {
      mu_hat(j) = fmax(mu_hat(j), minmu);
    }
    ridge = diagmat(lambda);
    dev = 0.0;
    dev_old = 0.0;
    if (useQR) {
      // make an orthonormal design matrix including
      // the ridge penalty
      for (int t = 0; t < maxit; t++) {
	iter(i)++;
	alpha_mu = mu_hat;
	for (unsigned int c = 0; c < alpha_hat.n_cols; c++) {
		alpha_mu[c] *= alpha_hat.row(i)[c];
	}
	w = diagmat(mu_hat/(1.0 + alpha_mu));
	// prepare matrices
	weighted_x_ridge = join_cols(sqrt(w) * x, sqrt(ridge));
	qr(q, r, weighted_x_ridge);
	big_w_diag = arma::ones(y_m + x_p);
	big_w_diag(arma::span(0, y_m - 1)) = diagvec(w);
	big_w = diagmat(big_w_diag);
	big_z = arma::zeros(y_m + x_p);
	z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
	big_z(arma::span(0,y_m - 1)) = z;
	// IRLS with Q matrix for X    
	gamma_hat = q.t() * sqrt(big_w) * big_z;
	solve(beta_hat, r, gamma_hat);
	if (sum(abs(beta_hat) > large) > 0) {
	  iter(i) = maxit;
	  break;
	}
	mu_hat = nfrow % exp(x * beta_hat);
	for (int j = 0; j < y_m; j++) {
	  mu_hat(j) = fmax(mu_hat(j), minmu);
	}
	dev = 0.0;
	for (int j = 0; j < y_m; j++) {
	  // note the order for Rf_dnbinom_mu: x, sz, mu, lg
	  dev = dev + -2.0 * Rf_dnbinom_mu(yrow[j], 1.0/alpha_hat.row(i)[j], mu_hat[j], 1);
	}
	conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
	if (std::isnan(conv_test)) {
	  iter(i) = maxit;
	  break;
	}
	if ((t > 0) & (conv_test < tol)) {
	  break;
	}
	dev_old = dev;
      }
    } else {
      // use the standard design matrix x
      // and matrix inversion
      for (int t = 0; t < maxit; t++) {
	iter(i)++;
	alpha_mu = mu_hat;
	for (unsigned int c=0; c<alpha_hat.n_cols; c++){
		alpha_mu[c] *= alpha_hat.row(i)[c];
	}

	w = diagmat(mu_hat/(1.0 + alpha_mu));
	z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
	solve(beta_hat, x.t() * w * x + ridge, x.t() * w * z);
	if (sum(abs(beta_hat) > large) > 0) {
	  iter(i) = maxit;
	  break;
	}
	mu_hat = nfrow % exp(x * beta_hat);
	for (int j = 0; j < y_m; j++) {
	  mu_hat(j) = fmax(mu_hat(j), minmu);
	}
	dev = 0.0;
	for (int j = 0; j < y_m; j++) {
	  // note the order for Rf_dnbinom_mu: x, sz, mu, lg
	  dev = dev + -2.0 * Rf_dnbinom_mu(yrow[j], 1.0/alpha_hat.row(i)[j], mu_hat[j], 1);
	}
	conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
	if (std::isnan(conv_test)) {
	  iter(i) = maxit;
	  break;
	}
	if ((t > 0) & (conv_test < tol)) {
	  break;
	}
	dev_old = dev;
      }
    }
    deviance(i) = dev;
    beta_mat.row(i) = beta_hat.t();
    // recalculate w so that this is identical if we start with beta_hat
	alpha_mu = mu_hat;
	for (unsigned int c=0; c<alpha_hat.n_cols; c++){
		alpha_mu[c] *= alpha_hat.row(i)[c];
	}
    w = diagmat(mu_hat/(1.0 + alpha_mu));
    hat_matrix = sqrt(w) * x * (x.t() * w * x + ridge).i(true) * x.t() * sqrt(w);
    hat_diagonals.row(i) = diagvec(hat_matrix).t();
    // sigma is the covariance matrix for the betas
    sigma = (x.t() * w * x + ridge).i(true) * x.t() * w * x * (x.t() * w * x + ridge).i(true);
    contrast_num.row(i) = contrast.t() * beta_hat;
    contrast_denom.row(i) = sqrt(contrast.t() * sigma * contrast);
    beta_var_mat.row(i) = diagvec(sigma).t();
  }

  return Rcpp::List::create(Rcpp::Named("beta_mat",beta_mat),
			    Rcpp::Named("beta_var_mat",beta_var_mat),
			    Rcpp::Named("iter",iter),
			    Rcpp::Named("hat_diagonals",hat_diagonals),
			    Rcpp::Named("contrast_num",contrast_num),
			    Rcpp::Named("contrast_denom",contrast_denom),
			    Rcpp::Named("deviance",deviance));
}
 
