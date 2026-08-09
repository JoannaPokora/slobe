#include <RcppEigen.h>
#include <slobe/slobe_admm.hpp>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List slobe_admm_cpp(
    const Eigen::VectorXd& start,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Xinit,
    const Eigen::VectorXd& y,
    double a_prior,
    double b_prior,
    const Eigen::MatrixXd& covariance,
    double sigma,
    double fdr,
    double tolerance,
    bool known_sigma,
    int max_iterations,
    bool verbose,
    bool bh,
    bool known_covariance
) {
  slobe::SLOBEResult result =
    slobe::slobe_admm(
      start,
      X,
      Xinit,
      y,
      a_prior,
      b_prior,
      covariance,
      sigma,
      fdr,
      tolerance,
      known_sigma,
      max_iterations,
      verbose,
      bh,
      known_covariance
    );

  return Rcpp::List::create(
    Rcpp::Named("coefficients") = result.coefficients,
    Rcpp::Named("sigma")        = result.sigma,
    Rcpp::Named("theta")        = result.theta,
    Rcpp::Named("c")            = result.c,
    Rcpp::Named("w")            = result.w,
    Rcpp::Named("converged")    = result.converged,
    Rcpp::Named("X")            = result.X,
    Rcpp::Named("Sigma")        = result.Sigma,
    Rcpp::Named("mu")           = result.mu,
    Rcpp::Named("lambda")       = result.lambda
  );
}
