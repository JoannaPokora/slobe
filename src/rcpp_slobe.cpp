#include <slobe/slobe_fit.hpp>

#include <RcppEigen.h>
#include <string>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List slobe_fit_cpp(
    const Eigen::VectorXd& start,
    const Eigen::MatrixXd& Xinit,
    const Eigen::VectorXd& y,
    Rcpp::Nullable<Rcpp::NumericMatrix> Xmis,
    const std::string& family,
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
  Eigen::MatrixXd Xmis_cpp;
  const Eigen::MatrixXd* Xmis_ptr = nullptr;
  
  if (Xmis.isNotNull()) {
    Rcpp::NumericMatrix tmp(Xmis);
    Xmis_cpp = Rcpp::as<Eigen::MatrixXd>(tmp);
    Xmis_ptr = &Xmis_cpp;
  }
  
  slobe::SLOBEResult result =
    slobe::slobe_fit(
      start,
      Xinit,
      y,
      Xmis_ptr,
      family,
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
    Rcpp::Named("sigma") = result.sigma,
    Rcpp::Named("theta") = result.theta,
    Rcpp::Named("c") = result.c,
    Rcpp::Named("w") = result.w,
    Rcpp::Named("converged") = result.converged,
    Rcpp::Named("X") = result.X,
    Rcpp::Named("Sigma") = result.Sigma,
    Rcpp::Named("mu") = result.mu,
    Rcpp::Named("lambda") = result.lambda
  );
}
