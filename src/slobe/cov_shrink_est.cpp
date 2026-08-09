#include "internal/cov_shrink_est.hpp"

namespace slobe {

// Ledoit-Wolf linear shrinkage covariance
void linshrink_cov(const Eigen::MatrixXd& X,
                          Eigen::MatrixXd& S,
                          int n, int p) {
    // means per column
    Eigen::RowVectorXd means = X.colwise().sum() / (double)n;
  
    // S = X'X then centered
    S = X.transpose() * X;
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
            S(i, j) -= n * means(i) * means(j);
            S(i, j) /= (double)n;
        }
    }
  
    double m = S.trace() / (double)p;
    double d2 = S.norm();              // Frobenius
    double b_bar2 = d2 * d2;
  
    d2 = (b_bar2 - p * m * m) / (double)p;
    b_bar2 *= n;
  
    double prod;
    double sum_prod;
  
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
            sum_prod = 0.0;
            for (int k = 0; k < n; ++k) {
                prod = (X(k, i) - means(i)) * (X(k, j) - means(j));
                b_bar2 += prod * prod;
                sum_prod += prod;
            }
            b_bar2 -= 2.0 * S(i, j) * sum_prod;
        }
    }
  
    b_bar2 /= p * (n - 1.0) * (n - 1.0);
    double b2 = (b_bar2 < d2) ? b_bar2 : d2;
    double a2 = d2 - b2;
  
    S *= (a2 / d2);
    m *= (b2 / d2);
    for (int i = 0; i < p; ++i) S(i, i) += m;
}

}