#pragma once

#include <Eigen/Core>

namespace slobe {

void create_lambda(
    Eigen::VectorXd& lam,
    int p,
    double FDR,
    bool BH
);

void div_X_by_w(
    Eigen::MatrixXd& X_div_w,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& w,
    int n, int p
);

Eigen::VectorXd slope_libslope_fit(
    Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& lambda_sorted_decreasing
);

}
