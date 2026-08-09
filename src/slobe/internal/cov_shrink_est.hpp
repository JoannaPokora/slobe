#pragma once

#include <Eigen/Core>

namespace slobe {

void linshrink_cov(
    const Eigen::MatrixXd& X,
    Eigen::MatrixXd& S,
    int n,
    int p
);

}