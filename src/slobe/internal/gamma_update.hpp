#pragma once

#include <Eigen/Core>

namespace slobe {
    
void gamma_mean_update(
    const Eigen::VectorXd& abs_beta_ord,
    const Eigen::VectorXd& lambda,
    int p,
    Eigen::VectorXd& gamma_h
);

}
