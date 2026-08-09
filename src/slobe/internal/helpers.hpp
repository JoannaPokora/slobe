#pragma once

#include <Eigen/Core>
#include <vector>

namespace slobe {
    
Eigen::VectorXd sign(const Eigen::VectorXd& x);

void argsort_desc(
    const Eigen::VectorXd& w,
    std::vector<int>& ord
);

bool is_finite(double x);

double sum_nz(const Eigen::VectorXd& x);

double sum_abs(const Eigen::VectorXd& x);

}