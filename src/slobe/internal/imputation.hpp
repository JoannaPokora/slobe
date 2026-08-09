#pragma once

#include <Eigen/Core>
#include <vector>

namespace slobe {
    
void impute_advance(
    const Eigen::VectorXd& beta,
    Eigen::MatrixXd& X,
    const Eigen::VectorXd& Y,
    const Eigen::MatrixXd& S,
    double sigma_sq,
    int n, int p,
    const Eigen::VectorXd& mu,
    const std::vector<std::vector<bool>>& XisFin,
    const std::vector<int>& anyNanXrows,
    const std::vector<std::vector<int>>& nanIndInRows
);

}