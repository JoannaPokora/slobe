#pragma once

#include <Eigen/Core>

namespace slobe {

    int connection_test();

    struct SLOBEResult {
        Eigen::VectorXd coefficients;
        double sigma;
        double theta;
        double c;
        Eigen::VectorXd w;
        bool converged;
        Eigen::MatrixXd X;
        Eigen::MatrixXd Sigma;
        Eigen::VectorXd mu;
        Eigen::VectorXd lambda;
    };

    SLOBEResult slobe_admm(
        const Eigen::VectorXd& start,
        const Eigen::MatrixXd& Xmis_r,
        const Eigen::MatrixXd& Xinit,
        const Eigen::VectorXd& Y_r,
        double a_prior,
        double b_prior,
        const Eigen::MatrixXd& Covmat_r,
        double sigma = 1.0,
        double FDR = 0.05,
        double tol = 1e-4,
        bool known_sigma = false,
        int max_iter = 100,
        bool verbose = false,
        bool BH = true,
        bool known_cov = false
    );

}