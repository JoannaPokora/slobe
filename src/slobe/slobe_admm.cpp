#include <slobe/slobe_admm.hpp>

#include "internal/center_and_scale.hpp"
#include "internal/cov_shrink_est.hpp"
#include "internal/gamma_dist.hpp"
#include "internal/gamma_update.hpp"
#include "internal/helpers.hpp"
#include "internal/imputation.hpp"
#include "internal/slope_fit.hpp"

#include <iostream>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <slope.h>
#include <stats.hpp>

namespace slobe {

slobe::SLOBEResult slobe_admm(
    const Eigen::VectorXd& start,
    const Eigen::MatrixXd& Xinit,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd* Xmis,
    double a_prior,
    double b_prior,
    const Eigen::MatrixXd& Covmat,
    double sigma,
    double FDR,
    double tol,
    bool known_sigma,
    int max_iter,
    bool verbose,
    bool BH,
    bool known_cov
) {
    const int p = static_cast<int>(start.size());
    const int n = static_cast<int>(y.size());

    Eigen::VectorXd beta = start;
    Eigen::VectorXd beta_new(p);
    Eigen::VectorXd beta_e = beta;
    Eigen::VectorXd w = Eigen::VectorXd::Ones(p);
    Eigen::VectorXd w_e = Eigen::VectorXd::Ones(p);

    Eigen::VectorXd wbeta(p);
    Eigen::VectorXd gamma(p);
    Eigen::VectorXd gamma_h(p);
    Eigen::VectorXd lambda_sigma(p);

    std::vector<int> order(p);
    Eigen::MatrixXd X_div_w(n, p);

    double error = 0.0;
    double swlambda = 0.0;
    double RSS = 0.0;
    double sigma_sq = 1.0;
    
    const bool has_missing = (Xmis != nullptr);
    
    std::vector<std::vector<bool>> XisFin(n, std::vector<bool>(p));
    std::vector<int> anyNanXrows;
    std::vector<std::vector<int>> nanIndicesInRow;

    if (has_missing) {
      for (int i = 0; i < n; ++i) {
        bool anyNan = false;
        std::vector<int> nanInd;
        for (int j = 0; j < p; ++j) {
          XisFin[i][j] = is_finite((*Xmis)(i, j));
          if (!XisFin[i][j]) {
            nanInd.push_back(j);
            anyNan = true;
          }
        }
        if (anyNan) {
          anyNanXrows.push_back(i);
          nanIndicesInRow.push_back(nanInd);
        }
      } 
    }

    // X starts from Xinit, then center/scale
    Eigen::MatrixXd X = Xinit;
    center_and_scale(X);

    // Sigma + precision S
    Eigen::MatrixXd Sigma(p, p);
    if (!known_cov) {
        linshrink_cov(X, Sigma, n, p);
    } else {
        Sigma = Covmat;
    }

    Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("Cholesky failed for Sigma.");
    }
    Eigen::MatrixXd S = llt.solve(Eigen::MatrixXd::Identity(p, p));

    // mu (column means)
    Eigen::VectorXd mu = X.colwise().mean().transpose();

    // lambda
    Eigen::VectorXd lambda(p);
    create_lambda(lambda, p, FDR, BH);

    // init sigma, c, theta, gamma
    double sstart = sum_nz(start);
    double pstart = (sstart > 0.0) ? sstart : 1.0;
    
    std::cout << "n=" << n << ", pstart=" << pstart;

    if (!known_sigma) {
        sigma = std::sqrt((X * beta_e - y).squaredNorm() / (n - pstart));
    }

    lambda_sigma = lambda * sigma;

    argsort_desc(beta, order);

    double c = 0.0;
    if (sstart > 0.0) {
        double h = (sstart + 1.0) / sum_abs(beta) * (sigma / lambda(p - 1));
        c = (h < 1.0) ? h : 1.0;
    } else {
        c = 1.0;
    }

    double theta = (sstart + a_prior) / (a_prior + b_prior + p);

    // compute gamma_h on abs vals of ord coeffs
    argsort_desc(beta, order);
    Eigen::VectorXd beta_abs_ordered(p);
    for (int i = 0; i < p; ++i) {
        beta_abs_ordered(i) = std::abs(beta(order[i]));
    }
    gamma_mean_update(
        beta_abs_ordered,
        lambda,
        p,
        gamma_h
    );
    gamma_h = gamma_h * (c - 1.0) / sigma;
    gamma_h =
        (theta * c) /
        (theta * c +
        (1.0 - theta) * gamma_h.array().exp()).array();
    for (int i = 0; i < p; ++i) {
        gamma(order[i]) = gamma_h(i);
    }


    bool converged = false;
    int iter = 0;

    while (iter < max_iter) {
        if (verbose) {
            std::cout << "Iteration: " << iter << "/" << max_iter << std::endl;
        }

        div_X_by_w(X_div_w, X, w_e, n, p);

        Eigen::VectorXd lambda_sorted = lambda_sigma;
        std::sort(lambda_sorted.data(), lambda_sorted.data() + p, std::greater<double>());

        Eigen::VectorXd just_lambda_sorted = lambda;
        std::sort(just_lambda_sorted.data(), just_lambda_sorted.data() + p, std::greater<double>());

        Eigen::VectorXd beta_hat = slope_libslope_fit(X_div_w, y, lambda_sorted);

        wbeta = beta_hat.cwiseAbs();
        argsort_desc(wbeta, order);

        // unweight
        for (int i = 0; i < p; ++i) beta_hat(i) /= w_e(i);
        beta_new = beta_hat;

        if (!known_sigma) {
            RSS = (X * beta_hat - y).squaredNorm();
            Eigen::VectorXd wbeta_sorted = wbeta;
            std::sort(wbeta_sorted.data(), wbeta_sorted.data() + p, std::greater<double>());
            swlambda = (wbeta_sorted.array() * lambda.array()).sum();
            sigma = (swlambda + std::sqrt(swlambda * swlambda + 4.0 * (n + 2.0) * RSS))
                    / (2.0 * (n + 2.0));
        }

        lambda_sigma = lambda * sigma;
        sigma_sq     = sigma * sigma;
        
        if (has_missing) {
          if (!known_cov) {
            linshrink_cov(X, Sigma, n, p);
            Eigen::LLT<Eigen::MatrixXd> llt2(Sigma);
            if (llt2.info() != Eigen::Success) {
              throw std::runtime_error("Cholesky failed for Sigma (iter).");
            }
            S = llt2.solve(Eigen::MatrixXd::Identity(p, p));
          }
          
          mu = X.colwise().mean().transpose();
          
          impute_advance(beta_hat, X, y, S, sigma_sq, n, p, mu, XisFin,
                         anyNanXrows, nanIndicesInRow);
          center_and_scale(X);  
        }

        Eigen::VectorXd beta_new_abs_ordered(p);
        for (int i = 0; i < p; ++i) {
            beta_new_abs_ordered(i) =
                std::abs(beta_new(order[i]));
        }
        gamma_mean_update(
            beta_new_abs_ordered,
            lambda,
            p,
            gamma_h
        );
        gamma_h = gamma_h * (c - 1.0) / sigma;
        gamma_h = (theta * c) / (theta * c +
            (1.0 - theta) * gamma_h.array().exp()).array();
        for (int i = 0; i < p; ++i) {
            gamma(order[i]) = gamma_h(i);
        }

        double sum_gamma = gamma.sum();
        double b_sum = 0.0;
        for (int k = 0; k < p; ++k) {
            int idx = order[k];
            b_sum += gamma[idx] * std::fabs(beta_new[idx]) * lambda[k];
        }
        b_sum /= sigma;

        if (sum_gamma > 0.0) {
            if (b_sum > 0.0) {
                c = EX_trunc_gamma(sum_gamma, b_sum);
            } else {
                c = (sum_gamma + 1.0) / (sum_gamma + 2.0);
            }
        } else {
            c = 0.5;
        }

        theta = (sum_gamma + a_prior) / (p + a_prior + b_prior);

        error = (beta - beta_new).array().abs().sum();
        if (error < tol) {
            converged = true;
            break;
        }

        if (verbose) {
            std::cout << "Error = " << error << " sigma = " << sigma
                      << " theta = " << theta << " c = " << c << std::endl;
        }

        beta = beta_new;
        ++iter;
    }

    slobe::SLOBEResult result;
    result.coefficients = beta_new;
    result.sigma        = sigma;
    result.theta        = theta;
    result.c            = c;
    result.w            = w;
    result.converged    = converged;
    result.X            = X;
    result.Sigma        = Sigma;
    result.mu           = mu;
    result.lambda       = lambda;

    return result;
}

}
