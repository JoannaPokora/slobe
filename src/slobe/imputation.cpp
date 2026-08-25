#include "internal/imputation.hpp"

#include <Eigen/LU>

namespace slobe {

namespace {

void impute_row_advance(
    const Eigen::VectorXd& beta,
    Eigen::MatrixXd& X,
    const Eigen::VectorXd& Y,
    const Eigen::MatrixXd& S,
    double sigma_sq,
    const std::vector<std::vector<bool>>& XisFin,
    int n, int p,
    int row,
    const std::vector<int>& nanCols,
    const Eigen::VectorXd& m,
    const Eigen::VectorXd& tau_sq
)
{
    const int l = static_cast<int>(nanCols.size());
    if (l == 0) return;

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(l, l);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(l);

    double r = Y(row);
    int u_ind = 0;

    for (int i = 0; i < p; ++i) {
        if (!XisFin[row][i]) {
            for (int j = 0; j < p; ++j) {
                if (XisFin[row][j]) {
                    u(u_ind) += X(row, j) * S(j, i);
                }
            }
            ++u_ind;
        } else {
            r -= beta(i) * X(row, i);
        }
    }

    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < l; ++j) {
            if (i == j) {
                A(i, j) = 1.0;
            } else {
                int s = nanCols[i];
                int t = nanCols[j];
                A(i, j) = (beta(s) * beta(t) / sigma_sq + S(s, t)) / tau_sq(s);
            }
        }
    }

    Eigen::VectorXd b = Eigen::VectorXd::Zero(l);
    for (int i = 0; i < l; ++i) {
        int t = nanCols[i];
        b(i) = ((r * beta(t)) / sigma_sq + m(t) - u(i)) / tau_sq(t);
    }

    Eigen::VectorXd sol = A.fullPivLu().solve(b);

    for (int i = 0; i < l; ++i) {
        int t = nanCols[i];
        X(row, t) = sol(i);
    }
}

}

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
) {
    Eigen::VectorXd tau_sq = (beta.array().square() / sigma_sq).matrix();
    for (int i = 0; i < p; ++i) tau_sq(i) = tau_sq(i) + S(i, i);
    
    Eigen::VectorXd m = Eigen::VectorXd::Zero(p);
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
            m(i) += mu(j) * S(i, j);
        }
    }
    
    for (int idx = 0; idx < static_cast<int>(anyNanXrows.size()); ++idx) {
        int row = anyNanXrows[idx];
        impute_row_advance(beta, X, Y, S, sigma_sq, XisFin, n, p, row,
                           nanIndInRows[idx], m, tau_sq);
    }
}

}