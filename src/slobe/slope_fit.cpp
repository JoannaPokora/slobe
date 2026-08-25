#include "internal/slope_fit.hpp"

#include <stats.hpp>
#include <slope.h>

#include <set>

namespace slobe {

// Creates a vector of weights (lambda) for SLOPE
void create_lambda(Eigen::VectorXd& lam, int p, double FDR, bool BH) {
    Eigen::VectorXd h(p);
    if (BH) {
        for (int i = 0; i < p; ++i) {
            h(i) = 1.0 - (FDR * (i + 1) / (2.0 * p));
        }
    } else {
        double const_val = FDR * 1.0 / (2.0 * p);
        h.setConstant(1.0 - const_val);
    }
    for (int i = 0; i < p; ++i) {
        lam(i) = stats::qnorm(h(i));
    }
}

// Scaling matrix X by weight vector w (divide each column by w[j])
void div_X_by_w(Eigen::MatrixXd& X_div_w,
                       const Eigen::MatrixXd& X,
                       const Eigen::VectorXd& w,
                       int n, int p) {
    X_div_w.resize(n, p);
    for (int j = 0; j < p; ++j) {
        X_div_w.col(j) = X.col(j) / w(j);
    }
}

Eigen::VectorXd slope_libslope_fit(
    Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& lambda,
    const std::string family
) {
    const int n = X.rows();
    const int p = X.cols();

    slope::Slope model;

    std::set<double> unique_values(y.data(), y.data() + y.size());
    int n_unique = static_cast<int>(unique_values.size());

    if (family == "binomial") {
      model.setLoss("logistic");
    } else {
      model.setLoss("quadratic");
    }
    
    model.setScaling("none");
    model.setCentering("none");
    model.setIntercept(false);

    Eigen::MatrixXd ymat(n, 1);
    ymat.col(0) = y;

    Eigen::ArrayXd lam = lambda.array();

    auto fit = model.fit(
        X,
        ymat,
        1.0 / static_cast<double>(n),
        lam
    );

    auto B = fit.getCoefs(false);

    Eigen::VectorXd beta =
        Eigen::VectorXd::Zero(p);

    for (int k = 0; k < B.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, k);
             it;
             ++it) {
            beta(it.row()) = it.value();
        }
    }

    return beta;
}

}
