#include "internal/helpers.hpp"

#include <cmath>
#include <numeric>

namespace slobe {

// argsort decreasing by values (for Eigen::VectorXd)
void argsort_desc(const Eigen::VectorXd& w, std::vector<int>& ord) {
    ord.resize(w.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&w](int a, int b) {
        return w(a) > w(b);
    });
}

// sign
Eigen::VectorXd sign(const Eigen::VectorXd& x) {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
    result = (x.array() > 0.0).select(1.0, result);
    result = (x.array() < 0.0).select(-1.0, result);
    return result;
}

// finite
bool is_finite(double x) {
  return std::isfinite(x);
}

// sum(x != 0)
double sum_nz(const Eigen::VectorXd& x) {
    double s = 0.0;
    for (int i = 0; i < x.size(); ++i)
        s += (x(i) != 0.0) ? 1.0 : 0.0;
    return s;
}

// sum(abs(x))
double sum_abs(const Eigen::VectorXd& x) {
    return x.array().abs().sum();
}

}
