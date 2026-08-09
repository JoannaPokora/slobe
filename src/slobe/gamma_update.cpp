#include "internal/gamma_update.hpp"

namespace slobe {

// gamma mean update
void gamma_mean_update(
    const Eigen::VectorXd& abs_beta_ord,
    const Eigen::VectorXd& lambda,
    int p,
    Eigen::VectorXd& gamma_h
) {
    double bl_sum = lambda(0) * abs_beta_ord(0);
    double bl_mean;
    int equals = 1;

    for (int i = 1; i < p; ++i) {
        if (abs_beta_ord(i) == abs_beta_ord(i - 1)) {
            bl_sum += lambda(i) * abs_beta_ord(i);
            equals += 1;
        } else {
            bl_mean = bl_sum / equals;
            for (int j = i - equals; j < i; ++j) gamma_h(j) = bl_mean;
            bl_sum = lambda(i) * abs_beta_ord(i);
            equals = 1;
        }
    }

    bl_mean = bl_sum / equals;
    for (int j = p - equals; j < p; ++j) gamma_h(j) = bl_mean;
}

}