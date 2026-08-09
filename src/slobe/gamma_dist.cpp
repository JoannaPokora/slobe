#include "internal/gamma_dist.hpp"

#include <stats.hpp>
#include <cmath>

namespace slobe {

namespace {

// Expectation of truncated gamma distribution
double pgamma(double x, double shape, double scale, int lower_tail, int log_p) {
    double val = stats::pgamma(x, shape, scale);
    
    if (!lower_tail) val = 1.0 - val;
    if (log_p) return std::log(val);
    
    return val;
}

}

double EX_trunc_gamma(double a, double b) {
    // E[Gamma(a,b) | X <= 1] = a/b * P(X<=1 | shape=a+1) / P(X<=1 | shape=a)
    double log_p1 = pgamma(1.0, a + 1.0, 1.0 / b, 1, 1);
    double log_p2 = pgamma(1.0, a,     1.0 / b, 1, 1);
    double c = std::exp(log_p1 - log_p2);
    c /= b;
    c *= a;
    return c;
}

}