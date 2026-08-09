#include "internal/center_and_scale.hpp"

namespace slobe {

void center_and_scale(Eigen::MatrixXd& X) {
    const int p = X.cols();
    
    // center
    Eigen::RowVectorXd means = X.colwise().mean();
    X.rowwise() -= means;
    
    // normalise columns (l2 norm)
    for (int j = 0; j < p; ++j) {
        double norm = X.col(j).norm();
        if (norm > 0) X.col(j) /= norm;
    }
}

}