#include "simplex/preprocess.hpp"

LinearProgram make_b_nonnegative(const LinearProgram& lp)
{
    int m = static_cast<int>(lp.b.size());
    bool any_neg = false;
    for (int i = 0; i < m; ++i) {
        if (lp.b[i] < 0.0) { any_neg = true; break; }
    }
    if (!any_neg) return lp;

    Eigen::MatrixXd A2 = lp.A;
    Eigen::VectorXd b2 = lp.b;
    for (int i = 0; i < m; ++i) {
        if (b2[i] < 0.0) {
            A2.row(i) *= -1.0;
            b2[i]     *= -1.0;
        }
    }
    return LinearProgram(A2, b2, lp.c);
}
