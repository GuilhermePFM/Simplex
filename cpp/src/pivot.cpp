#include "simplex/pivot.hpp"
#include <cmath>
#include <vector>

std::optional<int> LargestCoefficient::entering_index(const Eigen::VectorXd& cr,
                                                        double tol) const
{
    int    best_j = -1;
    double best_v = tol;
    for (int j = 0; j < static_cast<int>(cr.size()); ++j) {
        if (cr[j] > best_v) {
            best_v = cr[j];
            best_j = j;
        }
    }
    if (best_j < 0) return std::nullopt;
    return best_j;
}

std::optional<int> Bland::entering_index(const Eigen::VectorXd& cr,
                                          double tol) const
{
    for (int j = 0; j < static_cast<int>(cr.size()); ++j) {
        if (cr[j] > tol) return j;
    }
    return std::nullopt;
}

std::optional<int> leaving_index(const Eigen::VectorXd& xb,
                                  const Eigen::VectorXd& d,
                                  const Eigen::MatrixXd& BinvN,
                                  double tol)
{
    int m = static_cast<int>(d.size());

    std::vector<int> candidates;
    candidates.reserve(m);
    for (int i = 0; i < m; ++i) {
        if (d[i] > tol) candidates.push_back(i);
    }
    if (candidates.empty()) return std::nullopt;

    int best = candidates[0];
    for (std::size_t k = 1; k < candidates.size(); ++k) {
        int i = candidates[k];
        double ratio_i    = xb[i]    / d[i];
        double ratio_best = xb[best] / d[best];

        if (ratio_i < ratio_best - tol) {
            best = i;
        } else if (std::abs(ratio_i - ratio_best) <= tol) {
            // Lexicographic tie-break: compare scaled rows of BinvN exactly,
            // mirroring Julia's Tuple(row_i) < Tuple(row_best).
            Eigen::VectorXd row_i    = BinvN.row(i).transpose()   / d[i];
            Eigen::VectorXd row_best = BinvN.row(best).transpose() / d[best];
            for (int c = 0; c < static_cast<int>(row_i.size()); ++c) {
                if (row_i[c] < row_best[c]) { best = i; break; }
                if (row_i[c] > row_best[c]) {           break; }
            }
        }
    }
    return best;
}
