#include "simplex/phase1.hpp"
#include "simplex/preprocess.hpp"
#include <Eigen/LU>
#include <stdexcept>

// Forward declaration for the helper
static void fix_artificials_in_basis(Eigen::MatrixXd& Aw,
                                      std::vector<int>& basic,
                                      std::vector<int>& nonbasic,
                                      int art_start);

std::tuple<BasisState, SolveStatus>
phase1(const LinearProgram& lp_in, const PivotRule& rule,
       const SimplexLogger& logger, int maxit, double tol)
{
    LinearProgram lp = make_b_nonnegative(lp_in);
    int m = static_cast<int>(lp.A.rows());
    int n = static_cast<int>(lp.A.cols());

    // Aw = [A | I_m], artificials start at column index n (0-based)
    Eigen::MatrixXd Aw(m, n + m);
    Aw.leftCols(n)  = lp.A;
    Aw.rightCols(m) = Eigen::MatrixXd::Identity(m, m);

    Eigen::VectorXd cw(n + m);
    cw.head(n).setZero();
    cw.tail(m).setConstant(-1.0);

    int art_start = n; // 0-based

    std::vector<int> basic(m), nonbasic(n);
    for (int i = 0; i < m; ++i) basic[i]    = art_start + i;
    for (int i = 0; i < n; ++i) nonbasic[i] = i;

    logger.log_phase(1);

    for (int it = 1; it <= maxit; ++it) {
        // Build B and N from current basis
        Eigen::MatrixXd B(m, m), N(m, n);
        for (int j = 0; j < m; ++j) B.col(j) = Aw.col(basic[j]);
        for (int j = 0; j < n; ++j) N.col(j) = Aw.col(nonbasic[j]);

        Eigen::PartialPivLU<Eigen::MatrixXd> BF(B);
        Eigen::VectorXd xb    = BF.solve(lp.b);
        Eigen::MatrixXd BinvN = BF.solve(N);

        Eigen::VectorXd cb(m);
        for (int i = 0; i < m; ++i) cb[i] = cw[basic[i]];

        Eigen::VectorXd cr_nb(n);
        for (int j = 0; j < n; ++j) cr_nb[j] = cw[nonbasic[j]];
        cr_nb -= BinvN.transpose() * cb;

        double z = cb.dot(xb);

        logger.log_iteration(it, BasisState{basic, nonbasic}, xb, z);

        auto opt_j = rule.entering_index(cr_nb, tol);

        if (!opt_j.has_value()) {
            if (z < -tol) {
                return { BasisState{{}, {}}, SolveStatus::INFEASIBLE };
            }
            fix_artificials_in_basis(Aw, basic, nonbasic, art_start);

            std::vector<int> orig_basic, orig_nonbasic;
            for (int idx : basic)    if (idx < n) orig_basic.push_back(idx);
            for (int idx : nonbasic) if (idx < n) orig_nonbasic.push_back(idx);
            return { BasisState{orig_basic, orig_nonbasic}, SolveStatus::OPTIMAL };
        }

        int j = *opt_j;
        Eigen::VectorXd d = BinvN.col(j);
        auto opt_i = leaving_index(xb, d, BinvN, tol);

        if (!opt_i.has_value())
            throw std::runtime_error("Phase 1 is unexpectedly unbounded");

        int i = *opt_i;
        std::swap(basic[i], nonbasic[j]);
    }

    throw std::runtime_error("Phase 1 did not converge in " +
                              std::to_string(maxit) + " iterations");
}

static void fix_artificials_in_basis(Eigen::MatrixXd& Aw,
                                      std::vector<int>& basic,
                                      std::vector<int>& nonbasic,
                                      int art_start)
{
    int m = static_cast<int>(basic.size());
    bool changed = true;
    while (changed) {
        changed = false;
        for (int pos = 0; pos < m; ++pos) {
            if (basic[pos] < art_start) continue;

            Eigen::MatrixXd B(Aw.rows(), m);
            for (int j = 0; j < m; ++j) B.col(j) = Aw.col(basic[j]);
            Eigen::PartialPivLU<Eigen::MatrixXd> BF(B);
            Eigen::MatrixXd BinvAw = BF.solve(Aw);
            Eigen::VectorXd row = BinvAw.row(pos).transpose();

            int swap_j = -1;
            for (int j = 0; j < static_cast<int>(nonbasic.size()); ++j) {
                int ncol = nonbasic[j];
                if (ncol < art_start && std::abs(row[ncol]) > 1e-10) {
                    swap_j = j;
                    break;
                }
            }

            if (swap_j < 0) {
                Aw.row(pos).setZero();
            } else {
                std::swap(basic[pos], nonbasic[swap_j]);
                changed = true;
                break;
            }
        }
    }
}
