#include "simplex/phase2.hpp"
#include <Eigen/LU>
#include <limits>
#include <stdexcept>

SimplexResult
phase2(const LinearProgram& lp, const BasisState& init_basis,
       const PivotRule& rule, const SimplexLogger& logger,
       int maxit, double tol)
{
    const Eigen::MatrixXd& A = lp.A;
    const Eigen::VectorXd& b = lp.b;
    const Eigen::VectorXd& c = lp.c;
    int n = static_cast<int>(A.cols());

    std::vector<int> basic    = init_basis.basic;
    std::vector<int> nonbasic = init_basis.nonbasic;
    int nb = static_cast<int>(nonbasic.size());
    int m  = static_cast<int>(basic.size());

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    {
        Eigen::MatrixXd B(m, m);
        for (int j = 0; j < m; ++j) B.col(j) = A.col(basic[j]);
        Eigen::VectorXd xb_init = B.lu().solve(b);
        for (int j = 0; j < m; ++j) x[basic[j]] = xb_init[j];
    }

    logger.log_phase(2);

    for (int it = 1; it <= maxit; ++it) {
        Eigen::MatrixXd B(m, m), N(m, nb);
        for (int j = 0; j < m;  ++j) B.col(j) = A.col(basic[j]);
        for (int j = 0; j < nb; ++j) N.col(j) = A.col(nonbasic[j]);

        Eigen::PartialPivLU<Eigen::MatrixXd> BF(B);
        Eigen::VectorXd xb    = BF.solve(b);
        Eigen::MatrixXd BinvN = BF.solve(N);

        Eigen::VectorXd cb(m);
        for (int i = 0; i < m; ++i) cb[i] = c[basic[i]];

        Eigen::VectorXd cr(nb);
        for (int j = 0; j < nb; ++j) cr[j] = c[nonbasic[j]];
        cr -= BinvN.transpose() * cb;

        double z = cb.dot(xb);

        logger.log_iteration(it, BasisState{basic, nonbasic}, xb, z);

        auto opt_j = rule.entering_index(cr, tol);
        if (!opt_j.has_value()) {
            x.setZero();
            for (int i = 0; i < m; ++i) x[basic[i]] = xb[i];
            return SimplexResult{ x, z, SolveStatus::OPTIMAL, it };
        }

        int j = *opt_j;
        Eigen::VectorXd d = BinvN.col(j);
        auto opt_i = leaving_index(xb, d, BinvN, tol);

        if (!opt_i.has_value()) {
            Eigen::VectorXd ray = Eigen::VectorXd::Zero(n);
            ray[nonbasic[j]] = 1.0;
            for (int i = 0; i < m; ++i) ray[basic[i]] = -d[i];
            return SimplexResult{ ray, std::numeric_limits<double>::infinity(),
                                  SolveStatus::UNBOUNDED, it };
        }

        int i = *opt_i;
        double theta = xb[i] / d[i];

        x.setZero();
        for (int k = 0; k < m;  ++k) x[basic[k]]    = xb[k] - theta * d[k];
        x[nonbasic[j]] = theta;
        x[basic[i]]    = 0.0;

        std::swap(basic[i], nonbasic[j]);
    }

    throw std::runtime_error("Phase 2 did not converge in " +
                              std::to_string(maxit) +
                              " iterations \xe2\x80\x94 possible cycling");
}
