#include "simplex/simplex.hpp"
#include <iostream>

SimplexResult solve(const LinearProgram& lp,
                    const PivotRule&     rule,
                    LogLevel             verbosity)
{
    SimplexLogger logger(verbosity, std::cout);
    logger.log_problem(lp);

    auto [basis, p1_status] = phase1(lp, rule, logger);

    if (p1_status == SolveStatus::INFEASIBLE) {
        return SimplexResult{
            Eigen::VectorXd::Zero(lp.A.cols()),
            0.0,
            SolveStatus::INFEASIBLE,
            0
        };
    }

    SimplexResult result = phase2(lp, basis, rule, logger);
    logger.log_result(result);
    return result;
}

SimplexResult solve_phase2(const LinearProgram& lp,
                            const BasisState&    init_basis,
                            const PivotRule&     rule,
                            LogLevel             verbosity)
{
    SimplexLogger logger(verbosity, std::cout);
    SimplexResult result = phase2(lp, init_basis, rule, logger);
    logger.log_result(result);
    return result;
}
