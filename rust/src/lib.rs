pub mod logger;
pub mod phase1;
pub mod phase2;
pub mod pivot;
pub mod preprocess;
pub mod types;

use logger::SimplexLogger;
use phase1::phase1 as run_phase1;
use phase2::phase2 as run_phase2;

pub use logger::LogLevel;
pub use pivot::{Bland, LargestCoefficient, PivotRule};
pub use types::{BasisState, LinearProgram, SimplexResult, SolveStatus};

/// Solve a linear program using the two-phase Revised Simplex Method.
pub fn solve(
    lp: LinearProgram,
    rule: &dyn PivotRule,
    verbosity: LogLevel,
    logfile: Option<&str>,
) -> SimplexResult {
    let mut logger = match logfile {
        Some(path) => SimplexLogger::to_file(path, LogLevel::Debug)
            .unwrap_or_else(|_| SimplexLogger::silent()),
        None => {
            if verbosity == LogLevel::Silent {
                SimplexLogger::silent()
            } else {
                SimplexLogger::new(verbosity)
            }
        }
    };

    logger.log_problem(&lp);
    let n = lp.ncols();
    // Phase 2 uses the original LP (same feasible set; preprocess only affects phase 1)
    let lp_for_p2 = lp.clone();

    let (basis, p1_status) = run_phase1(lp, rule, &mut logger, 1000, 1e-8);

    if p1_status == SolveStatus::Infeasible {
        return SimplexResult {
            x: nalgebra::DVector::zeros(n),
            z: 0.0,
            status: SolveStatus::Infeasible,
            iterations: 0,
        };
    }

    let result = run_phase2(&lp_for_p2, &basis, rule, &mut logger, 1000, 1e-10);
    logger.log_result(&result);
    result
}

/// Run Phase 2 only, starting from a known basic feasible solution.
pub fn solve_phase2(
    lp: &LinearProgram,
    init_basis: &BasisState,
    rule: &dyn PivotRule,
    verbosity: LogLevel,
) -> SimplexResult {
    let mut logger = if verbosity == LogLevel::Silent {
        SimplexLogger::silent()
    } else {
        SimplexLogger::new(verbosity)
    };
    let result = run_phase2(lp, init_basis, rule, &mut logger, 1000, 1e-10);
    logger.log_result(&result);
    result
}
