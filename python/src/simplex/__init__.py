from __future__ import annotations

import numpy as np

from .logger import (
    DEBUG,
    INFO,
    SILENT,
    LogLevel,
    SimplexLogger,
)
from .phase1 import phase1
from .phase2 import phase2
from .pivot import Bland, LargestCoefficient, PivotRule
from .types import BasisState, LinearProgram, SimplexResult, SolveStatus

__all__ = [
    "LinearProgram",
    "BasisState",
    "SimplexResult",
    "SolveStatus",
    "PivotRule",
    "LargestCoefficient",
    "Bland",
    "LogLevel",
    "SILENT",
    "INFO",
    "DEBUG",
    "SimplexLogger",
    "solve",
    "solve_phase2",
]


def solve(
    lp: LinearProgram,
    *,
    rule: PivotRule = None,
    verbosity: LogLevel = SILENT,
    logfile: str | None = None,
) -> SimplexResult:
    """Solve a linear program using the two-phase Revised Simplex Method."""
    if rule is None:
        rule = LargestCoefficient()

    if logfile is not None:
        logger = SimplexLogger.from_path(logfile)
    else:
        logger = SimplexLogger(verbosity)

    logger.log_problem(lp)

    basis, p1_status = phase1(lp, rule, logger)

    if p1_status == SolveStatus.INFEASIBLE:
        logger.close()
        return SimplexResult(np.zeros(lp.A.shape[1]), 0.0, SolveStatus.INFEASIBLE, 0)

    result = phase2(lp, basis, rule, logger)
    logger.log_result(result)
    logger.close()
    return result


def solve_phase2(
    lp: LinearProgram,
    init_basis: BasisState,
    *,
    rule: PivotRule = None,
    verbosity: LogLevel = SILENT,
) -> SimplexResult:
    """Run Phase 2 only, starting from a known basic feasible solution."""
    if rule is None:
        rule = LargestCoefficient()

    logger = SimplexLogger(verbosity)
    result = phase2(lp, init_basis, rule, logger)
    logger.log_result(result)
    return result
