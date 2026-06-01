from __future__ import annotations

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from .logger import SimplexLogger
from .pivot import PivotRule, entering_index, leaving_index
from .types import BasisState, LinearProgram, SimplexResult, SolveStatus


def phase2(lp: LinearProgram, init_basis: BasisState, rule: PivotRule,
           logger: SimplexLogger, maxit: int = 1000,
           tol: float = 1e-10) -> SimplexResult:

    A, b, c = lp.A, lp.b, lp.c
    n = A.shape[1]

    basic    = list(init_basis.basic)
    nonbasic = list(init_basis.nonbasic)

    x = np.zeros(n)
    x[basic] = np.linalg.solve(A[:, basic], b)

    logger.log_phase(2)

    for it in range(1, maxit + 1):
        B     = A[:, basic]
        N     = A[:, nonbasic]
        BF    = lu_factor(B)
        xb    = lu_solve(BF, b)
        BinvN = lu_solve(BF, N)

        cb = c[basic]
        cr = c[nonbasic] - BinvN.T @ cb
        z  = float(cb @ xb)

        logger.log_iteration(it, BasisState(list(basic), list(nonbasic)), xb, z)

        j = entering_index(rule, cr, tol)
        if j is None:
            x          = np.zeros(n)
            x[basic]   = xb
            return SimplexResult(x, z, SolveStatus.OPTIMAL, it)

        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN, tol)

        if i is None:
            ray = np.zeros(n)
            ray[nonbasic[j]] = 1.0
            ray[basic]       = -d
            return SimplexResult(ray, float("inf"), SolveStatus.UNBOUNDED, it)

        theta = xb[i] / d[i]

        x          = np.zeros(n)
        x[basic]   = xb - theta * d
        x[nonbasic[j]] = theta
        x[basic[i]]    = 0.0  # numerical cleanup for the leaving variable

        basic[i], nonbasic[j] = nonbasic[j], basic[i]

    raise RuntimeError(f"Phase 2 did not converge in {maxit} iterations — possible cycling")
