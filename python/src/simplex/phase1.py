from __future__ import annotations

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from .logger import SimplexLogger
from .pivot import PivotRule, entering_index, leaving_index
from .preprocess import make_b_nonnegative
from .types import BasisState, LinearProgram, SolveStatus


def phase1(lp: LinearProgram, rule: PivotRule, logger: SimplexLogger,
           maxit: int = 1000, tol: float = 1e-8) -> tuple[BasisState, SolveStatus]:

    lp   = make_b_nonnegative(lp)
    m, n = lp.A.shape

    Aw = np.hstack([lp.A, np.eye(m)])
    cw = np.concatenate([np.zeros(n), -np.ones(m)])
    art_start = n  # 0-based: artificials are columns n..n+m-1

    basic    = list(range(n, n + m))
    nonbasic = list(range(n))

    logger.log_phase(1)

    for it in range(1, maxit + 1):
        B     = Aw[:, basic]
        N     = Aw[:, nonbasic]
        BF    = lu_factor(B)
        xb    = lu_solve(BF, lp.b)
        BinvN = lu_solve(BF, N)

        cb = cw[basic]
        cr = cw[nonbasic] - BinvN.T @ cb
        z  = float(cb @ xb)

        logger.log_iteration(it, BasisState(list(basic), list(nonbasic)), xb, z)

        j = entering_index(rule, cr, tol)

        if j is None:
            if z < -tol:
                return BasisState([], []), SolveStatus.INFEASIBLE

            fix_artificials_in_basis(Aw, basic, nonbasic, art_start)
            orig_basic    = [i for i in basic    if i < n]
            orig_nonbasic = [i for i in nonbasic if i < n]
            return BasisState(orig_basic, orig_nonbasic), SolveStatus.OPTIMAL

        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN, tol)

        if i is None:
            raise RuntimeError("Phase 1 is unexpectedly unbounded")

        basic[i], nonbasic[j] = nonbasic[j], basic[i]

    raise RuntimeError(f"Phase 1 did not converge in {maxit} iterations")


def fix_artificials_in_basis(Aw: np.ndarray, basic: list[int],
                              nonbasic: list[int], art_start: int) -> None:
    changed = True
    while changed:
        changed = False
        for pos, col in enumerate(basic):
            if col < art_start:
                continue

            BF     = lu_factor(Aw[:, basic])
            BinvAw = lu_solve(BF, Aw)
            row    = BinvAw[pos, :]

            swap_j = None
            for j, ncol in enumerate(nonbasic):
                if ncol < art_start and abs(row[ncol]) > 1e-10:
                    swap_j = j
                    break

            if swap_j is None:
                Aw[pos, :] = 0.0
            else:
                basic[pos], nonbasic[swap_j] = nonbasic[swap_j], basic[pos]
                changed = True
                break
