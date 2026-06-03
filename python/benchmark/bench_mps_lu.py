# Usage: python bench_mps_lu.py <data.json> [N=10]
#
# Two-phase Revised Simplex with hand-rolled LU factorization.
# No scipy.linalg calls for the basis solve — numpy arrays only for storage.

from __future__ import annotations

import json
import sys
import time
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# LU factorization with partial pivoting (Doolittle, in-place)
# ---------------------------------------------------------------------------

def lu_factor(B: np.ndarray) -> tuple[np.ndarray, list[int]]:
    m = B.shape[0]
    M = B.copy().astype(float, copy=False)
    perm = list(range(m))
    for k in range(m):
        r = int(np.argmax(np.abs(M[k:, k]))) + k
        if r != k:
            M[[k, r]] = M[[r, k]]
            perm[k], perm[r] = perm[r], perm[k]
        if abs(M[k, k]) < 1e-14:
            continue
        M[k + 1:, k] /= M[k, k]
        M[k + 1:, k + 1:] -= np.outer(M[k + 1:, k], M[k, k + 1:])
    return M, perm


def lu_solve_vec(M: np.ndarray, perm: list[int], b: np.ndarray) -> np.ndarray:
    y = b[perm].astype(float)
    m = len(y)
    for i in range(1, m):
        y[i] -= M[i, :i] @ y[:i]
    x = y.copy()
    for i in range(m - 1, -1, -1):
        if M[i, i] != 0:
            x[i] = (x[i] - M[i, i + 1:] @ x[i + 1:]) / M[i, i]
    return x


def lu_solve_mat(M: np.ndarray, perm: list[int], N: np.ndarray) -> np.ndarray:
    result = np.empty_like(N, dtype=float)
    for j in range(N.shape[1]):
        result[:, j] = lu_solve_vec(M, perm, N[:, j])
    return result


# ---------------------------------------------------------------------------
# Pivot helpers — same logic as simplex package
# ---------------------------------------------------------------------------

def entering_index(cr: np.ndarray, tol: float = 1e-10) -> Optional[int]:
    j = int(np.argmax(cr))
    return j if cr[j] > tol else None


def leaving_index(xb: np.ndarray, d: np.ndarray, BinvN: np.ndarray,
                  tol: float = 1e-10) -> Optional[int]:
    candidates = [i for i in range(len(d)) if d[i] > tol]
    if not candidates:
        return None
    best = candidates[0]
    for i in candidates[1:]:
        ri = xb[i] / d[i]
        rb = xb[best] / d[best]
        if ri < rb - tol:
            best = i
        elif abs(ri - rb) <= tol:
            if tuple(BinvN[i, :] / d[i]) < tuple(BinvN[best, :] / d[best]):
                best = i
    return best


# ---------------------------------------------------------------------------
# Phase 1 — artificial variable method
# ---------------------------------------------------------------------------

def phase1(A: np.ndarray, b: np.ndarray, n: int,
           maxit: int = 10_000, tol: float = 1e-8
           ) -> tuple[Optional[list[int]], Optional[list[int]], str]:
    mask = b < 0
    if np.any(mask):
        A = A.copy(); b = b.copy()
        A[mask] *= -1; b[mask] *= -1

    m = A.shape[0]
    Aw = np.hstack([A, np.eye(m)])
    cw = np.concatenate([np.zeros(n), -np.ones(m)])
    art = n

    basic = list(range(n, n + m))
    nonbasic = list(range(n))

    for _ in range(maxit):
        B = Aw[:, basic]; N = Aw[:, nonbasic]
        LU, perm = lu_factor(B)
        xb = lu_solve_vec(LU, perm, b)
        BinvN = lu_solve_mat(LU, perm, N)
        cb = cw[basic]
        cr = cw[nonbasic] - BinvN.T @ cb
        z = float(cb @ xb)
        j = entering_index(cr, tol)
        if j is None:
            if z < -tol:
                return None, None, "infeasible"
            _fix_arts(Aw, basic, nonbasic, art)
            return ([i for i in basic if i < n],
                    [i for i in nonbasic if i < n], "optimal")
        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN, tol)
        if i is None:
            raise RuntimeError("Phase 1 unexpectedly unbounded")
        basic[i], nonbasic[j] = nonbasic[j], basic[i]

    raise RuntimeError("Phase 1 did not converge")


def _fix_arts(Aw, basic, nonbasic, art):
    changed = True
    while changed:
        changed = False
        for pos, col in enumerate(basic):
            if col < art:
                continue
            LU, perm = lu_factor(Aw[:, basic])
            BinvAw = lu_solve_mat(LU, perm, Aw)
            row = BinvAw[pos, :]
            swap_j = next(
                (j for j, nc in enumerate(nonbasic) if nc < art and abs(row[nc]) > 1e-10),
                None)
            if swap_j is None:
                Aw[pos, :] = 0.0
            else:
                basic[pos], nonbasic[swap_j] = nonbasic[swap_j], basic[pos]
                changed = True
                break


# ---------------------------------------------------------------------------
# Phase 2 — optimise
# ---------------------------------------------------------------------------

def phase2(A: np.ndarray, b: np.ndarray, c: np.ndarray,
           basic: list[int], nonbasic: list[int],
           maxit: int = 10_000, tol: float = 1e-10
           ) -> tuple[Optional[float], int, str]:
    basic = list(basic); nonbasic = list(nonbasic)
    for it in range(1, maxit + 1):
        B = A[:, basic]; N = A[:, nonbasic]
        LU, perm = lu_factor(B)
        xb = lu_solve_vec(LU, perm, b)
        BinvN = lu_solve_mat(LU, perm, N)
        cb = c[basic]; cr = c[nonbasic] - BinvN.T @ cb
        z = float(cb @ xb)
        j = entering_index(cr, tol)
        if j is None:
            return z, it, "optimal"
        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN, tol)
        if i is None:
            return float("inf"), it, "unbounded"
        basic[i], nonbasic[j] = nonbasic[j], basic[i]
    raise RuntimeError("Phase 2 did not converge")


def solve(A: np.ndarray, b: np.ndarray, c: np.ndarray) -> tuple[Optional[float], int]:
    n = A.shape[1]
    basic, nonbasic, st = phase1(A, b, n)
    if st != "optimal":
        return None, 0
    all_cols = set(range(n))
    nonbasic = nonbasic + sorted(all_cols - set(basic) - set(nonbasic))
    z, iters, st2 = phase2(A, b, c, basic, nonbasic)
    return (z if st2 == "optimal" else None), iters


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: bench_mps_lu.py <data.json> [N=10]", file=sys.stderr)
        sys.exit(1)
    data_path = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) >= 3 else 10

    with open(data_path) as f:
        data = json.load(f)
    A = np.array(data["A"], dtype=float)
    b = np.array(data["b"], dtype=float)
    c = np.array(data["c"], dtype=float)

    times_ms: list[float] = []
    z_result = None; iters_result = 0
    for _ in range(N):
        t0 = time.perf_counter()
        z_result, iters_result = solve(A, b, c)
        times_ms.append((time.perf_counter() - t0) * 1e3)

    print(json.dumps({"z": z_result, "iterations": iters_result, "times_ms": times_ms}))


if __name__ == "__main__":
    main()
