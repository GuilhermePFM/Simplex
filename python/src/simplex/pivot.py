from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class PivotRule(ABC):
    @abstractmethod
    def entering_index(self, cr: np.ndarray, tol: float = 1e-10) -> int | None:
        ...


class LargestCoefficient(PivotRule):
    """Pick the non-basic variable with the largest (most positive) reduced cost."""

    def entering_index(self, cr: np.ndarray, tol: float = 1e-10) -> int | None:
        j = int(np.argmax(cr))
        return j if cr[j] > tol else None


class Bland(PivotRule):
    """Bland's rule: always pick the lowest-index positive reduced cost."""

    def entering_index(self, cr: np.ndarray, tol: float = 1e-10) -> int | None:
        for j, v in enumerate(cr):
            if v > tol:
                return j
        return None


def entering_index(rule: PivotRule, cr: np.ndarray, tol: float = 1e-10) -> int | None:
    return rule.entering_index(cr, tol)


def leaving_index(xb: np.ndarray, d: np.ndarray, BinvN: np.ndarray,
                  tol: float = 1e-10) -> int | None:
    """Lexicographic ratio test — selects the leaving variable to prevent cycling."""
    candidates = [i for i in range(len(d)) if d[i] > tol]
    if not candidates:
        return None

    best = candidates[0]
    for i in candidates[1:]:
        ratio_i    = xb[i]    / d[i]
        ratio_best = xb[best] / d[best]

        if ratio_i < ratio_best - tol:
            best = i
        elif abs(ratio_i - ratio_best) <= tol:
            row_i    = BinvN[i,    :] / d[i]
            row_best = BinvN[best, :] / d[best]
            if tuple(row_i) < tuple(row_best):
                best = i

    return best
