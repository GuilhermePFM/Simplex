from __future__ import annotations

import numpy as np

from .types import LinearProgram


def make_b_nonnegative(lp: LinearProgram) -> LinearProgram:
    """Return a copy of lp with all right-hand side values non-negative."""
    mask = lp.b < 0
    if not np.any(mask):
        return lp
    A2 = lp.A.copy()
    b2 = lp.b.copy()
    A2[mask, :] *= -1
    b2[mask]    *= -1
    return LinearProgram(A2, b2, lp.c)
