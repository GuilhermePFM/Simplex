from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np


class SolveStatus(IntEnum):
    OPTIMAL    =  1
    UNBOUNDED  = -1
    INFEASIBLE = -2


@dataclass
class LinearProgram:
    """Linear program in standard equality form: maximize c'x s.t. Ax=b, x>=0."""
    A: np.ndarray
    b: np.ndarray
    c: np.ndarray

    def __post_init__(self) -> None:
        self.A = np.array(self.A, dtype=float)
        self.b = np.array(self.b, dtype=float)
        self.c = np.array(self.c, dtype=float)
        m, n = self.A.shape
        if self.b.shape != (m,):
            raise ValueError(f"b has {self.b.shape[0]} rows; A has {m}")
        if self.c.shape != (n,):
            raise ValueError(f"c has {self.c.shape[0]} entries; A has {n} columns")


@dataclass
class BasisState:
    """Index partition of variables into basic and non-basic (0-based)."""
    basic:    list[int]
    nonbasic: list[int]


@dataclass
class SimplexResult:
    """Complete output of solve or solve_phase2."""
    x:          np.ndarray
    z:          float
    status:     SolveStatus
    iterations: int
