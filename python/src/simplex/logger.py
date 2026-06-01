from __future__ import annotations

import sys
from enum import Enum
from typing import IO

from .types import BasisState, LinearProgram, SimplexResult


class LogLevel(Enum):
    SILENT = "SILENT"
    INFO   = "INFO"
    DEBUG  = "DEBUG"


SILENT = LogLevel.SILENT
INFO   = LogLevel.INFO
DEBUG  = LogLevel.DEBUG


class SimplexLogger:
    def __init__(self, level: LogLevel = LogLevel.SILENT, stream: IO | None = None) -> None:
        self.level  = level
        self._stream = stream
        self._owns_stream = False

    @classmethod
    def from_path(cls, path: str, level: LogLevel = DEBUG) -> "SimplexLogger":
        logger = cls(level, open(path, "w"))
        logger._owns_stream = True
        return logger

    def _write(self, s: str) -> None:
        stream = self._stream if self._stream is not None else sys.stdout
        print(s, file=stream)

    def log_problem(self, lp: LinearProgram) -> None:
        if self.level == LogLevel.SILENT:
            return
        self._write("=======================")
        self._write(" Revised Simplex Solver")
        self._write("=======================")
        self._write(f"A =\n{lp.A}")
        self._write(f"b = {lp.b}")
        self._write(f"c = {lp.c}")
        self._write("")

    def log_phase(self, phase: int) -> None:
        if self.level == LogLevel.SILENT:
            return
        header = ("Phase 1 — Finding initial BFS" if phase == 1
                  else "Phase 2 — Optimizing")
        self._write(header)
        self._write("─" * len(header))

    def log_iteration(self, it: int, basis: BasisState,
                      xb: "np.ndarray", z: float) -> None:
        if self.level != LogLevel.DEBUG:
            return
        self._write(f"  iter {it}:")
        self._write(f"    xb       = {xb}")
        self._write(f"    basic    = {basis.basic}")
        self._write(f"    nonbasic = {basis.nonbasic}")
        self._write(f"    z        = {z}")
        self._write("")

    def log_result(self, result: SimplexResult) -> None:
        if self.level == LogLevel.SILENT:
            return
        self._write("")
        self._write(f"Result: {result.status.name}")
        self._write(f"  x          = {result.x}")
        self._write(f"  z          = {result.z}")
        self._write(f"  iterations = {result.iterations}")

    def close(self) -> None:
        if self._owns_stream and self._stream is not None:
            self._stream.close()
