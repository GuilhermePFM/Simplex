# Usage: python bench_mps.py <data.json> [N=10]

from __future__ import annotations

import json
import sys
import time

import numpy as np

from simplex import LinearProgram, SolveStatus, solve


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: bench_mps.py <data.json> [N=10]", file=sys.stderr)
        sys.exit(1)

    data_path = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) >= 3 else 10

    with open(data_path) as f:
        data = json.load(f)

    A = np.array(data["A"], dtype=float)
    b = np.array(data["b"], dtype=float)
    c = np.array(data["c"], dtype=float)

    lp = LinearProgram(A, b, c)

    times_ms: list[float] = []
    result = None
    for _ in range(N):
        t0 = time.perf_counter()
        result = solve(lp)
        times_ms.append((time.perf_counter() - t0) * 1e3)

    z = result.z if result is not None and result.status == SolveStatus.OPTIMAL else None
    iters = result.iterations if result is not None else 0

    print(json.dumps({"z": z, "iterations": iters, "times_ms": times_ms}))


if __name__ == "__main__":
    main()
