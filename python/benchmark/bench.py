from __future__ import annotations

import time

import numpy as np

from simplex import BasisState, LinearProgram, solve, solve_phase2

N_RUNS = 10


def random_lp(m: int, n: int, seed: int = 42) -> LinearProgram:
    rng   = np.random.default_rng(seed)
    A_raw = rng.random((m, n - m))
    A     = np.hstack([A_raw, np.eye(m)])
    b     = np.abs(rng.random(m)) + 0.1
    c     = np.concatenate([rng.random(n - m), np.zeros(m)])
    return LinearProgram(A, b, c)


def benchmark(name: str, fn, n_runs: int = N_RUNS) -> None:
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    median_ms = np.median(times) * 1e3
    print(f"{name}: median={median_ms:.3f} ms")


def main() -> None:
    lp_small    = LinearProgram(
        np.array([[2, 1, 1, 0], [1, 2, 0, 1]], dtype=float),
        np.array([4, 4], dtype=float),
        np.array([4, 3, 0, 0], dtype=float),
    )
    basis_small = BasisState(basic=[2, 3], nonbasic=[0, 1])

    benchmark("small/phase2_only", lambda: solve_phase2(lp_small, basis_small))
    benchmark("small/two_phase",   lambda: solve(lp_small))

    for m, n in [(10, 20), (50, 100), (100, 200)]:
        lp  = random_lp(m, n)
        tag = f"m={m}_n={n}"
        benchmark(f"random/{tag}", lambda lp=lp: solve(lp))


if __name__ == "__main__":
    main()
