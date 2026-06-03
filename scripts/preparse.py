"""
Generate reproducible random LP instances and save to data/{name}.json.

Each LP has the identity block appended as slack columns, so the initial basis
{n-m, ..., n-1} is immediately feasible — Phase 1 drives out artificials quickly.

The reference optimal_z is computed by running Python's own solver.

Usage:
    python scripts/preparse.py
"""

from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

sys.path.insert(0, str(ROOT / "python" / "src"))
from simplex import solve, LinearProgram

PROBLEMS = {
    "small":   (10,  20),
    "medium":  (50,  100),
    "large":   (100, 200),
    "xlarge":  (200, 400),
}

SEED = 42


def random_lp(m: int, n: int, seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    A_raw = rng.random((m, n - m))
    A = np.hstack([A_raw, np.eye(m)])
    b = np.abs(rng.random(m)) + 0.1
    c = np.concatenate([rng.random(n - m), np.zeros(m)])
    return {"A": A, "b": b, "c": c}


def main() -> None:
    for name, (m, n) in PROBLEMS.items():
        print(f"[{name}]  m={m}, n={n} ...", end=" ", flush=True)

        lp_arrays = random_lp(m, n)
        lp = LinearProgram(lp_arrays["A"], lp_arrays["b"], lp_arrays["c"])

        result = solve(lp)
        if result.status.name != "OPTIMAL":
            print(f"WARNING: solver returned {result.status.name}")
            optimal_z = None
        else:
            optimal_z = float(result.z)

        data = {
            "name":      name,
            "m":         m,
            "n":         n,
            "A":         lp_arrays["A"].tolist(),
            "b":         lp_arrays["b"].tolist(),
            "c":         lp_arrays["c"].tolist(),
            "optimal_z": optimal_z,
        }

        out = DATA / f"{name}.json"
        with open(out, "w") as f:
            json.dump(data, f)

        print(f"optimal_z={optimal_z:.6f}  saved {out.name}")


if __name__ == "__main__":
    main()
