"""
Benchmark orchestrator: runs each language's bench_mps entry point on all
Netlib LP problems, collects timing + peak RSS, writes results/results.csv.

Usage:
    python scripts/run_benchmark.py [--problems afiro blend sc205 ship04l]
                                    [--langs julia python typescript go rust cpp
                                             julia_lu python_lu go_lu rust_lu cpp_lu]
                                    [--runs 10]
"""

from __future__ import annotations
import argparse
import csv
import json
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

# ── Language registry ──────────────────────────────────────────────────────────
# Each entry: (display_name, command_builder)
# command_builder(data_json: Path, n_runs: int) -> list[str]

def _julia_cmd(data: Path, n: int) -> list[str]:
    julia = Path.home() / ".juliaup" / "bin" / "julia"
    script = ROOT / "julia" / "benchmark" / "bench_mps.jl"
    return [str(julia), "--project=" + str(ROOT / "julia"), str(script), str(data), str(n)]

def _python_cmd(data: Path, n: int) -> list[str]:
    script = ROOT / "python" / "benchmark" / "bench_mps.py"
    uv = shutil.which("uv") or "uv"
    return [uv, "run", "--project", str(ROOT / "python"), "python", str(script), str(data), str(n)]

def _typescript_cmd(data: Path, n: int) -> list[str]:
    script = ROOT / "typescript" / "benchmark" / "bench_mps.ts"
    tsx = str(ROOT / "typescript" / "node_modules" / ".bin" / "tsx")
    return [tsx, str(script), str(data), str(n)]

def _go_cmd(data: Path, n: int) -> list[str]:
    binary = ROOT / "go" / "cmd" / "bench_mps" / "bench_mps"
    if not binary.exists():
        raise FileNotFoundError(f"Go binary not found: {binary}\nRun: cd go && go build -o cmd/bench_mps/bench_mps ./cmd/bench_mps/")
    return [str(binary), str(data), str(n)]

def _rust_cmd(data: Path, n: int) -> list[str]:
    binary = ROOT / "rust" / "target" / "release" / "bench_mps"
    if not binary.exists():
        raise FileNotFoundError(f"Rust binary not found: {binary}\nRun: cd rust && cargo build --release --bin bench_mps")
    return [str(binary), str(data), str(n)]

def _cpp_cmd(data: Path, n: int) -> list[str]:
    binary = ROOT / "cpp" / "build" / "bench_mps"
    if not binary.exists():
        raise FileNotFoundError(f"C++ binary not found: {binary}\nRun: cd cpp && cmake --build build")
    return [str(binary), str(data), str(n)]

# ── LU (hand-rolled, no external LA lib) variants ──────────────────────────

def _julia_lu_cmd(data: Path, n: int) -> list[str]:
    julia = Path.home() / ".juliaup" / "bin" / "julia"
    script = ROOT / "julia" / "benchmark" / "bench_mps_lu.jl"
    return [str(julia), "--project=" + str(ROOT / "julia"), str(script), str(data), str(n)]

def _python_lu_cmd(data: Path, n: int) -> list[str]:
    script = ROOT / "python" / "benchmark" / "bench_mps_lu.py"
    uv = shutil.which("uv") or "uv"
    return [uv, "run", "--project", str(ROOT / "python"), "python", str(script), str(data), str(n)]

def _go_lu_cmd(data: Path, n: int) -> list[str]:
    binary = ROOT / "go" / "cmd" / "bench_mps_lu" / "bench_mps_lu"
    if not binary.exists():
        raise FileNotFoundError(
            f"Go LU binary not found: {binary}\n"
            "Run: cd go && go build -o cmd/bench_mps_lu/bench_mps_lu ./cmd/bench_mps_lu/")
    return [str(binary), str(data), str(n)]

def _rust_lu_cmd(data: Path, n: int) -> list[str]:
    binary = ROOT / "rust" / "target" / "release" / "bench_mps_lu"
    if not binary.exists():
        raise FileNotFoundError(
            f"Rust LU binary not found: {binary}\n"
            "Run: cd rust && cargo build --release --bin bench_mps_lu")
    return [str(binary), str(data), str(n)]

def _cpp_lu_cmd(data: Path, n: int) -> list[str]:
    binary = ROOT / "cpp" / "build" / "bench_mps_lu"
    if not binary.exists():
        raise FileNotFoundError(
            f"C++ LU binary not found: {binary}\nRun: cd cpp && cmake --build build")
    return [str(binary), str(data), str(n)]


LANGUAGES: dict[str, tuple[str, object]] = {
    # External-library LA backends
    "julia":      ("Julia",         _julia_cmd),
    "python":     ("Python",        _python_cmd),
    "typescript": ("TypeScript",    _typescript_cmd),
    "go":         ("Go",            _go_cmd),
    "rust":       ("Rust",          _rust_cmd),
    "cpp":        ("C++",           _cpp_cmd),
    # Hand-rolled LU variants (no external LA library)
    "julia_lu":   ("Julia-LU",     _julia_lu_cmd),
    "python_lu":  ("Python-LU",    _python_lu_cmd),
    "go_lu":      ("Go-LU",        _go_lu_cmd),
    "rust_lu":    ("Rust-LU",      _rust_lu_cmd),
    "cpp_lu":     ("C++-LU",       _cpp_lu_cmd),
}

NETLIB_PROBLEMS  = ["afiro", "blend", "sc205", "ship04l"]
ALL_PROBLEMS     = NETLIB_PROBLEMS


# ── Helpers ────────────────────────────────────────────────────────────────────

def parse_time_v_rss(stderr: str) -> float | None:
    """Extract peak RSS (KB) from /usr/bin/time -v stderr output."""
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
    return float(m.group(1)) / 1024 if m else None  # return MB


def run_one(cmd: list[str], timeout: int = 300) -> tuple[dict, float | None]:
    """
    Run bench_mps under /usr/bin/time -v.
    Returns (stdout_json, peak_rss_mb).
    """
    time_bin = "/usr/bin/time"
    full_cmd = [time_bin, "-v"] + cmd

    proc = subprocess.run(
        full_cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {proc.returncode}):\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr[-1000:]}"
        )

    result = json.loads(proc.stdout)
    rss_mb = parse_time_v_rss(proc.stderr)
    return result, rss_mb


def load_problem_meta(name: str) -> dict:
    p = DATA / f"{name}.json"
    if not p.exists():
        raise FileNotFoundError(f"{p} not found — run: python scripts/preparse.py")
    with open(p) as f:
        d = json.load(f)
    return d


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--problems", nargs="+", default=ALL_PROBLEMS)
    ap.add_argument("--langs", nargs="+", default=list(LANGUAGES))
    ap.add_argument("--runs", type=int, default=10)
    ap.add_argument("--tol", type=float, default=1e-4,
                    help="Relative tolerance for objective verification")
    args = ap.parse_args()

    rows: list[dict] = []

    for lang_key in args.langs:
        if lang_key not in LANGUAGES:
            print(f"Unknown language: {lang_key}", file=sys.stderr)
            continue
        display, cmd_fn = LANGUAGES[lang_key]

        for prob in args.problems:
            meta = load_problem_meta(prob)
            optimal_z = meta["optimal_z"]
            data_path = DATA / f"{prob}.json"

            print(f"[{display:12s} / {prob:10s}] ", end="", flush=True)

            try:
                cmd = cmd_fn(data_path, args.runs)
                result, rss_mb = run_one(cmd)
            except Exception as e:
                print(f"ERROR: {e}")
                rows.append({"language": display, "problem": prob,
                             "status": "ERROR", "error": str(e)[:80]})
                continue

            times_ms: list[float] = result["times_ms"]
            z: float = result["z"]
            iters: int = result.get("iterations", -1)
            med_ms = statistics.median(times_ms)
            std_ms = statistics.stdev(times_ms) if len(times_ms) > 1 else 0.0

            # Verify objective
            rel_err = abs(z - optimal_z) / max(1.0, abs(optimal_z))
            ok = rel_err <= args.tol
            status = "OK" if ok else f"WRONG_Z(got={z:.6g},exp={optimal_z:.6g})"

            print(f"median={med_ms:8.2f} ms  rss={rss_mb or 0:.1f} MB  "
                  f"iters={iters}  {status}")

            rows.append({
                "language":    display,
                "problem":     prob,
                "status":      status,
                "median_ms":   round(med_ms, 3),
                "std_ms":      round(std_ms, 3),
                "peak_rss_mb": round(rss_mb or 0, 2),
                "iterations":  iters,
                "z":           round(z, 6),
                "optimal_z":   round(optimal_z, 6),
            })

    # Write CSV
    csv_path = RESULTS / "results.csv"
    fieldnames = ["language", "problem", "status", "median_ms", "std_ms",
                  "peak_rss_mb", "iterations", "z", "optimal_z"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nResults written to {csv_path}")

    # Write markdown summary
    _write_markdown(rows, RESULTS / "results.md")
    print(f"Summary written to {RESULTS / 'results.md'}")


def _write_markdown(rows: list[dict], path: Path) -> None:
    problems = list(dict.fromkeys(r["problem"] for r in rows))
    languages = list(dict.fromkeys(r["language"] for r in rows))

    # Build lookup: (lang, prob) -> row
    by = {(r["language"], r["problem"]): r for r in rows}

    with open(path, "w") as f:
        for prob in problems:
            f.write(f"## {prob}\n\n")
            f.write("| Language | Status | Median (ms) | Std (ms) | Peak RSS (MB) | Iterations |\n")
            f.write("|---|---|---|---|---|---|\n")
            for lang in languages:
                r = by.get((lang, prob), {})
                if not r:
                    continue
                med = r.get('median_ms')
                std = r.get('std_ms')
                f.write(
                    f"| {lang} | {r.get('status','')} | "
                    f"{'%.3f' % med if med is not None else 'ERR'} | "
                    f"{'%.3f' % std if std is not None else 'ERR'} | "
                    f"{r.get('peak_rss_mb','')} | {r.get('iterations','')} |\n"
                )
            f.write("\n")


if __name__ == "__main__":
    main()
