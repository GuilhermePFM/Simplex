# Simplex Benchmark

Cross-language implementation and benchmark of the **two-phase Revised Simplex Method**, comparing R, Python, Julia, C++, Go, and TypeScript.

---

## Goal

Each language implements the same algorithm against the same problem instances. The benchmark measures:

- **Wall-clock time** — median solve time over 10 runs, excluding I/O and JIT warm-up where applicable
- **Memory usage** — peak resident set size (RSS) during the solve

All implementations are allowed to use each language's idiomatic linear-algebra library for the basis solve step (`B \ N`): NumPy, Eigen, BLAS/LAPACK bindings, etc. The simplex logic itself — pivot selection, ratio test, anti-cycling — must be implemented from scratch.

---

## Languages

| Language   | LA library       | Entry point                  |
|------------|------------------|------------------------------|
| Julia      | `LinearAlgebra`  | `julia/Simplex.jl`           |
| Python     | NumPy            | `python/simplex.py`          |
| C++        | Eigen            | `cpp/simplex.cpp`            |
| Go         | Gonum `gonum/mat` (pure Go) | `go/simplex/`           |
| Rust          | base `solve()`   | `rust/simplex.R`                |
| TypeScript | `mathjs`         | `ts/simplex.ts`              |

---

## Dataset — Netlib LP Test Suite

The benchmark uses problems from the [Netlib LP library](https://www.netlib.org/lp/data/), the canonical test suite for LP solvers. Problems are distributed as `.mps` files with known optimal values, making them ideal for both performance measurement and correctness validation.

### Selected problems

| Problem    | Rows  | Columns | Optimal value  | Difficulty  |
|------------|-------|---------|----------------|-------------|
| `afiro`    | 27    | 51      | −464.753       | warm-up     |
| `blend`    | 74    | 114     | −30.812        | small       |
| `sc205`    | 205   | 317     | −52.202        | medium      |
| `ship04l`  | 402   | 2,118   | 1,793,091.56   | large       |

`afiro` and `blend` serve as correctness checks (compare final objective against the known optimum to within 1e-6). `sc205` and `ship04l` are the primary timing targets.

### Obtaining the data

```bash
mkdir -p data
for prob in afiro blend sc205 ship04l; do
  curl -o data/${prob}.mps.gz \
    https://www.netlib.org/lp/data/${prob}
done
```

All four problems are small enough to store in this repository under `data/` if preferred.

---

## Benchmark Methodology

1. **Parse** the `.mps` file into `(A, b, c)` — a shared parser is provided in `scripts/parse_mps.py` so each implementation reads identical matrices.
2. **Warm up** — for JIT-compiled languages (Julia, TypeScript/V8), run one cold solve before timing begins.
3. **Time** — run 10 solves; record median and standard deviation.
4. **Memory** — measure peak RSS using `/usr/bin/time -v` (Linux) or `time.proc_time` / `valgrind --tool=massif` as appropriate.
5. **Verify** — check that the returned objective matches the known optimum to within 1e-6.

Results are written to `results/results.csv` by `scripts/run_benchmark.sh`.

---

## Results

> Results will be populated once all implementations are complete.

| Language   | `sc205` median (ms) | `sc205` peak RSS (MB) | `ship04l` median (ms) | `ship04l` peak RSS (MB) |
|------------|---------------------|-----------------------|-----------------------|-------------------------|
| C++        |                     |                       |                       |                         |
| Julia      |                     |                       |                       |                         |
| Go         |                     |                       |                       |                         |
| Python     |                     |                       |                       |                         |
| TypeScript |                     |                       |                       |                         |
| Rust          |                     |                       |                       |                         |

---

## How the Algorithm Works

All six implementations follow the same specification.

### Core Idea

The feasible region of an LP is a convex polytope. Optimal solutions always occur at a **vertex** (basic feasible solution). The algorithm moves from vertex to vertex along edges, strictly improving the objective, until it reaches the optimum.

![2D feasible polytope: simplex path along edges from a start vertex to the optimum](docs/simplex1.svg)
*Figure: Optimal solutions lie at vertices; each pivot moves along an edge to a neighboring vertex with higher objective (until none improve).*

A **basis** B is a set of m linearly independent columns of A that fully determine a vertex: the basic variables are `x_B = B⁻¹b` and all non-basic variables are zero.

### Phase 1 — Finding a Feasible Starting Point

Before optimizing, we need a basic feasible solution. Phase 1 constructs an auxiliary problem by adding **artificial variables** w:

```
maximize  -sum(w)
subject to [A | I] [x; w] = b,  x,w >= 0
```

Any row with `b[i] < 0` is negated first so that `b ≥ 0` holds, making the all-artificials basis `w = b` immediately feasible. Phase 1 then runs the simplex iterations. If it terminates with objective value 0, the artificials are driven out and we have a BFS for the original problem. Otherwise the problem is **infeasible**.

### Phase 2 — Optimizing

Given the BFS from Phase 1, Phase 2 iterates:

1. **Compute reduced costs** `c_r = c_N' - c_B' B⁻¹ N`
2. **Optimality check** — if all `c_r <= 0`, the current vertex is optimal; stop.
3. **Entering variable** — pick the non-basic variable j with the largest (most positive) reduced cost (`argmax(c_r)`).
4. **Ratio test (leaving variable)** — compute `θ = min { x_B[i] / (B⁻¹ A_j)[i] : (B⁻¹ A_j)[i] > 0 }`. The basic variable achieving the minimum leaves the basis.
   - If no ratio is finite, the problem is **unbounded** in direction d.
5. **Pivot** — swap the entering and leaving variables; update x and the index sets.

### Anti-Cycling: Lexicographic Rule

When the ratio test produces a tie (degenerate pivot), the algorithm could cycle forever. This implementation breaks ties using the **lexicographic smallest** rule: among tied rows, the one that is lexicographically smallest in the scaled non-basic matrix `N / (B⁻¹ A_j)[i]` is chosen, guaranteeing finite termination.

### Return Values

| Value    | Meaning |
|----------|---------|
| `x`      | Optimal point (or extreme ray if unbounded) |
| `z`      | Objective value at optimum |
| `status` | `1` = optimal, `-1` = unbounded, `-2` = iteration limit reached |
| `it`     | Number of pivot iterations |

---

## Complexity

### Worst Case: Exponential

Simplex is **not** a polynomial-time algorithm. The Klee-Minty cube (1972) is a family of LP instances where the standard simplex method (with the largest-coefficient pivot rule) visits all 2ⁿ vertices before finding the optimum. So worst-case complexity is **O(2ⁿ)**.

No pivot rule has ever been proven to make simplex polynomial in the worst case, and finding one remains an open problem.

### LP is in P

Even though simplex itself is not polynomial, **Linear Programming as a problem is in P**:

- **Ellipsoid Method** (Khachiyan, 1979) — first polynomial-time LP algorithm, theoretically important but slow in practice.
- **Interior Point Methods** (Karmarkar, 1984) — polynomial and competitive with simplex on large dense problems.

LP is therefore **not NP-complete**, and not believed to be NP-hard.

### Integer LP is NP-Complete

When variables are restricted to integers (**ILP**), the problem becomes NP-complete. In practice, ILP solvers use simplex heavily — they solve LP *relaxations* (dropping the integrality constraints) inside a branch-and-bound search tree, making fast LP solving critical.

### Why Simplex Dominates in Practice: Smoothed Complexity

The gap between exponential worst case and fast practical performance was formally explained by Spielman & Teng (2004) with the concept of **smoothed complexity**: if you add an arbitrarily small random perturbation to any problem instance, the expected number of pivot steps becomes polynomial. Real-world data is never adversarially constructed, so simplex behaves well on it.

### Open Problem: Strongly Polynomial LP

Whether LP can be solved in time polynomial in m and n alone — independent of the magnitude of the coefficients — is still **unknown**. This is one of Smale's Millennium problems for mathematics and theoretical CS.
