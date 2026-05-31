# Simplex

Revised Simplex Method implementation in Julia.

---

## What is the Simplex Algorithm?

The Simplex Algorithm is the classical method for solving **Linear Programming (LP)** problems — optimization problems of the form:

```
maximize   c'x
subject to Ax = b
           x >= 0
```

It is used wherever a linear objective must be optimized under linear constraints, including:

- **Production planning** — maximize profit given resource limits
- **Transportation & logistics** — minimize shipping cost across routes
- **Network flow** — route commodities optimally through a graph
- **Finance** — portfolio allocation under budget and risk constraints
- **Scheduling** — assign jobs to machines to minimize makespan

---

## Why is it Good?

Despite its worst-case exponential complexity, the Simplex Method is the workhorse of practical LP solvers because:

- **Empirically fast** — on real-world problems it takes O(m) to O(2m) pivot steps for an m-constraint problem, making it polynomial in practice.
- **Exact solutions** — it walks along vertices of the feasible polytope and terminates at a certifiably optimal vertex, unlike gradient-based methods that only converge asymptotically.
- **Rich termination info** — it correctly identifies three distinct outcomes: *optimal*, *unbounded*, and *infeasible*, and returns a certificate for each.
- **Warm-starting** — an existing basis can be reused when problem data changes slightly, which is critical in branch-and-bound solvers for integer programming.
- **Numerical stability** — the revised form re-solves the basis system `B x = N` at each pivot using an LU factorization, rather than carrying the full explicit tableau, keeping the working matrix small and better-conditioned.

---

## How it Works

This implementation uses the **two-phase Revised Simplex Method** with **lexicographic anti-cycling**.

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

The Simplex Algorithm has one of the most interesting complexity stories in computer science.

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

---

## Usage

```julia
include("Simplex.jl")

A = Float64[2 1 1 0;
             1 2 0 1]
b = Float64[4; 4]
c = Float64[4; 3; 0; 0]

x, z, status, it = Simplex(A, b, c)
```

If the initial basis is already known to be feasible, `SimplexFase2` can be called directly to skip Phase 1.
