"""
    Simplex

Two-phase Revised Simplex Method for Linear Programming.

## Quick start
```julia
using Simplex

lp = LinearProgram(
    Float64[2 1 1 0;
            1 2 0 1],
    Float64[4; 4],
    Float64[4; 3; 0; 0]
)

result = solve(lp)
# SimplexResult(x=[4/3, 4/3, 0, 0], z=9.33, status=OPTIMAL, iterations=2)
```

## API surface
- `LinearProgram(A, b, c)` — construct an LP in standard equality form
- `solve(lp; rule, verbosity, logfile)` — two-phase solve (Phase 1 + Phase 2)
- `solve_phase2(lp, basis; rule, verbosity)` — Phase 2 only (needs a known BFS)
- `SimplexResult` — returned by both solve functions
- `SolveStatus`: `OPTIMAL`, `UNBOUNDED`, `INFEASIBLE`
- `LargestCoefficient`, `Bland` — pivot rules
- `LogLevel`: `SILENT`, `INFO`, `DEBUG`
"""
module Simplex

using LinearAlgebra

include("types.jl")
include("logging.jl")
include("pivot.jl")
include("preprocess.jl")
include("phase1.jl")
include("phase2.jl")

export LinearProgram, BasisState, SimplexResult
export SolveStatus, OPTIMAL, UNBOUNDED, INFEASIBLE
export PivotRule, LargestCoefficient, Bland
export LogLevel, SILENT, INFO, DEBUG
export SimplexLogger
export solve, solve_phase2

"""
    solve(lp; rule, verbosity, logfile) -> SimplexResult

Solve a linear program using the two-phase Revised Simplex Method.

## Keyword arguments
| Argument     | Default               | Meaning                                  |
|--------------|-----------------------|------------------------------------------|
| `rule`       | `LargestCoefficient()`| Pivot rule: `LargestCoefficient`, `Bland`|
| `verbosity`  | `SILENT`              | `SILENT`, `INFO`, or `DEBUG` to stdout   |
| `logfile`    | `nothing`             | Path for a full debug log file           |
"""
function solve(lp::LinearProgram;
               rule::PivotRule       = LargestCoefficient(),
               verbosity::LogLevel   = SILENT,
               logfile               = nothing) :: SimplexResult

    logger = isnothing(logfile) ? SimplexLogger(verbosity) : SimplexLogger(logfile)

    log_problem(logger, lp)

    basis, p1_status = phase1(lp, rule, logger)

    if p1_status == INFEASIBLE
        close_log(logger)
        return SimplexResult(zeros(size(lp.A, 2)), 0.0, INFEASIBLE, 0)
    end

    result = phase2(lp, basis, rule, logger)
    log_result(logger, result)
    close_log(logger)
    return result
end

"""
    solve_phase2(lp, basis; rule, verbosity) -> SimplexResult

Run Phase 2 only, starting from a known basic feasible solution `basis`.
Useful when you already have a BFS (e.g., when the identity columns of A
form a valid starting basis).
"""
function solve_phase2(lp::LinearProgram, init_basis::BasisState;
                       rule::PivotRule     = LargestCoefficient(),
                       verbosity::LogLevel = SILENT) :: SimplexResult
    logger = SimplexLogger(verbosity)
    result = phase2(lp, init_basis, rule, logger)
    log_result(logger, result)
    return result
end

end # module Simplex
