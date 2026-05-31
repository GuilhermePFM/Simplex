"""
A linear program in standard equality form:
    maximize  c'x
    subject to Ax = b,  x >= 0
"""
struct LinearProgram
    A :: Matrix{Float64}
    b :: Vector{Float64}
    c :: Vector{Float64}

    function LinearProgram(A::AbstractMatrix, b::AbstractVector, c::AbstractVector)
        m, n = size(A)
        length(b) == m || throw(DimensionMismatch("b has $(length(b)) rows; A has $m"))
        length(c) == n || throw(DimensionMismatch("c has $(length(c)) entries; A has $n columns"))
        new(Float64.(A), Float64.(b), Float64.(c))
    end
end

"""
Index partition of variables into basic (in basis B) and non-basic (fixed at zero).
Both vectors hold column indices into the LP's A matrix.
"""
struct BasisState
    basic    :: Vector{Int}
    nonbasic :: Vector{Int}
end

"""
Termination status of a Simplex solve.
- `OPTIMAL`   : a finite optimal solution was found
- `UNBOUNDED` : the objective is unbounded above; `x` is an extreme ray
- `INFEASIBLE`: the feasible region is empty
"""
@enum SolveStatus begin
    OPTIMAL    =  1
    UNBOUNDED  = -1
    INFEASIBLE = -2
end

"""
Complete output of `solve` or `solve_phase2`.

| Field        | Meaning                                        |
|--------------|------------------------------------------------|
| `x`          | Optimal point, or extreme ray if `UNBOUNDED`   |
| `z`          | Objective value (`Inf` if `UNBOUNDED`)         |
| `status`     | `OPTIMAL`, `UNBOUNDED`, or `INFEASIBLE`        |
| `iterations` | Total pivot count                              |
"""
struct SimplexResult
    x          :: Vector{Float64}
    z          :: Float64
    status     :: SolveStatus
    iterations :: Int
end
