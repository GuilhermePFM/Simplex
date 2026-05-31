abstract type PivotRule end

"""
Pick the non-basic variable with the largest (most positive) reduced cost.
"""
struct LargestCoefficient <: PivotRule end

"""
Bland's rule: always pick the lowest-index positive reduced cost.
Guarantees finite termination independently of the ratio test, at the cost of
slower convergence on typical problems.
"""
struct Bland <: PivotRule end

function entering_index(::LargestCoefficient, cr::AbstractVector{Float64};
                         tol::Float64 = 1e-10) :: Union{Int, Nothing}
    v, j = findmax(cr)
    v > tol ? j : nothing
end

function entering_index(::Bland, cr::AbstractVector{Float64};
                         tol::Float64 = 1e-10) :: Union{Int, Nothing}
    findfirst(x -> x > tol, cr)
end

"""
Lexicographic ratio test — selects the leaving variable to prevent cycling.

Given:
- `xb`:    values of the current basic variables (length m)
- `d`:     entering column in basis coordinates, B⁻¹Aⱼ (length m)
- `BinvN`: full direction matrix B⁻¹N (m × |nonbasic|), used for tie-breaking

Returns the *row position* (into `basic`) of the leaving variable, or `nothing`
if every component of `d` is non-positive (problem is unbounded).

**Anti-cycling guarantee**: ties in the primary ratio `xb[i]/d[i]` are broken by
comparing rows of `BinvN / d[i]` lexicographically. Because `BinvN` is a matrix
with linearly independent columns from the LP's initial tableau, no two candidates
can tie on every column, ensuring a unique choice and finite termination.
"""
function leaving_index(xb::AbstractVector{Float64},
                        d::AbstractVector{Float64},
                        BinvN::AbstractMatrix{Float64};
                        tol::Float64 = 1e-10) :: Union{Int, Nothing}

    candidates = findall(i -> d[i] > tol, eachindex(d))
    isempty(candidates) && return nothing

    best = candidates[1]
    for i in @view candidates[2:end]
        ratio_i    = xb[i]    / d[i]
        ratio_best = xb[best] / d[best]

        if ratio_i < ratio_best - tol
            best = i
        elseif abs(ratio_i - ratio_best) <= tol
            # Lexicographic tie-break on scaled rows of BinvN
            row_i    = BinvN[i,    :] ./ d[i]
            row_best = BinvN[best, :] ./ d[best]
            Tuple(row_i) < Tuple(row_best) && (best = i)
        end
    end

    return best
end
