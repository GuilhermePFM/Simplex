# Usage: julia --project=../julia bench_mps_lu.jl <data.json> [N=10]
#
# Two-phase Revised Simplex with hand-rolled LU factorization.
# No LinearAlgebra.lu or \ — plain Julia array operations only.

using JSON3

# ---------------------------------------------------------------------------
# LU factorization with partial pivoting (Doolittle, in-place)
# ---------------------------------------------------------------------------

function lu_factor!(M::Matrix{Float64})
    m = size(M, 1)
    perm = collect(1:m)
    for k in 1:m
        col_k = view(M, k:m, k)
        r = argmax(abs.(col_k)) + k - 1
        if r != k
            M[[k, r], :] = M[[r, k], :]
            perm[k], perm[r] = perm[r], perm[k]
        end
        abs(M[k, k]) < 1e-14 && continue
        M[k+1:m, k] ./= M[k, k]
        M[k+1:m, k+1:m] .-= M[k+1:m, k] * M[k, k+1:m]'
    end
    return perm
end

function lu_solve_vec(M::Matrix{Float64}, perm::Vector{Int}, b::Vector{Float64})
    m = length(b)
    y = b[perm]
    for i in 2:m
        y[i] -= dot(M[i, 1:i-1], y[1:i-1])
    end
    x = copy(y)
    for i in m:-1:1
        if abs(M[i, i]) > 1e-14
            x[i] = (x[i] - dot(M[i, i+1:m], x[i+1:m])) / M[i, i]
        end
    end
    return x
end

function lu_solve_mat(M::Matrix{Float64}, perm::Vector{Int}, N::Matrix{Float64})
    result = similar(N)
    for j in 1:size(N, 2)
        result[:, j] = lu_solve_vec(M, perm, N[:, j])
    end
    return result
end

# ---------------------------------------------------------------------------
# Pivot helpers
# ---------------------------------------------------------------------------

function entering_index(cr::Vector{Float64}; tol=1e-10)
    j = argmax(cr)
    return cr[j] > tol ? j : nothing
end

function leaving_index(xb, d, BinvN; tol=1e-10)
    m = length(d)
    cands = [i for i in 1:m if d[i] > tol]
    isempty(cands) && return nothing
    best = cands[1]
    for i in cands[2:end]
        ri = xb[i] / d[i]; rb = xb[best] / d[best]
        if ri < rb - tol
            best = i
        elseif abs(ri - rb) <= tol
            row_i = BinvN[i, :] ./ d[i]
            row_b = BinvN[best, :] ./ d[best]
            if row_i < row_b
                best = i
            end
        end
    end
    return best
end

# ---------------------------------------------------------------------------
# Phase 1
# ---------------------------------------------------------------------------

function phase1(A::Matrix{Float64}, b::Vector{Float64}, n::Int;
                maxit=10000, tol=1e-8)
    mask = b .< 0
    if any(mask)
        A = copy(A); b = copy(b)
        A[mask, :] .*= -1; b[mask] .*= -1
    end
    m = size(A, 1)
    Aw = hcat(A, Matrix{Float64}(I, m, m))
    cw = vcat(zeros(n), -ones(m))
    art = n + 1  # 1-based: artificials are columns n+1..n+m

    basic = collect(n+1:n+m)
    nonbasic = collect(1:n)

    for _ in 1:maxit
        B = Aw[:, basic]; N = Aw[:, nonbasic]
        LU = copy(B); perm = lu_factor!(LU)
        xb = lu_solve_vec(LU, perm, b)
        BinvN = lu_solve_mat(LU, perm, N)
        cb = cw[basic]; cr = cw[nonbasic] .- BinvN' * cb
        z = dot(cb, xb)
        j = entering_index(cr; tol)
        if isnothing(j)
            z < -tol && return nothing, nothing, :infeasible
            fix_arts!(Aw, basic, nonbasic, art)
            return filter(i -> i < art, basic), filter(i -> i < art, nonbasic), :optimal
        end
        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN; tol)
        isnothing(i) && error("Phase 1 unexpectedly unbounded")
        basic[i], nonbasic[j] = nonbasic[j], basic[i]
    end
    error("Phase 1 did not converge")
end

function fix_arts!(Aw, basic, nonbasic, art)
    changed = true
    while changed
        changed = false
        for (pos, col) in enumerate(basic)
            col < art && continue
            LU = copy(Aw[:, basic]); perm = lu_factor!(LU)
            BinvAw = lu_solve_mat(LU, perm, Aw)
            row = BinvAw[pos, :]
            swap_j = findfirst(j -> nonbasic[j] < art && abs(row[nonbasic[j]]) > 1e-10,
                               1:length(nonbasic))
            if isnothing(swap_j)
                Aw[pos, :] .= 0.0
            else
                basic[pos], nonbasic[swap_j] = nonbasic[swap_j], basic[pos]
                changed = true; break
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Phase 2
# ---------------------------------------------------------------------------

function phase2(A, b, c, basic, nonbasic; maxit=10000, tol=1e-10)
    basic = copy(basic); nonbasic = copy(nonbasic)
    n = size(A, 2)
    for it in 1:maxit
        B = A[:, basic]; N = A[:, nonbasic]
        LU = copy(B); perm = lu_factor!(LU)
        xb = lu_solve_vec(LU, perm, b)
        BinvN = lu_solve_mat(LU, perm, N)
        cb = c[basic]; cr = c[nonbasic] .- BinvN' * cb
        z = dot(cb, xb)
        j = entering_index(cr; tol)
        isnothing(j) && return z, it, :optimal
        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN; tol)
        isnothing(i) && return Inf, it, :unbounded
        basic[i], nonbasic[j] = nonbasic[j], basic[i]
    end
    error("Phase 2 did not converge")
end

function simplex_solve(A, b, c)
    n = size(A, 2)
    bas, nonbas, st = phase1(A, b, n)
    st == :infeasible && return nothing, 0
    all_cols = Set(1:n)
    nonbas = vcat(nonbas, sort(collect(setdiff(all_cols, Set(bas), Set(nonbas)))))
    z, iters, st2 = phase2(A, b, c, bas, nonbas)
    return st2 == :optimal ? z : nothing, iters
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

function main()
    length(ARGS) >= 1 || error("Usage: bench_mps_lu.jl <data.json> [N=10]")
    data_path = ARGS[1]
    N = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10

    data = JSON3.read(read(data_path, String))
    A = Matrix{Float64}(reduce(hcat, [Float64.(collect(row)) for row in data["A"]])')
    b = Float64.(collect(data["b"]))
    c = Float64.(collect(data["c"]))

    # Warmup (JIT)
    simplex_solve(A, b, c); simplex_solve(A, b, c)

    times_ms = Vector{Float64}(undef, N)
    local z_result, iters_result
    for i in 1:N
        t0 = time_ns()
        z_result, iters_result = simplex_solve(A, b, c)
        times_ms[i] = (time_ns() - t0) / 1e6
    end

    times_str = join(times_ms, ", ")
    if isnothing(z_result)
        print("{\"z\": null, \"iterations\": $iters_result, \"times_ms\": [$times_str]}")
    else
        print("{\"z\": $z_result, \"iterations\": $iters_result, \"times_ms\": [$times_str]}")
    end
end

# Need I for identity matrix
using LinearAlgebra: I, dot
main()
