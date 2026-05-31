function phase1(lp::LinearProgram, rule::PivotRule, logger::SimplexLogger;
                maxit::Int = 1000, tol::Float64 = 1e-8) :: Tuple{BasisState, SolveStatus}

    lp  = make_b_nonnegative(lp)
    m, n = size(lp.A)

    Aw = [lp.A  Matrix{Float64}(I, m, m)]
    cw = [zeros(n); -ones(m)]
    art_start = n + 1

    basic    = collect(art_start : n + m)
    nonbasic = collect(1:n)

    log_phase(logger, 1)

    for it in 1:maxit
        B     = Aw[:, basic]
        N     = Aw[:, nonbasic]
        BF    = lu(B)
        xb    = BF \ lp.b
        BinvN = BF \ N

        cb = cw[basic]
        cr = cw[nonbasic] .- BinvN' * cb
        z  = dot(cb, xb)

        log_iteration(logger, it, BasisState(copy(basic), copy(nonbasic)), xb, z)

        j = entering_index(rule, cr; tol)

        if isnothing(j)
            z >= -tol || return BasisState(Int[], Int[]), INFEASIBLE

            fix_artificials_in_basis!(Aw, basic, nonbasic, art_start)
            orig_basic    = filter(i -> i <= n, basic)
            orig_nonbasic = filter(i -> i <= n, nonbasic)
            return BasisState(orig_basic, orig_nonbasic), OPTIMAL
        end

        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN; tol)

        # Phase 1 auxiliary problem is always bounded (artificials provide a floor)
        isnothing(i) && error("Phase 1 is unexpectedly unbounded")

        basic[i], nonbasic[j] = nonbasic[j], basic[i]
    end

    error("Phase 1 did not converge in $maxit iterations")
end


function fix_artificials_in_basis!(Aw::Matrix{Float64},
                                    basic::Vector{Int},
                                    nonbasic::Vector{Int},
                                    art_start::Int)
    changed = true
    while changed
        changed = false
        for (pos, col) in enumerate(basic)
            col < art_start && continue

            BinvAw = lu(Aw[:, basic]) \ Aw
            row    = BinvAw[pos, :]

            swap_j = findfirst(1:length(nonbasic)) do j
                ncol = nonbasic[j]
                ncol < art_start && abs(row[ncol]) > 1e-10
            end

            if isnothing(swap_j)
                Aw[pos, :] .= 0.0
            else
                basic[pos], nonbasic[swap_j] = nonbasic[swap_j], basic[pos]
                changed = true
                break
            end
        end
    end
end
