function phase2(lp::LinearProgram, init_basis::BasisState, rule::PivotRule,
                logger::SimplexLogger;
                maxit::Int = 1000, tol::Float64 = 1e-10) :: SimplexResult

    A, b, c = lp.A, lp.b, lp.c
    n       = size(A, 2)

    basic    = copy(init_basis.basic)
    nonbasic = copy(init_basis.nonbasic)

    x         = zeros(n)
    x[basic]  = A[:, basic] \ b

    log_phase(logger, 2)

    for it in 1:maxit
        B     = A[:, basic]
        N     = A[:, nonbasic]
        BF    = lu(B)
        xb    = BF \ b
        BinvN = BF \ N

        cb = c[basic]
        cr = c[nonbasic] .- BinvN' * cb
        z  = dot(cb, xb)

        log_iteration(logger, it, BasisState(copy(basic), copy(nonbasic)), xb, z)

        j = entering_index(rule, cr; tol)
        if isnothing(j)
            x          .= 0.0
            x[basic]    = xb
            return SimplexResult(x, z, OPTIMAL, it)
        end

        d = BinvN[:, j]
        i = leaving_index(xb, d, BinvN; tol)

        if isnothing(i)
            ray          = zeros(n)
            ray[nonbasic[j]] = 1.0
            ray[basic]       = -d
            return SimplexResult(ray, Inf, UNBOUNDED, it)
        end

        theta = xb[i] / d[i]

        x          .= 0.0
        x[basic]    = xb .- theta .* d
        x[nonbasic[j]] = theta
        x[basic[i]]    = 0.0  # numerical cleanup for the leaving variable

        basic[i], nonbasic[j] = nonbasic[j], basic[i]
    end

    error("Phase 2 did not converge in $maxit iterations — possible cycling")
end
