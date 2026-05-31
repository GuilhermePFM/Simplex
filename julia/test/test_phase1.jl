@testset "Phase 1 — full solve" begin

    @testset "feasible: production problem with Phase 1" begin
        # Same as the direct Phase 2 test, but now going through solve() so Phase 1 runs
        lp = LinearProgram(
            Float64[2 1 1 0; 1 2 0 1],
            Float64[4; 4],
            Float64[4; 3; 0; 0])

        r = solve(lp)

        @test r.status == OPTIMAL
        @test r.z      ≈ 28/3  atol=1e-8
        @test r.x[1]   ≈ 4/3   atol=1e-8
        @test r.x[2]   ≈ 4/3   atol=1e-8
    end

    @testset "feasible: problem with a negative-b row" begin
        # Constraint row 3 has b=-1; make_b_nonnegative flips it before Phase 1
        # Original: 2x1+x2+x3=4, x1+2x2+x4=4, -x1-x2+x5=-1 (i.e., x1+x2-x5=1)
        lp = LinearProgram(
            Float64[2 1 1 0 0;
                    1 2 0 1 0;
                   -1 -1 0 0 1],
            Float64[4; 4; -1],
            Float64[4; 3; 0; 0; 0])

        r = solve(lp)

        @test r.status == OPTIMAL
        @test r.z > 0
        # Feasibility check: Ax ≈ b
        @test lp.A * r.x ≈ lp.b  atol=1e-6
        @test all(r.x .>= -1e-10)
    end

    @testset "infeasible: contradictory constraints" begin
        # x1 + x2 = 1  AND  x1 + x2 = 0  — impossible
        lp = LinearProgram(
            Float64[1 1; 1 1],
            Float64[1; 0],
            Float64[1; 1])

        r = solve(lp)

        @test r.status == INFEASIBLE
    end

    @testset "infeasible: conflicting bounds x1<=1 and x1>=2" begin
        # x1 + s1 = 1  (x1 <= 1)
        # x1 - s2 = 2  (x1 >= 2, written as -x1 + s2 = -2, negated by preprocess)
        # No x1 >= 0 can satisfy both simultaneously.
        lp = LinearProgram(
            Float64[1 1 0; -1 0 1],
            Float64[1; -2],
            Float64[1; 0; 0])

        r = solve(lp)

        @test r.status == INFEASIBLE
    end

end
