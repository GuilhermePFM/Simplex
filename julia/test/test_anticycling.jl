@testset "Anti-cycling — degenerate tie-breaking" begin

    # maximize x1 + x2
    # s.t.  x1 + x2 <= 2
    #       x1      <= 1
    #            x2 <= 1
    #
    # Standard form: A = [[1,1,1,0,0],[1,0,0,1,0],[0,1,0,0,1]], b=[2,1,1]
    # The optimal vertex x1=x2=1 lies at the intersection of ALL THREE constraints
    # (s1=s2=s3=0), making it a highly degenerate vertex.  The ratio test ties
    # at iteration 2 (both s1 and s3 have ratio xb/d = 1), and the lexicographic
    # rule must break the tie to guarantee termination.

    A = Float64[1 1 1 0 0;
                1 0 0 1 0;
                0 1 0 0 1]
    b = Float64[2; 1; 1]
    c = Float64[1; 1; 0; 0; 0]

    lp    = LinearProgram(A, b, c)
    basis = BasisState([3, 4, 5], [1, 2])

    @testset "LargestCoefficient — terminates at tied ratio" begin
        r = solve_phase2(lp, basis; rule=LargestCoefficient())

        @test r.status     == OPTIMAL
        @test r.z          ≈ 2.0  atol=1e-8
        @test r.x[1]       ≈ 1.0  atol=1e-8
        @test r.x[2]       ≈ 1.0  atol=1e-8
        @test r.iterations == 3           # 2 pivots + final optimality check
        # Feasibility certificate
        @test A * r.x ≈ b      atol=1e-8
        @test all(r.x .>= -1e-10)
    end

    @testset "Bland — same result" begin
        r = solve_phase2(lp, basis; rule=Bland())

        @test r.status == OPTIMAL
        @test r.z      ≈ 2.0  atol=1e-8
        @test r.x[1]   ≈ 1.0  atol=1e-8
        @test r.x[2]   ≈ 1.0  atol=1e-8
    end

end

@testset "Anti-cycling — degenerate start (Kotiah-Steinberg)" begin

    # maximize 10x1 - 57x2 - 9x3 - 24x4
    # s.t.
    #   -0.5x1 + 5.5x2 + 2.5x3 - 9x4 + x5 = 0
    #   -0.5x1 + 2.5x2 + 0.5x3 -  x4 + x6 = 0
    #    x1                              + x7 = 1
    #
    # Initial BFS: x5=x6=0, x7=1 — degenerate (two basic variables at zero).
    # Optimal: x1=1, z=10.

    A = Float64[
        -0.5   5.5   2.5  -9   1  0  0;
        -0.5   2.5   0.5  -1   0  1  0;
         1.0   0.0   0.0   0   0  0  1
    ]
    b = Float64[0; 0; 1]
    c = Float64[10; -57; -9; -24; 0; 0; 0]

    lp    = LinearProgram(A, b, c)
    basis = BasisState([5, 6, 7], [1, 2, 3, 4])

    r = solve_phase2(lp, basis; rule=LargestCoefficient())

    @test r.status == OPTIMAL
    @test r.z      ≈ 10.0  atol=1e-8
    @test r.x[1]   ≈ 1.0   atol=1e-8
    @test A * r.x  ≈ b     atol=1e-8
end
