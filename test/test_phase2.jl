@testset "Phase 2 — solve_phase2" begin

    @testset "optimal: production problem" begin
        # maximize 4x1 + 3x2  s.t. 2x1+x2<=4, x1+2x2<=4, x1,x2>=0
        # Standard form (slacks x3, x4): basis starts at {x3, x4}
        lp    = LinearProgram(
            Float64[2 1 1 0; 1 2 0 1],
            Float64[4; 4],
            Float64[4; 3; 0; 0])
        basis = BasisState([3, 4], [1, 2])

        r = solve_phase2(lp, basis)

        @test r.status == OPTIMAL
        @test r.z      ≈ 28/3  atol=1e-8
        @test r.x[1]   ≈ 4/3   atol=1e-8
        @test r.x[2]   ≈ 4/3   atol=1e-8
        @test r.x[3]   ≈ 0.0   atol=1e-8
        @test r.x[4]   ≈ 0.0   atol=1e-8
        @test r.iterations == 3   # 2 pivots + 1 final optimality check
    end

    @testset "optimal: already at vertex" begin
        # x1=1, x2=1 is the only feasible point; objective is c'x=5
        lp    = LinearProgram(
            Float64[1 0 1 0; 0 1 0 1],
            Float64[1; 1],
            Float64[3; 2; 0; 0])
        basis = BasisState([1, 2], [3, 4])

        r = solve_phase2(lp, basis)

        @test r.status == OPTIMAL
        @test r.z      ≈ 5.0  atol=1e-8
        @test r.x[1]   ≈ 1.0  atol=1e-8
        @test r.x[2]   ≈ 1.0  atol=1e-8
    end

    @testset "unbounded" begin
        # maximize x1+x2  s.t.  0.5x1 - x2 + x3 = 0.5,  -4x1 + x2 + x4 = 1
        # Feasible but unbounded: increasing x1 and x2 together is unconstrained
        lp    = LinearProgram(
            Float64[0.5 -1 1 0; -4 1 0 1],
            Float64[0.5; 1],
            Float64[1; 1; 0; 0])
        basis = BasisState([3, 4], [1, 2])

        r = solve_phase2(lp, basis)

        @test r.status == UNBOUNDED
        @test r.z      == Inf
    end

    @testset "bland pivot rule gives same optimal" begin
        lp    = LinearProgram(
            Float64[2 1 1 0; 1 2 0 1],
            Float64[4; 4],
            Float64[4; 3; 0; 0])
        basis = BasisState([3, 4], [1, 2])

        r = solve_phase2(lp, basis; rule=Bland())

        @test r.status == OPTIMAL
        @test r.z      ≈ 28/3  atol=1e-8
    end

end
