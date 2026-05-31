using BenchmarkTools

include("../src/Simplex.jl")
using .Simplex

"""Generate a random LP with m constraints and n total variables (n > m)."""
function random_lp(m::Int, n::Int; seed::Int=42)
    rng = MersenneTwister(seed)
    A_raw = rand(rng, m, n - m)
    A = [A_raw Matrix{Float64}(I, m, m)]   # append slack columns
    b = abs.(rand(rng, m)) .+ 0.1            # strictly positive RHS
    c = [rand(rng, n - m); zeros(m)]         # only original vars in objective
    LinearProgram(A, b, c)
end

"""Build benchmark random suite"""
function build_benchmark_suite()
    suite = BenchmarkGroup()

    # Small canonical problem (used in tests)
    lp_small = LinearProgram(
        Float64[2 1 1 0; 1 2 0 1],
        Float64[4; 4],
        Float64[4; 3; 0; 0])
    basis_small = BasisState([3, 4], [1, 2])

    suite["small/phase2_only"] = @benchmarkable solve_phase2($lp_small, $basis_small)
    suite["small/two_phase"] = @benchmarkable solve($lp_small)

    # Medium random LPs
    for (m, n) in [(10, 20), (50, 100), (100, 200)]
        lp = random_lp(m, n)
        tag = "m=$(m)_n=$(n)"
        suite["random/$tag"] = @benchmarkable solve($lp)
    end
    return suite
end

function run_benchmark_suite(suite)
    println("Tuning benchmarks…")
    tune!(suite)
    println("Running benchmarks…")
    results = run(suite, verbose=true)

    println("\n=== Results ===")
    for (name, trial) in leaves(results)
        println("$name: median=$(BenchmarkTools.median(trial))")
    end
end

run_benchmark_suite(build_benchmark_suite())