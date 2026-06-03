# Usage: julia --project=../julia bench_mps.jl <data.json> [N=10]

include("../src/Simplex.jl")
using .Simplex
using JSON3

function main()
    length(ARGS) >= 1 || error("Usage: bench_mps.jl <data.json> [N=10]")
    data_path = ARGS[1]
    N = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10

    data = JSON3.read(read(data_path, String))

    A = reduce(hcat, [Float64.(collect(row)) for row in data["A"]])'  # m×n
    b = Float64.(collect(data["b"]))
    c = Float64.(collect(data["c"]))

    lp = LinearProgram(A, b, c)

    # Warmup solves (Julia JIT)
    solve(lp)
    solve(lp)

    times_ms = Vector{Float64}(undef, N)
    local result
    for i in 1:N
        t0 = time_ns()
        result = solve(lp)
        times_ms[i] = (time_ns() - t0) / 1e6
    end

    z = result.status == OPTIMAL ? result.z : nothing
    iters = result.iterations

    times_str = join(times_ms, ", ")
    if isnothing(z)
        print("{\"z\": null, \"iterations\": $iters, \"times_ms\": [$times_str]}")
    else
        print("{\"z\": $z, \"iterations\": $iters, \"times_ms\": [$times_str]}")
    end
end

main()
