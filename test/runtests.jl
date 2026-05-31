using Test

include("../src/Simplex.jl")
using .Simplex

@testset "Simplex" begin
    include("test_phase2.jl")
    include("test_phase1.jl")
    include("test_anticycling.jl")
end
