using Test
using CircStats

@testset "All tests" begin
    include("parameters.jl")
    include("hypothesistests.jl")
end
