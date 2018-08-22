# Parameter estimation

using CircStats
using Compat.Test

@testset "Distance" begin
    @test cdist(0, 0) == 0
    @test cdist(0, 1) == 1
    @test cdist(-1, 0) == 1
    @test cdist(0, -1) == -1
    @test cdist(2.1, 2.0) ≈ -0.1
    @test cdist(340.1, 0.1, true) ≈ 20.0
    let (a, b) = rand(2)
        @test cdist(a, b, false) == cdist(a, b)
        @test cdist(a, b) ≈ deg2rad(cdist(rad2deg(a), rad2deg(b), true))
        @test cdist(a, b, true) ≈ rad2deg(cdist(deg2rad(a), deg2rad(b)))
    end
end

@testset "Parameter est" begin
    # Example 1.1 of Mardia & Jupp (2000)
    let data = [43, 45, 52, 61, 75, 88, 88, 279, 357], μ = 51, R̄ = 0.711, med = 52
        @test cmean(data, true) ≈ μ atol=0.1
        @test cresultant(data, true) ≈ R̄ atol=0.001
        @test cmedian(data, true) ≈ med
        @test cvariance(data, true) ≈ 1 - R̄ atol=0.001
    end
end