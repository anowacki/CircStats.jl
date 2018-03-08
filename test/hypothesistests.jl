# Significance tests
using CircStats
using Base.Test

@testset "Tests" begin
    # Watson's U²n test for goodness-of-fit to a distribution
    # Test 96 of Kanji (2006) 100 Statistical Tests
    let θ = [20, 135, 145, 165, 170, 200, 300, 325, 335, 350, 350, 350, 355]
        fit_is_significant, U², U²crit = watson_U2n(θ, x->x/2π, 0.05, true)
        @test !fit_is_significant
        @test U² ≈ 0.1361 atol=0.0001
        @test U²crit ≈ 0.184 atol=0.001
    end
    
    # V-test (modified Rayleigh) to test if a distribution has a preferred direction
    let θ = []
    end
end
