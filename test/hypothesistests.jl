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
    # Test 95 of Kanji (2006)
    let θ = [250, 275, 285, 285, 290, 290, 295, 300, 305, 310, 315, 320, 330, 330, 5],
            θ₀ = 265, α = 0.01, V = 3.884, Vcrit = 2.302, is_not_random = true
        is_not_random′, V′, Vcrit′ = V_test(θ, θ₀, α, true)
        @test is_not_random′ == is_not_random
        @test V′ ≈ V atol=0.1
        @test Vcrit′ ≈ Vcrit atol=0.01
    end
end
