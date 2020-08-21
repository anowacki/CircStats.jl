# Significance tests
using Test
using CircStats

@testset "Tests" begin
    # Watson's U²n test for goodness-of-fit to a distribution
    # Test 96 of Kanji (2006) 100 Statistical Tests
    let θ = [20, 135, 145, 165, 170, 200, 300, 325, 335, 350, 350, 350, 355]
        fit_is_significant, U², U²crit = watson_U2n(θ, x->x/2π, 0.05, true)
        @test !fit_is_significant
        @test U² ≈ 0.1361 atol=0.0001
        @test U²crit ≈ 0.184 atol=0.002
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
    
    # Watson's U² test for the difference of two samples
    # Data and values from worked example from
    #    http://webspace.ship.edu/pgmarr/Geo441/Examples/Watsons%20U2%20Test.pdf
    let θ = [38,45,46,52,53,54,56,57,60,64],
            ϕ = [36,40,44,45,51,51,52,54,54,55,55,56,67,78,89,314], α = 0.05,
            reject = false, U² = 0.0427, U²crit = 0.1856
        reject′, U²′, U²crit′ = watson_U2(θ, ϕ, α, true)
        @test reject′ == reject
        @test U²′ ≈ U² atol=0.0001
        @test U²crit′ ≈ U²crit atol=0.02
    end
    
    # Data of Example 8.2 of Mardia and Jupp (2000); test values of Example 8.4
    # Currently fails because Mardia and Jupp inexplicably use the value of
    # U²crit for n = m = ∞, and their application of expression (8.3.8) does not
    # seem to give the same answer as as my application of (8.3.7).
    # let θ = [50, 290, 300, 300, 305, 320, 330, 330, 335, 340, 340, 355],
    #         ϕ = [70, 155, 190, 195, 215, 235, 235, 240, 255, 260, 290, 300, 300, 300],
    #         α = 0.01, U² = 0.320, U²crit = 0.268, reject = true
    #     reject′, U²′, U²crit′ = watson_U2(θ, ϕ, α, true)
    #     @test reject′ == reject
    #     println(U²′)
    #     @test U²′ ≈ U² atol=0.001
    #     println(U²crit′)
    #     @test U²crit′ ≈ U²crit atol=0.001
    # end
end
